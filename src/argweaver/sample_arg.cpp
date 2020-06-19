//=============================================================================
// sample full ARGs
//

// c++ includes
#include <vector>

// arghmm includes
#include "common.h"
#include "local_tree.h"
#include "logging.h"
#include "model.h"
#include "sample_arg.h"
#include "sample_thread.h"
#include "sequences.h"
#include "total_prob.h"


namespace argweaver {

using namespace std;


// sequentially sample an ARG from scratch
// sequences are sampled in the order given unless random is true
void sample_arg_seq(const ArgModel *model, Sequences *sequences,
                    LocalTrees *trees, bool random, int num_buildup)
{
    const int nseqs = sequences->get_num_seqs();
    const int seqlen = sequences->length();

    vector<int> seqids;
    for (int i=0; i<nseqs; i++)
        seqids.push_back(i);
    if (random)
        shuffle(&seqids[0], seqids.size());

    if (trees->get_num_leaves() == 0) {
        // initialize ARG as trunk
        const int capacity = 2 * sequences->get_num_seqs() - 1;
        int start = trees->start_coord;
        int end = trees->end_coord;
        if (end != seqlen) {
            start = 0;
            end = seqlen;
        }
        int pop_path = ( model->pop_tree == NULL ? 0 :
                         model->pop_tree->most_likely_path(sequences->pops[seqids[0]]) );
        trees->make_trunk(start, end, seqids[0], pop_path, capacity);
    }
    assert_trees(trees, model->pop_tree);

    // record which sequences are already in the tree
    bool has_sequence[nseqs];
    
    fill(has_sequence, has_sequence+nseqs, false);
    for (int i=0; i<trees->get_num_leaves(); i++)
        has_sequence[trees->seqids[i]] = true;

    // add more chromosomes one by one
    for (int i=0; i<nseqs; i++) {
        int new_chrom = seqids[i];
        if (!has_sequence[new_chrom]) {
            printLog(LOG_LOW, "add sequence %d of %d (%s)\n",
                     trees->get_num_leaves() + 1, nseqs,
                     sequences->names[new_chrom].c_str());
            sample_arg_thread(model, sequences, trees, new_chrom);
            assert_trees(trees, model->pop_tree);
            printLog(LOG_LOW, "\n");
	    for (int buildup=1; buildup < num_buildup; buildup++) {
		printLog(LOG_LOW, "buildup rep %d of %d\n", buildup, num_buildup);
		resample_arg_random_leaf(model, sequences, trees);
		printLog(LOG_LOW, "\n");
	    }
        }
    }
}


// resample the threading of all the chromosomes
void resample_arg(const ArgModel *model, Sequences *sequences,
                  LocalTrees *trees)
{
    const int nleaves = trees->get_num_leaves();

    // cycle through chromosomes

    for (int chrom=0; chrom<nleaves; chrom++)
	resample_arg_leaf(model, sequences, trees, chrom);
}


// resample the threading of an internal branch
void resample_arg_all(const ArgModel *model, Sequences *sequences,
                      LocalTrees *trees, double prob_path_switch=.1)
{
    const int maxtime = model->get_removed_root_time();
    int *removal_path = new int [trees->get_num_trees()];

    // ramdomly choose a removal path
    int node = irand(trees->nnodes);
    int pos = irand(trees->start_coord, trees->end_coord);
    sample_arg_removal_path(trees, node, pos, removal_path, prob_path_switch);

    remove_arg_thread_path(trees, removal_path, maxtime, model->pop_tree);
    sample_arg_thread_internal(model, sequences, trees);

    delete [] removal_path;
}


// resample the threading of a leaf of an ARG
void resample_arg_leaf(const ArgModel *model, Sequences *sequences,
                       LocalTrees *trees, int node)
{
    const int maxtime = model->get_removed_root_time();
    int *removal_path = new int [trees->get_num_trees()];
    assert_trees(trees, model->pop_tree);

    sample_arg_removal_leaf_path(trees, node, removal_path);

    remove_arg_thread_path(trees, removal_path, maxtime, model->pop_tree);

    int mintime = sequences->ages[trees->seqids[node]];
    if (mintime > 0) {
        for (LocalTrees::iterator it=trees->begin(); it != trees->end(); ++it) {
            LocalTree *tree = it->tree;
            tree->nodes[node].age = mintime;
            assert(it->spr.recomb_node != node);
            assert(it->spr.coal_node != node);
        }
    }

    PhaseProbs *phase_pr = NULL;
    if (model->unphased)
        phase_pr = new PhaseProbs(trees->seqids[node], node,
                                  sequences, trees, model);
    sample_arg_thread_internal(model, sequences, trees, mintime, phase_pr);
    if (model->unphased)
        delete phase_pr;

    delete [] removal_path;
}


void resample_arg_random_leaf(const ArgModel *model, Sequences *sequences,
			      LocalTrees *trees)
{
    int node = irand(trees->get_num_leaves());
    resample_arg_leaf(model, sequences, trees, node);
}



// resample the threading of an internal branch using MCMC
bool resample_arg_mcmc(const ArgModel *model, Sequences *sequences,
                       LocalTrees *trees)
{
    const int maxtime = model->get_removed_root_time();
    int *removal_path = new int [trees->get_num_trees()];

    // save a copy of the local trees
    LocalTrees trees2;
    trees2.copy(*trees);

    // ramdomly choose a removal path
    double npaths = sample_arg_removal_path_uniform(trees, removal_path);
    remove_arg_thread_path(trees, removal_path, maxtime, model->pop_tree);
    sample_arg_thread_internal(model, sequences, trees);
    double npaths2 = count_total_arg_removal_paths(trees);

    // perform reject if needed
    double accept_prob = exp(npaths - npaths2);
    bool accept = (frand() < accept_prob);
    if (!accept)
        trees->copy(trees2);

    // logging
    printLog(LOG_LOW, "accept_prob = exp(%lf - %lf) = %f, accept = %d\n",
             npaths, npaths2, accept_prob, (int) accept);


    // clean up
    delete [] removal_path;

    return accept;
}

// resample the threading of an internal branch using MCMC
// Also sometimes resample leaves specifically
void resample_arg_mcmc_all(const ArgModel *model, Sequences *sequences,
                           LocalTrees *trees, bool do_leaf,
                           int window, int niters, double heat,
                           bool no_resample_mig)
{
    if (do_leaf) {
        resample_arg_random_leaf(model, sequences, trees);
        printLog(LOG_LOW, "resample_arg_leaf: accept=%f\n", 1.0);
    } else {
        // if there are migration events in tree, then choose a time interval
        // for focused resampling from the time intervals where migrations occur.
        // choose uniformly between each of these intervals and normal
        // subtree pruning.
        int time_interval = -1;
        int hap=-1;
        if (model->pop_tree != NULL && !no_resample_mig && frand() < 0.5) {
            vector<int> mig_times;
            vector<int> possible_haps;
            for (int i=0; i < model->ntimes-1; i++) {
                if (model->pop_tree->has_migration(i)) {
                    for (int hap=0; hap < trees->get_num_leaves(); hap++) {
                        int start_pop = sequences->get_pop(hap);
                        if (model->pop_tree->num_sub_path[0][i][start_pop] <
                            model->pop_tree->num_sub_path[0][i+1][start_pop] &&
                            trees->front().tree->nodes[hap].age <= i) {
                            mig_times.push_back(i);
                            possible_haps.push_back(hap);
                        }
                    }
                }
            }
            if (mig_times.size() > 0) {
                int val = irand(mig_times.size());
                time_interval = mig_times[val];
                hap = possible_haps[val];
            }
        }
        if (time_interval >= 0) {
            int num_break = resample_arg_by_time_and_hap(model, sequences,
                 trees, time_interval, hap);
            printLog(LOG_LOW, "resample_arg_by_hap (%i %s numbreak=%i): accept=1.0\n",
                     time_interval, sequences->names[hap].c_str(), num_break);
        } else {
            double accept_rate = resample_arg_regions(
              model, sequences, trees, window, niters, heat);
            printLog(LOG_LOW, "resample_arg_regions: accept=%f\n", accept_rate);
        }
    }
}



// resample the threading of an internal branch with preference for recombs
void resample_arg_recomb(const ArgModel *model, Sequences *sequences,
                         LocalTrees *trees, double recomb_preference)
{
    const int maxtime = model->get_removed_root_time();
    int *removal_path = new int [trees->get_num_trees()];

    // ramdomly choose a removal path weighted by recombinations
    sample_arg_removal_path_recomb(trees, recomb_preference, removal_path);
    remove_arg_thread_path(trees, removal_path, maxtime, model->pop_tree);
    sample_arg_thread_internal(model, sequences, trees);

    delete [] removal_path;
}



// resample an ARG heuristically and aggressively to high joint probability
void resample_arg_climb(const ArgModel *model, Sequences *sequences,
                        LocalTrees *trees, double recomb_preference)
{
    resample_arg_recomb(model, sequences, trees, recomb_preference);
}




/*
// sample an ARG with both sequential and gibbs iterations
void sample_arg_seq_gibbs(const ArgModel *model, const Sequences *sequences,
                          LocalTrees *trees, int seqiters, int gibbsiters)
{
    const int nseqs = sequences->get_num_seqs();
    const int seqlen = sequences->length();
    const int minseqs = 3;

    // initialize ARG as trunk
    const int capacity = 2 * nseqs - 1;
    trees->make_trunk(0, seqlen, capacity);

    int nleaves = 1;
    Timer time;
    while (nleaves < nseqs) {
        time.start();

        // sequential stage
        int nleaves2 = max(min(nleaves + seqiters, nseqs), minseqs);

        // add more chromosomes one by one
        for (int nchroms=nleaves+1; nchroms<=nleaves2; nchroms++) {
            // use first nchroms sequences
            Sequences sequences2(sequences, nchroms);
            int new_chrom = nchroms - 1;
            sample_arg_thread(model, &sequences2, trees, new_chrom);
        }
        nleaves = nleaves2;


        // gibbs stage
        // randomly choose gibbsiters chromosomes
        int chroms[nleaves];
        for (int i=0; i<nleaves; i++)
            chroms[i] = i;
        shuffle(chroms, nleaves);

        for (int i=0; i<min(gibbsiters, nleaves); i++)
            resample_arg_all(model, sequences, trees);

        printTimerLog(time, LOG_QUIET,
                      "seq_gibbs stage (%3d leaves):       ", nleaves);
    }
}
*/


//=============================================================================
// sub-region resampling


 // not updated for pop paths, not currently called
State find_state_sub_tree(
    const LocalTree *full_tree, const vector<int> &full_seqids,
    const LocalTree *partial_tree, const vector<int> &partial_seqids,
    int new_chrom)
{
    // reconcile full tree to partial tree
    int recon[full_tree->nnodes];
    map_congruent_trees(full_tree, &full_seqids[0],
                        partial_tree, &partial_seqids[0], recon);

    // find new chrom in full_tree
    int ptr = find_array(&full_seqids[0], full_seqids.size(), new_chrom);
    assert(ptr != -1);

    // walk up from new chrom until we hit a reconciled portion of full tree
    while (recon[ptr] == -1)
        ptr = full_tree->nodes[ptr].parent;

    return State(recon[ptr], full_tree->nodes[ptr].age);
}



// sequentially sample an ARG from scratch
// sequences are sampled in the order given
/*void cond_sample_arg_seq(const ArgModel *model, Sequences *sequences,
                         LocalTrees *trees,
                         LocalTree *start_tree, LocalTree *end_tree,
                         const vector<int> &full_seqids)
{
    // initialize ARG as trunk
    const int capacity = 2 * sequences->get_num_seqs() - 1;
    trees->make_trunk(trees->start_coord, trees->end_coord, capacity);

    // add more chromosomes one by one
    for (int nchroms=2; nchroms<=sequences->get_num_seqs(); nchroms++) {
        // use first nchroms sequences
        Sequences sequences2(sequences, nchroms);
        int new_chrom = nchroms - 1;

        // determine start and end states from given trees
        LocalTree *first_tree = trees->front().tree;
        LocalTree *last_tree = trees->back().tree;
        State start_state = find_state_sub_tree(
            start_tree, full_seqids, first_tree, trees->seqids, new_chrom);
        State end_state = find_state_sub_tree(
            end_tree, full_seqids, last_tree, trees->seqids, new_chrom);

        cond_sample_arg_thread(model, &sequences2, trees, new_chrom,
                               start_state, end_state);

        assert_trees(trees);
    }
    }*/



// sequentially sample an ARG only for a given region
// sequences are sampled in the order given
/*void sample_arg_seq_region(const ArgModel *model, Sequences *sequences,
                           LocalTrees *trees, int region_start, int region_end)
{
    // ensure region is within ARG
    assert(region_start > trees->start_coord);
    assert(region_end < trees->end_coord);
    assert(region_start < region_end);

    // partion trees into three segments
    LocalTrees *trees2 = partition_local_trees(trees, region_start);
    LocalTrees *trees3 = partition_local_trees(trees2, region_end);
    assert(trees2->length() == region_end - region_start);

    // resample region conditioning on starting and ending trees
    cond_sample_arg_seq(model, sequences,
                        trees2, trees->back().tree, trees3->front().tree,
                        trees->seqids);

    // rejoin trees
    append_local_trees(trees, trees2);
    append_local_trees(trees, trees3);

    // clean up
    delete trees2;
    delete trees3;
    }*/


State find_state_sub_tree_internal(const ArgModel *model,
        const LocalTree *full_tree, const LocalTree *partial_tree, int maxtime)
{
    if (partial_tree->nodes[partial_tree->root].age < maxtime) {
        // fully specified tree
        return State(-1, -1, -1);
    }

    // NOTE: do not assume internal nodes have same naming scheme between
    // trees

    int subtree_root = partial_tree->nodes[partial_tree->root].child[0];

    // identify node by path length from left most leaf
    int count = 0;
    int leaf = subtree_root;
    while (!partial_tree->nodes[leaf].is_leaf()) {
        leaf = partial_tree->nodes[leaf].child[0];
        count++;
    }

    // find equivalent node in full tree
    int ptr = leaf;
    while (count > 0) {
        ptr = full_tree->nodes[ptr].parent;
        count--;
    }

    // find sibling and age of coalescence
    int sib = full_tree->get_sibling(ptr);
    assert(sib != -1);
    int parent = full_tree->nodes[ptr].parent;
    assert(parent != -1);
    int coal_time = full_tree->nodes[parent].age;
    int pop_path = 0;
    if (model->pop_tree != NULL) {
        pop_path = model->pop_tree->path_to_root(full_tree->nodes,
                                                 ptr);
    }

    // identify sibling by leaf and path length
    count = 0;
    leaf = sib;
    while (!full_tree->nodes[leaf].is_leaf()) {
        leaf = full_tree->nodes[leaf].child[0];
        count++;
    }

    // map sib back to partial tree
    ptr = leaf;
    while (count > 0) {
        ptr = partial_tree->nodes[ptr].parent;
        count--;
    }

    return State(ptr, coal_time, pop_path);
}



int resample_arg_by_time_and_hap(
    const ArgModel *model, Sequences *sequences,
    LocalTrees *trees, int time_interval, int hap)
{
    const int maxtime = model->get_removed_root_time();
    static int count=0;
    const bool open_ended=true;
    LocalTrees orig_trees;
    decLogLevel();
    orig_trees.copy(*trees);

    int orig_numtree = trees->get_num_trees();
    int *removal_path = new int[orig_numtree];
    vector<int> break_coords;
    assert(time_interval >= 0 && time_interval < model->ntimes - 1);
    assert(hap >= 0 && hap < trees->get_num_leaves());

    get_arg_removal_path_by_ind_and_time(trees, time_interval, hap,
                                         removal_path,
                                         break_coords);
    int num_break = (int)break_coords.size();
    for (int i=0; i <= num_break; i++) {
        count++;
        int region_start, region_end;
        if (i == 0) {
            region_start = trees->start_coord;
        } else {
            region_start = break_coords[i-1]-1;
        }
        if (i == num_break) {
            region_end = trees->end_coord;
        } else region_end = break_coords[i]+1;

        LocalTrees otrees;
        otrees.copy(*trees);

        // partion trees into three segments
        LocalTrees *trees2 = partition_local_trees(trees, region_start, true);
        LocalTrees otrees2;
        otrees2.copy(*trees2);
        LocalTrees *trees3 = partition_local_trees(trees2, region_end, true);
        Spr stub_spr;
        int *stub_mapping=NULL;
        stub_spr.set_null();
        if (i != num_break) {
            if (trees2->trees.back().blocklen == 0) {
                stub_spr = trees2->back().spr;
                stub_mapping = new int[trees2->nnodes];
                for (int j=0; j < trees2->nnodes; j++)
                    stub_mapping[j] = trees2->back().mapping[j];
                trees2->trees.pop_back();
            }
            assert(trees2->trees.back().blocklen == 1);
        }
        vector<int> temp1, temp2;
        int curr_numtree = trees2->get_num_trees();
        if (curr_numtree > 2) {
            int *curr_removal_path = new int[curr_numtree];
            get_arg_removal_path_by_ind_and_time(trees2, time_interval, hap,
                                                 curr_removal_path, temp1,
       // do not use the first site to decide removal path; it may have changed from previous threading
                                                 i==0 || trees2->front().blocklen > 1);

            if (i != num_break) {
                const LocalTrees::iterator it = orig_trees.get_block(region_end-1);
                const LocalTrees::iterator it2 = orig_trees.get_block(region_end-2);
                int next_nodes[2];
                get_next_removal_nodes(it2->tree, it->tree, it->spr, it->mapping,
                                       curr_removal_path[curr_numtree-2],
                                       next_nodes);
                // we just want to choose a legal node to remove here
                // it doesn't matter what it is because we are going to condition
                // on end state
                curr_removal_path[curr_numtree-1] = next_nodes[0];
                assert(curr_removal_path[curr_numtree-1] != -1);
            }



            if (i != 0) {
                // same as above (except for beginning state)
                // the removal node at first tree does not matter
                int prev_nodes[2];
                list<LocalTreeSpr>::iterator it = trees2->begin();
                ++it;
                get_prev_removal_nodes(trees2->front().tree,
                                       it->tree, it->spr,
                                       it->mapping, curr_removal_path[1],
                                       prev_nodes);
                curr_removal_path[0] = prev_nodes[0];
                assert(curr_removal_path[0] != -1);
            }

            printLog(LOG_LOW, "region sample: iter=%d, region=(%d, %d)\n",
                     i, region_start, region_end);

            // get starting and ending trees
            LocalTree start_tree(*trees2->front().tree);
            LocalTree end_tree(*trees2->back().tree);

            // remove internal branch from trees2
            remove_arg_thread_path(trees2, curr_removal_path,
                                   maxtime, model->pop_tree);
            assert_trees(trees2, model->pop_tree, true);

            // determine start and end states from start and end trees
            LocalTree *start_tree_partial = trees2->front().tree;
            LocalTree *end_tree_partial = trees2->back().tree;
            State start_state = find_state_sub_tree_internal(model,
                  &start_tree, start_tree_partial, maxtime);
            State end_state = find_state_sub_tree_internal(model,
                  &end_tree, end_tree_partial, maxtime);

            // set start/end state to null if open ended is requested
            if (open_ended) {
                if (region_start == trees->start_coord)
                    start_state.set_null();
                if (region_end == trees3->end_coord)
                    end_state.set_null();
            }
            // sample new ARG conditional on start and end states
            decLogLevel();
            cond_sample_arg_thread_internal(model, sequences, trees2,
                                            start_state, end_state);
            incLogLevel();
            assert_trees(trees2, model->pop_tree);

            delete [] curr_removal_path;
        }

        // rejoin trees
        assert_trees(trees, model->pop_tree);
        assert_trees(trees2, model->pop_tree);
        assert_trees(trees3, model->pop_tree);
        append_local_trees(trees, trees2, true, model->pop_tree);
        if (trees3->get_num_trees() > 0) {
            trees3->front().spr = stub_spr;
            assert(trees3->front().mapping == NULL);
            trees3->front().mapping = stub_mapping;
        }

        append_local_trees(trees, trees3, true, model->pop_tree);
        assert_trees(trees, model->pop_tree);

        // clean up
        delete trees2;
        delete trees3;
    }
    delete [] removal_path;
    incLogLevel();
    return num_break;
}


// resample an ARG only for a given region
// all branches are possible to resample
// open_ended -- If true and region touches start or end of local trees do not
//               conditioned on state.
double resample_arg_region(
    const ArgModel *model, Sequences *sequences,
    LocalTrees *trees, int region_start, int region_end, int niters,
    bool open_ended, double heat)
{
    const int maxtime = model->get_removed_root_time();
    static int count=0;

    // special case: zero length region
    if (region_start == region_end)
        return 1.0;

    // assert region is within trees
    assert(region_start >= trees->start_coord);
    assert(region_end <= trees->end_coord);
    assert(region_start < region_end);

    // partion trees into three segments
    LocalTrees *trees2 = partition_local_trees(trees, region_start);
    LocalTrees *trees3 = partition_local_trees(trees2, region_end);
    assert(trees2->length() == region_end - region_start);

    // TODO: refactor
    // extend stub (zero length block) if it happens to exist
    bool stub = (trees2->trees.back().blocklen == 0);
    if (stub) {
        trees2->trees.back().blocklen += 1;
        trees2->end_coord++;
    }


    // perform several iterations of resampling
    int accepts = 0;
    for (int i=0; i<niters; i++) {
        count++;
        printLog(LOG_LOW, "region sample: iter=%d, region=(%d, %d)\n",
                 i, region_start, region_end);

        // save a copy of the local trees
        LocalTrees old_trees2;
        old_trees2.copy(*trees2);

        // get starting and ending trees
        LocalTree start_tree(*trees2->front().tree);
        LocalTree end_tree(*trees2->back().tree);

        // remove internal branch from trees2
        int *removal_path = new int [trees2->get_num_trees()];
        double npaths;
        npaths = sample_arg_removal_path_uniform(trees2, removal_path);
        remove_arg_thread_path(trees2, removal_path, maxtime, model->pop_tree);
        delete [] removal_path;
        assert_trees(trees2, model->pop_tree, true);

        // determine start and end states from start and end trees
        LocalTree *start_tree_partial = trees2->front().tree;
        LocalTree *end_tree_partial = trees2->back().tree;
        State start_state = find_state_sub_tree_internal(model,
            &start_tree, start_tree_partial, maxtime);
        State end_state = find_state_sub_tree_internal(model,
            &end_tree, end_tree_partial, maxtime);

        // set start/end state to null if open ended is requested
        if (open_ended) {
            if (region_start == trees->start_coord)
                start_state.set_null();
            if (region_end == trees3->end_coord)
                end_state.set_null();
        }
        // sample new ARG conditional on start and end states
        decLogLevel();
        cond_sample_arg_thread_internal(model, sequences, trees2,
                                        start_state, end_state);
        incLogLevel();
        assert_trees(trees2, model->pop_tree);

        double npaths2 = count_total_arg_removal_paths(trees2);

            // perform reject if needed
        double accept_prob = exp(heat*(npaths - npaths2));
        bool accept = (frand() < accept_prob);

        if (!accept) {
            trees2->copy(old_trees2);
        } else {
            accepts++;
        }

        // logging
        printLog(LOG_LOW, "accept_prob = exp(%lf - %lf) = %f, accept = %d\n",
                 npaths, npaths2, accept_prob, (int) accept);
        /*        printf("%i accept_prob = exp(%lf - %lf) = %f, accept = %d\n",
                  i, npaths, npaths2, accept_prob, (int) accept);*/
    }

    // remove stub if it exists
    if (stub) {
        trees2->trees.back().blocklen -= 1;
        trees2->end_coord--;
    }

    // rejoin trees
    append_local_trees(trees, trees2, true, model->pop_tree);
    append_local_trees(trees, trees3, true, model->pop_tree);

    // clean up
    delete trees2;
    delete trees3;

    return accepts / double(niters);
}




// resample an ARG a region at a time in a sliding window
double resample_arg_regions(
    const ArgModel *model, Sequences *sequences,
    LocalTrees *trees, int window, int niters, double heat)
{
    decLogLevel();
    double accept_rate = 0.0;
    int nwindows = 0;
    int currwindow = irand(window - window/4, window + window/4);
    int currstep = (int)currwindow/2+1;
    for (int start=trees->start_coord;
         start == trees->start_coord || start+currwindow/2 <trees->end_coord;
         start+=currstep)
    {
        nwindows++;
        int end = min(start + currwindow, trees->end_coord);
        accept_rate += resample_arg_region(
             model, sequences, trees, start, end, niters, true, heat);
    }
    incLogLevel();

    accept_rate /= nwindows;
    return accept_rate;
}

void resample_migrates(ArgModel *model,
                       const LocalTrees *trees,
                       vector<Spr> &invisible_recombs) {
    if (model->pop_tree == NULL) return;
    for (unsigned int i=0; i < model->pop_tree->mig_params.size(); i++) {
        MigParam mp = model->pop_tree->mig_params[i];
        int from_pop = mp.from_pop;
        int to_pop = mp.to_pop;
        int time_idx = mp.time_idx;
        int count, total;
        count_mig_events(from_pop, to_pop, time_idx,
                         model, trees, &invisible_recombs, &count, &total);
        double alpha = (double)count+mp.alpha;
        double beta = (double)(total-count)+mp.beta;
        double new_migrate=1.0;
        while (new_migrate > 0.5)
            new_migrate = rand_beta(alpha, beta);
        double curr_migrate = model->pop_tree->mig_matrix[time_idx].get(from_pop, to_pop);
        double diff = new_migrate - curr_migrate;
        double curr_self_rate = model->pop_tree->mig_matrix[time_idx].get(from_pop, from_pop);
        model->pop_tree->mig_matrix[mp.time_idx].set(mp.from_pop,
                                                     mp.to_pop,
                                                     new_migrate);
        model->pop_tree->mig_matrix[mp.time_idx].set(mp.from_pop,
                                                     mp.from_pop,
                                                     curr_self_rate - diff);
        model->pop_tree->update_population_probs();
    }
}

/*
// cut a branch in the ARG and resample branch
double resample_arg_cut(
    const ArgModel *model, const Sequences *sequences, LocalTrees *trees,
    int window_start, int window_end)
{
    const int maxtime = model->get_removed_root_time();

    // cut branch and remove thread
    int cuttime;
    int *removal_path = new int [trees->get_num_trees()];
    int region_start, region_end;
    sample_arg_removal_path_cut(trees, model->ntimes, removal_path, &cuttime,
                                &region_start, &region_end,
                                window_start, window_end);

    // partion trees into three segments
    // trees2 will be resampled
    // do not perform SPR trimming since we will manually adjust SPRs
    const bool trim = false;
    LocalTrees *trees2 = partition_local_trees(trees, region_start, trim);
    LocalTrees *trees3 = partition_local_trees(trees2, region_end, trim);
    assert(trees2->length() == region_end - region_start);
    assert(trees2->trees.back().blocklen > 0);

    // suppress first SPR and mapping of trees2
    // since the resampling code excepts no SPR on first tree
    LocalTrees::iterator it = trees2->begin();
    int *mapping = NULL;
    if (it->mapping) {
        mapping = it->mapping;
        it->mapping = NULL;
    }
    Spr spr = it->spr;
    it->spr.set_null();

    // ensure remove path is compatiable with trees2
    int ntrees = trees2->get_num_trees();
    int *removal_path2 = removal_path;
    while (removal_path2[0] == -1) removal_path2++;
    assert(removal_path2[ntrees-1] != -1);
    assert(trees3->length() == 0 || removal_path2[ntrees] == -1);

    // remove and resample region
    printLog(LOG_LOW, "branch cut: time=%d, region=[%d,%d]\n",
             cuttime, region_start, region_end);
    remove_arg_thread_path(trees2, removal_path2, maxtime, model->pop_tree);
    sample_arg_thread_internal(model, sequences, trees2, cuttime);
    //assert_trees(trees2);

    // restore previously suppressed spr and mapping
    if (trees->length() > 0) {
        it = trees2->begin();
        if (mapping)
            it->mapping = mapping;
        it->spr = spr;
    }

    // rejoin trees
    append_local_trees(trees, trees2);
    append_local_trees(trees, trees3);
    assert_trees(trees);

    // clean up
    delete [] removal_path;
    delete trees2;
    delete trees3;

    return 1;
}


double resample_arg_cut(
    const ArgModel *model, const Sequences *sequences, LocalTrees *trees,
    int window, int step, int niters)
{
    decLogLevel();
    for (int start=trees->start_coord;
         start == trees->start_coord || start+window/2 <trees->end_coord;
         start+=step)
    {
        int end = min(start + window, trees->end_coord);
        for (int i=0; i<niters; i++)
            resample_arg_cut(model, sequences, trees, start, end);
    }
    incLogLevel();

    return 1;
}
*/


//=============================================================================
// C interface
extern "C" {

// sequentially sample until all chromosomes are present
LocalTrees *arghmm_complete_arg(
    LocalTrees *trees, ArgModel *model, Sequences *sequences)
{
    const int nseqs = sequences->get_num_seqs();
    for (int new_chrom=trees->get_num_leaves(); new_chrom<nseqs; new_chrom++)
        sample_arg_thread(model, sequences, trees, new_chrom);
    return trees;
}


// sequentially sample an ARG
LocalTrees *arghmm_sample_arg_seq(
    double *times, int ntimes,
    double **popsizes, double rho, double mu,
    char **seqs, int nseqs, int seqlen)
{
    // setup model, local trees, sequences
    ArgModel model(ntimes, times, popsizes, rho, mu);
    Sequences sequences(seqs, nseqs, seqlen);
    LocalTrees *trees = new LocalTrees();

    sample_arg_seq(&model, &sequences, trees);

    return trees;
}


// resample an ARG with gibbs
LocalTrees *arghmm_resample_arg(
    LocalTrees *trees, double *times, int ntimes,
    double **popsizes, double rho, double mu,
    char **seqs, int nseqs, int seqlen, int niters, int nremove)
{
    // setup model, local trees, sequences
    ArgModel model(ntimes, times, popsizes, rho, mu);
    Sequences sequences(seqs, nseqs, seqlen);

    // sequentially sample until all chromosomes are present
    arghmm_complete_arg(trees, &model, &sequences);

    // gibbs sample
    for (int i=0; i<niters; i++)
        resample_arg(&model, &sequences, trees);

    return trees;
}


// resample all branches in an ARG with gibbs
LocalTrees *arghmm_resample_all_arg(
    LocalTrees *trees, double *times, int ntimes,
    double **popsizes, double rho, double mu,
    char **seqs, int nseqs, int seqlen, int niters, double prob_path_switch)
{
    // setup model, local trees, sequences
    ArgModel model(ntimes, times, popsizes, rho, mu);
    Sequences sequences(seqs, nseqs, seqlen);

    // sequentially sample until all chromosomes are present
    arghmm_complete_arg(trees, &model, &sequences);

    // gibbs sample
    for (int i=0; i<niters; i++)
        resample_arg_all(&model, &sequences, trees, prob_path_switch);

    return trees;
}


// resample all branches in an ARG with mcmc
LocalTrees *arghmm_resample_mcmc_arg(
    LocalTrees *trees, double *times, int ntimes,
    double **popsizes, double rho, double mu,
    char **seqs, int nseqs, int seqlen, int niters, int niters2, int window)
{
    // setup model, local trees, sequences
    double frac_leaf = 0.5;

    ArgModel model(ntimes, times, popsizes, rho, mu);
    Sequences sequences(seqs, nseqs, seqlen);

    // sequentially sample until all chromosomes are present
    arghmm_complete_arg(trees, &model, &sequences);

    // gibbs sample
    for (int i=0; i<niters; i++) {
        printLog(LOG_LOW, "sample %d\n", i);
        resample_arg_mcmc_all(&model, &sequences, trees, frand() < frac_leaf,
                              window, niters2);
    }

    return trees;
}


LocalTrees *arghmm_resample_arg_leaf(
    LocalTrees *trees, double *times, int ntimes,
    double **popsizes, double rho, double mu,
    char **seqs, int nseqs, int seqlen, int niters)
{
    // setup model, local trees, sequences
    ArgModel model(ntimes, times, popsizes, rho, mu);
    Sequences sequences(seqs, nseqs, seqlen);

    // sequentially sample until all chromosomes are present
    arghmm_complete_arg(trees, &model, &sequences);

    // gibbs sample
    for (int i=0; i<niters; i++)
        resample_arg_random_leaf(&model, &sequences, trees);

    return trees;
}


/*
LocalTrees *arghmm_resample_arg_cut(
    LocalTrees *trees, double *times, int ntimes,
    double **popsizes, double rho, double mu,
    char **seqs, int nseqs, int seqlen, int niters)
{
    // setup model, local trees, sequences
    ArgModel model(ntimes, times, popsizes, rho, mu);
    Sequences sequences(seqs, nseqs, seqlen);

    // gibbs sample
    for (int i=0; i<niters; i++)
        resample_arg_cut(&model, &sequences, trees);

    return trees;
}
*/


// resample ARG focused on recombinations
LocalTrees *arghmm_resample_climb_arg(
    LocalTrees *trees, double *times, int ntimes,
    double **popsizes, double rho, double mu,
    char **seqs, int nseqs, int seqlen, int niters, double recomb_preference)
{
    // setup model, local trees, sequences
    ArgModel model(ntimes, times, popsizes, rho, mu);
    Sequences sequences(seqs, nseqs, seqlen);

    // sequentially sample until all chromosomes are present
    arghmm_complete_arg(trees, &model, &sequences);

    // gibbs sample
    for (int i=0; i<niters; i++)
        resample_arg_climb(&model, &sequences, trees, recomb_preference);

    return trees;
}


LocalTrees *arghmm_resample_arg_region(
    LocalTrees *trees, double *times, int ntimes,
    double **popsizes, double rho, double mu,
    char **seqs, int nseqs, int seqlen,
    int region_start, int region_end, int niters)
{
    // setup model, local trees, sequences
    ArgModel model(ntimes, times, popsizes, rho, mu);
    Sequences sequences(seqs, nseqs, seqlen);

    resample_arg_region(&model, &sequences, trees,
                        region_start, region_end, niters);

    return trees;
}

} // extern C

} // namespace argweaver
