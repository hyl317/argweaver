
// C/C++ includes
#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "err.h"
#include <iostream>

// argweaver includes
#include "compress.h"
#include "common.h"
#include "local_tree.h"
#include "logging.h"
#include "parsing.h"
#include "pop_model.h"

// tskit includes
#include <tskit.h>
#include "rspr.h"
#include "lgt.h"

#define check_tsk_error(val)                                                            \
    if (val < 0) {                                                                      \
        errx(EXIT_FAILURE, "line %d: %s", __LINE__, tsk_strerror(val));                 \
    }


namespace argweaver {

struct edge
{
    tsk_id_t parent;
    tsk_id_t child;
    int start;
};

//=============================================================================
// tree methods


LocalNode null_node;


// Counts the number of lineages in a tree for each time segment
//
// NOTE: Nodes in the tree are not allowed to exist at the top time point
// point (ntimes - 1).
//
// tree      -- local tree to count
// ntimes    -- number of time segments
// nbranches -- number of branches that exists between time i and i+1
// nrecombs  -- number of possible recombination points at time i
// ncoals    -- number of possible coalescing points at time i
void count_lineages(const LocalTree *tree, int ntimes,
                    int *nbranches, int *nrecombs,
                    int **nbranches_pop, int **ncoals_pop,
                    const PopulationTree *pop_tree)
{
    const LocalNode *nodes = tree->nodes;
    int npop = ( pop_tree == NULL ? 1 : pop_tree->npop );

    // initialize counts
    for (int i=0; i < ntimes; i++) {
        nbranches[i] = 0;
        nrecombs[i] = 0;
    }
    for (int i=0; i < npop; i++) {
        for (int j=0; j < 2*ntimes; j++)
            nbranches_pop[i][j] = 0;
	for (int j=0; j<ntimes; j++)
	    ncoals_pop[i][j] = 0;
    }

    // iterate over the branches of the tree
    for (int i=0; i<tree->nnodes; i++) {
        assert(nodes[i].age < ntimes - 1);
        const int parent = nodes[i].parent;
        const int parent_age = ((parent == -1) ? ntimes - 2 :
                                nodes[parent].age);

        // add counts for every segment along branch
        for (int j=nodes[i].age; j<parent_age; j++) {
            int pop = nodes[i].get_pop(j, pop_tree);
            nbranches[j]++;
            nrecombs[j]++;
            nbranches_pop[pop][2*j]++;
            ncoals_pop[pop][j]++;
            pop = nodes[i].get_pop(j+1, pop_tree);
            nbranches_pop[pop][2*j+1]++;
        }


        int pop = nodes[i].get_pop(parent_age, pop_tree);
        // recomb and coal are also allowed at the top of a branch
        nrecombs[parent_age]++;
        ncoals_pop[pop][parent_age]++;
        if (parent == -1) {
            nbranches[parent_age]++;
            nbranches_pop[pop][2*parent_age]++;
            pop = nodes[i].get_pop(parent_age+1, pop_tree);
            nbranches_pop[pop][2*parent_age+1]++;
        }
    }

    // ensure last time segment always has one branch
    nbranches[ntimes - 1] = 1;
    int final_pop = (pop_tree == NULL ? 0 : pop_tree->final_pop() );
    for (int i=0; i < npop; i++) {
	nbranches_pop[i][2*(ntimes - 1)] =
            nbranches_pop[i][2*ntimes - 1] =
            (i == final_pop ? 1 : 0);
        ncoals_pop[i][ntimes-1] = (i == final_pop ? 1 : 0);
    }
}


// Counts the number of lineages in a tree for each time segment
//
// NOTE: Nodes in the tree are not allowed to exist at the top time point
// point (ntimes - 1).
//
// tree      -- local tree to count
// ntimes    -- number of time segments
// nbranches -- number of branches that exists between time i and i+1
// nrecombs  -- number of possible recombination points at time i
// ncoals    -- number of possible coalescing points at time i
void count_lineages_internal(const LocalTree *tree, int ntimes,
                             int *nbranches, int *nrecombs,
                             int **nbranches_pop, int **ncoals_pop,
                             const PopulationTree *pop_tree)
{
    const LocalNode *nodes = tree->nodes;
    const int subtree_root = nodes[tree->root].child[0];
    //    const int minage = nodes[subtree_root].age;
    int npop = ( pop_tree == NULL ? 1 : pop_tree->npop );

    // initialize counts
    for (int i=0; i<ntimes; i++) {
        nbranches[i]=0;
        nrecombs[i]=0;
    }
    for (int i=0; i < npop; i++) {
        for (int j=0; j < 2*ntimes; j++)
            nbranches_pop[i][j] = 0;
	for (int j=0; j<ntimes; j++)
	    ncoals_pop[i][j] = 0;
    }

    // iterate over the branches of the tree
    for (int i=0; i<tree->nnodes; i++) {
        // skip virtual branches
        if (i == subtree_root || i == tree->root)
            continue;

        assert(nodes[i].age < ntimes - 1);
        const int parent = nodes[i].parent;
        const int parent_age = ((parent == tree->root) ? ntimes - 2 :
                                nodes[parent].age);

        // add counts for every segment along branch
        for (int j=nodes[i].age; j<parent_age; j++) {
            int pop = nodes[i].get_pop(j, pop_tree);
            nbranches[j]++;
            nrecombs[j]++;
            nbranches_pop[pop][2*j]++;
            ncoals_pop[pop][j]++;
            pop = nodes[i].get_pop(j+1, pop_tree);
            nbranches_pop[pop][2*j+1]++;
            assert(j < ntimes-1);
        }

        // recomb and coal are also allowed at the top of a branch
        int pop = nodes[i].get_pop(parent_age, pop_tree);
        nrecombs[parent_age]++;
        ncoals_pop[pop][parent_age]++;
        if (parent == tree->root) {
            nbranches[parent_age]++;
            nbranches_pop[pop][2*parent_age]++;
            pop = nodes[i].get_pop(parent_age+1, pop_tree);
            nbranches_pop[pop][2*parent_age+1]++;
        }
    }

    // ensure last time segment always has one branch
    nbranches[ntimes-1]=1;
    int final_pop = ( pop_tree == NULL ? 0 : pop_tree->final_pop() );
    for (int i=0; i < npop; i++) {
        if (i == final_pop) {
            nbranches_pop[i][2*ntimes - 2] = 1;
            nbranches_pop[i][2*ntimes - 1] = 1;
            ncoals_pop[i][ntimes-1] = 1;
        } else assert(nbranches_pop[i][2*ntimes - 1] == 0 &&
                      nbranches_pop[i][2*ntimes - 2] == 0 &&
                      ncoals_pop[i][ntimes-1] == 0);
    }
}



// Calculate tree length according to ArgHmm rules
double get_treelen(const LocalTree *tree, const double *times, int ntimes,
                   bool use_basal)
{
    double treelen = 0.0;
    const LocalNode *nodes = tree->nodes;

    for (int i=0; i<tree->nnodes; i++) {
        int parent = nodes[i].parent;
        int age = nodes[i].age;
        if (parent == -1) {
            // add basal stub
            if (use_basal)
                treelen += times[age+1] - times[age];
        } else {
            treelen += times[nodes[parent].age] - times[age];
        }
    }

    return treelen;
}


double get_treelen_internal(const LocalTree *tree, const double *times,
                            int ntimes)
{
    double treelen = 0.0;
    const LocalNode *nodes = tree->nodes;

    for (int i=0; i<tree->nnodes; i++) {
        int parent = nodes[i].parent;
        int age = nodes[i].age;
        if (parent == tree->root || parent == -1) {
            // skip virtual branches
        } else {
            treelen += times[nodes[parent].age] - times[age];
            assert(!isnan(treelen));
        }
    }

    return treelen;
}


double get_treelen_branch(const LocalTree *tree, const double *times,
                          int ntimes, int node, int time,
                          double treelen, bool use_basal)
{
    double root_time;
    int rooti = tree->nodes[tree->root].age;

    if (treelen < 0.0)
        treelen = get_treelen(tree, times, ntimes, false);

    double blen = times[time];
    double treelen2 = treelen + blen;
    if (node == tree->root) {
        treelen2 += blen - times[tree->nodes[tree->root].age];
        root_time = times[time+1] - times[time];
    } else {
        rooti = tree->nodes[tree->root].age;
        root_time = times[rooti+1] - times[rooti];
    }

    if (use_basal)
        return treelen2 + root_time;
    else
        return treelen2;
}


double get_basal_branch(const LocalTree *tree, const double *times, int ntimes,
                        int node, int time)
{
    double root_time;

    if (node == tree->root) {
        root_time = times[time+1] - times[time];
    } else {
        int rooti = tree->nodes[tree->root].age;
        root_time = times[rooti+1] - times[rooti];
    }

    return root_time;
}


// time_idx2 is based on half-time intervals and should be odd,
// since migrations occur between time intervals
void count_mig_events(int from_pop, int to_pop,
                      int time_idx2, const ArgModel *model,
                      const LocalTrees *trees,
                      const vector<Spr> *invisible_recombs,
                      int *count, int *total) {
    assert(time_idx2 % 2 == 1);
    int lower_time = time_idx2 / 2;
    int upper_time = lower_time+1;
    *count = *total = 0;
    LocalTree *tree = trees->front().tree;
    for (int i=0; i < tree->nnodes; i++) {
        if (tree->nodes[i].age <= lower_time &&
            ( i == tree->root ||
              tree->nodes[tree->nodes[i].parent].age >= upper_time)) {
            if (model->get_pop(tree->nodes[i].pop_path, lower_time) == from_pop) {
                (*total)++;
                if (model->get_pop(tree->nodes[i].pop_path, upper_time) == to_pop) {
                    (*count)++;
                }
            }
        }
    }
    for (LocalTrees::const_iterator it=trees->begin(); it != trees->end(); ++it) {
        const Spr *spr = &it->spr;
        if (spr->is_null()) continue;
        if (spr->recomb_time > lower_time || spr->coal_time < upper_time)
            continue;
        if (model->get_pop(spr->pop_path, lower_time) != from_pop)
            continue;
        (*total)++;
        if (model->get_pop(spr->pop_path, upper_time) == to_pop)
            (*count)++;
    }
    if (invisible_recombs == NULL) return;
    for (unsigned int i=0; i < invisible_recombs->size(); i++) {
        const Spr spr = (*invisible_recombs)[i];
        if (spr.is_null()) continue;
        if (spr.recomb_time > lower_time ||
            spr.coal_time < upper_time) continue;
        if (model->get_pop(spr.pop_path, lower_time) != from_pop)
            continue;
        (*total)++;
        if (model->get_pop(spr.pop_path, upper_time) == to_pop)
            (*count)++;
    }
}


// modify a local tree by Subtree Pruning and Regrafting
void apply_spr(LocalTree *tree, const Spr &spr,
               const PopulationTree *pop_tree)
{
    // before SPR:
    //       bp          cp
    //      / \           \       .
    //     rc              c
    //    / \                     .
    //   r   rs

    // after SPR:
    //    bp         cp
    //   /  \         \           .
    //  rs             rc
    //                /  \        .
    //               r    c

    // key:
    // r = recomb branch
    // rs = sibling of recomb branch
    // rc = recoal node (broken node)
    // bp = parent of broken node
    // c = coal branch
    // cp = parent of coal branch

    // updating population paths:
    // bp, cp, c paths are unchanged
    // r->pop_path becomes path consistent with r->path and spr->path
    // rs->pop_path becomes path consistent with rc->path and original rs->path
    // rc->pop_path becomes c->pop_path

    LocalNode *nodes = tree->nodes;

    // trival case
    if (spr.recomb_node == tree->root) {
        assert(0); // Melissa added: how can this happen?
        assert(spr.coal_node == tree->root);
        nodes[tree->root].pop_path = spr.pop_path; // not sure this is correct. check if assert above ever fails.
        return;
    }

    if (spr.recomb_node == spr.coal_node) {
        assert(pop_tree != NULL);
        int path1 = tree->nodes[spr.recomb_node].pop_path;
        int path2 = spr.pop_path;
        assert(! pop_tree->paths_equal(path1, path2,
                                       spr.recomb_time, spr.coal_time));
        int path3 =
            pop_tree->consistent_path(path1, path2,
                                      tree->nodes[spr.recomb_node].age,
                                      spr.recomb_time, spr.coal_time);
        tree->nodes[spr.recomb_node].pop_path =
            pop_tree->consistent_path(path3, path1,
                                      tree->nodes[spr.recomb_node].age,
                                      spr.coal_time, -1);
        return;
    }

    // recoal is also the node we are breaking
    int recoal = nodes[spr.recomb_node].parent;

    // find recomb node sibling and broke node parent
    int *c = nodes[recoal].child;
    int other = (c[0] == spr.recomb_node ? 1 : 0);
    int recomb_sib = c[other];
    int broke_parent =  nodes[recoal].parent;
    if (pop_tree != NULL)
        nodes[recomb_sib].pop_path = pop_tree->path_to_root(nodes, recomb_sib);

    // fix recomb sib pointer
    nodes[recomb_sib].parent = broke_parent;

    // fix parent of broken node
    int x = 0;
    if (broke_parent != -1) {
        c = nodes[broke_parent].child;
        x = (c[0] == recoal ? 0 : 1);
        nodes[broke_parent].child[x] = recomb_sib;
    }

    // reuse node as recoal
    if (spr.coal_node == recoal) {
        // we just broke coal_node, so use recomb_sib
        nodes[recoal].child[other] = recomb_sib;
        nodes[recoal].parent = nodes[recomb_sib].parent;
        nodes[recomb_sib].parent = recoal;
        if (broke_parent != -1)
            nodes[broke_parent].child[x] = recoal;
        if (pop_tree != NULL)
            nodes[recoal].pop_path = nodes[recomb_sib].pop_path;
    } else {
        nodes[recoal].child[other] = spr.coal_node;
        nodes[recoal].parent = nodes[spr.coal_node].parent;
        nodes[recoal].pop_path = nodes[spr.coal_node].pop_path;
        nodes[spr.coal_node].parent = recoal;

        // fix coal_node parent
        int parent = nodes[recoal].parent;
        if (parent != -1) {
            c = nodes[parent].child;
            if (c[0] == spr.coal_node)
                c[0] = recoal;
            else
                c[1] = recoal;
        }
    }
    if (pop_tree != NULL) {
        int path1 = pop_tree->consistent_path(nodes[spr.recomb_node].pop_path,
                                              spr.pop_path,
                                              nodes[spr.recomb_node].age,
                                              spr.recomb_time,
                                              spr.coal_time);
        nodes[spr.recomb_node].pop_path =
            pop_tree->consistent_path(path1, nodes[spr.coal_node].pop_path,
                                      nodes[spr.recomb_node].age, spr.coal_time, -1);
    }
    nodes[recoal].age = spr.coal_time;

    // set new root
    int root;
    if (spr.coal_node == tree->root)
        root = recoal;
    else if (recoal == tree->root) {
        if (spr.coal_node == recomb_sib)
            root = recoal;
        else
            root = recomb_sib;
    } else {
        root = tree->root;
    }
    tree->root = root;
}


//=============================================================================
// local trees methods


LocalTrees::LocalTrees(int **ptrees, int**ages, int **isprs, int *blocklens,
                       int ntrees, int nnodes, int capacity, int start) :
    chrom("chr"),
    start_coord(start),
    nnodes(nnodes)
{
    if (capacity < nnodes)
        capacity = nnodes;

    // copy data
    int pos = start;
    for (int i=0; i<ntrees; i++) {
        end_coord = pos + blocklens[i];

        // make mapping
        int *mapping = NULL;
        if (i > 0) {
            mapping = new int [nnodes];
            make_node_mapping(ptrees[i-1], nnodes, isprs[i][0], mapping);
        }

        trees.push_back(LocalTreeSpr(new LocalTree(ptrees[i], nnodes, ages[i],
                                                   NULL, capacity),
                                     isprs[i], blocklens[i], mapping));

        pos = end_coord;
    }

    set_default_seqids();
}


// Copy tree structure from another tree
void LocalTrees::copy(const LocalTrees &other)
{
    // clear previous data
    clear();

    // copy over information
    chrom = other.chrom;
    start_coord = other.start_coord;
    end_coord = other.end_coord;
    nnodes = other.nnodes;
    seqids = other.seqids;

    // copy local trees
    for (const_iterator it=other.begin(); it != other.end(); ++it) {
        const int nnodes = it->tree->nnodes;
        LocalTree *tree2 = new LocalTree();
        tree2->copy(*it->tree);

        int *mapping = it->mapping;
        int *mapping2 = NULL;
        if (mapping) {
            mapping2 = new int [nnodes];
            for (int i=0; i<nnodes; i++)
                mapping2[i] = mapping[i];
        }

        trees.push_back(LocalTreeSpr(tree2, it->spr, it->blocklen, mapping2));
    }
}


// get total ARG length
double get_arglen(const LocalTrees *trees, const double *times)
{
    double arglen = 0.0;

    for (LocalTrees::const_iterator it=trees->begin(); it!=trees->end(); ++it) {
        const LocalNode *nodes = it->tree->nodes;
        const int nnodes = it->tree->nnodes;

        double treelen = 0.0;
        for (int i=0; i<nnodes; i++) {
            int parent = nodes[i].parent;
            if (parent != -1)
                treelen += times[nodes[parent].age] - times[nodes[i].age];
        }

        arglen += treelen * it->blocklen;
    }

    return arglen;
}


// removes a null SPR from one local tree
bool remove_null_spr(LocalTrees *trees, LocalTrees::iterator it,
                     const PopulationTree *pop_tree)
{
    // look one tree ahead
    LocalTrees::iterator it2 = it;
    ++it2;
    if (it2 == trees->end())
        return false;

    // get spr from next tree, skip it if it is not null
    Spr *spr2 = &it2->spr;
    if (!spr2->is_null())
        return false;

    int nnodes = it2->tree->nnodes;

    int subtree_root = it->tree->nodes[it->tree->root].child[0];
    for (int i=0; i < it2->tree->nnodes; i++) {
        assert(it->tree->nodes[i].age == it2->tree->nodes[it2->mapping[i]].age);
        if (i != it->tree->root) {
            assert(it->tree->nodes[it->tree->nodes[i].parent].age ==
                   it2->tree->nodes[it2->tree->nodes[it2->mapping[i]].parent].age);
        }
        assert(i==subtree_root ||
               ( pop_tree == NULL ||
                 pop_tree->paths_equal(it->tree->nodes[i].pop_path,
                                       it2->tree->nodes[it2->mapping[i]].pop_path,
                                       it->tree->nodes[i].age,
                                       i == it->tree->root ? -1 :
                                       it->tree->nodes[it->tree->nodes[i].parent].age)));
    }


    if (it->mapping == NULL) {
        // it2 will become first tree and therefore does not need a mapping
        delete [] it2->mapping;
        it2->mapping = NULL;
    } else {
        // compute transitive mapping
        int *M1 = it->mapping;
        int *M2 = it2->mapping;
        int mapping[nnodes];
        for (int i=0; i<nnodes; i++) {
            if (M1[i] != -1)
                mapping[i] = M2[M1[i]];
            else
                mapping[i] = -1;
        }

        // set mapping
        for (int i=0; i<nnodes; i++)
            M2[i] = mapping[i];

        // copy over non-null spr
        *spr2 = it->spr;
        assert(!spr2->is_null());
    }


    // delete this tree
    it2->blocklen += it->blocklen;
    it->clear();
    trees->trees.erase(it);

    return true;
}



// Removes trees with null SPRs from the local trees
void remove_null_sprs(LocalTrees *trees, const PopulationTree *pop_tree)
{
    for (LocalTrees::iterator it=trees->begin(); it != trees->end();) {
        LocalTrees::iterator it2 = it;
        ++it2;
        remove_null_spr(trees, it, pop_tree);
        it = it2;
    }
}


// find recoal node, it is the node with no inward mappings
int get_recoal_node(const LocalTree *tree, const Spr &spr, const int *mapping)
{
    const int nnodes = tree->nnodes;
    bool mapped[nnodes];
    fill(mapped, mapped + nnodes, false);

    for (int i=0; i<nnodes; i++)
        if (mapping[i] != -1)
            mapped[mapping[i]] = true;

    for (int i=0; i<nnodes; i++)
        if (!mapped[i])
            return i;

    // this can happen for self-recombinations, though spr.recomb_node may
    // not equal spr.coal_node if it has already been renamed in
    // remove_arg_thread_path
    return spr.coal_node;
    assert(false);
    return -1;
}


void get_inverse_mapping(const int *mapping, int size, int *inv_mapping)
{
    // make inverse mapping
    fill(inv_mapping, inv_mapping + size, -1);
    for (int i=0; i<size; i++)
        if (mapping[i] != -1)
            inv_mapping[mapping[i]] = i;
}


LocalTrees *partition_local_trees(LocalTrees *trees, int pos,
                                  LocalTrees::iterator it, int it_start,
                                  bool trim)
{
    // create new local trees
    LocalTrees *trees2 = new LocalTrees(pos, trees->end_coord, trees->nnodes);
    trees2->chrom = trees->chrom;
    trees2->seqids.insert(trees2->seqids.end(), trees->seqids.begin(),
                          trees->seqids.end());

    // splice trees over
    trees2->trees.splice(trees2->begin(), trees->trees, it, trees->end());

    LocalTrees::iterator it2 = trees2->begin();
    if (trim) {
        // copy first tree back
            LocalTree *tree = it2->tree;
            LocalTree *last_tree = new LocalTree(tree->nnodes, tree->capacity);
            last_tree->copy(*tree);

            int *mapping = NULL;
            if (it2->mapping) {
                mapping = new int[trees->nnodes];
                for (int i=0; i<trees->nnodes; i++)
                    mapping[i] = it2->mapping[i];
            }

            trees->trees.push_back(
               LocalTreeSpr(last_tree, it2->spr, pos - it_start, mapping));

        // modify first tree of trees2
        if (it2->mapping)
            delete [] it2->mapping;
        it2->mapping = NULL;
        it2->spr.set_null();
    }

    trees->end_coord = pos;
    it2->blocklen -= pos - it_start;
    assert(it2->blocklen > 0);

    //assert_trees(trees);
    //assert_trees(trees2);

    return trees2;
}


// breaks a list of local trees into two separate trees
// Returns second list of local trees.
LocalTrees *partition_local_trees(LocalTrees *trees, int pos, bool trim)
{
    // special case (pos at beginning of local trees)
    if (pos == trees->start_coord) {
        LocalTrees *trees2 = new LocalTrees(pos, trees->end_coord,
                                            trees->nnodes);
        trees2->chrom = trees->chrom;
        trees2->seqids.insert(trees2->seqids.end(), trees->seqids.begin(),
                              trees->seqids.end());
        trees2->trees.splice(trees2->begin(), trees->trees,
                             trees->begin(), trees->end());
        trees->end_coord = pos;
        return trees2;
    }

    // special case (pos at end of local trees)
    if (pos == trees->end_coord) {
        LocalTrees *trees2 = new LocalTrees(pos, pos, trees->nnodes);
        trees2->chrom = trees->chrom;
        trees2->seqids.insert(trees2->seqids.end(), trees->seqids.begin(),
                              trees->seqids.end());
        trees2->seqids.insert(trees2->seqids.end(), trees->seqids.begin(),
                              trees->seqids.end());
        return trees2;
    }


    // find break point
    int start, end;
    LocalTrees::iterator it = trees->get_block(pos, start, end);
    if (it != trees->end())
        return partition_local_trees(trees, pos, it, start, trim);

    // break point was not found
    return NULL;
}


// Returns a mapping from nodes in tree1 to equivalent nodes in tree2
// If no equivalent is found, node maps to -1
void map_congruent_trees(const LocalTree *tree1, const int *seqids1,
                         const LocalTree *tree2, const int *seqids2,
                         int *mapping)
{
    const int nleaves1 = tree1->get_num_leaves();
    const int nleaves2 = tree2->get_num_leaves();

    for (int i=0; i<tree1->nnodes; i++)
        mapping[i] = -1;

    // reconcile leaves
    for (int i=0; i<nleaves1; i++) {
        const int seqid = seqids1[i];
        mapping[i] = -1;
        for (int j=0; j<nleaves2; j++) {
            if (seqids2[j] == seqid) {
                mapping[i] = j;
                break;
            }
        }
    }

    // postorder iterate over full tree to reconcile internal nodes
    int order[tree1->nnodes];
    tree1->get_postorder(order);
    LocalNode *nodes = tree1->nodes;
    for (int i=0; i<tree1->nnodes; i++) {
        const int j = order[i];
        const int *child = nodes[j].child;

        if (!nodes[j].is_leaf()) {
            if (mapping[child[0]] != -1) {
                if (mapping[child[1]] != -1) {
                    // both children mapping, so we map to their LCA
                    mapping[j] = tree2->nodes[mapping[child[0]]].parent;
                    assert(tree2->nodes[mapping[child[0]]].parent ==
                           tree2->nodes[mapping[child[1]]].parent);
                } else {
                    // single child maps, copy its mapping
                    mapping[j] = mapping[child[0]];
                }
            } else {
                if (mapping[child[1]] != -1) {
                    // single child maps, copy its mapping
                    mapping[j] = mapping[child[1]];
                } else {
                    // neither child maps, so neither do we
                    mapping[j] = -1;
                }
            }
        }
    }
}



// infer the mapping between two trees that differ by an SPR with known
// recombination node
void infer_mapping(const LocalTree *tree1, const LocalTree *tree2,
                   int recomb_node, int *mapping)
{
    const int nleaves1 = tree1->get_num_leaves();

    // map leaves
    for (int i=0; i<nleaves1; i++)
        mapping[i] = i;

    for (int i=nleaves1; i<tree1->nnodes; i++)
        mapping[i] = -1;

    // calculate mapping as much as possible
    int order[tree1->nnodes];
    tree1->get_postorder(order);
    LocalNode *nodes = tree1->nodes;
    for (int i=0; i<tree1->nnodes; i++) {
        const int j = order[i];
        const int *child = nodes[j].child;
        if (!nodes[j].is_leaf() && mapping[child[0]] != -1 &&
            mapping[child[1]] != -1) {
            // both children mapping, see if they share parent
            int a = tree2->nodes[mapping[child[0]]].parent;
            int b = tree2->nodes[mapping[child[1]]].parent;
            if (a == b)
                mapping[j] = a;
        }
    }

    // use mapping to find important nodes
    // at least the recombination node should be mapped
    int broken = tree1->nodes[recomb_node].parent;
    int other = tree1->get_sibling(recomb_node);
    int recomb = mapping[recomb_node];
    assert(recomb != -1);
    int recoal = tree2->nodes[recomb].parent;

    // map remaining nodes
    for (int i=0; i<tree1->nnodes; i++) {
        const int j = order[i];
        if (!nodes[j].is_leaf() && j != broken) {
            int a = nodes[j].child[0];
            int b = nodes[j].child[1];
            // skip over broken node
            if (a == broken) a = other;
            if (b == broken) b = other;
            int c = mapping[a];
            int d = mapping[b];
            c = tree2->nodes[c].parent;
            d = tree2->nodes[d].parent;
            // skip over recoal node
            if (c == recoal) c = tree2->nodes[c].parent;
            if (d == recoal) d = tree2->nodes[d].parent;
            assert(c == d);
            mapping[j] = c;
        }
    }

    // ensure broken node maps to nothing
    mapping[broken] = -1;
}


// Infer the SPR and mapping between two local trees.
// The local trees and recombination node and time must be correct.
// The population path should be correct as well.
// All other information is inferred.
void repair_spr(const LocalTree *last_tree, const LocalTree *tree, Spr &spr,
                int *mapping)
{
    // infer the mapping between local trees using the recombination node
    infer_mapping(last_tree, tree, spr.recomb_node, mapping);

    // determine the coal time
    int broken = last_tree->nodes[spr.recomb_node].parent;
    int recomb = mapping[spr.recomb_node];
    assert(recomb != -1);
    int recoal = tree->nodes[recomb].parent;
    spr.coal_time = tree->nodes[recoal].age;

    // determine the coal node
    int other = tree->get_sibling(recomb);
    int inv_mapping[tree->nnodes];
    get_inverse_mapping(mapping, tree->nnodes, inv_mapping);
    spr.coal_node = inv_mapping[other];

    // adjust coal node due to branch movement
    if (spr.coal_node == broken)
        spr.coal_node = last_tree->get_sibling(spr.recomb_node);
    int parent = last_tree->nodes[spr.coal_node].parent;
    if (parent != -1 && spr.coal_time > last_tree->nodes[parent].age)
        spr.coal_node = parent;
}



// appends the data in 'trees2' to 'trees'
// trees2 is then empty
// if merge is true, then merge identical neighboring local trees
void append_local_trees(LocalTrees *trees, LocalTrees *trees2, bool merge,
                        const PopulationTree *pop_tree)
{
    const int ntrees = trees->get_num_trees();
    const int ntrees2 = trees2->get_num_trees();

    // ensure seqids are the same
    for (unsigned int i=0; i<trees->seqids.size(); i++)
        assert(trees->seqids[i] == trees2->seqids[i]);
    assert(trees->nnodes == trees2->nnodes);

    // move trees2 onto end of trees
    LocalTrees::iterator it = trees->end();
    --it;
    trees->trees.splice(trees->end(),
                        trees2->trees, trees2->begin(), trees2->end());
    trees->end_coord = trees2->end_coord;
    trees2->end_coord = trees2->start_coord;

    // set the mapping the newly neighboring trees
    if (merge && ntrees > 0 && ntrees2 > 0) {
        LocalTrees::iterator it2 = it;
        ++it2;

        if (it2->spr.is_null()) {
            // there is no SPR between these trees
            // infer a congruent mapping and remove redunant local blocks
            if (it2->mapping == NULL)
                it2->mapping = new int [trees2->nnodes];
            map_congruent_trees(it->tree, &trees->seqids[0],
                                it2->tree, &trees2->seqids[0], it2->mapping);
            remove_null_spr(trees, it, pop_tree);
        } else {
            // there should be an SPR between these trees, repair it.
            repair_spr(it->tree, it2->tree, it2->spr, it2->mapping);
        }
    }

    //assert_trees(trees);
    //assert_trees(trees2);
}

void remove_population_paths(LocalTrees *trees) {
    for (LocalTrees::iterator it=trees->begin();
         it != trees->end(); ++it)
    {
        LocalTree *tree = it->tree;
        for (int i=0; i < tree->nnodes; i++)
            tree->nodes[i].pop_path = 0;
        Spr spr = it->spr;
        if (!spr.is_null())
            spr.pop_path = 0;
    }
}



//=============================================================================
// local tree alignment compression

void uncompress_local_trees(LocalTrees *trees,
                            const SitesMapping *sites_mapping)
{
    // get block lengths
    vector<int> blocklens;
    for (LocalTrees::iterator it=trees->begin(); it != trees->end(); ++it)
        blocklens.push_back(it->blocklen);

    // get uncompressed block lengths
    vector<int> blocklens2;
    sites_mapping->uncompress_blocks(blocklens, blocklens2);

    // apply block lengths to local trees
    int i = 0;
    for (LocalTrees::iterator it=trees->begin(); it != trees->end(); ++it, i++)
    {
        assert(blocklens2[i] > 0);
        it->blocklen = blocklens2[i];
    }

    trees->start_coord = sites_mapping->old_start;
    trees->end_coord = sites_mapping->old_end;

    //assert_trees(trees);
}


// compress local trees according to sites_mapping
void compress_local_trees(LocalTrees *trees, const SitesMapping *sites_mapping,
                          bool fuzzy)
{
    // get block lengths
    vector<int> blocklens;
    for (LocalTrees::iterator it=trees->begin(); it != trees->end(); ++it)
        blocklens.push_back(it->blocklen);

    // get compressed block lengths
    vector<int> blocklens2;
    sites_mapping->compress_blocks(blocklens, blocklens2);

    // enfore non-zero length blocks
    for (unsigned int i=0; i<blocklens2.size(); i++) {
        if (fuzzy && blocklens2[i] <= 0) {
            // shift block end to the right and compensate in next block
	    int diff = 1 - blocklens2[i];
	    blocklens2[i] = 1;
	    if (i < blocklens2.size() - 1) {
		blocklens2[i+1] -= diff;
	    } else {
		int j=i-1;
		for (; j >= 0; j--) {
		    if (blocklens2[j] > 1) {
			int remove = min(diff, blocklens2[j]-1);
			blocklens2[j] -= remove;
			diff -= remove;
			if (diff == 0) break;
		    }
		}
		if (j < 0)  {
		    fprintf(stderr, "Unable to compress local trees\n");
		    exit(1);
		}
	    }
        } else
            assert(blocklens2[i] > 0);
    }

    // apply new block lengths to local trees
    int i = 0;
    for (LocalTrees::iterator it=trees->begin(); it != trees->end(); ++it, ++i)
        it->blocklen = blocklens2[i];

    trees->start_coord = sites_mapping->new_start;
    trees->end_coord = sites_mapping->new_end;
}


// assert that local trees uncompress correctly
void assert_uncompress_local_trees(LocalTrees *trees,
                                   const SitesMapping *sites_mapping)
{
    vector<int> blocklens;

    for (LocalTrees::iterator it=trees->begin(); it != trees->end(); ++it)
        blocklens.push_back(it->blocklen);

    uncompress_local_trees(trees, sites_mapping);
    compress_local_trees(trees, sites_mapping);

    int i = 0;
    int pos = 0;
    for (LocalTrees::iterator it=trees->begin(); it != trees->end(); ++it, i++){
        int blocklen = it->blocklen;
        assert(blocklens[i] == blocklen);
        pos += blocklen;
    }
}



//=============================================================================
// local tree newick output

// write out the newick notation of a tree
void write_newick_node(FILE *out, const LocalTree *tree,
                       const char *const *names,
                       const double *times, int node, int depth, bool oneline,
                       bool pop_model)
{   
    //printLog(LOG_LOW, "accessing node %d\n", node);
    if (tree->nodes[node].is_leaf()) {
        if (!oneline)
            for (int i=0; i<depth; i++) fprintf(out, "  ");
        fprintf(out, "%s:%f[&&NHX:age=%f", names[node],
                tree->get_dist(node, times), times[tree->nodes[node].age]);
        if (pop_model)
            fprintf(out, ":pop_path=%i", tree->nodes[node].pop_path);
        fprintf(out, "]");
    } else {
        // indent
        if (oneline) {
            fprintf(out, "(");
        } else {
            for (int i=0; i<depth; i++) fprintf(out, "  ");
            fprintf(out, "(\n");
        }

        write_newick_node(out, tree, names, times,
                          tree->nodes[node].child[0], depth+1, oneline,
                          pop_model);
        if (oneline)
            fprintf(out, ",");
        else
            fprintf(out, ",\n");

        write_newick_node(out, tree, names, times,
                          tree->nodes[node].child[1], depth+1, oneline,
                          pop_model);
        if (!oneline) {
            fprintf(out, "\n");
            for (int i=0; i<depth; i++) fprintf(out, "  ");
        }
        fprintf(out, ")");

        if (depth > 0)
            fprintf(out, "%s:%f[&&NHX:age=%f",
                    names[node], tree->get_dist(node, times),
                    times[tree->nodes[node].age]);
        else
            fprintf(out, "%s[&&NHX:age=%f", names[node],
                    times[tree->nodes[node].age]);
        if (pop_model)
            fprintf(out, ":pop_path=%i", tree->nodes[node].pop_path);
        fprintf(out, "]");
    }
}


// write out the newick notation of a tree to a stream
void write_newick_tree(FILE *out, const LocalTree *tree,
                       const char *const *names,
                       const double *times, int depth, bool oneline,
                       bool pop_model)
{
    // setup default names
    char **names2 = (char **) names;
    char **default_names = NULL;
    if (names == NULL) {
        default_names = new char* [tree->nnodes];
        for (int i=0; i<tree->nnodes; i++) {
            default_names[i] = new char [16];
            snprintf(default_names[i], 15, "%d", i);
        }
        names2 = default_names;
    }


    write_newick_node(out, tree, names2, times, tree->root, 0, oneline, pop_model);
    if (oneline)
        fprintf(out, ";");
    else
        fprintf(out, ";\n");

    // clean up default names
    if (default_names) {
        for (int i=0; i<tree->nnodes; i++)
            delete [] default_names[i];
        delete [] default_names;
    }
}

string get_newick_rep_rSPR(const LocalTree *tree){
    string s;
    get_newick_rep_rSPR_helper(&s, tree, tree->root);
    s += ";";
    return s;
}

void get_newick_rep_rSPR_helper(string *s, const LocalTree *tree, int node){
    if (tree->nodes[node].is_leaf()){
        char str[50];
        sprintf(str, "%d", node);
        *s += str;
    }else{
        *s += "(";
        get_newick_rep_rSPR_helper(s, tree, tree->nodes[node].child[0]);
		*s += ",";
        get_newick_rep_rSPR_helper(s, tree, tree->nodes[node].child[1]);
		*s += ")";
	}
}


void write_newick_node_rSPR(FILE *out, const LocalTree *tree,
                        const char *const *names,
                       const double *times, int node, int depth)
{
    if (tree->nodes[node].is_leaf()) {
        fprintf(out, "%s", names[node]);
    } else {
        // for rSPR, don't print internal nodes
        // or do we want this?
        // looks like Chris's code assumes internal node's label comes from its pair of parenthesis,
        // which seems inconsistent with the convention, but whatever.
        
        fprintf(out, "(");
        write_newick_node_rSPR(out, tree, names, times,
                          tree->nodes[node].child[0], depth+1);
        fprintf(out, ",");
        write_newick_node_rSPR(out, tree, names, times,
                          tree->nodes[node].child[1], depth+1);
        fprintf(out, ")");
        //fprintf(out, "%s", names[node]);
    }
}

// write out newick notation of a local tree to a file stream for rSPR program to use
// this is a simplified newick format
void write_newick_tree_rSPR(FILE *out, const LocalTree *tree, const double *times)
{   
    // setup default names
    char **default_names;
    default_names = new char* [tree->nnodes];
    for (int i=0; i<tree->nnodes; i++) {
        default_names[i] = new char [16];
        snprintf(default_names[i], 15, "%d", i);
    }

    write_newick_node_rSPR(out, tree, default_names, times, tree->root, 0);
    fprintf(out, ";\n");
    // clean up default names
    if (default_names) {
        for (int i=0; i<tree->nnodes; i++)
            delete [] default_names[i];
        delete [] default_names;
    }


}

// write out the newick notation of a tree to a file
bool write_newick_tree(const char *filename, const LocalTree *tree,
                       const char *const *names, const double *times,
                       bool oneline, bool pop_model)
{
    FILE *out = NULL;

    if ((out = fopen(filename, "w")) == NULL) {
        printError("cannot write file '%s'\n", filename);
        return false;
    }

    write_newick_tree(out, tree, names, times, 0, oneline, pop_model);
    fclose(out);
    return true;
}


bool write_newick_tree_for_bedfile_recur(FILE *out, const LocalTree *tree,
                                         const char *const *names,
                                         const ArgModel *model,
                                         const Spr &spr, int node) {

    vector<string> nhx;
    const double *times = model->times;
    char tmpstr[1000];
    if (tree->nodes[node].is_leaf()) {
        fprintf(out, "%s", names[node]);
    } else {
        fprintf(out, "(");
        write_newick_tree_for_bedfile_recur(out, tree, names, model, spr,
                                            tree->nodes[node].child[0]);
        fprintf(out, ",");
        write_newick_tree_for_bedfile_recur(out, tree, names, model, spr,
                                            tree->nodes[node].child[1]);
        fprintf(out, ")");
    }
    if (node != tree->root) {
        int parent = tree->nodes[node].parent;
        fprintf(out, ":%.1f", times[tree->nodes[parent].age] -
                times[tree->nodes[node].age]);
    }

    if (model->pop_tree != NULL && tree->nodes[node].pop_path != 0) {
        sprintf(tmpstr, "pop_path=%i", tree->nodes[node].pop_path);
        nhx.push_back(string(tmpstr));
    }
    if (node == spr.recomb_node) {
        sprintf(tmpstr, "recomb_time=%.1f", times[spr.recomb_time]);
        nhx.push_back(string(tmpstr));
    }
    if (node == spr.coal_node) {
        sprintf(tmpstr, "coal_time=%.1f", times[spr.coal_time]);
        nhx.push_back(string(tmpstr));
    }
    if (node == spr.recomb_node && model->pop_tree != NULL && spr.pop_path != 0) {
        sprintf(tmpstr, "spr_pop_path=%i", spr.pop_path);
        nhx.push_back(string(tmpstr));
    }
    if (nhx.size() > 0) {
        fprintf(out, "[&&NHX:%s", nhx[0].c_str());
        for (int i=1; i < (int)nhx.size(); i++)
            fprintf(out, ",%s", nhx[i].c_str());
        fprintf(out, "]");
    }
    return true;
}


bool write_newick_tree_for_bedfile(FILE *out,
                                   const LocalTree *tree,
                                   const char *const *names,
                                   const ArgModel *model,
                                   const Spr &spr) {
    write_newick_tree_for_bedfile_recur(out, tree, names, model, spr,
                                        tree->root);
    fprintf(out, ";");
    return true;
}


//=============================================================================
// read local tree

// find closest time in times array
int find_time(double time, const double *times, int ntimes)
{
    double mindiff = INFINITY;
    int mini = -1;

    for (int i=0; i<ntimes; i++) {
        double diff = fabs(times[i] - time);
        if (diff < mindiff) {
            mindiff = diff;
            mini = i;
        }
    }
    assert(mini != -1);

    return mini;
}


// Iterates through the key-value pairs of a NHX comment
// NOTE: end is exclusive
// Example: "key1=value2:key2=value2"
bool iter_nhx_key_values(char *text, char *end,
                         char **key, char **key_end,
                         char **value, char **value_end)
{
    if (*key >= end)
        return false;

    // parse key
    *key_end = *key;
    while (**key_end != '=') {
        if (*key_end == end)
            return false;
        (*key_end)++;
    }

    // parse value
    *value = *key_end + 1;
    *value_end = *value;
    while (*value_end < end && **value_end != ':') (*value_end)++;

    return true;
}


// Parse the node age from a string 'text'
// NOTE: end is exclusive
// Example: "&&NHX:age=20"
bool parse_node_age(char* text, char *end, double *age)
{
    // ensure comment begins with "&&NHX:"
    if (strncmp(text, "&&NHX:", 6) != 0) {
        return false;
    }

    text += 6;

    char *key = text;
    char *key_end, *value, *value_end;
    while (iter_nhx_key_values(text, end, &key, &key_end, &value, &value_end)){
        if (strncmp(key, "age", 3) == 0 && key_end - key == 3) {
            if (sscanf(value, "%lf", age) != 1)
                return false;
            else {
                return true;
            }
        }

        key = value_end + 1;
    }

    return false;
}


// Parse the node age from a string 'text'
// NOTE: end is exclusive
// Example: "&&NHX:pop_path=1"
bool parse_node_pop_path(char* text, char *end, int *pop_path)
{
    *pop_path = 0; // default pop_path
    // ensure comment begins with "&&NHX:"
    if (strncmp(text, "&&NHX:", 6) != 0) {
        return false;
    }

    text += 6;

    char *key = text;
    char *key_end, *value, *value_end;
    while (iter_nhx_key_values(text, end, &key, &key_end, &value, &value_end)){
        if (strncmp(key, "pop_path", 8) == 0 && key_end - key == 8) {
            if (sscanf(value, "%d", pop_path) != 1)
                return false;
            else {
                return true;
            }
        }

        key = value_end + 1;
    }
    return false;
}


// Parses a local tree from a newick string
bool parse_local_tree(const char* newick, LocalTree *tree,
                      const double *times, int ntimes)
{
    const int len = strlen(newick);
    vector<int> ptree;
    vector<int> ages;
    vector<int> stack;
    vector<int> names;
    vector<int> pop_paths;

    // create root node
    ptree.push_back(-1);
    ages.push_back(-1);
    names.push_back(-1);
    pop_paths.push_back(0);
    int node = 0;

    for (int i=0; i<len; i++) {
        switch (newick[i]) {
        case '(': // new branchset
            ptree.push_back(node);
            ages.push_back(-1);
            pop_paths.push_back(0);
            names.push_back(-1);
            stack.push_back(node);
            node = ptree.size() - 1;
            break;

        case ',': // another branch
            ptree.push_back(stack.back());
            ages.push_back(-1);
            pop_paths.push_back(0);
            names.push_back(-1);
            node = ptree.size() - 1;
            break;

        case ')': // optional name next
            node = stack.back();
            stack.pop_back();
            break;

        case ':': // optional dist next
            break;

        case '[': { // comment next
            int j = i + 1;
            while (j<len && newick[j] != ']') j++;

            double age;
            if (newick[j] == ']') {
                // parse age field
                if (parse_node_age((char*) &newick[i+1],
                                   (char*) &newick[j], &age))
                    ages[node] = find_time(age, times, ntimes);
                parse_node_pop_path((char*) &newick[i+1],
                                    (char*) &newick[j], &(pop_paths[node]));
                i = j;
            } else {
                // error, quit early
                printError("bad newick: malformed NHX comment");
                i = len;
            }
            } break;


        case ';':
            break;

        default: {
            char last = newick[i-1];

            // skip leading whitespace
            while (newick[i] == ' ') i++;
            int j = i;
            // find end of token
            while (j < len && !inChars(newick[j], ")(,:;[")) j++;

            if (last == ')' || last == '(' || last == ',') {
                // name
                if (sscanf(&newick[i], "%d", &names[node]) != 1) {
                    // error, quit early
                    printError("bad newick: node name is not an integer");
                    i = len;
                }
            } else if (last == ':') {
                // ignore distance
            }

            i = j - 1;
        }
        }
    }

    if (stack.size() != 0)
        return false;


    // fill in local tree data structure
    int nnodes = ptree.size();
    tree->clear();

    // add nodes to tree
    int *order = &names[0];

    tree->ensure_capacity(nnodes);
    tree->nnodes = nnodes;

    for (int i=0; i<nnodes; i++) {
        int j = order[i];
        if (j == -1) {
            printError("unexpected error (%d)", i);
            return false;
        }
        if (ptree[i] != -1)
            tree->nodes[j].parent = order[ptree[i]];
        else {
            tree->nodes[j].parent = -1;
            tree->root = j;
        }
        tree->nodes[j].age = ages[i];
        tree->nodes[j].child[0] = -1;
        tree->nodes[j].child[1] = -1;
        tree->nodes[j].pop_path = pop_paths[i];
    }

    // set children
    for (int i=0; i<nnodes; i++) {
        if (ptree[i] != -1) {
            if (tree->add_child(order[ptree[i]], order[i]) == -1) {
                printError("local tree is not binary");
                return false;
            }
        }
    }

    // leaves default to age 0
    for (int i=0; i<nnodes; i++)
        if (tree->nodes[i].is_leaf() && tree->nodes[i].age == -1)
            tree->nodes[i].age = 0;

    // check for valid tree structure
    if (!assert_tree(tree))
        return false; 

    return true;
}


//=============================================================================
// output ARG as local trees


void write_local_tree(const LocalTree *tree) {
    for (int i=0; i < tree->nnodes; i++) {
        printf("node %i: parent=%i child=(%i,%i) age=%i path=%i\n", i, tree->nodes[i].parent,
               tree->nodes[i].child[0], tree->nodes[i].child[1],
               tree->nodes[i].age, tree->nodes[i].pop_path);
    }
}


void write_local_trees_as_bed(FILE *out, const LocalTrees *trees,
                              const vector<string> seqnames,
                              const ArgModel *model, int sample) {
    const int nnodes = trees->nnodes;
    char **nodeids = new char* [nnodes];
    int i = 0;
    Spr spr;

    for (i=0; i<trees->get_num_leaves(); i++) {
        nodeids[i] = new char [seqnames[trees->seqids[i]].length()+1];
        strcpy(nodeids[i], seqnames[trees->seqids[i]].c_str());
    }
    for (; i<nnodes; i++) {
        nodeids[i] = new char[1];
        nodeids[i][0]='\0';
    }

    int end = trees->start_coord;
    for (LocalTrees::const_iterator it=trees->begin();
         it != trees->end(); ++it)
    {
        int start = end;
        end += it->blocklen;
        assert(it->blocklen > 0);
        LocalTree *tree = it->tree;

        if (end - start > 0) {
            fprintf(out, "%s\t%i\t%i\t%i\t",
                    trees->chrom.c_str(), start, end, sample);

            LocalTrees::const_iterator it2 = it;
            ++it2;
            if (it2 != trees->end()) {
                spr = it2->spr;
            } else {
                spr.set_null();
            }

            write_newick_tree_for_bedfile(out, tree, nodeids, model, spr);
            fprintf(out, "\n");
        }
    }

    // cleanup
    for (int i=0; i<nnodes; i++)
        delete [] nodeids[i];
    delete [] nodeids;
}

void remove_edge(map<pair<tsk_id_t, tsk_id_t>, int> *edges, tsk_table_collection_t *tables, 
    tsk_id_t p, tsk_id_t c, int coord){
    if (tables->nodes.time[p] == tables->nodes.time[c]) {
        //printLog(LOG_LOW, "remove edge rejected\n");
        return;
    } // this is the zero-length branch case
    //remove the edge whose child node is tskit_id

    auto it = edges->find(make_pair(p, c));
    if (it == edges->end()){
        return;
        //printLog(LOG_LOW, "something wrong with the code: removing non-existing edge\n");
    }else{
        int ret = tsk_edge_table_add_row(&(tables->edges), it->second, coord, p, c);
        check_tsk_error(ret);
        //printLog(LOG_LOW, "removed successfully\n");
        edges->erase(make_pair(p, c));
    }
    
}

void insert_edge(map<pair<tsk_id_t, tsk_id_t>, int> *edges, tsk_id_t p, tsk_id_t c, 
        int coord, tsk_table_collection_t *tables){
    if (tables->nodes.time[p] == tables->nodes.time[c]) {
        //printLog(LOG_LOW, "insert edge rejected\n");    
        return;
    }
    auto it = edges->find(make_pair(p, c));
    if(it == edges->end()){
        edges->insert(make_pair(make_pair(p, c), coord));
    }else{
        return;
        //printLog(LOG_LOW, "something wrong with the code: duplicate insertion\n");
    }
    //printLog(LOG_LOW, "inserted successfully\n");
}

int init_nodes_mapping(const LocalTree *tree, int *nodes, tsk_table_collection_t *tables, 
        map<pair<tsk_id_t, tsk_id_t>, int> *edges, vector<tsk_id_t*> *node_maps, 
        const double *times, int start_coord){
    int num_samples = tree->get_num_leaves();
    int j, u, ret;
    bool visited[tree->nnodes];
    fill(visited, visited + tree->nnodes, false);

    for(int i = 0; i < num_samples; i++){
        ret = tsk_node_table_add_row(&(tables->nodes), 1, times[tree->get_node(i).age],
                TSK_NULL, TSK_NULL, NULL, 0);
        check_tsk_error(ret);
    }

    int counter = num_samples;
    for(j = 0; j < num_samples; j++){
        visited[j] = true;
        nodes[j] = j;
        u = j;
        while(u != -1){
            int p = tree->get_node(u).parent;
            if (p == -1) {break;}
            if (tree->get_node(p).age == tree->get_node(u).age){
                nodes[p] = nodes[u];
            }else{
                if (!visited[p]) {
                    nodes[p] = counter++;
                    ret = tsk_node_table_add_row(&(tables->nodes), 0, times[tree->get_node(p).age],
                            TSK_NULL, TSK_NULL, NULL, 0);
                    check_tsk_error(ret);
                }
                edges->insert(make_pair(make_pair(nodes[p], nodes[u]), start_coord));
            }
            if (visited[p]) {break;}
            visited[p] = true;
            u = p;
        }
    }

    tsk_id_t *tmp = new tsk_id_t[tree->nnodes];
    for(int i = 0; i < tree->nnodes; i++){
        tmp[i] = nodes[i];
    }
    node_maps->push_back(tmp);
    return counter;
}


// a simplistic version of writing tree sequence file
// for now, ignore population/sequence name/individuals, etc
void write_local_trees_ts(const char *filename, const LocalTrees *trees, 
        const Sequences *sequences, const SitesMapping *sitesmapping, const double *times){

#ifdef DEBUG
    printLog(LOG_LOW, "discrete time points\n");
    int i = 0;
    for(; i < 20; i++){
        printLog(LOG_LOW, "%lf\n", times[i]);
    }
#endif
    int ret;
    tsk_table_collection_t tables;
    ret = tsk_table_collection_init(&tables, 0);
    check_tsk_error(ret);
    // need to add start_coord otherwise tskit throws " Right coordinate > sequence length"
    // better if tskit can add an additional field for the table struct to indicate the starting position of a sequence
    tables.sequence_length = trees->length() + trees->start_coord;

    //add nodes from the first tree
    // this map is a bit complicated
    // the first int refers to the starting position of this edge
    // the second int refers to the multiplicity of this edge
    // in some corner cases, <parent, child> is not enough to fully identify an edge
    
    int nnodes = trees->nnodes;
    int nodes[nnodes];
    int coord = trees->start_coord;
    map<pair<tsk_id_t, tsk_id_t>, int> edges;
    vector<tsk_id_t*> node_maps;
    int id = init_nodes_mapping(trees->trees.front().tree, nodes, &tables, &edges, &node_maps, 
                                    times, trees->start_coord);

    LocalTree prev;
    int tree_id = 0;
    for (LocalTrees::const_iterator it=trees->begin(); it != trees->end(); ++it){
        //printLog(LOG_LOW, "at tree %d\n", tree_id);
        tree_id++;

        if(it == trees->begin()) {
            coord += it->blocklen;
            prev.copy(*(it->tree)); 
            continue;
        }
        //add the new node to the node table
        int new_node = get_recoal_node(it->tree, it->spr, it->mapping);
        int recomb_node = it->spr.recomb_node; // this is the id as in argweaver, not tskit
        int p = prev.get_node(recomb_node).parent;
        int sib = prev.get_sibling(recomb_node);
        // three cases: 2 edge removed, 2 edge added
        // or 3 edge removed, 3 edge added
        // or 4 edge removed, 4 edge added

        //printLog(LOG_LOW, "remove 1 attempted(%d->%d)\n", nodes[p], nodes[recomb_node]);
        remove_edge(&edges, &tables, nodes[p], nodes[recomb_node], coord);
        
        // recomb node is never the root node, so we always get a valid sibling
        //printLog(LOG_LOW, "revmoe 2 attempted(%d->%d)\n", nodes[p], nodes[sib]);
        remove_edge(&edges, &tables, nodes[p], nodes[sib], coord);
        
        if (prev.root != p){
            //printLog(LOG_LOW, "remove 3 attempted(%d->%d)\n", nodes[prev.get_node(p).parent], nodes[p]);
            remove_edge(&edges, &tables, nodes[prev.get_node(p).parent], nodes[p], coord);
        }
        if (it->spr.coal_node != sib && it->spr.coal_node != p && it->spr.coal_node != prev.root){
            //printLog(LOG_LOW, "remove 4 attempted(%d->%d)\n", nodes[prev.get_node(it->spr.coal_node).parent], 
            //        nodes[it->spr.coal_node]);
            remove_edge(&edges, &tables, nodes[prev.get_node(it->spr.coal_node).parent], 
                    nodes[it->spr.coal_node], coord);
        }

        //now added newly created edge in this local tree to the map
        int *child = it->tree->get_node(new_node).child;
        int parent_of_new_node = it->tree->get_node(new_node).parent;
        int new_node_age = it->tree->get_node(new_node).age;

        int tmp[nnodes];
        memcpy(tmp, nodes, nnodes*sizeof(int));
        for (int i = 0; i < nnodes; i++){
            if(it->mapping[i] != -1){
                nodes[it->mapping[i]] = tmp[i];
            }
        }

        int child1_age = it->tree->get_node(child[0]).age;
        int child2_age = it->tree->get_node(child[1]).age;
        if(new_node_age == child1_age){
            nodes[new_node] = nodes[child[0]];
        }else if(new_node_age == child2_age){
            nodes[new_node] = nodes[child[1]];
        }else if(parent_of_new_node != -1 && new_node_age == it->tree->get_node(parent_of_new_node).age){
            nodes[new_node] = nodes[parent_of_new_node];
        }else{
                ret = tsk_node_table_add_row(&(tables.nodes), 0, times[it->tree->get_node(new_node).age], 
                        TSK_NULL, TSK_NULL, NULL, 0);
                check_tsk_error(ret);
                nodes[new_node] = id++; 
        }

        //printLog(LOG_LOW, "insert 1 attempted(%d->%d)\n", nodes[new_node], 
        //        nodes[it->tree->get_sibling(it->mapping[recomb_node])]);
        insert_edge(&edges, nodes[new_node], nodes[it->tree->get_sibling(it->mapping[recomb_node])], 
                    coord, &tables);
        //printLog(LOG_LOW, "insert 2 attempted(%d->%d)\n", nodes[new_node], nodes[it->mapping[recomb_node]]);
        insert_edge(&edges, nodes[new_node], nodes[it->mapping[recomb_node]], coord, &tables);

        if(it->tree->root != new_node){
            //printLog(LOG_LOW, "insert 3 attempted(%d->%d)\n", nodes[it->tree->get_node(new_node).parent], nodes[new_node]);
            insert_edge(&edges, nodes[it->tree->get_node(new_node).parent], nodes[new_node], coord, &tables);
        }

        if (it->spr.coal_node != sib && it->spr.coal_node != prev.get_node(sib).parent
                && prev.root != p){
            //printLog(LOG_LOW, "insert 4 attempted(%d->%d)\n", nodes[it->mapping[prev.get_node(p).parent]], 
            //        nodes[it->mapping[sib]]);
            insert_edge(&edges, nodes[it->mapping[prev.get_node(p).parent]], 
                    nodes[it->mapping[sib]], coord, &tables);
        }      
        coord += it->blocklen;//set up starting coordinate of next tree
        
        tsk_id_t* curr_map = new tsk_id_t[nnodes];
        for(int j=0; j < nnodes; j++){
            curr_map[j] = nodes[j];
        }
        node_maps.push_back(curr_map);
        prev.copy(*(it->tree));
    }

    // flush out edges that remain at the very end of the tree sequence
    for(auto it = edges.begin(); it != edges.end(); ++it){
        if (tables.nodes.time[it->first.first] != tables.nodes.time[it->first.second]){
            ret = tsk_edge_table_add_row(&(tables.edges), it->second, coord, 
                    it->first.first, it->first.second);
            check_tsk_error(ret);    
        }
    }

    ret = tsk_table_collection_sort(&tables, NULL, 0);
    check_tsk_error(ret);

    // add site and mutation information
    Sites sites;
    make_sites_from_sequences(sequences, &sites);
    uncompress_sites(&sites, sitesmapping);
    int nseqs = sites.get_num_seqs();
    if (sites.ref.size() != sites.get_num_sites() || sites.alt.size() != sites.get_num_sites()){
        printLog(LOG_LOW, "Can't output tree sequecne without ancestral allele info for every SNP site\n");
        exit(EXIT_FAILURE);
    }
    
    int end = trees->start_coord + trees->trees.front().blocklen; // this is the right coordinate of the first tree
    int curr_tree = 0;
    LocalTrees::const_iterator it = trees->begin();
    for(int i = 0; i < sites.get_num_sites(); i++){
        ret = tsk_site_table_add_row(&(tables.sites), sites.positions[i], 
                        &(sites.ref[i]), sizeof(char), NULL, NULL);
        check_tsk_error(ret);
        //printLog(LOG_LOW, "sites %d, ref: %c, alt: %c\n", i, sites.ref[i], sites.alt[i]);
        // add mutation
        int site_pos = sites.positions[i];
        while(end < site_pos && it != trees->end()){
            it++;
            curr_tree++;
            end += it->blocklen;
        }

        //printLog(LOG_LOW, "current tree: %d\n", curr_tree);
        // now we determine which branch should this mutation reside
        // may not be unique, may not be completely compatible ...
        tsk_id_t *node_map = node_maps[curr_tree];
        char *site = sites.cols[i];
        set<tsk_id_t> derived;
        for(int k = 0; k < nseqs; k++){
            if(site[k] == sites.alt[i]){derived.insert(k);}
        }
        
        map<int, set<tsk_id_t>*> descent_map;
        int postorder[nnodes];
        it->tree->get_postorder(postorder);
        bool mapped = false;
        for(int j = 0; j < nnodes; j++){
            //determine the set of leaf nodes that carry the derived alleles
            int node = postorder[j];
            //printLog(LOG_LOW, "postorder: %d\n", node);
            if (node < it->tree->get_num_leaves()){
                descent_map.insert(make_pair(node, new set<tsk_id_t>{node}));
            }else{
                int *child = it->tree->nodes[node].child;
                auto descent_union = new set<tsk_id_t>(*descent_map.at(child[0]));
                descent_union->insert(descent_map.at(child[1])->begin(), descent_map.at(child[1])->end());
                descent_map.insert(make_pair(node, descent_union));
            }

            if(*descent_map.at(node) == derived){
                mapped = true;
                ret = tsk_mutation_table_add_row(&(tables.mutations), i, node_map[node], 
                        -1, &(sites.alt[i]), sizeof(char), NULL, 0);
                check_tsk_error(ret);
                //printLog(LOG_LOW, "map mutation to %d at tree %d\n", node, curr_tree);
                break;
            }
        }
        if(!mapped){
            printLog(LOG_LOW, "can't unambiguous map mutation at site %d\n", sites.positions[i]);
        }

        // clean up
        for(int j = 0; j < nnodes; j++){
            delete descent_map[j];
        }

    }

    ret = tsk_table_collection_dump(&tables, filename, 0);
    check_tsk_error(ret);
    tsk_table_collection_free(&tables);
    for(int i = 0; i < trees->get_num_trees(); i++){
        delete [] node_maps[i];
    }

}


void write_local_trees(FILE *out, const LocalTrees *trees,
                       const char *const *names, const double *times,
                       bool pop_model, const vector<int> &self_recomb_pos,
                       const vector<Spr> &self_recombs)
{
    const int nnodes = trees->nnodes;
    const int nodeid_len = 10;

    assert(self_recomb_pos.size() == self_recombs.size());

    // print names
    if (names) {
        fprintf(out, "NAMES");
        for (int i=0; i<trees->get_num_leaves(); i++)
            fprintf(out, "\t%s", names[trees->seqids[i]]);
        fprintf(out, "\n");
    }

    // print region, convert to 1-index
    fprintf(out, "REGION\t%s\t%d\t%d\n",
            trees->chrom.c_str(), trees->start_coord + 1, trees->end_coord);

    int next_self_pos = self_recomb_pos.size() == 0 ?
        trees->end_coord + 1 : self_recomb_pos[0];
    int self_idx = 0;


    // setup nodeids
    char **nodeids = new char* [nnodes];
    int *total_mapping = new int [nnodes];
    int *tmp_mapping = new int [nnodes];
    for (int i=0; i<nnodes; i++) {
        nodeids[i] = new char [nodeid_len + 1];
        total_mapping[i] = i;
    }

    int end = trees->start_coord;

    for (LocalTrees::const_iterator it=trees->begin();
         it != trees->end(); ++it)
    {
        int start = end;
        end += it->blocklen;
        LocalTree *tree = it->tree;

        // compute nodeids
        for (int i=0; i<nnodes; i++)
            snprintf(nodeids[i], nodeid_len, "%d", total_mapping[i]);

        // write tree
        // convert to 1-index
        fprintf(out, "TREE\t%d\t%d\t", start+1, end);
        write_newick_tree(out, tree, nodeids, times, 0, true, pop_model);
        fprintf(out, "\n");

        while (next_self_pos < end) {
            assert(self_idx < (int)self_recomb_pos.size());
            fprintf(out, "SPR-INVIS\t%d\t%d\t%f\t%d\t%f\t%i\n",
                    next_self_pos + 1,
                    total_mapping[self_recombs[self_idx].recomb_node],
                    times[self_recombs[self_idx].recomb_time],
                    total_mapping[self_recombs[self_idx].recomb_node],
                    times[self_recombs[self_idx].coal_time],
                    self_recombs[self_idx].pop_path);
            self_idx++;
            if (self_idx < (int)self_recomb_pos.size())
                next_self_pos = self_recomb_pos[self_idx];
            else next_self_pos = trees->end_coord + 1;
        }

        LocalTrees::const_iterator it2 = it;
        ++it2;
        if (it2 != trees->end()) {
            // write SPR
            const Spr &spr = it2->spr;
            fprintf(out, "SPR\t%d\t%d\t%f\t%d\t%f", end,
                    total_mapping[spr.recomb_node], times[spr.recomb_time],
                    total_mapping[spr.coal_node], times[spr.coal_time]);
            if (pop_model)
                fprintf(out, "\t%i", spr.pop_path);
            fprintf(out, "\n");

            // update total mapping
            int *mapping = it2->mapping;
            for (int i=0; i<nnodes; i++)
                tmp_mapping[i] = total_mapping[i];
            for (int i=0; i<nnodes; i++) {
                if (mapping[i] != -1)
                    total_mapping[mapping[i]] = tmp_mapping[i];
                else {
                    int recoal = get_recoal_node(tree, spr, mapping);
                    total_mapping[recoal] = tmp_mapping[i];
                }
            }
        }
    }

    // clean nodeids
    for (int i=0; i<nnodes; i++)
        delete [] nodeids[i];
    delete [] nodeids;
    delete [] tmp_mapping;
    delete [] total_mapping;
}


bool write_local_trees(const char *filename, const LocalTrees *trees,
                       const char *const *names, const double *times,
                       bool pop_model, const vector<int> &self_recomb_pos,
                       const vector<Spr> &self_recombs)
{
    FILE *out = NULL;

    if ((out = fopen(filename, "w")) == NULL) {
        printError("cannot write file '%s'\n", filename);
        return false;
    }

    write_local_trees(out, trees, names, times, pop_model,
                      self_recomb_pos, self_recombs);
    fclose(out);
    return true;
}


void write_local_trees(FILE *out, const LocalTrees *trees,
                       const Sequences &seqs,
                       const double *times, bool pop_model,
                       const vector<int> &self_recomb_pos,
                       const vector<Spr> &self_recombs)
{
    // setup names
    char **names;
    const unsigned int nleaves = trees->get_num_leaves();
    names = new char* [nleaves];

    for (unsigned int i=0; i<nleaves; i++) {
        if (i < seqs.names.size()) {
            names[i] = new char [seqs.names[i].size()+1];
            strncpy(names[i], seqs.names[i].c_str(), seqs.names[i].size()+1);
        } else {
            // use ids
            names[i] = new char [11];
            snprintf(names[i], 10, "%d", i);
        }
    }

    write_local_trees(out, trees, names, times, pop_model,
                      self_recomb_pos, self_recombs);

    // clean up names
    for (unsigned int i=0; i<nleaves; i++)
        delete [] names[i];
    delete [] names;
}


bool write_local_trees(const char *filename, const LocalTrees *trees,
                       const Sequences &seqs, const double *times, bool pop_model)
{
    FILE *out = NULL;

    if ((out = fopen(filename, "w")) == NULL) {
        printError("cannot write file '%s'\n", filename);
        return false;
    }

    write_local_trees(out, trees, seqs, times, pop_model);
    fclose(out);
    return true;
}


//=============================================================================
// read local trees


bool read_local_trees(FILE *infile, const double *times, int ntimes,
                      LocalTrees *trees, vector<string> &seqnames,
                      vector<int> *invisible_recomb_pos,
                      vector<Spr> *invisible_recombs)
{
    const char *delim = "\t";
    char *line = NULL;

    assert((invisible_recomb_pos==NULL && invisible_recombs==NULL) ||
           (invisible_recomb_pos!=NULL && invisible_recombs!=NULL));

    // init tree
    seqnames.clear();
    trees->clear();
    LocalTree *last_tree = NULL;

    int nnodes = 0;
    Spr spr;
    spr.set_null();

    int lineno = 1;
    while ((line = fgetline(infile))) {
        chomp(line);

        if (strncmp(line, "NAMES", 5) == 0) {
            // parse names
            split(&line[6], delim, seqnames);
            nnodes = 2 * seqnames.size() - 1;

        } else if (strncmp(line, "RANGE", 5) == 0) {
            // parse range
            printError("deprecated RANGE line detected, use REGION instead (line %d)", lineno);
            delete [] line;
            return false;

        } else if (strncmp(line, "REGION\t", 7) == 0) {
            // parse range
            char chrom[51];
            if (sscanf(&line[7], "%50s\t%d\t%d",
                       chrom,
                       &trees->start_coord, &trees->end_coord) != 3) {
                printError("bad REGION line (line %d)", lineno);
                delete [] line;
                return false;
            }
            trees->chrom = chrom;
            trees->start_coord--; // convert start to 0-index

        } else if (strncmp(line, "TREE", 4) == 0) {
            // parse tree
            int start, end;
            if (sscanf(&line[5], "%d\t%d", &start, &end) != 2) {
                printError("bad TREE line (line %d)", lineno);
                delete [] line;
                return false;
            }

            // find newick
            char *newick_end = line + strlen(line);
            char *newick = find(line+5, newick_end, delim[0]) + 1;
            newick = find(newick, newick_end, delim[0]) + 1;

            LocalTree *tree = new LocalTree(nnodes);
            if (!parse_local_tree(newick, tree, times, ntimes)) {
                printError("bad newick format (line %d)", lineno);
                delete tree;
                delete [] line;
                return false;
            }

            // setup mapping
            int *mapping = NULL;
            if (!spr.is_null()) {
                mapping = new int [nnodes];
                for (int i=0; i<nnodes; i++)
                    mapping[i] = i;
                if (spr.recomb_node != spr.coal_node)
                    mapping[last_tree->nodes[spr.recomb_node].parent] = -1;
            }

            // convert start to 0-index
            int blocklen = end - start + 1;
            trees->trees.push_back(LocalTreeSpr(tree, spr, blocklen, mapping));

            last_tree = tree;
        } else if (strncmp(line, "SPR-INVIS", 9) == 0) {
            if (invisible_recombs != NULL) {
                int pos, val;
                double recomb_time, coal_time;
                Spr ispr;
                ispr.pop_path = 0;
                val = sscanf(&line[10],
                             "%d\t%d\t%lf\t%d\t%lf\t%i",
                             &pos, &ispr.recomb_node, &recomb_time,
                             &ispr.coal_node, &coal_time, &ispr.pop_path);
                if (val != 5 && val != 6) {
                    printError("bad SPR-INVIS line (line %d)", lineno);
                    delete [] line;
                    return false;
                }
                ispr.recomb_time = find_time(recomb_time, times, ntimes);
                ispr.coal_time = find_time(coal_time, times, ntimes);
                invisible_recombs->push_back(ispr);
                invisible_recomb_pos->push_back(pos);
            }
            // for now just ignore these; could add argument to read them
            // into a separate object
        } else if (strncmp(line, "SPR", 3) == 0) {
            // parse SPR

            int pos, val;
            double recomb_time, coal_time;

            spr.pop_path = 0;
            val = sscanf(&line[4], "%d\t%d\t%lf\t%d\t%lf\t%i",
                         &pos, &spr.recomb_node, &recomb_time,
                         &spr.coal_node, &coal_time, &spr.pop_path);
            if (val != 5 && val != 6) {
                printError("bad SPR line (line %d)", lineno);
                delete [] line;
                return false;
            }

            spr.recomb_time = find_time(recomb_time, times, ntimes);
            spr.coal_time = find_time(coal_time, times, ntimes);
        }


        lineno++;
        delete [] line;
    }

    // set trees info
    if (trees->get_num_trees() > 0) {
        trees->nnodes = trees->front().tree->nnodes;
        trees->set_default_seqids();
    }

    return true;
}


bool read_local_trees(const char *filename, const double *times,
                      int ntimes, LocalTrees *trees, vector<string> &seqnames)
{
    FILE *infile = NULL;

    if ((infile = fopen(filename, "r")) == NULL) {
        printError("cannot read file '%s'\n", filename);
        return false;
    }

    bool result = read_local_trees(infile, times, ntimes, trees, seqnames);

    fclose(infile);
    return result;
}

/////////////////////////////////// read from tsinfer ////////////////////////////////////

bool read_local_tree_from_tsinfer(tsk_tree_t *tree, int *ptree, int *ages,
                                    const double *times, int ntimes, map<tsk_id_t, int> *id_map){
    // maybe first traverse the tree upwards as you would normally do
    tsk_id_t *samples = tree->samples;
    tsk_size_t num_samples = tsk_treeseq_get_num_samples(tree->tree_sequence);
    int nnodes = 2*num_samples-1;
    double age_tmp[nnodes];
    fill(age_tmp, age_tmp+nnodes, -1);
    // map parent to all its children
    auto pcmap = map<int, set<int>>();
    auto visited = set<int>();
    auto id_mapping = map<int, int>(); //mapping between tsk_id_t and the id as in argweaver
    int index = num_samples; // keep track of how many non-sample nodes have been visited
    for(int j = 0; j < num_samples; j++){
        int u = samples[j];
        visited.insert(u);
        tsk_tree_get_time(tree, u, &age_tmp[u]);
        id_mapping.insert(pair<int, int>(u, j));
        tsk_tree_get_time(tree, u, &age_tmp[j]);
        while (u != TSK_NULL){
            tsk_id_t p = tree->parent[u];
            if (p == TSK_NULL) {
                ptree[id_mapping[u]] = -1;
                break;
            }
            if (visited.find(p) != visited.end()){
                // this parent node is already visited
                assert(pcmap.find(p) != pcmap.end());
                pcmap[p].insert(u);
                break;
            }else{
                pcmap.insert(pair<tsk_id_t, set<int>>(p, set<int>()));
                pcmap[p].insert(u);
                visited.insert(p);
                tsk_tree_get_time(tree, p, &age_tmp[index]);
                id_mapping.insert(pair<int, int>(p, index));
                index++;
                u = p;
            }
        }
    }


#ifdef DEBUG
    for(auto it = pcmap.begin(); it != pcmap.end(); it++){
        int p = it->first;
        auto children = it->second;
        printLog(LOG_LOW, "%d(->%d) has %d children\n", p, id_mapping[p], children.size());
        for(int child : children){
            printLog(LOG_LOW, "%d\n", child);
        }
    }

    //printLog(LOG_LOW, "index: %d\n", index);
    //for(int i = 0; i < index; i++){
    //    printLog(LOG_LOW, "age of node %d is %lf\n", i, age_tmp[i]);
    //}
#endif

    // now resolve polytomy
    // we want to resolve from bottom up; since C++ map is ordered, this should be fine
    // we make use of the node labeling rule in tskit: older node has bigger ids
    map<int, int> polynode_map;
    for(auto it = pcmap.begin(); it != pcmap.end(); ++it){
        int p_id = id_mapping[it->first];
        if (it->second.size() > 2){
            int prev_p = tree->parent[it->first] == TSK_NULL? -1:id_mapping[tree->parent[it->first]];
            int counter = 0;
            //printLog(LOG_LOW, "prev_p is %d\n", prev_p);
            // the C++ set is sorted in ascending order
            for(int child : it->second){
                int c_id = id_mapping[child];
                if (polynode_map.find(child) != polynode_map.end()) {c_id = polynode_map[child];}
                if (counter <= 1){ptree[c_id] = p_id;}
                else{
                    if (counter == 2) {
                        ptree[p_id] = index;
                    } else {
                        ptree[index-1] = index;
                    }
                    ptree[c_id] = index;
                    age_tmp[index] = age_tmp[p_id];
                    index++;
                }
                counter++;
            }
            ptree[index-1] = prev_p;
            polynode_map.insert(pair<int, int>(it->first, index-1));
        }else{
            assert(it->second.size() == 2);
            for(int child : it->second){
                if (pcmap.find(child) == pcmap.end() || pcmap[child].size() == 2){
                    ptree[id_mapping[child]] = p_id;
                }
            }
        }
    }

    assert(index == nnodes);

    //finally, discrete time
    for(int i = 0; i < nnodes; i++){
        assert(age_tmp[i] != -1);
        ages[i] = find_time(age_tmp[i], times, ntimes);
    }


}


bool read_local_trees_from_tsinfer(const char *ts_fileName, const double *times, int ntimes, 
                            LocalTrees *trees, vector<string> &seqnames,
                            int start_coord, int end_coord)
{
    // for testing purpose, for now this function will write a file containing the newick rep of
    // each of the local tree with polytomy resolved
    tsk_treeseq_t ts;
    tsk_tree_t tree;
    int ret = tsk_treeseq_load(&ts, ts_fileName, 0);
    check_tsk_error(ret);

    tsk_size_t numSamples = ts.num_samples;
    tsk_size_t numTrees = ts.num_trees;
    trees->clear();
    seqnames.clear();
    trees->start_coord = start_coord;
    trees->end_coord = end_coord;
    int nnodes = 2*numSamples -1;
    int iter;
    ret = tsk_tree_init(&tree, &ts, 0);
    check_tsk_error(ret);

    string s_prev;
    auto prev_map = map<tsk_id_t, int>(); // map from tsk_id_t to 
    Spr spr;
    spr.set_null();
    int *mapping = NULL;
    int id = 0;
    for(iter = tsk_tree_first(&tree); iter == 1; iter = tsk_tree_next(&tree)){
        int start = floor(tree.left);
        int end = floor(tree.right);
        // if (end < start_coord){
        //     printLog(LOG_LOW, "skipping tree %d\n", tsk_tree_get_index(&tree));
        //     continue;
        // }else if (start >= end_coord){
        //     printLog(LOG_LOW, "ignoring local tree from %d\n", tsk_tree_get_index(&tree));
        //     break;
        // }else if (start < start_coord){
        //     start = start_coord;
        // }else if(end > end_coord){
        //     end = end_coord;
        // }

        int ptree[nnodes];
        int ages[nnodes];
        auto curr_map = map<tsk_id_t, int>();
        read_local_tree_from_tsinfer(&tree, ptree, ages, times, ntimes, &curr_map);
        LocalTree *localtree = new LocalTree(ptree, nnodes, ages);
        string s_curr = get_newick_rep_rSPR(localtree);
        printLog(LOG_LOW, "parsing tree %d\n", id++);
        if (s_prev.empty()){
            s_prev = s_curr;
            trees->trees.push_back(LocalTreeSpr(localtree, spr, end - start, mapping));
            continue;
        }
        Node *prev = build_tree(s_prev);
        Node *curr = build_tree(s_curr);
        map<string, int> label_map= map<string, int>();
	    map<int, string> reverse_label_map = map<int, string>();
        prev->labels_to_numbers(&label_map, &reverse_label_map);
		curr->labels_to_numbers(&label_map, &reverse_label_map);
        queue<set<int>*> q1, q2;
        get_moves(prev, curr, &label_map, &reverse_label_map, &q1, &q2);
        assert(q1.size() == q2.size());
        cout << "SPR distance: " << q1.size() << endl;

        mapping = new int[nnodes];
        if (q1.empty()){
            // internal node has age changes
            continue;
        }else{
            map<int, set<int>*> descent_map;
            localtree->descent_leaf_map(&descent_map);
            while(!q1.empty()){
                set<int> *s1 = q1.front();
                set<int> *s2 = q2.front();
                // determine corresponding internal node with leaves s1 and s2 as its descendants
                delete s1;
                delete s2;
                q1.pop();
                q2.pop();
            }
            for(auto it = descent_map.begin(); it != descent_map.end(); ++it){
                cout << it->first << ": ";
                for(int descent : *(it->second)){
                    cout << descent << " ";
                }
                cout << endl;
                delete it->second;
            }
        }
        
        
        s_prev = s_curr;
        prev_map = map<tsk_id_t, int>(curr_map);
    }
    
    ret = tsk_tree_free(&tree);
    check_tsk_error(ret);
    ret = tsk_treeseq_free(&ts);
    check_tsk_error(ret);
    exit(EXIT_FAILURE);
    return true;
}



/////////////////////////////////////////////////////////////////////////////////////////////

void traverse_upwards(tsk_tree_t *tree, int* ptree, int* ages, map<tsk_id_t, int> *mapping,
                        int nnodes, const double *times, int ntimes)
{
    tsk_id_t *samples = tree->samples;
    tsk_size_t num_samples = tsk_treeseq_get_num_samples(tree->tree_sequence);
    tsk_size_t j;
    tsk_id_t u;

    auto visited = set<tsk_id_t>();
    double ages_tmp[nnodes];
    int order = 0; //denote the number of non-sample nodes encountered so far
    for (j = 0; j < num_samples; j++) {
        u = samples[j];
        visited.insert(u);
        mapping->insert(pair<tsk_id_t, int>(u, u));
        tsk_tree_get_time(tree, u, &ages_tmp[u]);
        while (u != TSK_NULL) {
            tsk_id_t p = tree->parent[u];
            if (visited.find(p) == visited.end()){
                 if (tsk_treeseq_is_sample(tree->tree_sequence, u)){
                    ptree[u] = p;
                    //printLog(LOG_HIGH, "%d gets parent %d\n", u, p);
                } else{
                    ptree[order + num_samples] = p;
                    mapping->insert(pair<tsk_id_t, int>(u, order + num_samples));
                    tsk_tree_get_time(tree, u, &ages_tmp[order + num_samples]);
                    //printLog(LOG_HIGH, "%d gets parent %d\n", order+num_samples, p);
                    order++;
                }
                visited.insert(p);
                u = p;
            }else{
                int tmp = u;
                if (!tsk_treeseq_is_sample(tree->tree_sequence, u)){
                    tmp = order + num_samples;
                    order++;
                }
                ptree[tmp] = p;
                mapping->insert(pair<tsk_id_t, int>(u, tmp));
                tsk_tree_get_time(tree, u, &ages_tmp[tmp]);
                //printLog(LOG_HIGH, "%d gets parent %d\n", u, p);
                break;
            }
           
        }
    }

    // convert tskit_id_t to the numbering for argweaver
    // also round the age of each node
    for (int i = 0; i < nnodes; i++){
        if (ptree[i] != -1){
            ptree[i] = mapping->at(ptree[i]);
        }
        ages[i] = find_time(ages_tmp[i], times, ntimes);
    }


    // for debugging only
#ifdef DEBUG
    // for(int i = 0; i < 2*num_samples-1; i++){
    //     printLog(LOG_LOW, "%d has parent %d\n", i, ptree[i]);\
    //     printLog(LOG_HIGH, "%d is %lf generations old in DSMC\n", i, times[ages[i]]);
    // }
    // for(auto iter = mapping->begin(); iter != mapping->end(); ++iter){
    //     printLog(LOG_LOW, "%d -> %d\n", iter->first, iter->second);
    // }
#endif

}

bool identify_1SPR(Spr *spr, int *mapping, const map<tsk_id_t, int> *prev, const map<tsk_id_t, int> *curr,
                    tsk_tree_t* prev_t, tsk_tree_t* curr_t, 
                    const double *times, int ntimes)
{
    //TODO: identify spr and fill in the mapping array
    //fill in mapping array
    int count = 0;
    tsk_id_t out = -1;
    for(auto it = prev->begin(); it != prev->end(); ++it){
        if (curr->count(it->first) == 0) {
            mapping[it->second] = -1;
            out = it->first;
            count++;
        }else{
            mapping[prev->at(it->first)] = curr->at(it->first);
        }
    }

    if (count == 0){
        // maybe we can also return true in this case? Just let spr remain its null state
        printLog(LOG_LOW, "Two consecutive trees are equivalent");
        exit(EXIT_FAILURE);
    }else if(count > 1){
        return false; // a single SPR is not sufficient
    }else{
        // find the "in" node_id
        tsk_id_t in = -1;
        for(auto it = curr->begin(); it != curr->end(); ++it){
            if (prev->count(it->first) == 0){
                in = it->first;
                break;
            }
        }
        //compare the topology around out and in node:
    #ifdef DEBUG
        printLog(LOG_LOW, "node out: %d, node in: %d\n", out, in);
    #endif

        tsk_id_t p_prev = prev_t->parent[out];
        tsk_id_t c1_prev = prev_t->left_child[out];
        tsk_id_t c2_prev = prev_t->right_child[out];
        tsk_id_t p_curr = curr_t->parent[in];
        tsk_id_t c1_curr = curr_t->left_child[in];
        tsk_id_t c2_curr = curr_t->right_child[in];

        double recomb_time;
        double coal_time;
        double age_out;
        tsk_id_t recomb_node;
        tsk_id_t coal_node;
        tsk_tree_get_time(prev_t, out, &age_out);
        tsk_tree_get_time(curr_t, in, &coal_time);
        
        if (p_prev == p_curr && ( (c1_prev == c1_curr && c2_prev == c2_curr) || (c1_prev == c2_curr && c2_prev == c1_curr) )){
            // this is the simplest case: the tree topology remains unchanged, only coalesce time changes
            //choose c2_prev is also equivalent
            recomb_node = c1_prev;
            double out_time;
            tsk_tree_get_time(prev_t, out, &out_time);
            coal_node = coal_time >= out_time? out : c2_prev;
            printLog(LOG_LOW, "case1\n");
        }else if (curr_t->parent[c1_prev] == in){
            recomb_node = c1_prev;
            coal_node = curr_t->left_child[in] != c1_prev? curr_t->left_child[in] : curr_t->right_child[in];
            printLog(LOG_LOW, "case2\n");
        }else if(curr_t->parent[c2_prev] == in){
            recomb_node = c2_prev;
            coal_node = curr_t->left_child[in] != c2_prev? curr_t->left_child[in] : curr_t->right_child[in];
            printLog(LOG_LOW, "case3\n");
        }else{
            printLog(LOG_LOW, "case4: single SPR is not enough\n");
            return false;
        }
        
        spr->recomb_node = prev->at(recomb_node);
        spr->coal_node = prev->at(coal_node);
        double lower_bound_recomb_time;
        tsk_tree_get_time(prev_t, recomb_node, &lower_bound_recomb_time);

        //sanity check
        if(age_out < lower_bound_recomb_time || coal_time < lower_bound_recomb_time){
            printLog(LOG_LOW, "coalescence time younger than lower bound, i.e, age of the recombination node");
            exit(EXIT_FAILURE);
        }

        recomb_time = rand() % (min((int)age_out, (int)coal_time) - (int)lower_bound_recomb_time) + (int)lower_bound_recomb_time;
        spr->recomb_time = find_time(recomb_time, times, ntimes);
        spr->coal_time = find_time(coal_time, times, ntimes);
    
    #ifdef DEBUG
        printLog(LOG_LOW, "recomb_node: %d\n", recomb_node);
        printLog(LOG_LOW, "recomb_time: %lf\n", recomb_time);
        printLog(LOG_LOW, "coal_node: %d\n", coal_node);
        printLog(LOG_LOW, "coal_time: %lf\n", coal_time);
    #endif

    }
    return true;
}

 bool read_local_trees_from_ts(const char *ts_fileName, const double *times, int ntimes, 
                                     LocalTrees *trees, vector<string> &seqnames,
                                     int start_coord, int end_coord)
{
    tsk_treeseq_t ts;
    tsk_tree_t tree;
    int ret = tsk_treeseq_load(&ts, ts_fileName, 0);
    check_tsk_error(ret);

    tsk_size_t numSamples = ts.num_samples;
    tsk_size_t numTrees = ts.num_trees;
    trees->clear();
    seqnames.clear();
    trees->start_coord = start_coord;
    trees->end_coord = end_coord;

    int iter;
    ret = tsk_tree_init(&tree, &ts, 0);
    check_tsk_error(ret);
    auto prev_map = map<tsk_id_t, int>();
    tsk_tree_t prev_tree;
    int nnodes = 2*numSamples-1;
    Spr spr;
    spr.set_null();
    for (iter = tsk_tree_first(&tree); iter == 1; iter = tsk_tree_next(&tree)) {
        int start = floor(tree.left);
        int end = floor(tree.right);
        if (end < start_coord){
            printLog(LOG_LOW, "skipping tree %d\n", tsk_tree_get_index(&tree));
            continue;
        }else if (start >= end_coord){
            printLog(LOG_LOW, "ignoring local tree from %d\n", tsk_tree_get_index(&tree));
            break;
        }else if (start < start_coord){
            start = start_coord;
        }else if(end > end_coord){
            end = end_coord;
        }

        printLog(LOG_LOW, "\ntree %d: %d - %d\n", tsk_tree_get_index(&tree), start, end);
        int ptree[nnodes];
        int ages[nnodes];
        auto curr_map = map<tsk_id_t, int>();
        traverse_upwards(&tree, ptree, ages, &curr_map, nnodes, times, ntimes);

        // make a local tree
        LocalTree *localtree = new LocalTree(ptree, nnodes, ages);
        int *mapping = NULL;
        if (prev_map.size() != 0) {
            //printLog(LOG_LOW, "try to identify SPR");
            mapping = new int[nnodes];
            // need to figure out how to do the mapping and identify SPR
            if(!identify_1SPR(&spr, mapping, &prev_map, &curr_map, &prev_tree, &tree, times, ntimes)){
               printLog(LOG_LOW, "consecutive trees are not reachable by one SPR");
               exit (EXIT_FAILURE);
            }
            tsk_tree_free(&prev_tree); //don't forget to free the copied tree
        }
        prev_map = map<tsk_id_t, int>(curr_map); // copy constructor
        tsk_tree_copy(&tree, &prev_tree, 0);
        trees->trees.push_back(LocalTreeSpr(localtree, spr, end - start, mapping));
    }

    check_tsk_error(iter);
    tsk_tree_free(&prev_tree); //don't forget to free the copied tree
    tsk_tree_free(&tree);
    tsk_treeseq_free(&ts);

    if (trees->get_num_trees() > 0) {
        trees->nnodes = trees->front().tree->nnodes;
        trees->set_default_seqids();
    }

    //don't forget to initilize seqnames!!!
    
    //for(int j = 0; j < numSamples; j++){
    //    seqnames.push_back(string());
    //}

    printLog(LOG_LOW, "number of samples in the tree sequences: %d\n", numSamples);
    printLog(LOG_LOW, "number of local trees in the tree sequences: %d\n", numTrees);
    printLog(LOG_LOW, "total number of local trees read: %d\n", trees->trees.size());
    printLog(LOG_LOW, "tree sequence starat at %d and ends at %d\n", trees->start_coord, trees->end_coord);


    return true;
 }

//=============================================================================
// debugging output

// dump raw information about a local tree
void print_local_tree(const LocalTree *tree, FILE *out)
{
    const LocalNode *nodes = tree->nodes;

    for (int i=0; i<tree->nnodes; i++) {
        fprintf(out, "%d: parent=%2d, child=(%2d, %2d), age=%d, path=%d\n",
                i, nodes[i].parent, nodes[i].child[0], nodes[i].child[1],
                nodes[i].age, nodes[i].pop_path);
    }
}




// dump raw information about a local tree
void draw_local_tree(const LocalTree *tree, FILE *out, int depth, int inode)
{
    const LocalNode &node = tree->nodes[inode];
    depth = tree->nodes[tree->root].age - tree->nodes[inode].age;
    for (int i=0; i<depth; i++)
        fprintf(out, " ");
    fprintf(out, "%d: age=%i\t(%i)%s", inode, node.age, node.pop_path,
            node.is_leaf() ? " (leaf)\n" : "\n");

    // recurse
    if (!node.is_leaf()) {
        draw_local_tree(tree, out, depth+2, node.child[0]);
        draw_local_tree(tree, out, depth+2, node.child[1]);
    }
}

void draw_local_tree(const LocalTree *tree, FILE *out, int depth)
{
    if (tree->root != -1)
        draw_local_tree(tree, out, depth, tree->root);
}


// dump raw information about a set of local trees
void print_local_trees(const LocalTrees *trees, FILE *out)
{
    int end = trees->start_coord;
    for (LocalTrees::const_iterator it=trees->begin();
         it != trees->end(); ++it)
    {
        int start = end;
        end += it->blocklen;
        LocalTree *tree = it->tree;

        fprintf(out, "%d-%d\n", start, end);
        print_local_tree(tree, out);

        LocalTrees::const_iterator it2 = it;
        ++it2;
        if (it2 != trees->end()) {
            const Spr &spr = it2->spr;
            fprintf(out, "spr: r=(%d, %d), c=(%d, %d) path=%i\n\n",
                    spr.recomb_node, spr.recomb_time,
                    spr.coal_node, spr.coal_time, spr.pop_path);
        }
    }
}


//=============================================================================
// assert functions

// Asserts that a postorder traversal is correct
bool assert_tree_postorder(const LocalTree *tree, const int *order)
{
    if (tree->root != order[tree->nnodes-1])
        return false;

    bool seen[tree->nnodes];
    fill(seen, seen + tree->nnodes, false);

    for (int i=0; i<tree->nnodes; i++) {
        int node = order[i];
        seen[node] = true;
        if (!tree->nodes[node].is_leaf()) {
            if (! seen[tree->nodes[node].child[0]] ||
                ! seen[tree->nodes[node].child[1]])
                return false;
        }
    }

    return true;
}


// Asserts structure of tree
bool assert_tree(const LocalTree *tree,
                 const PopulationTree *pop_tree)
{
    LocalNode *nodes = tree->nodes;
    int nnodes = tree->nnodes;

    for (int i=0; i<nnodes; i++) {
        int *c = nodes[i].child;

        // assert parent child links
        if (c[0] != -1) {
            if (c[0] < 0 || c[0] >= nnodes)
                return false;
            if (nodes[c[0]].parent != i)
                return false;
        }
        if (c[1] != -1) {
            if (c[1] < 0 || c[1] >= nnodes)
                return false;
            if (nodes[c[1]].parent != i)
                return false;
        }

        // check root
        if (nodes[i].parent == -1) {
            if (tree->root != i)
                return false;
        } else {
            if (nodes[i].parent < 0 || nodes[i].parent >= nnodes)
                return false;
        }

        if (pop_tree != NULL && nodes[i].parent != -1) {
            assert(pop_tree->get_pop(nodes[i].pop_path,
                                     nodes[nodes[i].parent].age) ==
                   pop_tree->get_pop(nodes[nodes[i].parent].pop_path,
                                     nodes[nodes[i].parent].age));
        }
    }

    // check root
    if (nodes[tree->root].parent != -1)
        return false;

    return true;
}


bool assert_spr(const LocalTree *last_tree, const LocalTree *tree,
                const Spr *spr, const int *mapping,
                const PopulationTree *pop_tree, bool pruned_internal)
{
    LocalNode *last_nodes = last_tree->nodes;
    LocalNode *nodes = tree->nodes;
    static int count=0;
    count++;

    if (spr->is_null()) {
        // just check that mapping is 1-to-1
        bool mapped[tree->nnodes];
        fill(mapped, mapped + tree->nnodes, false);

        for (int i=0; i<tree->nnodes; i++) {
            int i2 = mapping[i];
            assert(i2 != -1);
            assert(!mapped[i2]);
            mapped[i2] = true;
            // check for parental node
            assert((last_tree->nodes[i].parent == -1 &&
                    tree->nodes[i2].parent == -1) ||
                   (mapping[last_tree->nodes[i].parent] ==
                    tree->nodes[i2].parent));
            // check  for child node
            if (last_tree->nodes[i].child[0] == -1) {
                assert(last_tree->nodes[i].child[1] == -1);
                assert(tree->nodes[i2].child[0] == -1);
                assert(tree->nodes[i2].child[1] == -1);
            } else {
                assert((mapping[last_tree->nodes[i].child[0]] ==
                        tree->nodes[i2].child[0] &&
                        mapping[last_tree->nodes[i].child[1]] ==
                        tree->nodes[i2].child[1]) ||
                       (mapping[last_tree->nodes[i].child[0]] ==
                        tree->nodes[i2].child[1] &&
                        mapping[last_tree->nodes[i].child[1]] ==
                        tree->nodes[i2].child[0]));
            }
            // check for age
            assert(last_tree->nodes[i].age ==
                   tree->nodes[i2].age);
            assert(i == last_tree->nodes[last_tree->root].child[0] ||
                   pop_tree->paths_equal(last_tree->nodes[i].pop_path,
                                         tree->nodes[i2].pop_path,
                                         tree->nodes[i2].age,
                                         i2 == tree->root ? -1 :
                                         tree->nodes[tree->nodes[i2].parent].age));
        }
        return true;
    }

    if (pop_tree != NULL) {
        assert(pop_tree->get_pop(last_tree->nodes[spr->recomb_node].pop_path,
                                 spr->recomb_time) ==
               pop_tree->get_pop(spr->pop_path, spr->recomb_time));
        assert(pop_tree->get_pop(last_tree->nodes[spr->coal_node].pop_path,
                                 spr->coal_time) ==
               pop_tree->get_pop(spr->pop_path, spr->coal_time));
        assert(pop_tree->path_prob(spr->pop_path, spr->recomb_time, spr->coal_time) > 0);
    }

    if (spr->recomb_node == -1)
        assert(false);

    // coal time is older than recomb time
    // that is, we need recomb_time \leq coal_time
    if (spr->recomb_time > spr->coal_time)
        assert(false);

    // recomb cannot be on root branch
    if (pop_tree == NULL)
        assert(last_nodes[spr->recomb_node].parent != -1);

    // ensure recomb is within branch
    if ((last_nodes[spr->recomb_node].parent != -1 &&
         spr->recomb_time > last_nodes[last_nodes[spr->recomb_node].parent].age)
        || spr->recomb_time < last_nodes[spr->recomb_node].age)
        assert(false);

    // ensure coal is within branch
    if (spr->coal_time < last_nodes[spr->coal_node].age)
        assert(false);
    if (last_nodes[spr->coal_node].parent != -1) {
        if (spr->coal_time > last_nodes[last_nodes[spr->coal_node].parent].age)
            assert(false);
    }

    // recomb baring branch cannot be broken
    assert(mapping[spr->recomb_node] != -1);

    if (spr->recomb_node == spr->coal_node) {
        assert(pop_tree != NULL);
        assert(!pop_tree->paths_equal(last_tree->nodes[spr->recomb_node].pop_path,
                                      spr->pop_path,
                                      spr->recomb_time,
                                      spr->coal_time));
        assert(spr->recomb_time != spr->coal_time);
        assert(pop_tree->paths_equal(tree->nodes[mapping[spr->recomb_node]].pop_path,
                                     spr->pop_path,
                                     spr->recomb_time,
                                     spr->coal_time));
        assert(tree->nodes[mapping[spr->recomb_node]].age ==
               last_tree->nodes[spr->recomb_node].age);
        for (int i=0; i < last_tree->nnodes; i++) {
            assert(mapping[i] >= 0 && mapping[i] < last_tree->nnodes);
            int last_parent = last_tree->nodes[i].parent;
            int parent = tree->nodes[mapping[i]].parent;
            assert(last_tree->nodes[i].age == tree->nodes[mapping[i]].age);
            int parent_age;
            if (last_parent == -1) {
                assert(parent == -1);
                parent_age = -1;
            } else {
                assert(mapping[last_parent] == parent);
                parent_age = tree->nodes[parent].age;
                assert(parent_age == last_tree->nodes[last_parent].age);
            }
            if (i == spr->recomb_node) {
                assert(pop_tree->paths_equal(last_tree->nodes[i].pop_path,
                                             tree->nodes[mapping[i]].pop_path,
                                             last_tree->nodes[i].age,
                                             spr->recomb_time));
                assert(pop_tree->paths_equal(last_tree->nodes[i].pop_path,
                                             tree->nodes[mapping[i]].pop_path,
                                             spr->coal_time, parent_age));
            } else {
                assert((i == last_tree->nodes[last_tree->root].child[0] && pruned_internal) ||
                       pop_tree->paths_equal(last_tree->nodes[i].pop_path,
                                             tree->nodes[mapping[i]].pop_path,
                                             tree->nodes[mapping[i]].age,
                                             parent == -1 ? -1 :
                                             min(last_tree->nodes[last_parent].age,
                                                 tree->nodes[parent].age)));
            }
            if (last_tree->nodes[i].child[0] == -1) {
                assert(last_tree->nodes[i].child[1] == -1);
                assert(tree->nodes[mapping[i]].child[0] == -1);
                assert(tree->nodes[mapping[i]].child[1] == -1);
            } else {
                assert((mapping[last_tree->nodes[i].child[0]] == tree->nodes[mapping[i]].child[0] &&
                        mapping[last_tree->nodes[i].child[1]] == tree->nodes[mapping[i]].child[1]) ||
                       (mapping[last_tree->nodes[i].child[1]] == tree->nodes[mapping[i]].child[0] &&
                        mapping[last_tree->nodes[i].child[0]] == tree->nodes[mapping[i]].child[1]));
            }
        }
        return true;
    }

    // ensure spr matches the trees
    int recoal = tree->nodes[mapping[spr->recomb_node]].parent;
    assert(recoal != -1);
    int *c = tree->nodes[recoal].child;
    int other = (c[0] == mapping[spr->recomb_node] ? c[1] : c[0]);
    if (mapping[spr->coal_node] != -1) {
        // coal node is not broken, so it should map correctly
        assert(other == mapping[spr->coal_node]);
    } else {
        // coal node is broken
        int broken = last_tree->nodes[spr->recomb_node].parent;
        int *c = last_tree->nodes[broken].child;
        int last_other = (c[0] == spr->recomb_node ? c[1] : c[0]);
        assert(mapping[last_other] != -1);
        assert(tree->nodes[mapping[last_other]].parent == recoal);
    }

    // ensure mapped nodes don't change in age
    for (int i=0; i<last_tree->nnodes; i++) {
        int i2 = mapping[i];
        if (i2 != -1) assert(last_nodes[i].age == nodes[i2].age);
        if (pop_tree != NULL) {
            int subtree_root = nodes[tree->root].child[0];
            int last_subtree_root = last_nodes[last_tree->root].child[0];
            if (i2 == -1 && (i != last_subtree_root || !pruned_internal)) {
                int recomb_parent = last_nodes[spr->recomb_node].parent;
                assert(i = recomb_parent);
                if (last_nodes[spr->coal_node].parent == recomb_parent) {
                    assert(mapping[spr->recomb_node] != -1);
                    int mapped_node = nodes[mapping[spr->recomb_node]].parent;
                    int path1 = pop_tree->consistent_path(last_nodes[spr->coal_node].pop_path,
                                                          last_nodes[i].pop_path,
                                                          spr->coal_time,
                                                          last_nodes[i].age,
                                                          last_nodes[i].parent == -1 ? -1 :
                                                          last_nodes[last_nodes[i].parent].age);
                    int path2 = nodes[mapped_node].pop_path;
                    if (mapped_node != subtree_root || !pruned_internal) {
                        assert(pop_tree->paths_equal(path1, path2, nodes[mapped_node].age,
                                                     nodes[mapped_node].parent == -1 ? -1 :
                                                     nodes[nodes[mapped_node].parent].age));
                    }
                } else if (recomb_parent == spr->coal_node) {
                    int mapped_node = nodes[mapping[spr->recomb_node]].parent;
                    if (mapped_node != subtree_root || !pruned_internal) {
                        assert(pop_tree->paths_equal(last_nodes[spr->coal_node].pop_path,
                                                     nodes[mapped_node].pop_path,
                                                     spr->coal_time,
                                                     last_nodes[spr->coal_node].parent == -1
                                                     ? -1 : last_nodes[last_nodes[spr->coal_node].parent].age));
                    }
                } else {
                    int mapped_node = nodes[mapping[spr->coal_node]].parent;
                    if (mapped_node != subtree_root || !pruned_internal) {
                        assert(pop_tree->paths_equal(last_nodes[spr->coal_node].pop_path,
                                                     nodes[mapped_node].pop_path,
                                                     spr->coal_time, last_nodes[spr->coal_node].parent == -1 ? -1 :
                                                     last_nodes[last_nodes[spr->coal_node].parent].age));
                    }
                }
            } else if (i == spr->recomb_node) {
                int target_path = pop_tree->consistent_path(last_nodes[i].pop_path,
                                                            spr->pop_path,
                                                            last_nodes[i].age,
                                                            spr->recomb_time,
                                                            spr->coal_time);
                assert(pop_tree->paths_equal(nodes[i2].pop_path,
                                             target_path,
                                             nodes[i2].age,
                                             i2 == tree->root ? -1 :
                                             nodes[nodes[i2].parent].age));
            } else if ((i != last_tree->nodes[last_tree->root].child[0] || !pruned_internal) &&
                           i2 != tree->nodes[tree->root].child[0]) {
                int last_end = (i == last_tree->root ? -1 : last_nodes[last_nodes[i].parent].age);
                int end = (i2 == tree->root ? -1 : nodes[nodes[i2].parent].age);
                int end_time = (last_end == -1 && end == -1 ? -1 :
                                (last_end == -1 ? end :
                                 (end == -1 ? last_end :
                                  min(end, last_end))));
                assert(pop_tree->paths_equal(last_nodes[i].pop_path,
                                             nodes[i2].pop_path,
                                             last_nodes[i].age,
                                             end_time));
            }
        }
        if (last_nodes[i].is_leaf())
            assert(nodes[i2].is_leaf());
    }

    // test for bubbles
    assert(spr->recomb_node != spr->coal_node);

    return true;
}


// add a thread to an ARG
bool assert_trees(const LocalTrees *trees, const PopulationTree *pop_tree,
                  bool pruned_internal)
{
    LocalTree *last_tree = NULL;
    int seqlen = 0;

    // assert first tree has null mapping and spr
    if (trees->begin() != trees->end()) {
        assert(trees->begin()->spr.is_null());
        assert(!trees->begin()->mapping);
    }

    // loop through blocks
    for (LocalTrees::const_iterator it=trees->begin();
         it != trees->end(); ++it)
    {
        LocalTree *tree = it->tree;
        const Spr *spr = &it->spr;
        const int *mapping = it->mapping;
        seqlen += it->blocklen;

        assert(it->blocklen >= 0);
        assert(assert_tree(tree, pop_tree));

        if (last_tree)
            assert(assert_spr(last_tree, tree, spr, mapping, pop_tree, pruned_internal));
        last_tree = tree;
    }

    assert(seqlen == trees->length());

    return true;
}


//=============================================================================
// C inferface
extern "C" {


LocalTrees *arghmm_new_trees(
    int **ptrees, int **ages, int **sprs, int *blocklens,
    int ntrees, int nnodes, int start_coord)
{
    // setup model, local trees, sequences
    return new LocalTrees(ptrees, ages, sprs, blocklens, ntrees, nnodes,
                          -1, start_coord);
}


LocalTrees *arghmm_copy_trees(LocalTrees *trees)
{
    LocalTrees *trees2 = new LocalTrees();
    trees2->copy(*trees);
    return trees2;
}


int get_local_trees_ntrees(LocalTrees *trees)
{
    return trees->trees.size();
}


int get_local_trees_nnodes(LocalTrees *trees)
{
    return trees->nnodes;
}


void get_local_trees_ptrees(LocalTrees *trees, int **ptrees, int **ages,
                            int **sprs, int *blocklens)
{
    // setup permutation
    const int nleaves = trees->get_num_leaves();
    int perm[trees->nnodes];
    for (int i=0; i<nleaves; i++)
        perm[i] = trees->seqids[i];
    for (int i=nleaves; i<trees->nnodes; i++)
        perm[i] = i;

    // debug
    assert_trees(trees);

    // convert trees
    int i = 0;
    for (LocalTrees::iterator it=trees->begin(); it!=trees->end(); ++it, i++) {
        LocalTree *tree = it->tree;

        for (int j=0; j<tree->nnodes; j++) {
            int parent = tree->nodes[j].parent;
            if (parent != -1)
                parent = perm[parent];
            ptrees[i][perm[j]] = parent;
            ages[i][perm[j]] = tree->nodes[j].age;
        }
        blocklens[i] = it->blocklen;

        if (!it->spr.is_null()) {
            sprs[i][0] = perm[it->spr.recomb_node];
            sprs[i][1] = it->spr.recomb_time;
            sprs[i][2] = perm[it->spr.coal_node];
            sprs[i][3] = it->spr.coal_time;

            assert(it->spr.recomb_time >= ages[i-1][sprs[i][0]]);
            assert(it->spr.coal_time >= ages[i-1][sprs[i][2]]);

        } else {
            sprs[i][0] = it->spr.recomb_node;
            sprs[i][1] = it->spr.recomb_time;
            sprs[i][2] = it->spr.coal_node;
            sprs[i][3] = it->spr.coal_time;
        }

    }
}


void delete_local_trees(LocalTrees *trees)
{
    delete trees;
}


void write_local_trees(char *filename, LocalTrees *trees, char **names,
                       double *times, int ntimes, bool pop_model)
{
    write_local_trees(filename, trees, names, times, pop_model);
}


LocalTrees *read_local_trees(const char *filename, const double *times,
                             int ntimes)
{
    char ***names = NULL;
    LocalTrees *trees = new LocalTrees();
    vector<string> seqnames;

    CompressStream stream(filename, "r");
    if (stream.stream &&
        read_local_trees(stream.stream, times, ntimes, trees, seqnames))
    {
        if (names) {
            // copy names
            *names = new char* [seqnames.size()];
            for (unsigned int i=0; i<seqnames.size(); i++) {
                (*names)[i] = new char [seqnames[i].size()+1];
                strncpy((*names)[i], seqnames[i].c_str(), seqnames[i].size()+1);
            }
        }
    } else {
        delete trees;
        trees = NULL;
    }
    return trees;
}


void get_treelens(const LocalTrees *trees, const double *times, int ntimes,
                  double *treelens)
{
    const bool use_basal = false;
    int i = 0;
    for (LocalTrees::const_iterator it=trees->begin();
         it!=trees->end(); ++it, i++)
        treelens[i] = get_treelen(it->tree, times, ntimes, use_basal);
}


void get_local_trees_blocks(const LocalTrees *trees, int *starts, int *ends)
{
    int i = 0;
    int end = trees->start_coord;
    for (LocalTrees::const_iterator it=trees->begin();
         it!=trees->end(); ++it, i++)
    {
        int start = end;
        end += it->blocklen;
        starts[i] = start;
        ends[i] = end;
    }
}


} // extern C

} // namespace argweaver

