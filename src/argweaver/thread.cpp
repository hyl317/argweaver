
#include "thread.h"
#include "trans.h"
#include "model.h"

namespace argweaver {


// rename a node from src_node to dest_node while mantaining the
// structure of the tree
void rename_node(LocalTree *tree, int src_node, int dest_node)
{
    // special case
    if (src_node == dest_node)
        return;

    LocalNode *nodes = tree->nodes;

    // copy displaced node
    nodes[dest_node].copy(nodes[src_node]);

    // notify parent of displacement
    int parent = nodes[dest_node].parent;
    if (parent != -1) {
        int *c = nodes[parent].child;
        if (c[0] == src_node)
            c[0] = dest_node;
        else
            c[1] = dest_node;
    }

    // notify children of displacement
    int *c = nodes[dest_node].child;
    if (c[0] != -1)
        nodes[c[0]].parent = dest_node;
    if (c[1] != -1)
        nodes[c[1]].parent = dest_node;
}


// Add a single new leaf to a tree.
// Leaf connects to the point indicated by (node, time)
//
//   - New leaf will be named 'nleaves' (newleaf)
//   - Node currently named 'nleaves' will be displaced to 'nnodes' (displaced)
//   - New internal node will be named 'nnodes+1' (newcoal)
//v
void add_tree_branch(LocalTree *tree, int node, int time, int pop_path)
{
    // Ensure tree has capacity to add nodes.
    // This is the responsibility of the caller.
    assert(tree->capacity >= tree->nnodes + 2);

    // get tree info
    LocalNode *nodes = tree->nodes;
    int nleaves = tree->get_num_leaves();
    int nnodes = tree->nnodes;
    int nnodes2 = nnodes + 2;

    // get major node ids
    int newleaf = nleaves;
    int displaced = nnodes;
    int newcoal = nnodes + 1;

    // determine node displacement
    int node2 = (node != newleaf ? node : displaced);
    int parent = nodes[node].parent;
    int parent2 = (parent != newleaf ? parent : displaced);

    // displace node
    // Node currently named 'nleaves' will be displaced to 'nnodes' (displaced)
    if (newleaf < displaced)
        rename_node(tree, newleaf, displaced);

    // add new leaf
    nodes[newleaf].parent = newcoal;
    nodes[newleaf].child[0] = -1;
    nodes[newleaf].child[1] = -1;
    nodes[newleaf].age = 0;
    nodes[newleaf].pop_path = pop_path;

    // add new coal node
    nodes[newcoal].parent = parent2;
    nodes[newcoal].child[0] = newleaf;
    nodes[newcoal].child[1] = node2;
    nodes[newcoal].age = time;
    nodes[newcoal].pop_path = nodes[nodes[newcoal].child[1]].pop_path;

    // fix pointers
    nodes[node2].parent = newcoal;
    if (parent2 != -1) {
        int *c = nodes[parent2].child;
        if (c[0] == node2)
            c[0] = newcoal;
        else
            c[1] = newcoal;
    }

    // fix up tree data
    tree->nnodes = nnodes2;
    if (nodes[newcoal].parent == -1)
        tree->root = newcoal;
    else
        tree->root = (tree->root != newleaf ? tree->root : displaced);
}


// removes a leaf branch from a local tree
// any node displacements are recorded in the displace array
// If displace is NULL, displacements are not recorded
void remove_tree_branch(LocalTree *tree, int remove_leaf,
                        const ArgModel *model, int *displace)
{
    LocalNode *nodes = tree->nodes;
    int nnodes = tree->nnodes;
    int last_leaf = tree->get_num_leaves() - 1;

    // remove coal node
    int remove_coal = nodes[remove_leaf].parent;
    int *c = nodes[remove_coal].child;
    int coal_child = (c[0] == remove_leaf ? c[1] : c[0]);
    int coal_parent = nodes[remove_coal].parent;
    nodes[coal_child].parent = coal_parent;
    if (model->pop_tree != NULL) {
        nodes[coal_child].pop_path = model->pop_tree->consistent_path(
               nodes[coal_child].pop_path, nodes[remove_coal].pop_path,
               nodes[coal_child].age,
               nodes[remove_coal].age,
               coal_parent == -1 ? model->ntimes - 1 : nodes[coal_parent].age);

    }
    if (coal_parent != -1) {
        c = nodes[coal_parent].child;
        if (c[0] == remove_coal)
            c[0] = coal_child;
        else
            c[1] = coal_child;
    }

    // record displace nodes
    if (displace) {
        for (int i=0; i<nnodes; i++)
            displace[i] = i;
        displace[remove_leaf] = -1;
        displace[remove_coal] = -1;
    }

    // move last leaf into remove_leaf spot
    if (last_leaf != remove_leaf) {
        if (displace)
            displace[last_leaf] = remove_leaf;
        rename_node(tree, last_leaf, remove_leaf);
    }

    // move nodes in nnodes-2 and nnodes-1 into holes
    int hole = last_leaf;
    if (remove_coal != nnodes-2) {
        if (displace)
            displace[nnodes-2] = hole;
        rename_node(tree, nnodes-2, hole);
        hole = remove_coal;
    }
    if (remove_coal != nnodes-1) {
        if (displace)
            displace[nnodes-1] = hole;
        rename_node(tree, nnodes-1, hole);
    }

    // set tree data
    tree->nnodes -= 2;
    int root = tree->root;
    if (tree->root == remove_coal)
        root = coal_child;
    if (root == nnodes-2)
        root = last_leaf;
    if (root == nnodes-1)
        root = hole;
    tree->root = root;
}



// update an SPR and mapping after adding a new branch to two
// neighboring local trees
//
// newleaf -- new leaf node added to the tree
// displaced  -- node that used to be named 'newleaf'
// newcoal -- new node that has 'newleaf' as a child.
void add_spr_branch(const LocalTree *tree, const LocalTree *last_tree,
                    State state, State last_state,
                    Spr *spr, int *mapping,
                    int newleaf, int displaced, int newcoal)
{
    // get tree info
    const LocalNode *nodes = tree->nodes;
    const LocalNode *last_nodes = last_tree->nodes;
    bool fix_mapping = true;
#ifdef DEBUG
    const LocalTree orig_tree(*tree);
    const LocalTree orig_last_tree(*last_tree);
    const State orig_state(state);
    const State orig_last_state(last_state);
    const Spr orig_spr(*spr);
    static int count=0;
    count++;
#endif

    // determine node displacement
    int last_node2 = (last_state.node != newleaf ? last_state.node : displaced);
    int node2 = (state.node != newleaf ? state.node : displaced);

    // update mapping due to displacement
    mapping[displaced] = mapping[newleaf];
    mapping[newleaf] = newleaf;

    // set default new node mapping
    mapping[newcoal] = newcoal;

    for (int i=newleaf+1; i<tree->nnodes; i++) {
        if (mapping[i] == newleaf)
            mapping[i] = displaced;
    }


    // update spr due to displacement
    if (spr->recomb_node == newleaf)
        spr->recomb_node = displaced;
    if (spr->coal_node == newleaf)
        spr->coal_node = displaced;

    // parent of recomb node should be the recoal point
    // however, if it equals newcoal, then either (1) the recomb branch is
    // renamed, (2) there is mediation, or (3) new branch escapes
    int recoal = nodes[mapping[spr->recomb_node]].parent;
    if (spr->recomb_node == spr->coal_node) {
        if (last_state.node == spr->recomb_node) {
            if (last_state.time < spr->recomb_time ||
                (last_state.time == spr->recomb_time &&
                 state.time == spr->recomb_time)) {
                spr->recomb_node = last_nodes[spr->recomb_node].parent;
                spr->coal_node = spr->recomb_node;
                fix_mapping = false;
            } else if (last_state.time >= spr->recomb_time &&
                       last_state.time < spr->coal_time) {
                if (state.time <= last_state.time) {
                    assert(state.time >= spr->recomb_time);
                    spr->coal_node = last_tree->get_sibling(spr->recomb_node);
                    spr->coal_time = state.time;
                    mapping[newcoal] = -1;
                    fix_mapping = false;
                } else {
                    assert(state.time == spr->coal_time);
                    spr->coal_node = newcoal;
                    mapping[newcoal] = -1;
                    fix_mapping = false;
                }
            } else if (last_state.time == spr->coal_time) {
                spr->coal_time = state.time;
                spr->coal_node = last_tree->get_sibling(spr->recomb_node);
                mapping[newcoal] = -1;
                fix_mapping = false;
            } // if last_state.time > spr->coal_time do nothing

        }
    } else if (recoal == newcoal) {
        if (spr->recomb_node == spr->coal_node) {
            // (4) new branch breaks "bubble" recomb.
            if (last_state.time > spr->recomb_time &&
                last_state.time < spr->coal_time) {
                assert(state.time == spr->coal_time);
                spr->coal_node = recoal;
            } else if (last_state.time <= spr->recomb_time) {
                spr->recomb_node = recoal;
                spr->coal_node = recoal;
            }
        } else if (mapping[last_node2] == node2) {
            // (1) recomb is above coal state, we rename spr recomb node
            spr->recomb_node = newcoal;
        } else {
            // if this is a mediated coal, then state should equal recomb
            int state_node = (state.node != newleaf) ?
                state.node : displaced;
            if (state_node == mapping[spr->recomb_node]) {
                // (3) this is a mediated coal, rename coal node and time
                spr->coal_node = newleaf;
                // time should only become more recent, meaning pop_path
                // need not change
                assert(spr->coal_time >= state.time);
                spr->coal_time = state.time;
            } else {
                // (2) this is the new branch escaping
                // no other updates are necessary
            }
        }
    } else {
        // the other possibility is that newcoal is under recoal point
        // if newcoal is child of recoal, then coal is renamed
        const int *c = nodes[recoal].child;
        if (c[0] == newcoal || c[1] == newcoal) {
            // we either coal above the newcoal or our existing
            // node just broke and newcoal was underneath.

            // if newcoal was previously above spr->coal_node
            // then we rename the spr coal node
            if (last_nodes[spr->coal_node].parent == newcoal)
                spr->coal_node = newcoal;
        }
    }

    // determine if mapping of new node needs to be changed
    // newcoal was parent of recomb, it is broken
    if (fix_mapping) {
        if (last_nodes[spr->recomb_node].parent == newcoal) {
            if (spr->recomb_node != spr->coal_node) {
                mapping[newcoal] = -1;
                int p = last_nodes[newcoal].parent;
                if (p != -1)
                    mapping[p] = newcoal;
            }
        } else {
            // newcoal was not broken
            // find child without recomb or coal on it
            int x = newcoal;
            while (true) {
                int y = last_nodes[x].child[0];
                if (y == spr->coal_node || y == spr->recomb_node)
                    y = last_nodes[x].child[1];
                x = y;
                if (mapping[x] != -1)
                    break;
            }
            mapping[newcoal] = nodes[mapping[x]].parent;
        }
    }
}




// add a leaf thread to an ARG
void add_arg_thread(LocalTrees *trees, const StatesModel &states_model,
                    int ntimes, int *thread_path, int seqid,
                    vector<int> &recomb_pos, vector<Spr> &recombs,
		    const PopulationTree *pop_tree)
{
    unsigned int irecomb = 0;
    int nleaves = trees->get_num_leaves();
    int nnodes = trees->nnodes;
    int nnodes2 = nnodes + 2;

    // node names
    int newleaf = nleaves;
    int displaced = nnodes;
    int newcoal = nnodes + 1;

    States states;
    State last_state;
    LocalTree *last_tree = NULL;
#ifdef DEBUG
    static int count=0;
#endif

    assert_trees(trees, pop_tree, true);

    // update trees info
    trees->seqids.push_back(seqid);
    trees->nnodes = nnodes2;


    // loop through blocks
    int end = trees->start_coord;
    for (LocalTrees::iterator it=trees->begin(); it != trees->end(); ++it) {
#ifdef DEBUG
        count++;
#endif
        LocalTree *tree = it->tree;
        //Spr *spr = &(it->spr);
        int start = end;
        end += it->blocklen;
        states_model.get_coal_states(tree, states);

        // add new branch to local tree
        it->ensure_capacity(nnodes2);
        State state = states[thread_path[start]];
        add_tree_branch(tree, state.node, state.time, state.pop_path);

        // update mapping and spr
        int *mapping = it->mapping;
        if (mapping) {
            add_spr_branch(tree, last_tree, state, last_state,
                           &it->spr, mapping,
                           newleaf, displaced, newcoal);
            assert(assert_spr(last_tree, tree, &(it->spr), mapping, pop_tree));
        }

        // assert new branch is where it should be
        assert(tree->nodes[newcoal].age == states[thread_path[start]].time);

        // break this block for each new recomb within this block
        for (;irecomb < recombs.size() &&
              recomb_pos[irecomb] < end; irecomb++) {
            int pos = recomb_pos[irecomb];
            LocalNode *nodes = tree->nodes;
            state = states[thread_path[pos]];
            last_state = states[thread_path[pos-1]];

            // assert that thread time is still on track
            assert(tree->nodes[newcoal].age == last_state.time);

            // determine real name of recomb node
            // it may be different due to adding a new branch
            Spr spr2(recombs[irecomb]);
            if (spr2.recomb_node == newleaf)
                spr2.recomb_node = displaced;
            if (spr2.coal_node == newleaf)
                spr2.coal_node = displaced;
            assert(spr2.recomb_time <= tree->nodes[newcoal].age);

            // determine coal node and time
            //int istate = thread_path[pos];
            if (spr2.recomb_node == -1 && spr2.coal_node == -1) {
                spr2.recomb_node = spr2.coal_node = newleaf;
            } else if (spr2.recomb_node == -1) {
                // recomb on new branch, coal given thread
                spr2.recomb_node = newleaf;

                // fix coal node due to displacement
                if (spr2.coal_node == newleaf)
                    spr2.coal_node = displaced;

                // rename coal node due to newcoal underneath
                if (state.node == last_state.node &&
                    state.time > last_state.time)
                    spr2.coal_node = newcoal;

            } else {
                // recomb in ARG, coal on new branch
                if (state.time > last_state.time)
                    spr2.coal_node = nodes[newleaf].parent;
                else
                    spr2.coal_node = newleaf;
            }

            // determine mapping:
            // all nodes keep their name expect the broken node, which is the
            // parent of recomb
            int *mapping2 = new int [tree->capacity];
            for (int j=0; j<nnodes2; j++)
                mapping2[j] = j;
            if (spr2.recomb_node != spr2.coal_node)
                mapping2[nodes[spr2.recomb_node].parent] = -1;


            // make new local tree and apply SPR operation
            LocalTree *new_tree = new LocalTree(nnodes2, tree->capacity);
            new_tree->copy(*tree);
            apply_spr(new_tree, spr2, pop_tree);

            // calculate block end
            int block_end;
            if (irecomb < recombs.size() - 1)
                // use next recomb in this block to determine block end
                block_end = min(recomb_pos[irecomb+1], end);
            else
                // no more recombs in this block
                block_end = end;

            // insert new tree and spr into local trees list
            it->blocklen = pos - start;
            ++it;
            it = trees->trees.insert(it,
                LocalTreeSpr(new_tree, spr2, block_end - pos, mapping2));


            // assert tree and SPR
            //assert(assert_tree(new_tree));
            //assert(new_tree->nodes[newcoal].age == state.time);
            assert(assert_spr(tree, new_tree, &spr2, mapping2, pop_tree, false));

            // remember the previous tree for next iteration of loop
            tree = new_tree;
            nodes = tree->nodes;
            start = pos;
        }

        // remember the previous tree for next iteration of loop
        last_tree = tree;
        last_state = states[thread_path[end-1]];
        if (last_state.node == newleaf)
            last_state.node = displaced;
    }

    assert_trees(trees, pop_tree, false);
}



// Removes a leaf thread from an ARG
// NOTE: if remove_leaf is not last_leaf, nleaves - 1,
// last_leaf is renamed to remove_leaf
void remove_arg_thread(LocalTrees *trees, int remove_seqid,
                       const ArgModel *model)
{
    int nnodes = trees->nnodes;
    int nleaves = trees->get_num_leaves();
    int displace[nnodes];
    int last_leaf = nleaves - 1;

    // find leaf to remove from seqid
    int remove_leaf = -1;
    for (unsigned int i=0; i<trees->seqids.size(); i++) {
        if (trees->seqids[i] == remove_seqid) {
            remove_leaf = i;
            break;
        }
    }
    assert(remove_leaf != -1);

    // special case for trunk genealogy
    if (nnodes == 3) {
        assert(remove_leaf == 0 || remove_leaf == 1);
        int seqid = trees->seqids[1-remove_leaf];
        int pop_path = 0;
        // need to choose population for trunk genealogy; this should work
        // for now, but if we switch to admixture model this may need to be
        // changed. It might be better not to treat this as a special case.
        if (model->pop_tree != NULL) {
            const PopulationTree *pop_tree = model->pop_tree;
            LocalTree *tree = trees->begin()->tree;
            pop_path = tree->nodes[1-remove_leaf].pop_path;
            int start_pop = pop_tree->path_pop(pop_path, 0);
            pop_path = pop_tree->most_likely_path(start_pop);
        }
        trees->make_trunk(trees->start_coord, trees->end_coord, seqid, pop_path,
                          trees->begin()->tree->capacity);
        return;
    }


    // remove extra branch
    for (LocalTrees::iterator it=trees->begin(); it != trees->end(); ++it) {
        LocalTree *tree = it->tree;
        LocalNode *nodes = tree->nodes;

        // get information about removal
        int remove_coal = nodes[remove_leaf].parent;
        int coal_time = nodes[remove_coal].age;
        int *c = nodes[remove_coal].child;
        int coal_child = (c[0] == remove_leaf ? c[1] : c[0]);

        // remove branch from tree
        remove_tree_branch(tree, remove_leaf, model, displace);
        //assert_tree(tree);


        // fix this mapping due to displacement
        int *mapping = it->mapping;
        if (mapping) {
            for (int i=0; i<nnodes-2; i++)
                if (mapping[i] != -1)
                    mapping[i] = displace[mapping[i]];
                else
                    mapping[i] = -1;
        }

        // get next tree
        LocalTrees::iterator it2 = it;
        ++it2;
        if (it2 == trees->end())
            continue;

        // fix next mapping due to displacement
        mapping = it2->mapping;
        if (displace[last_leaf] != -1)
            mapping[displace[last_leaf]] = mapping[last_leaf];
        if (displace[nnodes-2] != -1)
            mapping[displace[nnodes-2]] = mapping[nnodes-2];
        if (displace[nnodes-1] != -1)
            mapping[displace[nnodes-1]] = mapping[nnodes-1];



        // fix SPR
        Spr *spr = &it2->spr;

        // get new name of coal_child
        coal_child = displace[coal_child];

        // if recomb is on branch removed, prune it
        if (spr->recomb_node == remove_leaf) {
            spr->set_null();
            continue;
        }

        // see if recomb node is renamed
        if (spr->recomb_node == remove_coal) {
            spr->recomb_node = coal_child;
        } else {
            // rename recomb_node due to displacement
            spr->recomb_node = displace[spr->recomb_node];
        }

        // if recomb is on root branch, prune it
        if (spr->recomb_node == coal_child && nodes[coal_child].parent == -1) {
            spr->set_null();
            continue;
        }

        // rename spr coal_node
        if (spr->coal_node == remove_leaf) {
            // mediated coal
            if (model->pop_tree != NULL)
                spr->pop_path = model->pop_tree->consistent_path(
                     spr->pop_path, nodes[coal_child].pop_path,
                     spr->recomb_time, spr->coal_time, coal_time);
            spr->coal_node = coal_child;
            spr->coal_time = coal_time;
        } else if (spr->coal_node == remove_coal) {
            // move coal down a branch
            spr->coal_node = coal_child;
        } else {
            // rename recomb_node due to displacement
            spr->coal_node = displace[spr->coal_node];
        }


        // check for bubbles
        if (spr->recomb_node == spr->coal_node &&
            (model->paths_equal(nodes[spr->recomb_node].pop_path,
                                spr->pop_path,
                                spr->recomb_time, spr->coal_time))) {
            spr->set_null();
            continue;
        }
    }


    // update trees info
    trees->seqids[remove_leaf] = trees->seqids[last_leaf];
    trees->seqids.resize(nleaves - 1);
    trees->nnodes -= 2;

    // remove extra trees
    remove_null_sprs(trees, model->pop_tree);

    assert_trees(trees, model->pop_tree, true);
}


//=============================================================================
// internal branch threading operations


// find the next possible branches in a removal path
void get_next_removal_nodes(const LocalTree *tree1, const LocalTree *tree2,
                            const Spr &spr2, const int *mapping2,
                            int node, int next_nodes[2], int recoal)
{

    if (spr2.coal_node == spr2.recomb_node) {
        next_nodes[0] = mapping2[node];
        next_nodes[1] = -1;
        assert(next_nodes[0] != -1);
        return;
    }

    if (recoal == -1)
        recoal = get_recoal_node(tree1, spr2, mapping2);

    // get passive transition
    next_nodes[0] = mapping2[node];
    if (next_nodes[0] == -1) {
        // node is broken by SPR
        // next node is then non-recomb child or recoal
        int sib = tree1->get_sibling(spr2.recomb_node);
        if (spr2.coal_node == sib)
            next_nodes[0] = recoal;
        else
            next_nodes[0] = mapping2[sib];
    }

    // get possible active transition
    // if recoal is on this branch (node) then there is a split in the path
    if (spr2.coal_node == node) {
        // find recoal node, its the node with no inward mappings
        next_nodes[1] = recoal;
    } else {
        // no second transition
        next_nodes[1] = -1;
    }
}


// find the next possible branches in all possible removal path
void get_all_next_removal_nodes(const LocalTree *tree1, const LocalTree *tree2,
                                const Spr &spr2, const int *mapping2,
                                int next_nodes[][2])
{
    const int recoal = get_recoal_node(tree1, spr2, mapping2);

    for (int node=0; node<tree1->nnodes; node++) {
        get_next_removal_nodes(tree1, tree2, spr2, mapping2,
                               node, next_nodes[node], recoal);
        assert(next_nodes[node][0] != next_nodes[node][1]);
    }
}


// find the previous possible branches in a removal path
void get_prev_removal_nodes(const LocalTree *tree1, const LocalTree *tree2,
                            const Spr &spr2, const int *mapping2,
                            int node, int prev_nodes[2], int *inv_mapping)
{
    const int nnodes = tree1->nnodes;

    if (spr2.coal_node == spr2.recomb_node) {
        prev_nodes[0] = prev_nodes[1] = -1;
        for (int i=0; i < nnodes; i++)
            if (mapping2[i] == node) {
                prev_nodes[0] = i;
                break;
            }
        assert(prev_nodes[0] != -1);
        return;
    }

    // make inverse mapping
    int inv_mapping_alloc[nnodes];
    if (!inv_mapping) {
        inv_mapping = inv_mapping_alloc;
        fill(inv_mapping, inv_mapping + nnodes, -1);
        for (int i=0; i<nnodes; i++)
            if (mapping2[i] != -1)
                inv_mapping[mapping2[i]] = i;
    }

    // get first transition
    prev_nodes[0] = inv_mapping[node];
    if (prev_nodes[0] == -1) {
        // there is no inv_mapping because node is recoal
        // the node spr.coal_node therefore is the previous node
        prev_nodes[0] = spr2.coal_node;

        // get optional second transition
        int sib = tree1->get_sibling(spr2.recomb_node);
        if (sib == spr2.coal_node) {
            prev_nodes[1] = tree1->nodes[sib].parent;
        } else
            prev_nodes[1] = -1;
    } else {
        // get optional second transition
        int sib = tree1->get_sibling(spr2.recomb_node);
        if (mapping2[sib] == node && sib != spr2.coal_node) {
            prev_nodes[1] = tree1->nodes[sib].parent;
        } else
            prev_nodes[1] = -1;
    }
    assert(prev_nodes[0] >= 0 || prev_nodes[1] >= 0);
}


// find the previous possible branches in a removal path
void get_all_prev_removal_nodes(const LocalTree *tree1, const LocalTree *tree2,
                                const Spr &spr2, const int *mapping2,
                                int prev_nodes[][2])
{
    const int nnodes = tree1->nnodes;

    // make inverse mapping
    int inv_mapping[nnodes];
    fill(inv_mapping, inv_mapping + nnodes, -1);
    for (int i=0; i<nnodes; i++)
        if (mapping2[i] != -1)
            inv_mapping[mapping2[i]] = i;

    for (int node=0; node<tree1->nnodes; node++) {
        get_prev_removal_nodes(tree1, tree2, spr2, mapping2, node,
                               prev_nodes[node], inv_mapping);
    }
}


// sample a removal path that only contains leaves
// TODO: this could be simplified, calculating leaf path are easy
void sample_arg_removal_leaf_path(const LocalTrees *trees, int node, int *path)
{
    int i = 0;
    path[i++] = node;
    LocalTree *last_tree = NULL;

    for (LocalTrees::const_iterator it=trees->begin();
         it != trees->end(); ++it)
    {
        LocalTree *tree = it->tree;

        if (last_tree) {
            int next_nodes[2];
            const Spr *spr = &it->spr;
            const int *mapping = it->mapping;
            get_next_removal_nodes(last_tree, tree, *spr, mapping,
                                   path[i-1], next_nodes);
            path[i++] = next_nodes[0];

            // ensure that a removal path re-enters the local tree correctly
            if (last_tree->root == path[i-2] && tree->root != path[i-1]) {
                assert(spr->coal_node == last_tree->root);
            }
        }

        last_tree = tree;
    }
}


// sample a removal path forward along an ARG
void sample_arg_removal_path_forward(
    const LocalTrees *trees, LocalTrees::const_iterator it,
    int node, int *path, int i, double prob_switch=.5)
{
    path[i++] = node;
    LocalTree *last_tree = NULL;

    for (; it != trees->end(); ++it) {
        LocalTree *tree = it->tree;

        if (last_tree) {
            int next_nodes[2];
            const Spr *spr = &it->spr;
            const int *mapping = it->mapping;
            get_next_removal_nodes(last_tree, tree, *spr, mapping,
                                   path[i-1], next_nodes);
            int j;
            if (next_nodes[1] == -1)
                j = 0;
            else
                j = int(rand() < prob_switch);
            path[i++] = next_nodes[j];

            // ensure that a removal path re-enters the local tree correctly
            if (last_tree->root == path[i-2] && tree->root != path[i-1]) {
                assert(spr->coal_node == last_tree->root);
            }
        }

        last_tree = tree;
    }
}


// sample a removal path backward along an ARG
void sample_arg_removal_path_backward(
    const LocalTrees *trees, LocalTrees::const_iterator it,
    int node, int *path, int i, double prob_switch=.5)
{
    path[i--] = node;

    const LocalTree *tree2 = it->tree;
    const Spr *spr2 = &it->spr;
    const int *mapping2 = it->mapping;
    --it;

    for (; it != trees->end(); --it) {
        int prev_nodes[2];
        LocalTree *tree1 = it->tree;
        assert(!spr2->is_null());

        get_prev_removal_nodes(tree1, tree2, *spr2, mapping2,
                               path[i+1], prev_nodes);
        int j;
        if (prev_nodes[1] == -1)
            j = 0;
        else
            j = int(rand() < prob_switch);
        path[i--] = prev_nodes[j];

        spr2 = &it->spr;
        mapping2 = it->mapping;
        tree2 = tree1;
    }
}


// sample a removal path that goes through a particular node and position
// in the ARG
void sample_arg_removal_path(const LocalTrees *trees, int node, int pos,
                             int *path, double prob_switch)
{
    // search for block with pos
    LocalTrees::const_iterator it = trees->begin();
    int end = trees->start_coord;
    int i = 0;
    for (; it != trees->end(); ++it, i++) {
        int start = end;
        end += it->blocklen;
        if (start <= pos && pos < end)
            break;
    }

    // search forward
    sample_arg_removal_path_forward(trees, it, node, path, i, prob_switch);
    sample_arg_removal_path_backward(trees, it, node, path, i, prob_switch);
}


// sample a removal path that starts at a particular node in the ARG
void sample_arg_removal_path(const LocalTrees *trees, int node, int *path)
{
    // search for block with pos
    LocalTrees::const_iterator it = trees->begin();
    sample_arg_removal_path_forward(trees, it, node, path, 0);
}


// sample a removal path that perfers recombination baring branches
void sample_arg_removal_path_recomb(
    const LocalTrees *trees, double recomb_preference, int *path)
{
    const int ntrees = trees->get_num_trees();
    const int nnodes = trees->nnodes;

    // build forward table for removal path sampling
    typedef int next_row[2];
    double **forward = new_matrix<double>(ntrees, nnodes);
    next_row **backptrs = new_matrix<next_row>(ntrees, nnodes);
    double **trans = new_matrix<double>(ntrees, nnodes);


    // calculate prior
    fill(forward[0], forward[0] + nnodes, 1.0 / nnodes);

    // compute forward table
    LocalTrees::const_iterator it= trees->begin();
    LocalTree const *last_tree = it->tree;
    int next_nodes[nnodes][2];
    ++it;
    for (int i=1; i<ntrees; i++, ++it) {
        LocalTree const *tree = it->tree;
        const int *mapping = it->mapping;
        next_row *prev_nodes = backptrs[i];

        // get next and previous transitions
        get_all_next_removal_nodes(last_tree, tree, it->spr, mapping,
                                   next_nodes);
        get_all_prev_removal_nodes(last_tree, tree, it->spr, mapping,
                                   prev_nodes);

        // assert transitions
        for (int j=0; j<nnodes; j++) {
            int k = next_nodes[j][0];
            assert(prev_nodes[k][0] == j || prev_nodes[k][1] == j);
            k = next_nodes[j][1];
            if (k != -1)
                assert(prev_nodes[k][0] == j || prev_nodes[k][1] == j);
        }

        // get next spr
        LocalTrees::const_iterator it2 = it;
        it2++;
        const Spr &spr2 = it->spr;

        // calc transition probs
        for (int j=0; j<nnodes; j++)
            trans[i-1][j] = (next_nodes[j][1] != -1 ? .5 : 1.0);

        // calc forward column
        double norm = 0.0;
        for (int j=0; j<nnodes; j++) {
            double sum = 0.0;
            for (int ki=0; ki<2; ki++) {
                int k = backptrs[i][j][ki];
                if (k == -1)
                    continue;
                sum += trans[i-1][j] * forward[i-1][k];
            }

            double emit = ((!spr2.is_null() && spr2.recomb_node == j) ?
                           recomb_preference : 1.0 - recomb_preference);
            forward[i][j] = sum * emit;

            assert(!isnan(forward[i][j]));

            norm += forward[i][j];
        }

        // normalize column for numerical stability
        for (int j=0; j<nnodes; j++)
            forward[i][j] /= norm;


        last_tree = tree;
    }


    // choose last branch
    int i = ntrees-1;
    path[i] = sample(forward[i], nnodes);

    // stochastic traceback
    int j = path[i];
    i--;
    for (; i>=0; i--) {
        if (backptrs[i+1][j][1] == -1) {
            // only one path
            j = path[i] = backptrs[i+1][j][0];
        } else {
            // fork, sample path
            double probs[2];
            probs[0] = forward[i][backptrs[i+1][j][0]] * trans[i][j];
            probs[1] = forward[i][backptrs[i+1][j][1]] * trans[i][j];
            int ji = sample(probs, 2);
            j = path[i] = backptrs[i+1][j][ji];
        }
    }


    // clean up
    delete_matrix<double>(forward, ntrees);
    delete_matrix<next_row>(backptrs, ntrees);
    delete_matrix<double>(trans, ntrees);
}


//=============================================================================
// sample removal paths uniformly

// count number of removal paths
void count_arg_removal_paths(const LocalTrees *trees,
                             RemovalPaths &removal_paths)
{
    const int ntrees = trees->get_num_trees();
    const int nnodes = trees->nnodes;
    double **counts = removal_paths.counts;
    RemovalPaths::next_row **backptrs = removal_paths.backptrs;

    // compute forward table
    LocalTrees::const_iterator it= trees->begin();
    LocalTree const *last_tree = it->tree;

    // calculate first column
    fill(counts[0], counts[0] + nnodes, 0.0);

    ++it;
    for (int i=1; i<ntrees; i++, ++it) {
        LocalTree const *tree = it->tree;
        const int *mapping = it->mapping;

        // get back pointers
        get_all_prev_removal_nodes(last_tree, tree, it->spr, mapping,
                                   backptrs[i]);

        // calc counts column
        for (int j=0; j<nnodes; j++) {
            int *ptrs = backptrs[i][j];
            if (ptrs[0] < 0 && ptrs[1] < 0)
                counts[i][j] = -INFINITY;
            else if (ptrs[1] == -1)
                counts[i][j] = counts[i-1][ptrs[0]];
            else if (ptrs[0] == -1)
                counts[i][j] = counts[i-1][ptrs[1]];
            else
                counts[i][j] = logadd(counts[i-1][ptrs[0]],
                                      counts[i-1][ptrs[1]]);
        }
        last_tree = tree;
    }
}


// count total number of removal paths
double count_total_arg_removal_paths(const RemovalPaths &removal_paths)
{
    // count total number of paths
    return logsum(removal_paths.counts[removal_paths.ntrees - 1],
                  removal_paths.nnodes);
}


//sets path to the series of nodes that relate to hap's ancestry
// at time between time_interval and time_interval+1.
// returns integer vector of coordinates where this path is broken
// due to recomb on node of interest at or below time_interval
void get_arg_removal_path_by_ind_and_time(const LocalTrees *trees,
                                          int time_interval,
                                          int hap,
                                          int *path,
                                          vector<int> &break_coords,
                                          bool do_first) {
    int i=0;
    int coord = trees->start_coord;
    LocalTree *last_tree = NULL;
    int next_nodes[2];
    //    int prev_pop_path;

    for (LocalTrees::const_iterator it=trees->begin();
         it != trees->end(); ++it) {
        LocalTree *tree = it->tree;
        int node = hap;
        assert(tree->nodes[node].age <= time_interval);
        while (node != tree->root) {
            assert(tree->nodes[node].age <= time_interval);
            int parent = tree->nodes[node].parent;
            if (tree->nodes[parent].age > time_interval)
                break;
            node = parent;
        }

        if (last_tree) {
            get_next_removal_nodes(last_tree, tree, it->spr, it->mapping,
                                   path[i-1], next_nodes);
            int possible[2] = {-1,-1};
            for (int j=0; j < 2; j++) {
                if (next_nodes[j] == -1) continue;
                if (tree->nodes[next_nodes[j]].age <= time_interval &&
                    (next_nodes[j] == tree->root ||
                     tree->nodes[tree->nodes[next_nodes[j]].parent].age > time_interval))
                    possible[j] = next_nodes[j];
            }
            if (possible[0] == -1 && possible[1] == -1) {
                break_coords.push_back(coord);
            } else if (possible[0] >= 0 && possible[1] == -1) {
                node = possible[0];
            } else if (possible[1] >= 0 && possible[0] == -1) {
                node = possible[1];
            } else {
                // Need to decide what to do in this situation
                // want to follow the same pop_path but may need additional criteria
                //look at some examples
                assert(0);
            }
        }
        coord += it->blocklen;
        path[i] = node;
        if (do_first || i > 0)
            last_tree = tree;
        i++;
    }
}


// sample a removal path uniformly from all paths and return total path count
double sample_arg_removal_path_uniform(const LocalTrees *trees, int *path)
{
    // compute path counts table
    RemovalPaths removal_paths(trees);
    count_arg_removal_paths(trees, removal_paths);

    // convenience variables
    const int ntrees = trees->get_num_trees();
    const int nnodes = trees->nnodes;
    double **counts = removal_paths.counts;
    RemovalPaths::next_row **backptrs = removal_paths.backptrs;

    // sample last branch of path first weighted by path counts
    double weights[nnodes];
    double norm = logsum(counts[ntrees - 1], nnodes);
    for (int j=0; j<nnodes; j++)
        weights[j] = exp(counts[ntrees - 1][j] - norm);
    path[ntrees - 1] = sample(weights, nnodes);

    for (int i=ntrees-1; i>0; i--) {
        const int *ptrs = backptrs[i][path[i]];
        if (ptrs[1] == -1) {
            // single trace back
            path[i-1] = ptrs[0];
        } else if (ptrs[0] == -1) {
            path[i-1] = ptrs[1];
        } else {
            // sample traceback
            const double p1 = counts[i-1][ptrs[0]];
            const double p2 = counts[i-1][ptrs[1]];

            if (log(frand()) < (p1 - logadd(p1, p2)))
                path[i-1] = ptrs[0];
            else
                path[i-1] = ptrs[1];
        }
        assert(path[i-1] >= 0);
    }

    // count total number of paths
    return count_total_arg_removal_paths(removal_paths);
}


// count total number of removal paths
double count_total_arg_removal_paths(const LocalTrees *trees)
{
    // compute path counts table
    RemovalPaths removal_paths(trees);
    count_arg_removal_paths(trees, removal_paths);

    // count total number of paths
    return count_total_arg_removal_paths(removal_paths);
}


/*
Under development.

//=============================================================================
// branch cut sampling


void sample_arg_cut(const LocalTrees *trees, int ntimes,
                    LocalTrees::const_iterator &it, int *time, int *branch,
                    int window_start, int window_end)
{
    // sample site
    if (window_start == -1)
        window_start = trees->start_coord;
    if (window_end == -1)
        window_end = trees->end_coord;
    int site = irand(window_start, window_end);

    // find branches at site and time
    it = trees->get_block(site);
    assert(it != trees->end());
    const LocalTree *tree = it->tree;

    // sample time
    double weights[ntimes-1];
    for (int i=0; i<ntimes-1; i++)
        weights[i] = exp(-i / double(ntimes-2));
    *time = sample(weights, ntimes-1);

    // find branches that enter the time point
    vector<int> branches;
    for (int i=0; i<tree->nnodes; i++) {
        if ((tree->nodes[i].is_leaf() || tree->nodes[i].age < *time) &&
            (tree->nodes[i].parent == -1 ||
             *time <= tree->nodes[tree->nodes[i].parent].age))
            branches.push_back(i);
    }

    // sample branch
    *branch = branches[irand(branches.size())];
}

// sample a removal path using the branch cut method
void sample_arg_removal_path_cut(const LocalTrees *trees, int ntimes,
                                 int *path, int *cuttime,
                                 int *region_start, int *region_end,
                                 int window_start, int window_end)
{
    // sample branch to cut
    LocalTrees::const_iterator center_it;
    int branch;
    sample_arg_cut(trees, ntimes, center_it, cuttime, &branch,
                   window_start, window_end);

    // find region where branch exists

    // find equivalent branches to the right of cut site
    LocalTrees::const_iterator it = center_it;
    LocalTrees::const_iterator last_it = center_it;
    LocalTrees::const_iterator right_it = center_it;
    ++it;
    int branch2 = branch;
    vector<int> right_branches;
    right_branches.push_back(branch2);
    while (it != trees->end()) {
        // stop growing region if recombination occurs right below cut site
        if (it->spr.recomb_node == branch2 && it->spr.recomb_time < *cuttime)
            break;

        // find equivalent branch in next tree
        int next_nodes[2];
        get_next_removal_nodes(last_it->tree, it->tree, it->spr, it->mapping,
                               branch2, next_nodes);
        LocalNode *nodes = it->tree->nodes;
        int parent = nodes[next_nodes[0]].parent;
        if (parent != -1 && nodes[parent].age < *cuttime)
            branch2 = next_nodes[1];
        else
            branch2 = next_nodes[0];
        assert(branch2 != -1);

        right_branches.push_back(branch2);

        // advance local block iterators
        right_it = it;
        last_it = it;
        ++it;
    }


    // search left for recombination the removes branch below cuttime
    it = center_it;
    last_it = center_it;
    LocalTrees::const_iterator left_it = center_it;
    --last_it;
    vector<int> left_branches;
    branch2 = branch;
    while (last_it != trees->end()) {
        // find equivalent branch in previous tree
        int prev_nodes[2];
        LocalTree *last_tree = last_it->tree;
        LocalTree *tree = it->tree;
        int *mapping = it->mapping;
        const Spr &spr = it->spr;
        get_prev_removal_nodes(last_tree, tree, spr, mapping,
                               branch2, prev_nodes);
        LocalNode *nodes = last_it->tree->nodes;
        int parent = nodes[prev_nodes[0]].parent;
        if (parent != -1 && nodes[parent].age < *cuttime && prev_nodes[1] != -1)
            branch2 = prev_nodes[1];
        else
            branch2 = prev_nodes[0];
        assert(branch2 != -1);

        // stop growing region if recombination occurs right below cut site
        if (it->spr.recomb_node == branch2 && it->spr.recomb_time < *cuttime)
            break;

        left_branches.push_back(branch2);

        // advance local block iterators
        left_it = last_it;
        it = last_it;
        --last_it;
    }


    // get removal path
    int blocki = 0;
    int starti = 0;
    *region_start = trees->start_coord;
    for (it=trees->begin(); it!=left_it; ++it) {
        path[blocki++] = -1;
        *region_start += it->blocklen;
    }
    starti = blocki;
    *region_end = *region_start;
    for (unsigned int i=0; i<left_branches.size(); i++, ++it) {
        path[blocki++] = left_branches[left_branches.size()-1-i];
        *region_end += it->blocklen;
    }
    assert(it == center_it);
    for (unsigned int i=0; i<right_branches.size(); i++, ++it) {
        path[blocki++] = right_branches[i];
        *region_end += it->blocklen;
    }
    int endi = blocki;
    --it;
    assert(it == right_it);
    it = right_it; ++it;
    int len = 0;
    for (; it!=trees->end(); ++it) {
        path[blocki++] = -1;
        len += it->blocklen;
    }
    assert(*region_end + len == trees->end_coord);


    // assert remove path
    LocalTrees::const_iterator end = right_it;
    right_it++;
    int i = starti;
    int x = path[i];
    int next_nodes[2];
    for (it=left_it; it!=end; ++it, i++) {
        assert(x == path[i]);
        LocalTrees::const_iterator it2 = it;
        it2++;
        if (i < endi - 1) {
            assert(it2 != trees->end());
            get_next_removal_nodes(it->tree, it2->tree, it2->spr, it2->mapping,
                                   x, next_nodes);
            if (next_nodes[0] == path[i+1])
                x = next_nodes[0];
            else
                x = next_nodes[1];
        }
    }
}

*/


//=============================================================================
// internal branch adding and removing


// update an SPR and mapping after adding a new internal branch
void add_spr_branch(const LocalTree *tree, const LocalTree *last_tree,
                    State state, State last_state,
                    Spr *spr, int *mapping,
                    int subtree_root, int last_subtree_root,
                    const PopulationTree *pop_tree)
{
    // get tree info
    const LocalNode *nodes = tree->nodes;
    const LocalNode *last_nodes = last_tree->nodes;
    int node2 = state.node;
    int last_newcoal = last_nodes[last_subtree_root].parent;
    static int count=0;
    bool fix_mapping=true;
    count++;

#ifdef DEBUG
    Spr orig_spr(*spr);
    int orig_mapping[tree->nnodes];
    for (int i=0; i < tree->nnodes; i++)
        orig_mapping[i] = mapping[i];
#endif

    // determine newcoal
    int newcoal;
    if (state.node != -1) {
        newcoal = nodes[subtree_root].parent;
    } else {
        // fully specified tree
        if (mapping[last_subtree_root] != -1)
            newcoal = nodes[mapping[last_subtree_root]].parent;
        else {
            int sib = last_tree->get_sibling(spr->recomb_node);
            assert(mapping[sib] != -1);
            newcoal = nodes[mapping[sib]].parent;
        }
    }

    // set default new node mapping
    mapping[last_newcoal] = newcoal;

    // parent of recomb node should be the recoal point
    // however, if it equals newcoal, then either (1) the recomb branch is
    // renamed, (2) there is mediation, or (3) new branch escapes
    int recoal = nodes[mapping[spr->recomb_node]].parent;
    if (spr->recomb_node == spr->coal_node) {
        if (last_state.node == spr->recomb_node) {
            if (last_state.time < spr->recomb_time ||
                (last_state.time == spr->recomb_time &&
                 state.time == spr->recomb_time)) {
                spr->recomb_node = last_nodes[spr->recomb_node].parent;
                spr->coal_node = spr->recomb_node;
                fix_mapping = false;
            } else if (last_state.time >= spr->recomb_time &&
                       last_state.time < spr->coal_time) {
                if (state.time <= last_state.time) {
                    assert(state.time >= spr->recomb_time);
                    spr->coal_node = last_tree->get_sibling(spr->recomb_node);
                    spr->coal_time = state.time;
                    mapping[last_newcoal] = -1;
                    fix_mapping = false;
                } else {
                    assert(state.time == spr->coal_time);
                    spr->coal_node = last_newcoal;
                    mapping[last_newcoal] = -1;
                    fix_mapping = false;
                }
            } else if (last_state.time == spr->coal_time) {
                spr->coal_time = state.time;
                spr->coal_node = last_tree->get_sibling(spr->recomb_node);
                mapping[last_newcoal] = -1;
                fix_mapping = false;
            } // if last_state.time > spr->coal_time do nothing
        }
    } else if (recoal == newcoal) {
        if (mapping[last_state.node] == node2) {
            // (1) recomb is above coal state, we rename spr recomb node
            spr->recomb_node = last_newcoal;
        } else {
            // if this is a mediated coal, then state should equal recomb
            int state_node = state.node;
            if (spr->coal_time == last_nodes[last_newcoal].age &&
                state_node == mapping[spr->recomb_node]) {
                // (3) this is a mediated coal, rename coal node and time
                if (state.time < last_nodes[last_subtree_root].age) {
                    spr->coal_node = last_tree->get_sibling(spr->recomb_node);
                    assert(spr->coal_time >= state.time);
                    spr->coal_time = state.time;
                } else{
                    spr->coal_node = last_subtree_root;
                    assert(spr->coal_time >= state.time);
                    spr->coal_time = state.time;
                }
                assert(spr->coal_time >= last_nodes[spr->coal_node].age);
            } else {
                // (2) this is the new branch escaping
                // no oUther updates are necessary
            }
        }
    } else {
        // the other possibility is that newcoal is under recoal point
        // if newcoal is child of recoal, then coal is renamed
        const int *c = nodes[recoal].child;
        if (c[0] == newcoal || c[1] == newcoal) {
            // we either coal above the newcoal or our existing
            // node just broke and newcoal was underneath.

            // if newcoal was previously above spr->coal_node
            // then we rename the spr coal node
            if (last_nodes[spr->coal_node].parent == last_newcoal)
                spr->coal_node = last_newcoal;
            assert(spr->coal_time >= last_nodes[spr->coal_node].age);
            int p = last_nodes[spr->coal_node].parent;
            if (p != -1)
                assert(spr->coal_time <= last_nodes[p].age);
        }
    }

    // determine if mapping of new node needs to be changed
    // newcoal was parent of recomb, it is broken
    if (fix_mapping) {
        if (last_nodes[spr->recomb_node].parent == last_newcoal) {
            if (spr->recomb_node != spr->coal_node) {
                mapping[last_newcoal] = -1;
                int p = last_nodes[last_newcoal].parent;
                if (fix_mapping && p != -1 && spr->coal_node != last_newcoal)
                    mapping[p] = newcoal;
            }
        } else {
            // newcoal was not broken
            // find child without recomb or coal on it
            int x = last_newcoal;
            int y = last_nodes[x].child[0];
            if (y == spr->coal_node)
                y = last_nodes[x].child[1];
            if (mapping[y] == -1)
                y = last_tree->get_sibling(spr->recomb_node);
            if (y == spr->coal_node)
                y = last_nodes[x].child[1];
            mapping[last_newcoal] = nodes[mapping[y]].parent;
        }
    }
    assert(assert_spr(last_tree, tree, spr, mapping, pop_tree, false));
}



// Add a branch to a partial ARG
void add_arg_thread_path(LocalTrees *trees, const StatesModel &states_model,
                         int ntimes, const int *thread_path,
                         vector<int> &recomb_pos, vector<Spr> &recombs,
			 const PopulationTree *pop_tree)
{
    States states;
    LocalTree *last_tree = NULL;
    State last_state;
    int last_subtree_root = -1;
    unsigned int irecomb = 0;
    int end = trees->start_coord;
#ifdef DEBUG
    static int count=0;
    LocalTree orig_tree;
    LocalTree orig_last_tree;
    Spr orig_spr;
#endif

    for (LocalTrees::iterator it=trees->begin(); it != trees->end(); ++it)
    {
        LocalTree *tree = it->tree;
        LocalNode *nodes = tree->nodes;
        Spr *spr = &(it->spr);  // this is spr to left of tree
        State state;
        int start = end;
        end += it->blocklen;
        const int subtree_root = nodes[tree->root].child[0];
        states_model.get_coal_states(tree, states);
        int nstates = states.size();
 #ifdef DEBUG
        orig_tree.copy(*tree);
        orig_spr.copy(*spr);
        count++;
 #endif

        // detect whether local tree is partial
        if (nodes[tree->root].age > ntimes) {
            assert(nstates > 0);
            // modify local tree according to thread path
            assert(thread_path[start] >= 0 &&
                   thread_path[start] < (int)states.size());
            state = states[thread_path[start]];
            // this should only be true if we are re-threading an ancient lineage;
            // there will be a mismatch between the subtree age (which will be zero)
            // and the minage; need to make sure the paths are consistent
            if (pop_tree != NULL && states_model.minage > nodes[subtree_root].age) {
                state.pop_path = pop_tree->consistent_path(nodes[subtree_root].pop_path,
                                                           state.pop_path,
                                                           nodes[subtree_root].age,
                                                           states_model.minage,
                                                           state.time);
            }
            Spr add_spr(subtree_root, nodes[subtree_root].age,
                        state.node, state.time, state.pop_path);
            apply_spr(tree, add_spr, pop_tree);
        } else {
            // set null state
            state.set_null();
        }

        // fix spr
        // update mapping and spr
        int *mapping = it->mapping;
        if (mapping) {
            if (last_state.node != -1) {
                add_spr_branch(tree, last_tree, state, last_state,
                               spr, mapping, subtree_root, last_subtree_root,
                               pop_tree);
            }
        }

        // break this block for each new recomb within this block
        for (;irecomb < recombs.size() &&
              recomb_pos[irecomb] < end; irecomb++) {
            int pos = recomb_pos[irecomb];
            LocalNode *nodes = tree->nodes;
            assert(nstates > 0);

            state = states[thread_path[pos]];
            last_state = states[thread_path[pos-1]];
            int newcoal = nodes[subtree_root].parent;

            // assert that thread time is still on track
            assert(tree->nodes[newcoal].age == last_state.time);

            // determine real name of recomb node
            // it may be different due to adding a new branch
            Spr spr2(recombs[irecomb]);
            assert(spr2.recomb_time <= tree->nodes[newcoal].age);

            // determine coal node and time
            if (spr2.recomb_node == subtree_root && spr2.coal_node == subtree_root) {
                // only path of new branch changes; do nothing here
            } if (spr2.recomb_node == subtree_root) {
                // recomb on new branch, coal given thread
                // rename coal node due to newcoal underneath
                if (state.node == last_state.node &&
                    state.time > last_state.time)
                    spr2.coal_node = newcoal;
            } else {
                // recomb in maintree, coal on new branch
                if (state.time > last_state.time)
                    spr2.coal_node = nodes[subtree_root].parent;
            }

            // determine mapping:
            // all nodes keep their name except the broken node, which is the
            // parent of recomb
            int *mapping2 = new int [tree->capacity];
            for (int j=0; j<tree->nnodes; j++)
                mapping2[j] = j;
            if (spr2.recomb_node != spr2.coal_node)
                mapping2[nodes[spr2.recomb_node].parent] = -1;


            // make new local tree and apply SPR operation
            LocalTree *new_tree = new LocalTree(tree->nnodes, tree->capacity);
            new_tree->copy(*tree);
            apply_spr(new_tree, spr2, pop_tree);

            // calculate block end
            int block_end;
            if (irecomb < recombs.size() - 1)
                // use next recomb in this block to determine block end
                block_end = min(recomb_pos[irecomb+1], end);
            else
                // no more recombs in this block
                block_end = end;

            // insert new tree and spr into local trees list
            it->blocklen = pos - start;
            ++it;
            it = trees->trees.insert(it,
                LocalTreeSpr(new_tree, spr2, block_end - pos, mapping2));

            // remember the previous tree for next iteration of loop
            tree = new_tree;
            nodes = tree->nodes;
            start = pos;
        }

        // record previous local tree information
#ifdef DEBUG
        orig_last_tree.copy(orig_tree);
#endif
        last_tree = tree;
        last_state = state;
        last_subtree_root = subtree_root;
    }

    assert_trees(trees, pop_tree, false);
}



// get the new names for nodes due to collapsing null SPRs
LocalTree *get_actual_nodes(LocalTrees *trees, LocalTrees::iterator it,
                            int *nodes)
{
    LocalTree *tree = it->tree;
    const int nnodes = tree->nnodes;

    // get current nodes
    for (int i=0; i<nnodes; i++)
        nodes[i] = i;

    LocalTrees::iterator it2 = it;
    ++it2;
    for (;it2 != trees->end() && it2->spr.is_null(); ++it2) {
        // compute transitive mapping
        int *mapping = it2->mapping;
        for (int i=0; i<nnodes; i++) {
            if (mapping[i] != -1)
                nodes[i] = mapping[nodes[i]];
            else
                nodes[i] = -1;
        }
    }

    --it2;
    return it2->tree;
}


// Removes a thread path from an ARG and returns a partial ARG
void remove_arg_thread_path(LocalTrees *trees, const int *removal_path,
                            int maxtime, const PopulationTree *pop_tree,
                            int *original_thread)
{
    LocalTree *tree = NULL;
    State *original_states = NULL;
#ifdef DEBUG
    static int count=0;
    static int call_count=0;
    call_count++;
    assert_trees(trees, pop_tree, false);
    LocalTree orig_last_tree;
    LocalTree orig_tree;
    Spr orig_spr;
    Spr last_orig_spr;
#endif

    // prepare original thread array if requested
    if (original_thread) {
        original_states = new State [trees->length()];
    }

    int i = 0;
    int end = trees->start_coord;
    for (LocalTrees::iterator it=trees->begin(); it != trees->end(); ++it, i++)
    {
        LocalTree *last_tree = tree;
        tree = it->tree;
        LocalNode *nodes = tree->nodes;
        int start = end;
        end += it->blocklen;

#ifdef DEBUG
        count++;
        if (tree != NULL)
            orig_last_tree.copy(orig_tree);
        orig_tree.copy(*tree);
#endif

        int removal_node = removal_path[i];
        if (removal_node == tree->root) {
            // fix previous mapping
            if (it->mapping && removal_path[i-1] != last_tree->root)
                it->mapping[last_tree->root] = -1;

            // record thread
            if (original_states) {
                for (int j=start; j<end; j++) {
                    original_states[j-trees->start_coord].set(-1, -1, -1) ;
                }
           }

            // removal path has "fallen off the top" there is nothing to edit
            continue;
        }

        // modify local into subtree-maintree format
        int broken_node = nodes[removal_node].parent;
        int coal_time = nodes[broken_node].age;
        int broken_child = tree->get_sibling(removal_node);
        int orig_path = nodes[removal_node].pop_path;
        int pop_path = ( pop_tree == NULL ? 0 :
                         pop_tree->path_to_root(nodes, removal_node) );
        Spr removal_spr(removal_node, nodes[removal_node].age,
                        tree->root, maxtime, pop_path);

        apply_spr(tree, removal_spr, pop_tree);


        // determine subtree and maintree roots
        int subtree_root = removal_node;
        int maintree_root = tree->get_sibling(subtree_root);

        // ensure subtree is the first child of the root
        int *c = nodes[tree->root].child;
        if (c[0] == maintree_root) {
            c[0] = subtree_root;
            c[1] = maintree_root;
        }

        // fix previous mapping
        if (it->mapping) {
            assert(last_tree);
            if (removal_path[i-1] != last_tree->root)
                it->mapping[last_tree->root] = tree->root;
        }

        // record thread
        if (original_states) {
            for (int j=start; j<end; j++)
                original_states[j-trees->start_coord].set(
                    broken_child, coal_time, orig_path);
        }


        // get next tree
        LocalTrees::iterator it2 = it;
        ++it2;
        if (it2 == trees->end())
            continue;


        // fix SPR
        Spr *spr = &it2->spr;
        int *mapping = it2->mapping;

#ifdef DEBUG
        if (last_tree != NULL) last_orig_spr.copy(orig_spr);
        orig_spr.copy(*spr);
        int orig_mapping[tree->nnodes];
        for (int j=0; j < tree->nnodes; j++) orig_mapping[j] = mapping[j];
#endif

        // if recomb is on branch removed, prune it
        if (spr->recomb_node == removal_node) {
            int p = nodes[spr->recomb_node].parent;
            assert(mapping[p] != -1 || p == tree->root);
            spr->set_null();
        } else {

        // see if recomb node is renamed
        if (spr->recomb_node == broken_node) {
            spr->recomb_node = broken_child;
        }

        if (spr->coal_node == removal_node) {
            // detect branch path splits
            int next_nodes[2];
            get_next_removal_nodes(tree, it2->tree, *spr, mapping,
                                   removal_path[i], next_nodes);
            if (removal_path[i+1] == next_nodes[0]) {
                // removal path chooses lower path
                // note: sister_node = broken_child

                if (spr->recomb_node == broken_child) {
                    // spr is now bubble, prune it
                    if (pop_tree == NULL ||
                        (pop_tree->paths_equal(tree->nodes[spr->recomb_node].pop_path,
                                               spr->pop_path,
                                               spr->recomb_time,
                                               spr->coal_time) &&
                         pop_tree->paths_equal(tree->nodes[spr->coal_node].pop_path,
                                               tree->nodes[spr->recomb_node].pop_path,
                                               spr->coal_time, coal_time))) {
                        int p = nodes[spr->recomb_node].parent;
                        assert(mapping[p] != -1 || p == tree->root);
                        spr->set_null();
                    } else {
                        spr->coal_node = spr->recomb_node;
                        spr->pop_path = pop_tree->consistent_path(spr->pop_path,
                            tree->nodes[removal_node].pop_path,
                            spr->recomb_time, spr->coal_time, -1);
                        spr->coal_time = coal_time;
                    }
                } else {
                    // recomb is on non-sister branch, therefore it is a
                    // mediated coalescence
                    if (pop_tree != NULL)
                        spr->pop_path = pop_tree->consistent_path(
                             spr->pop_path, tree->nodes[removal_node].pop_path,
                             spr->recomb_time, spr->coal_time,
                             coal_time);
                    spr->coal_node = broken_child;
                    spr->coal_time = coal_time;
                }
            } else if (removal_path[i+1] == next_nodes[1]) {
                // removal path chooses upper path
                // keep spr recoal where it is
                // assert that upper path is the new recoal node
                // nobody should map to the new recoal node
                for (int j=0; j<tree->nnodes; j++)
                    assert(mapping[j] != removal_path[i+1]);
            } else {
                assert(0);
            }
        } else if (spr->coal_node == broken_node) {
            // rename spr recoal
            spr->coal_node = broken_child;
        }
        }

        // check for bubbles
        if (spr->is_null()) {
        } else if (spr->recomb_node == spr->coal_node) {
            if (pop_tree == NULL ||
                pop_tree->paths_equal(tree->nodes[spr->recomb_node].pop_path,
                                      spr->pop_path,
                                      spr->recomb_time,
                                      spr->coal_time)) {
                int p = nodes[spr->recomb_node].parent;
                assert(mapping[p] != -1 || p == tree->root);
                for (int j=0; j < tree->nnodes; j++) {
                    if (j != removal_node && j != tree->root) {
                        int maxt = min(maxtime-2, tree->nodes[tree->nodes[j].parent].age);
                        int parent2 = it2->tree->nodes[it2->mapping[j]].parent;
                        if (parent2 != -1)
                            maxt = min(maxt, it2->tree->nodes[parent2].age);
                        assert(pop_tree->paths_equal(tree->nodes[j].pop_path,
                                                     it2->tree->nodes[it2->mapping[j]].pop_path,
                                                     tree->nodes[j].age, maxt));
                    }
                    assert(mapping[j] != -1 || j== tree->root);
                }
                spr->set_null();
                continue;
            }
        } else {
            // ensure broken node maps to -1
            int spr_broken_node = nodes[spr->recomb_node].parent;
            mapping[spr_broken_node] = -1;
        }
#ifdef DEBUG
        // assert spr
        if (last_tree)
            assert_spr(last_tree, tree, &it->spr, it->mapping, pop_tree, true);
#endif
    }

    // record original thread
    if (original_states) {
        // TODO: may not want to assume ntimes = maxtime - 1
        const int ntimes = maxtime - 1;
        const int nnodes = trees->nnodes;
        States states;

        int end = trees->start_coord;
        for (LocalTrees::iterator it=trees->begin(); it!=trees->end(); ++it) {
            int start = end;
            end += it->blocklen;

            int nodes_lookup[nnodes];
            LocalTree *tree2 = get_actual_nodes(trees, it, nodes_lookup);

            get_coal_states_internal(tree2, ntimes, states);
            int nstates = states.size();
            NodeStateLookup lookup(states, 0, pop_tree);
            for (int i=start; i<end; i++) {
                if (nstates == 0) {
                    original_thread[i-trees->start_coord] = 0;
                } else {
                    int statei = lookup.lookup(
                        nodes_lookup[original_states[i-trees->start_coord].node],
                        original_states[i-trees->start_coord].time,
                        original_states[i-trees->start_coord].pop_path);
                    assert(statei != -1);
                    original_thread[i-trees->start_coord] = statei;
                }
            }
        }

        delete [] original_states;
    }

#ifdef DEBUG
    assert_trees(trees, pop_tree, true);
#endif

    // remove extra trees
    remove_null_sprs(trees, pop_tree);

    assert_trees(trees, pop_tree, true);
}



//=============================================================================
// C interface
/*
extern "C" {

void arghmm_sample_arg_removal_path(LocalTrees *trees, int node, int *path)
{
    sample_arg_removal_path(trees, node, path);
}

void arghmm_sample_arg_removal_path2(LocalTrees *trees, int node,
                                     int pos, int *path)
{
    sample_arg_removal_path(trees, node, pos, path);
}

void arghmm_sample_arg_removal_leaf_path(LocalTrees *trees, int node,
                                           int *path)
{
    sample_arg_removal_leaf_path(trees, node, path);
}

void arghmm_sample_arg_removal_path_recomb(LocalTrees *trees,
                                           double recomb_preference, int *path)
{
    sample_arg_removal_path_recomb(trees, recomb_preference, path);
}


void arghmm_remove_arg_thread_path(LocalTrees *trees, int *removal_path,
                                   int maxtime)
{
    remove_arg_thread_path(trees, removal_path, maxtime);
}

void arghmm_remove_arg_thread_path2(LocalTrees *trees, int *removal_path,
                                    int maxtime, int *original_thread)
{
    remove_arg_thread_path(trees, removal_path, maxtime, original_thread);
}


void arghmm_get_thread_times(LocalTrees *trees, int ntimes, int *path,
                             int *path_times)
{
    States states;

    int end = trees->start_coord;
    for (LocalTrees::iterator it=trees->begin(); it != trees->end(); ++it) {
        int start = end;
        end += it->blocklen;
        LocalTree *tree = it->tree;
        get_coal_states_internal(tree, ntimes, states);

        for (int i=start; i<end; i++)
            path_times[i] = states[path[i]].time;
    }
}



} // extern C
*/
} // namespace argweaver
