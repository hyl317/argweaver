//=============================================================================
// Local trees

#ifndef ARGWEAVER_LOCAL_TREES_H
#define ARGWEAVER_LOCAL_TREES_H

// c++ includes
#include <assert.h>
#include <list>
#include <vector>
#include <string.h>
#include <stdio.h>
#include <iostream>
#include <memory>

// arghmm includes
#include "sequences.h"
#include "model.h"
//#include "rspr.h"

namespace argweaver {

using namespace std;

//class PopulationTree;



// A block within a sequence alignment
class Block
{
public:
    Block() {}
    Block(int start, int end) :
        start(start), end(end) {}
    int start;
    int end;

    int length() {
        return end - start;
    }
};


// A Subtree Pruning and Regrafting operation
class Spr
{
public:
    Spr() {}
    Spr(int recomb_node, int recomb_time,
        int coal_node, int coal_time, int pop_path=-1) :
            recomb_node(recomb_node), recomb_time(recomb_time),
            coal_node(coal_node), coal_time(coal_time),
            pop_path(pop_path) {}
    Spr(const Spr &other) {
        copy(other);
    }
    void copy(const Spr &other) {
        recomb_node = other.recomb_node;
        recomb_time = other.recomb_time;
        coal_node = other.coal_node;
        coal_time = other.coal_time;
        pop_path = other.pop_path;
    }

    // sets the SPR to a null value
    void set_null()
    {
        recomb_node = -1;
        recomb_time = -1;
        coal_node = -1;
        coal_time = -1;
        pop_path = -1;
    }

    bool is_null() const
    {
        return recomb_node == -1;
    }

    void write() const {
        printf("rn=%i rt=%i cn=%i ct=%i pp=%i\n", recomb_node,
               recomb_time, coal_node, coal_time, pop_path);
    }

    int recomb_node;
    int recomb_time;
    int coal_node;
    int coal_time;
    int pop_path;
};



// A node in a local tree
class LocalNode
{
public:
    LocalNode()
    {}

    LocalNode(int parent, int left_child, int right_child, int age=-1,
              int pop_path=0) :
        parent(parent), age(age), pop_path(pop_path)
    {
        child[0] = left_child;
        child[1] = right_child;
    }
    LocalNode(const LocalNode &other) :
        parent(other.parent),
        age(other.age)
    {
        child[0] = other.child[0];
        child[1] = other.child[1];
        pop_path = other.pop_path;
    }

    ~LocalNode() {
    }

    inline bool is_leaf() const
    {
        return child[0] == -1;
    }

    inline int add_child(int child_node)
    {
        if (child[0] == -1) {
            child[0] = child_node;
            return 0;
        } else if (child[1] == -1) {
            child[1] = child_node;
            return 1;
        }

        // already have two children
        return -1;
    }

    inline void copy(const LocalNode &other)
    {
        parent = other.parent;
        age = other.age;
        child[0] = other.child[0];
        child[1] = other.child[1];
        pop_path = other.pop_path;
    }

    void set_pop_path(int path) {
        pop_path = path;
    }

    inline int get_pop(int time, const PopulationTree *pop_tree) const {
        if (pop_tree == NULL) return 0;
        return pop_tree->get_pop(pop_path, time);
    }

    int parent;
    int child[2];
    int age;
    int pop_path;
};

extern LocalNode null_node;


// A local tree in a set of local trees
//
//   Leaves are always listed first in nodes array
//
class LocalTree
{
public:
    LocalTree() :
        nnodes(0),
        capacity(0),
        root(-1),
        nodes(NULL)
    {}

    LocalTree(int nnodes, int capacity=0) :
        nnodes(nnodes),
        capacity(capacity),
        root(-1)
    {
        if (capacity < nnodes)
            capacity = nnodes;
        nodes = new LocalNode [capacity];
    }


    LocalTree(int *ptree, int nnodes, int *ages=NULL, int *paths=NULL,
           int capacity=-1) :
        nnodes(nnodes),
        capacity(0),
        root(-1),
        nodes(NULL)
    {
        set_ptree(ptree, nnodes, ages, paths, capacity);
    }

    LocalTree(const LocalTree &other) :
        nnodes(0),
        capacity(0),
        root(-1),
        nodes(NULL)
    {
        copy(other);
    }


    ~LocalTree() {
        if (nodes) {
            delete [] nodes;
            nodes = NULL;
        }
    }

    // initialize a local tree by on a parent array
    void set_ptree(int *ptree, int _nnodes, int *ages=NULL, int *paths=NULL,
                   int _capacity=-1)
    {
        // delete existing nodes if they exist
        if (nodes) {
	    delete [] nodes;
	}

        nnodes = _nnodes;
        if (_capacity >= 0)
            capacity = _capacity;
        if (capacity < nnodes)
            capacity = nnodes;

        nodes = new LocalNode [capacity];

        // populate parent pointers
        for (int i=0; i<nnodes; i++) {
            nodes[i].parent = ptree[i];
            nodes[i].child[0] = -1;
            nodes[i].child[1] = -1;
        }

        // initialize age if given
        if (ages)
            for (int i=0; i<nnodes; i++)
                nodes[i].age = ages[i];

        if (paths){
            for (int i=0; i<nnodes; i++)
                nodes[i].pop_path = paths[i];
        }else{
            for(int i=0; i<nnodes; i++){
                nodes[i].pop_path = 0; //need to set a default/dummy path, otherwise we get "use of uninitialized error"
            }
        }

        // populate children pointers
        for (int i=0; i<nnodes; i++) {
            const int parent = ptree[i];

            if (parent != -1) {
                int *child = nodes[parent].child;
                if (child[0] == -1)
                    child[0] = i;
                else
                    child[1] = i;
            } else {
                root = i;
            }
        }
    }


    // Sets the root of the tree by finding node without a parent
    void set_root()
    {
        for (int j=0; j<nnodes; j++) {
            if (nodes[j].parent == -1) {
                root = j;
                break;
            }
        }
    }


    // Sets a new capacity for the allocated data structures
    void set_capacity(int _capacity)
    {
        if (_capacity == capacity)
            return;

        LocalNode *tmp = new LocalNode[_capacity];
	assert(tmp);

        std::copy(nodes, nodes + capacity, tmp);
        delete [] nodes;

        nodes = tmp;
        capacity = _capacity;
    }


    // Ensures that we have a certain capacity.  If capacity is not >=
    // than _capacity, increase capacity
    void ensure_capacity(int _capacity)
    {
        if (capacity < _capacity)
            set_capacity(_capacity);
    }



    // Returns the postorder traversal of the nodes
    void get_postorder(int *order) const
    {
        char visit[nnodes];
        int i;

        // initialize array of number of visits to a node
        for (i=0; i<nnodes; i++)
            visit[i] = 0;

        // list leaves first
        for (i=0; i<nnodes; i++) {
            if (!nodes[i].is_leaf())
                break;
            order[i] = i;
        }

        // add the remaining nodes
        int end = i;
        for (i=0; i<nnodes; i++) {
            int parent = nodes[order[i]].parent;
            if (parent != -1) {
                visit[parent]++;

                // add parent to queue if both children have been seen
                if (visit[parent] == 2)
                    order[end++] = parent;
            }
        }
    }


    void get_preorder(int node, int *order, int &norder) const
    {
        // start queue
        int queue[nnodes];
        int queuei = 0;
        queue[queuei++] = node;
        norder = 0;

        while (queuei > 0) {
            int node2 = queue[--queuei];
            order[norder++] = node2;

            if (!nodes[node2].is_leaf()) {
                queue[queuei++] = nodes[node2].child[0];
                queue[queuei++] = nodes[node2].child[1];
            }
        }
    }


    // Returns the number of leaves in the tree
    inline int get_num_leaves() const
    {
        return (nnodes + 1) / 2;
    }


    // Convenience method for accessing nodes
    inline LocalNode &operator[](int name) const
    {
        return nodes[name];
    }


    inline LocalNode &get_node(int name) const
    {
        return nodes[name];
    }

    inline LocalNode &get_root() const
    {
        if (root == -1)
            return null_node;
        else
            return nodes[root];
    }


    inline double get_dist(int node, const double *times) const
    {
        int parent = nodes[node].parent;
        if (parent != -1)
            return times[nodes[parent].age] - times[nodes[node].age];
        else
            return 0.0;
    }


    // add a node to the tree and return its name
    inline int add_node()
    {
        nnodes++;
        if (nnodes > capacity)
            ensure_capacity(2*nnodes);
        return nnodes - 1;
    }


    // clear nodes
    inline void clear()
    {
        nnodes = 0;
        root = -1;
    }


    // Copy tree structure from another tree
    inline void copy(const LocalTree &other)
    {
        // copy tree info
        nnodes = other.nnodes;
        ensure_capacity(other.capacity);
        root = other.root;

        // copy node info
        for (int i=0; i<nnodes; i++)
            nodes[i].copy(other.nodes[i]);
    }


    // get the sibling of a node
    inline int get_sibling(int node) const {
        int parent = nodes[node].parent;
        if (parent == -1)
            return -1;
        int *c = nodes[parent].child;
        if (c[0] == node)
            return c[1];
        else
            return c[0];
    }


    // add a child to a node in the tree
    inline int add_child(int parent, int child) {
        int childi = nodes[parent].add_child(child);
        if (childi == -1)
            return -1;

        nodes[child].parent = parent;

        return childi;
    }

    // inline map<set<int>, int> descent_leaf_map(){
    //     int order[nnodes];
    //     get_postorder(order);
    //     int num_leaves = get_num_leaves();
    //     map<set<int>, int> descent_map;
    //     map<int, set<int>> reverse_descent_map;
    //     for(int i = 0; i < nnodes; i++){
    //         int j = order[i];
    //         set<int> tmp;
    //         if (j < num_leaves){
    //             tmp.insert(j);
    //             descent_map.insert(pair<set<int>, int>(tmp, j));
    //             reverse_descent_map.insert(pair<int, set<int>>(j, tmp));
    //         }else{
    //             int *child = nodes[j].child;
    //             set<int> leaf_set1 = reverse_descent_map.find(child[0])->second;
    //             set<int> leaf_set2 = reverse_descent_map.find(child[1])->second;
    //             //set_union(leaf_set1->begin(), leaf_set1->end(),
    //             //            leaf_set2->begin(), leaf_set2->end(),
    //             //            insert_iterator<set<int>>(*tmp, tmp->begin()));
    //             tmp.insert(leaf_set1.begin(), leaf_set1.end());
    //             tmp.insert(leaf_set2.begin(), leaf_set2.end());
    //             descent_map.insert(pair<set<int>, int>(tmp, j));
    //             reverse_descent_map.insert(pair<int, set<int>>(j, tmp));
    //         }
    //     }
    //     return descent_map;
    // }

    set<int> get_descent_leaves(int node) const{
        set<int> leaves;
        if (node < get_num_leaves()){
            leaves.insert(node);
        }else{
            set<int> leaves1 = get_descent_leaves(nodes[node].child[0]);
            set<int> leaves2 = get_descent_leaves(nodes[node].child[1]);
            set_union(leaves1.begin(), leaves1.end(),
                        leaves2.begin(), leaves2.end(),
                        insert_iterator<set<int>>(leaves, leaves.begin()));
        }
        return leaves;
    }

    int find_mrca(shared_ptr<set<int>> leaves) const{

        //cout << "find mrca for: " << endl;
        //for(int leaf : *leaves){
        //    cout << leaf << " ";
        //}
        //cout << endl;

        // find mrca of leaves contained in the given set
        if (leaves->size() == get_num_leaves()){
            return root;
        }else{
            set<int> leaves_so_far;
            int curr = *(leaves->begin());
            leaves_so_far.insert(curr);
            while (!includes(leaves_so_far.begin(), leaves_so_far.end(), leaves->begin(), leaves->end())){
                int p = nodes[curr].parent;
                assert(p != -1);
                int *child = nodes[p].child;
                int other = (child[1] == curr ? child[0] : child[1]);
                set<int> leaves_to_insert = get_descent_leaves(other);
                leaves_so_far.insert(leaves_to_insert.begin(), leaves_to_insert.end());
                curr = p;
            }
            return curr;
        }
    }

    int find_mrca(set<int> *leaves) const{

        //cout << "find mrca for: " << endl;
        //for(int leaf : *leaves){
        //    cout << leaf << " ";
        //}
        //cout << endl;

        // find mrca of leaves contained in the given set
        if (leaves->size() == get_num_leaves()){
            return root;
        }else{
            set<int> leaves_so_far;
            int curr = *(leaves->begin());
            leaves_so_far.insert(curr);
            while (!includes(leaves_so_far.begin(), leaves_so_far.end(), leaves->begin(), leaves->end())){
                cout << "curr: " << curr << endl;
                int p = nodes[curr].parent;
                //cout << "p: " << p;
                assert(p != -1);
                int *child = nodes[p].child;
                int other = (child[1] == curr ? child[0] : child[1]);
                //cout << ", other: " << other << endl;
                set<int> leaves_to_insert = get_descent_leaves(other);
                leaves_so_far.insert(leaves_to_insert.begin(), leaves_to_insert.end());
                curr = p;
            }
            return curr;
        }
    }


    int nnodes;        // number of nodes in tree
    int capacity;      // capacity of nodes array
    int root;          // id of root node
    LocalNode *nodes;  // nodes array
};




// A tree within a set of local trees
//
// Specifically this structure describes the block over which the
// local exists and the SPR operation to the left of the local
// tree.
class LocalTreeSpr
{
public:
     LocalTreeSpr(LocalTree *tree, int *ispr, int blocklen, int *mapping=NULL) :
        tree(tree),
        spr(ispr[0], ispr[1], ispr[2], ispr[3]),
        mapping(mapping),
        blocklen(blocklen)
    {}

     LocalTreeSpr(LocalTree *tree, Spr spr, int blocklen, int *mapping=NULL) :
        tree(tree),
        spr(spr),
        mapping(mapping),
        blocklen(blocklen)
    {}

    // deallocate associated data
    void clear() {
        if (tree) {
            delete tree;
            tree = NULL;
        }

        if (mapping) {
            delete [] mapping;
            mapping = NULL;
        }
    }

    // set allocation capacity of underlying tree and node mapping
    void set_capacity(int _capacity)
    {
        if (tree->capacity == _capacity)
            return;

        tree->set_capacity(_capacity);

        // ensure capacity of mapping
        if (mapping) {
            int *tmp = new int [_capacity];
            assert(tmp);

            std::copy(mapping, mapping + _capacity, tmp);
            delete [] mapping;

            mapping = tmp;
        }
    }

    // ensure capacity of tree and mapping
    void ensure_capacity(int _capacity)
    {
        if (tree->capacity < _capacity)
            set_capacity(_capacity);
    }


    LocalTree *tree;  // local tree
    Spr spr;          // SPR operation to the left of local tree
    int *mapping;     // node mapping between previous tree and this tree
    int blocklen;     // length of sequence block
};



// A set of local trees that together specify an ARG
//
// This structure specifies a list of local trees, their blocks, and
// SPR operations.  Together this specifies an ancestral recombination
// graph (ARG), which because of our assumptions is an SMC-style ARG.
class LocalTrees
{
public:
    LocalTrees(int start_coord=0, int end_coord=0, int nnodes=0) :
        chrom("chr"),
        start_coord(start_coord),
        end_coord(end_coord),
        nnodes(nnodes) {}
    LocalTrees(int **ptrees, int**ages, int **isprs, int *blocklens,
               int ntrees, int nnodes, int capacity=-1, int start=0);
    ~LocalTrees()
    {
        clear();
    }

    // iterators for the local trees
    typedef list<LocalTreeSpr>::iterator iterator;
    typedef list<LocalTreeSpr>::reverse_iterator reverse_iterator;
    typedef list<LocalTreeSpr>::const_iterator const_iterator;
    typedef list<LocalTreeSpr>::const_reverse_iterator const_reverse_iterator;


    // Returns iterator for first local tree
    iterator begin() { return trees.begin(); }
    reverse_iterator rbegin() { return trees.rbegin(); }
    LocalTreeSpr &front() { return trees.front(); }

    const_iterator begin() const { return trees.begin(); }
    const_reverse_iterator rbegin() const { return trees.rbegin(); }
    const LocalTreeSpr &front() const { return trees.front(); }


    // Returns the ending iterator
    iterator end() { return trees.end(); }
    reverse_iterator rend() { return trees.rend(); }
    LocalTreeSpr &back() { return trees.back(); }

    const_iterator end() const { return trees.end(); }
    const_reverse_iterator rend() const { return trees.rend(); }
    const LocalTreeSpr &back() const { return trees.back(); }


    // Returns number of leaves
    inline int get_num_leaves() const
    {
        return (nnodes + 1) / 2;
    }

    // Returns sequence length
    inline int length() const
    {
        return end_coord - start_coord;
    }

    // Returns number of local trees
    inline int get_num_trees() const
    {
        return trees.size();
    }

    // Copy trees from another set of local trees
    void copy(const LocalTrees &other);

    // deallocate local trees
    void clear()
    {
        for (iterator it=begin(); it!=end(); it++)
            it->clear();
        trees.clear();
    }

    // make trunk genealogy
    void make_trunk(int start, int end, int seqid, int pop_path,
                    int capacity=-1)
    {
        clear();

        start_coord = start;
        end_coord = end;
        nnodes = 1;

        int ptree[] = {-1};
        int ages[] = {0};
        LocalTree *tree = new LocalTree(ptree, 1, ages, NULL, capacity);
        tree->nodes[0].pop_path = pop_path;
        trees.push_back(
        LocalTreeSpr(tree, Spr(-1, -1, -1, -1, -1), end - start, NULL));
        seqids.clear();
        seqids.push_back(seqid);
    }

    // set default sequence IDs
    void set_default_seqids()
    {
        const int nleaves = get_num_leaves();
        seqids.clear();
        for (int i=0; i<nleaves; i++)
            seqids.push_back(i);
    }


    // set sequence IDs according to a permutation of sequence names
    bool set_seqids(const vector<string> &names,
                    const vector<string> &new_order)
    {
        const int nnames = names.size();

        for (int i=0; i<nnames; i++) {
            seqids[i] = -1;
            for (unsigned int j=0; j<new_order.size(); j++) {
                if (names[i] == new_order[j]) {
                    seqids[i] = j;
                    break;
                }
            }
            // check for unmapped seqids
            if (seqids[i] == -1)
                return false;
        }

        return true;
    }


    // return local block containing site
    const_iterator get_block(int site, int &start, int &end) const
    {
        start = end = start_coord;
        for (const_iterator it = begin(); it != this->end(); ++it) {
            start = end;
            end += it->blocklen;
            if (start <= site && site < end)
                return it;
        }
        return this->end();
    }

    // return local block containing site
    const_iterator get_block(int site) const
    {
        int start, end;
        return get_block(site, start, end);
    }

    // return local block containing site
    iterator get_block(int site, int &start, int &end)
    {
        start = end = start_coord;
        for (iterator it = begin(); it != this->end(); ++it) {
            start = end;
            end += it->blocklen;
            if (start <= site && site < end)
                return it;
        }
        return this->end();
    }

    // return local block containing site
    iterator get_block(int site)
    {
        int start, end;
        return get_block(site, start, end);
    }



    string chrom;              // chromosome name of region
    int start_coord;           // start coordinate of whole tree list
                               // 0-based coordinate system
    int end_coord;             // end coordinate of whole tree list
    int nnodes;                // number of nodes in each tree
    list<LocalTreeSpr> trees;  // linked list of local trees

    vector<int> seqids;        // mapping from tree leaves to sequence ids
};


// count the lineages in a tree
void count_lineages(const LocalTree *tree, int ntimes,
                    int *nbranches, int *nrecombs,
                    int **nbranches_pop, int **ncoals_pop,
                    const PopulationTree *pop_tree);
void count_lineages_internal(const LocalTree *tree, int ntimes,
                             int *nbranches, int *nrecombs,
                             int **nbranches_pop, int **ncoals_pop,
                             const PopulationTree *pop_tree);
 void remove_population_paths(LocalTrees *trees);


// A structure that stores the number of lineages within each time segment
class LineageCounts
{
public:
   LineageCounts(int ntimes, int npops) :
    ntimes(ntimes), npops(npops)
    {
        nbranches = new int [ntimes];
        nrecombs = new int [ntimes];
	nbranches_pop = new int* [npops];
        ncoals_pop = new int* [npops];
	for (int i=0; i < npops; i++) {
	    nbranches_pop[i] = new int [2*ntimes];
	    ncoals_pop[i] = new int [ntimes];
	}
    }

    ~LineageCounts()
    {
	for (int i=0; i < npops; i++) {
	    delete [] nbranches_pop[i];
	    delete [] ncoals_pop[i];
	}
        delete [] nbranches;
        delete [] nbranches_pop;
        delete [] nrecombs;
        delete [] ncoals_pop;

    }

    // Counts the number of lineages for a tree
    inline void count(const LocalTree *tree, const PopulationTree *pop_tree,
                      bool internal=false) {
        if (internal)
            count_lineages_internal(tree, ntimes, nbranches, nrecombs,
                                    nbranches_pop, ncoals_pop, pop_tree);
        else
            count_lineages(tree, ntimes, nbranches, nrecombs,
                           nbranches_pop, ncoals_pop, pop_tree);
    }

    int ntimes;       // number of time points
    int npops;        // number of populations
    int *nbranches;  // number of branches per time slice
    int *nrecombs;   // number of recombination points per time slice
    int **ncoals_pop;     // number of coalescing points per time slice
    int **nbranches_pop;
};


//=============================================================================
// tree functions

void count_mig_events(int from_pop, int to_pop, int time_idx2,
                      const ArgModel *model,
                      const LocalTrees *trees,
                      const vector<Spr> *invisible_recombs,
                      int *count, int *total);

void apply_spr(LocalTree *tree, const Spr &spr,
               const PopulationTree *pop_tree=NULL);
               
void set_up_spr(Spr *spr, int coal_node, int recomb_node, int recomb_time_upper_bound,
                int recomb_time_lower_bound, int recoal_time, const double *times);
// return a new local tree that results from applying spr to prev_tree
LocalTree* apply_spr_new(LocalTree *prev_tree, const Spr &spr, int *mapping);

double get_treelen(const LocalTree *tree, const double *times, int ntimes,
                    bool use_basal=true);
double get_treelen_internal(const LocalTree *tree, const double *times,
                            int ntimes);
double get_treelen_branch(const LocalTree *tree, const double *times,
                          int ntimes, int node, int time, double treelen=-1.0,
                          bool use_basal=true);
double get_basal_branch(const LocalTree *tree, const double *times,
                        int ntimes, int node, int time);

//=============================================================================
// trees functions

double get_arglen(const LocalTrees *trees, const double *times);

void map_congruent_trees(const LocalTree *tree1, const int *seqids1,
                         const LocalTree *tree2, const int *seqids2,
                         int *mapping);
void infer_mapping(const LocalTree *tree1, const LocalTree *tree2,
                   int recomb_node, int *mapping);
void repair_spr(const LocalTree *last_tree, const LocalTree *tree, Spr &spr,
                int *mapping);

bool remove_null_spr(LocalTrees *trees, LocalTrees::iterator it,
                     const PopulationTree *pop_tree);
void remove_null_sprs(LocalTrees *trees, const PopulationTree *pop_tree);
void get_inverse_mapping(const int *mapping, int size, int *inv_mapping);

LocalTrees *partition_local_trees(LocalTrees *trees, int pos,
                                  LocalTrees::iterator it, int it_start,
                                  bool trim=true);
LocalTrees *partition_local_trees(LocalTrees *trees, int pos, bool trim=true);
void append_local_trees(LocalTrees *trees, LocalTrees *trees2, bool merge=true,
                        const PopulationTree *pop_tree=NULL);

void uncompress_local_trees(LocalTrees *trees,
                            const SitesMapping *sites_mapping);
void compress_local_trees(LocalTrees *trees, const SitesMapping *sites_mapping,
                          bool fuzzy=false);
void assert_uncompress_local_trees(LocalTrees *trees,
                                   const SitesMapping *sites_mapping);
int get_recoal_node(const LocalTree *tree, const Spr &spr, const int *mapping);


// Make a mapping for nodes between two local trees in SMC
// Assume every node maps to the same node name except the broken node
// (the parent of the recomb node).
inline void make_node_mapping(const int *ptree, int nnodes, int recomb_node,
                              int *mapping)
{
    for (int j=0; j<nnodes; j++)
        mapping[j] = j;

    // parent of recomb is broken and therefore does not map
    const int broken = ptree[recomb_node];
    mapping[broken] = -1;
}

void draw_local_tree(const LocalTree *tree, FILE *out, int depth=0);
void print_local_tree(const LocalTree *tree, FILE *out);
void print_local_trees(const LocalTrees *trees, FILE *out=stdout);

//=============================================================================
// input and output

void write_local_tree(const LocalTree *tree);
void write_local_trees_ts(const char *filename, const LocalTrees *trees, 
                            const Sequences *sequences, const SitesMapping *sites_mapping, 
                            const double *time);
void write_newick_tree(FILE *out, const LocalTree *tree,
                       const char *const *names,
                       const double *times, int depth, bool oneline,
                       bool pop_model=false);
bool write_newick_tree(const char *filename, const LocalTree *tree,
                       const char *const *names, const double *times,
                       bool oneline, bool pop_model=false);
void write_local_trees(FILE *out, const LocalTrees *trees,
                       const char *const *names, const double *times,
                       bool pop_model=false);
bool write_local_trees(const char *filename, const LocalTrees *trees,
                       const char *const *names, const double *times,
                       bool pop_model=false,
                       const vector<int> &self_recomb_pos=vector<int>(),
                       const vector<Spr> &self_recombs=vector<Spr>());
void write_local_trees(FILE *out, const LocalTrees *trees,
                       const Sequences &seqs, const double *times,
                       bool pop_model=false,
                       const vector<int> &self_recomb_pos=vector<int>(),
                       const vector<Spr> &self_recombs=vector<Spr>());
bool write_local_trees(const char *filename, const LocalTrees *trees,
                       const Sequences &seqs, const double *times,
                       bool pop_model=false,
                       const vector<int> &self_recomb_pos=vector<int>(),
                       const vector<Spr> &self_recombs=vector<Spr>());
bool parse_local_tree(const char* newick, LocalTree *tree,
                      const double *times, int ntimes);
bool read_local_trees(FILE *infile, const double *times, int ntimes,
                      LocalTrees *trees, vector<string> &seqnames,
                      vector<int> *invisible_recomb_pos=NULL,
                      vector<Spr> *invisible_recombs=NULL);
bool read_local_trees(const char *filename, const double *times, int ntimes,
                      LocalTrees *trees, vector<string> &seqnames);
bool read_local_trees_from_ts(const char *ts_fileName, const double *times, int ntimes,
                      LocalTrees *trees, vector<string> &seqnames, int start_coord, int end_coord);
void node_mapping(const LocalTree *lasttree, const LocalTree *localtree, int *mapping);
void display_localtree(const LocalTree *localtree);
bool read_local_trees_from_tsinfer(const char *ts_fileName, const double *times, int ntimes,
                      LocalTrees *trees, vector<string> &seqnames, int start_coord, int end_coord,
                      int maxIter = 10);
void write_newick_Tree_for_bedfile(FILE *out, const LocalTree *tree,
                                   const char *const *names,
                                   const ArgModel *model,
                                   const Spr &spr);
void write_local_trees_as_bed(FILE *out, const LocalTrees *trees,
                              const vector<string> seqnames,
                              const ArgModel *model, int sample);
string get_newick_rep_rSPR(const LocalTree *tree);
void get_newick_rep_rSPR_helper(string *s, const LocalTree *tree, int node);

//=============================================================================
// assert functions

bool assert_tree_postorder(const LocalTree *tree, const int *order);
bool assert_tree(const LocalTree *tree, const PopulationTree *pop_tree=NULL);
bool assert_spr(const LocalTree *last_tree, const LocalTree *tree,
                const Spr *spr, const int *mapping,
                const PopulationTree *pop_tree=NULL,
                bool pruned_internal=false);
 bool assert_trees(const LocalTrees *trees, const PopulationTree *pop_tree=NULL,
                   bool pruned_internal=false);




} // namespace argweaver

#endif // ARGWEAVER_LOCAL_TREES
