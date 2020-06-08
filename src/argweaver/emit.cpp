
#include "common.h"
#include "emit.h"
#include "seq.h"
#include "thread.h"

namespace argweaver {


//=============================================================================
// invariant sites
// TODO: make this mask aware


// Returns true if position 'pos' is invariant
static inline bool is_invariant_site(const char *const *seqs,
                                     const int nseqs, const int pos)
{
    const char c = seqs[0][pos];
    for (int j=1; j<nseqs; j++) {
        if (seqs[j][pos] != c) {
            return false;
        }
    }
    return true;
}

// Returns true if position 'pos' is invariant
static inline bool is_invariant_site(const char *const *seqs,
                                     const int nseqs, const int pos,
                                     const vector<vector<BaseProbs> > &base_probs)
{
    if (base_probs.size() == 0)
        return is_invariant_site(seqs, nseqs, pos);

    const char c = seqs[0][pos];
    if (c != 'N' && !base_probs[0][pos].is_certain())
        return false;
    for (int j=0; j<nseqs; j++) {
        if (seqs[j][pos] != c || ! base_probs[j][pos].is_certain() ) {
            return false;
        }

    }
    return true;
}



// Populates array 'variant' with false for invariant sites and true otherwise
void find_variant_sites(const char *const *seqs, int nseqs, int seqlen,
                        bool *variant,
                        const vector<vector<BaseProbs> > &base_probs)
{
    bool have_base_probs = ( base_probs.size() > 0 );
    // find invariant sites
    for (int i=0; i<seqlen; i++) {
        char c = seqs[0][i];
        if (c != 'N' && have_base_probs && ! base_probs[0][i].is_certain())
            variant[i] = true;
        else {
            variant[i] = false;
            for (int j=1; j<nseqs; j++) {
                if (seqs[j][i] != c ||
                    (have_base_probs && ! base_probs[j][i].is_certain())) {
                    variant[i] = true;
                    break;
                }
            }
        }
    }
}


void find_masked_sites(const char *const *seqs, int nseqs, int seqlen,
                       bool *masked, bool *variant)
{
    if (variant) {
        for (int i=0; i<seqlen; i++)
            masked[i] = (seqs[0][i] == 'N' && !variant[i]);
    } else {
        for (int i=0; i<seqlen; i++)
            masked[i] = (seqs[0][i] == 'N' && is_invariant_site(seqs, nseqs, i));
    }
}


int count_alleles(const char *const *seqs,
                  const int nseqs, const int pos)
{
    int counts[4] = {0, 0, 0, 0};

    for (int j=0; j<nseqs; j++) {
        int c = dna2int[(int) seqs[j][pos]];
        if (c != -1)
            counts[c]++;
    }
    int alleles = int(counts[0] > 0) +
        int(counts[1] > 0) +
        int(counts[2] > 0) +
        int(counts[3] > 0);

    return alleles;
}


//=============================================================================
// calculate mutation probabilities


// Juke-Cantor
// probability of mutation over time 't' with mutation rate 'mu'
static inline double prob_branch(double t, double mu, bool mut)
{
    const double f = 4. / 3.;
    if (!mut)
        return .25 * (1.0 + 3. * exp(-f*mu*t));
    else
        return .25 * (1.0 - exp(-f*mu*t));
}


// get mutation probabilities
void prob_tree_mutation(const LocalTree *tree, const ArgModel *model,
                        double *muts, double *nomuts)
{
    const double *times = model->times;
    const int nnodes = tree->nnodes;
    const LocalNode *nodes = tree->nodes;

    for (int i=0; i<nnodes; i++) {
        if (i == tree->root)
            continue;
        int parent_age = nodes[nodes[i].parent].age;
        if (parent_age == model->get_removed_root_time())
            continue;

        double t = ( parent_age == nodes[i].age ? model->get_mintime(parent_age)
                     : times[parent_age] - times[nodes[i].age] );
        muts[i] = prob_branch(t, model->mu, true);
        nomuts[i] = prob_branch(t, model->mu, false);
    }
}


//============================================================================
// emissions

// table of partial likelihood values
typedef double lk_row[4];

class LikelihoodTable
{
public:
    LikelihoodTable(int seqlen, int nnodes) :
        seqlen(seqlen),
        nnodes(nnodes)
    {
        data = new lk_row* [seqlen];
        for (int i=0; i<seqlen; i++)
            data[i] = new double [nnodes][4];
    }

    ~LikelihoodTable()
    {
        for (int i=0; i<seqlen; i++)
            delete [] data[i];
        delete [] data;
    }

    int seqlen;
    int nnodes;
    lk_row **data;
};


// calculate inner partial likelihood for one node and site
inline void likelihood_site_node_inner(
    const LocalTree *tree, const int node,
    const char *const *seqs,
    const vector<vector<BaseProbs> > &base_probs,
    const int pos,
    const double *muts, const double *nomuts, lk_row* inner)
{
    const LocalNode* nodes = tree->nodes;
    const int j = node;

    if (nodes[j].is_leaf()) {
        // leaf case
        const char c = seqs[j][pos];
        if (c == 'N') {
            inner[j][0] = 1.0;
            inner[j][1] = 1.0;
            inner[j][2] = 1.0;
            inner[j][3] = 1.0;
        } else if (base_probs.size() > 0) {
            for (int k=0; k < 4; k++)
                inner[j][k] = base_probs[j][pos].prob[k];
        } else {
            inner[j][0] = 0.0;
            inner[j][1] = 0.0;
            inner[j][2] = 0.0;
            inner[j][3] = 0.0;
            inner[j][dna2int[(int) c]] = 1.0;
        }
    } else {
        // internal node case
        int c1 = nodes[j].child[0];
        int c2 = nodes[j].child[1];

        for (int a=0; a<4; a++) {
            double p1 = 0.0;
            double p2 = 0.0;

            for (int b=0; b<4; b++) {
                if (a == b) {
                    p1 += inner[c1][b] * nomuts[c1];
                    p2 += inner[c2][b] * nomuts[c2];
                } else {
                    p1 += inner[c1][b] * muts[c1];
                    p2 += inner[c2][b] * muts[c2];
                }
            }

            inner[j][a] = p1 * p2;
        }
    }
}


// calculate outer partial likelihood for one node and site
inline void likelihood_site_node_outer(
    const LocalTree *tree, const int root, const int node,
    const double *muts, const double *nomuts,
    lk_row* outer, lk_row* inner)
{
    const LocalNode* nodes = tree->nodes;
    const int j = node;

    if (node == root) {
        // root case
        outer[j][0] = 1.0;
        outer[j][1] = 1.0;
        outer[j][2] = 1.0;
        outer[j][3] = 1.0;
    } else {
        // non-root case
        int sib = tree->get_sibling(j);
        int parent = nodes[j].parent;

        if (parent != root) {
            for (int a=0; a<4; a++) {
                double p1 = 0.0;
                double p2 = 0.0;
                for (int b=0; b<4; b++) {
                    if (a == b) {
                        p1 += inner[sib][b] * nomuts[sib];
                        p2 += outer[parent][b] * nomuts[parent];
                    } else {
                        p1 += inner[sib][b] * muts[sib];
                        p2 += outer[parent][b] * muts[parent];
                    }
                }
                outer[j][a] = p1 * p2;
            }
        } else {
            for (int a=0; a<4; a++) {
                double p1 = 0.0;
                for (int b=0; b<4; b++) {
                    if (a == b)
                        p1 += inner[sib][b] * nomuts[sib];
                    else
                        p1 += inner[sib][b] * muts[sib];
                }
                outer[j][a] = p1;
            }
        }
    }
}


// calculate entire inner partial likelihood table
double likelihood_site_inner(
    const LocalTree *tree, const char *const *seqs,
    const vector<vector<BaseProbs> > &base_probs,
    const int pos, const int *order, const int norder,
    const double *muts, const double *nomuts, lk_row* inner)
{
    // iterate postorder through nodes
    for (int i=0; i<norder; i++)
        likelihood_site_node_inner(
            tree, order[i], seqs, base_probs, pos, muts, nomuts, inner);

    // sum over root node
    double p = 0.0;
    int root = tree->root;
    for (int a=0; a<4; a++)
        p += inner[root][a]  * .25;

    return p;
}


// calculate eniter outer partial likelihood table
void likelihood_site_outer(
    const LocalTree *tree,
    const double *muts, const double *nomuts, bool internal,
    lk_row *inner, lk_row *outer)
{
    int queue[tree->nnodes];
    int top = 0;

    // process in preorder
    int maintree_root = internal ? tree->nodes[tree->root].child[1] :
        tree->root;
    queue[top++] = maintree_root;
    while (top > 0) {
        int node = queue[--top];
        likelihood_site_node_outer(tree, maintree_root, node,
                                   muts, nomuts, outer, inner);

        // recurse
        if (!tree->nodes[node].is_leaf()) {
            queue[top++] = tree->nodes[node].child[0];
            queue[top++] = tree->nodes[node].child[1];
        }
    }
}


void calc_inner_outer(const LocalTree *tree, const ArgModel *model,
                      const char *const *seqs,
                      const vector<vector<BaseProbs> > &base_probs,
                      const int seqlen,
                      const bool *variant, bool internal,
                      lk_row **inner, lk_row **outer)
{
    // get postorder
    int norder = tree->nnodes;
    int order[tree->nnodes];
    tree->get_postorder(order);


    // get mutation probabilities and treelen
    double muts[tree->nnodes];
    double nomuts[tree->nnodes];
    prob_tree_mutation(tree, model, muts, nomuts);

    // calculate emissions for tree at each site
    for (int i=0; i<seqlen; i++) {
        if (variant[i]) {
            likelihood_site_inner(tree, seqs, base_probs, i, order, norder,
                                  muts, nomuts, inner[i]);
            likelihood_site_outer(tree, muts, nomuts, internal,
                                  inner[i], outer[i]);
        }
    }
}



void likelihood_sites(const LocalTree *tree, const ArgModel *model,
                      const char *const *seqs,
                      const vector<vector<BaseProbs> > &base_probs,
                      const int seqlen, const int statei,
                      const bool *variant,
                      double **emit, lk_row **table,
                      const int prev_node=-1, const int new_node=-1)
{
    const double *times = model->times;
    const LocalNode *nodes = tree->nodes;
    const double mintime = model->get_mintime();

    // get postorder
    int order[tree->nnodes];
    int norder;

    if (prev_node == -1) {
        tree->get_postorder(order);
        norder = tree->nnodes;
    } else {
        // find partial postorder

        // find dirty entries
        bool dirty[tree->nnodes];
        fill(dirty, dirty+tree->nnodes, false);
        for (int j=prev_node; j!=-1; j=nodes[j].parent)
            dirty[j] = true;

        // walk up root path
        norder = 0;
        for (int j=new_node; !dirty[j]; j=nodes[j].parent)
            order[norder++] = j;
        for (int j=prev_node; j!=-1; j=nodes[j].parent)
            order[norder++] = j;
    }


    // get mutation probabilities and treelen
    double muts[tree->nnodes];
    double nomuts[tree->nnodes];
    double treelen = 0.0;
    for (int i=0; i<tree->nnodes; i++) {
        if (i != tree->root) {
            double t = ( nodes[nodes[i].parent].age == nodes[i].age ?
                         model->get_mintime(nodes[i].age) :
                         times[nodes[nodes[i].parent].age] - times[nodes[i].age]);
            muts[i] = prob_branch(t, model->mu, true);
            nomuts[i] = prob_branch(t, model->mu, false);
            treelen += t;
        }
    }


    // calculate invariant_lk
    double invariant_lk = .25 * exp(- model->mu * max(treelen, mintime));

    // calculate emissions for tree at each site
    for (int i=0; i<seqlen; i++) {
        if (variant && !variant[i]) {
            // use precommuted invariant site likelihood
            if (seqs[0][i] == 'N')
                emit[i][statei] = 1.0; // masked case
            else
                emit[i][statei] = invariant_lk;
        } else {
            emit[i][statei] = likelihood_site_inner(
                tree, seqs, base_probs, i, order, norder, muts, nomuts, table[i]);
        }
    }
}



double likelihood_tree(const LocalTree *tree, const ArgModel *model,
                       const char *const *seqs,
                       const vector<vector<BaseProbs> > &base_probs,
                       const int nseqs,
                       const int start, const int end)
{
    const double *times = model->times;
    const int nnodes = tree->nnodes;
    const LocalNode *nodes = tree->nodes;
    double invariant_lk = -1;
    lk_row table[nnodes];

    // get postorder
    int order[tree->nnodes];
    tree->get_postorder(order);

    // get mutation probabilities
    double muts[tree->nnodes];
    double nomuts[tree->nnodes];
    for (int i=0; i<tree->nnodes; i++) {
        if (i != tree->root) {
            double t = ( nodes[nodes[i].parent].age == nodes[i].age ?
                  model->get_mintime(nodes[i].age) :
                  times[nodes[nodes[i].parent].age] - times[nodes[i].age] );
            muts[i] = prob_branch(t, model->mu, true);
            nomuts[i] = prob_branch(t, model->mu, false);
        }
    }


    // calculate emissions for tree at each site
    double lnl = 0.0;
    for (int i=start; i<end; i++) {
        double lk;
        bool invariant = is_invariant_site(seqs, nseqs, i, base_probs);
        if (invariant && seqs[0][i] == 'N')
            continue;

        if (invariant && invariant_lk > 0)
            // use precommuted invariant site likelihood
            lk = invariant_lk;
        else {
            lk = likelihood_site_inner(tree, seqs, base_probs, i, order,
                                       tree->nnodes, muts, nomuts, table);

            // save invariant likelihood
            if (invariant)
                invariant_lk = lk;
        }

        lnl += log(lk);
    }

    return lnl;
}



//=============================================================================
// emission calculation


void get_infinite_sites_states(const States &states, const LocalTree *tree,
                               const char *const *seqs, int nseqs, int seqlen,
                               bool *variant,
                               bool internal, bool **valid_states)
{
    const int nstates = states.size();
    const int nnodes = tree->nnodes;

    if (internal) {
        // internal branch case
        const int maintree_root = tree->nodes[tree->root].child[1];
        const int subtree_root = tree->nodes[tree->root].child[0];

        // get nodes within trees
        int nsubnodes = 0, subnodes[nnodes];
        tree->get_preorder(subtree_root, subnodes, nsubnodes);
        reverse(subnodes, subnodes + nsubnodes);
        int nmainnodes = 0, mainnodes[nnodes];
        tree->get_preorder(maintree_root, mainnodes, nmainnodes);
        reverse(mainnodes, mainnodes + nmainnodes);

        for (int i=0; i<seqlen; i++) {
            if (!variant[i])
                continue;

            // get set of all bases in the subtree
            char subset = 0;
            for (int k=0; k<nsubnodes; k++) {
                int j = subnodes[k];
                if (tree->nodes[j].is_leaf() && seqs[j][i] != 'N')
                    subset |= 1 << dna2int[(int) seqs[j][i]];
            }

            // get set of all bases in the maintree
            char mainset = 0;
            for (int k=0; k<nmainnodes; k++) {
                int j = mainnodes[k];
                if (tree->nodes[j].is_leaf() && seqs[j][i] != 'N')
                    mainset |= 1 << dna2int[(int) seqs[j][i]];
            }

            // if subtree and main have distinct bases then all states are valid
            if (!(subset & mainset)) {
                for (int j=0; j<nstates; j++)
                    valid_states[i][j] = true;
                continue;
            }


            // infer ancestral reconstruction
            char bases[nnodes];
            parsimony_ancestral_set(
                tree, seqs, i, subnodes, nsubnodes, bases);
            int cset = bases[subtree_root];
            parsimony_ancestral_set(
                tree, seqs, i, mainnodes, nmainnodes, bases);

            bool valid_nodes[nnodes];
            for (int k=0; k<nmainnodes; k++) {
                int j = mainnodes[k];
                int parent = tree->nodes[j].parent;
                valid_nodes[j] = bool(cset & bases[j]) ||
                    (j != maintree_root && (cset & bases[parent]));
            }
            for (int j=0; j<nstates; j++)
                valid_states[i][j] = valid_nodes[states[j].node];

            // check for at least one valid state
            bool valid = false;
            for (int j=0; j<nstates; j++) {
                if (valid_states[i][j]) {
                    valid = true;
                    break;
                }
            }
            if (!valid) {
                printLog(LOG_LOW, "unable to satisfy infinite sites assumption\n");
            }
        }

    } else {
        // external branch  case
        int newleaf = tree->get_num_leaves();
        int postorder[nnodes];
        tree->get_postorder(postorder);

        for (int i=0; i<seqlen; i++) {
            if (!variant[i])
                continue;

            // get set of all bases in the tree
            char set = 0;
            for (int j=0; j<newleaf; j++)
                if (seqs[j][i] != 'N')
                    set |= 1 << dna2int[(int) seqs[j][i]];

            char c = seqs[newleaf][i];
            char cset = ((c != 'N') ? 1 << dna2int[(int) c] : 0);
            if (!(cset & set)) {
                // newleaf has a new base, thus newleaf can go anywhere
                for (int j=0; j<nstates; j++)
                    valid_states[i][j] = true;
                continue;
            }

            // infer ancestral reconstruction
            char bases[nnodes];
            parsimony_ancestral_set(
                tree, seqs, i, postorder, nnodes, bases);

            bool valid_nodes[nnodes];
            for (int j=0; j<nnodes; j++) {
                int parent = tree->nodes[j].parent;
                valid_nodes[j] = bool(cset & bases[j]) ||
                    (parent != -1 && (cset & bases[parent]));
            }
            for (int j=0; j<nstates; j++)
                valid_states[i][j] = valid_nodes[states[j].node];

            // check for at least one valid state
            bool valid = false;
            for (int j=0; j<nstates; j++) {
                if (valid_states[i][j]) {
                    valid = true;
                    break;
                }
            }
            if (!valid) {
                printLog(LOG_LOW, "unable to satisfy infinite sites assumption");
            }
        }
    }
}


void get_infinite_sites_states(const States &states, const LocalTree *tree,
                               const char *const *seqs, int nseqs, int seqlen,
                               bool *variant,
                               bool internal, bool **valid_states,
                               PhaseProbs *phase_pr) {
    get_infinite_sites_states(states, tree, seqs, nseqs, seqlen, variant,
                              internal, valid_states);
    if (phase_pr != NULL &&
        phase_pr->treemap1 >= 0    && phase_pr->treemap2 >= 0 &&
        phase_pr->treemap1 < nseqs && phase_pr->treemap2 < nseqs) {
        int nstates = states.size();
        bool **valid_states2 = new_matrix<bool>(seqlen, nstates);
        const char *flipSeqs[nseqs];
        for (int i=0; i < nseqs; i++)
            flipSeqs[i] = seqs[i];
        flipSeqs[phase_pr->treemap1] = seqs[phase_pr->treemap2];
        flipSeqs[phase_pr->treemap2] = seqs[phase_pr->treemap1];
        get_infinite_sites_states(states, tree, flipSeqs, nseqs, seqlen,
                                  variant, internal, valid_states2);

        for (int i=0; i < seqlen; i++)
            for (int j=0; j < nstates; j++)
                valid_states[i][j] = (valid_states[i][j] || valid_states2[i][j]);
    }
}



double calc_emit(lk_row *in, lk_row *out, lk_row *in2,
		 int i, int node1, int node2, int maintree_root,
		 double *nomut, double *mut) {
    double emit=0.0;
    for (int a=0; a<4; a++) {
	double p1 = 0.0, p2 = 0.0, p3 = 0.0;
	for (int b=0; b<4; b++) {
	    if (a == b) {
		p1 += in2[node1][b] * nomut[0];
		p2 += in[node2][b] * nomut[1];
		p3 += out[node2][b] * nomut[2];
	    } else {
		p1 += in2[node1][b] * mut[0];
		p2 += in[node2][b] * mut[1];
		p3 += out[node2][b] * mut[2];
	    }
	}

	if (node2 != maintree_root) {
	    emit += p1 * p2 * p3 * .25;
	} else {
	    emit += p1 * p2 * .25;
	}
    }
    return emit;
}



// calculate emissions for external branch resampling
void calc_emissions(const States &states, const LocalTree *tree,
                    const char *const *seqs,
                    const vector<vector<BaseProbs> > &base_probs,
                    int nseqs, int seqlen,
                    const ArgModel *model, bool internal, double **emit,
		    PhaseProbs *phase_pr)
{
    const int nstates = states.size();
    const int newleaf = tree->get_num_leaves();
    const int maintree_root = internal ? tree->nodes[tree->root].child[1] :
        tree->root;
    const int subtree_root = internal ? tree->nodes[tree->root].child[0] :
        tree->root;


    // special case: ignore fully specified local tree
    if (internal && nstates == 0) {
        for (int i=0; i<seqlen; i++)
            emit[i][0] = 1.0;
        return;
    }


    // find invariant sites
    bool *variant = new bool [seqlen];
    bool *masked = new bool [seqlen];
    find_variant_sites(seqs, nseqs, seqlen, variant, base_probs);
    find_masked_sites(seqs, nseqs, seqlen, masked, variant);


    // compute inner and outer likelihood tables
    LikelihoodTable inner(seqlen, tree->nnodes);
    LikelihoodTable inner_subtree(seqlen, 1);
    LikelihoodTable outer(seqlen, tree->nnodes);
    calc_inner_outer(tree, model, seqs, base_probs, seqlen, variant, internal,
                     inner.data, outer.data);

    if (!internal) {
        // compute inner table for new leaf
        for (int i=0; i<seqlen; i++) {
            const char c = seqs[newleaf][i];
            if (c == 'N') {
                inner_subtree.data[i][0][0] = 1.0;
                inner_subtree.data[i][0][1] = 1.0;
                inner_subtree.data[i][0][2] = 1.0;
                inner_subtree.data[i][0][3] = 1.0;
            } else if (base_probs.size() > 0) {
                for (int j=0; j < 4; j++)
                    inner_subtree.data[i][0][j] = base_probs[newleaf][i].prob[j];
            } else {
                inner_subtree.data[i][0][0] = 0.0;
                inner_subtree.data[i][0][1] = 0.0;
                inner_subtree.data[i][0][2] = 0.0;
                inner_subtree.data[i][0][3] = 0.0;
                inner_subtree.data[i][0][dna2int[(int) c]] = 1.0;
            }
        }
    }

    LikelihoodTable inner2(seqlen, tree->nnodes);
    LikelihoodTable inner_subtree2(seqlen, 1);
    LikelihoodTable outer2(seqlen, tree->nnodes);
    bool *het = NULL;
    if (model->unphased && phase_pr != NULL &&
	phase_pr->treemap1 >= 0 && phase_pr->treemap1 < nseqs &&
	phase_pr->treemap2 >= 0 && phase_pr->treemap2 < nseqs) {
	const char *subseqs[nseqs];
	for (int i=0; i < nseqs; i++)
	    subseqs[i] = seqs[i];
	subseqs[phase_pr->treemap1] = seqs[phase_pr->treemap2];
	subseqs[phase_pr->treemap2] = seqs[phase_pr->treemap1];
	het = new bool[seqlen];
	for (int i=0; i < seqlen; i++) {
	    het[i] = (seqs[phase_pr->treemap1][i] != seqs[phase_pr->treemap2][i]);
            if (base_probs.size() > 0 && !het[i])
                het[i] = ! (base_probs[phase_pr->treemap1][i].is_equal(
                            base_probs[phase_pr->treemap2][i]));
        }
        vector<vector<BaseProbs> > base_probs2;
        if (base_probs.size() > 0) {
            for (int i=0; i < nseqs; i++) {
                if (i == phase_pr->treemap1)
                    base_probs2.push_back(base_probs[phase_pr->treemap2]);
                else if (i == phase_pr->treemap2)
                    base_probs2.push_back(base_probs[phase_pr->treemap1]);
                else base_probs2.push_back(base_probs[i]);
            }
        }

	calc_inner_outer(tree, model, subseqs, base_probs2, seqlen, het,
			 internal, inner2.data, outer2.data);

	if (!internal) {
	    // compute inner table for new leaf
	    for (int i=0; i<seqlen; i++) {
		const char c = subseqs[newleaf][i];
		if (c == 'N') {
		    inner_subtree2.data[i][0][0] = 1.0;
		    inner_subtree2.data[i][0][1] = 1.0;
		    inner_subtree2.data[i][0][2] = 1.0;
		    inner_subtree2.data[i][0][3] = 1.0;
                } else if (base_probs2.size() > 0) {
                    for (int j=0; j < 4; j++)
                        inner_subtree2.data[i][0][j] = base_probs2[newleaf][i].prob[j];
		} else {
		    inner_subtree2.data[i][0][0] = 0.0;
		    inner_subtree2.data[i][0][1] = 0.0;
		    inner_subtree2.data[i][0][2] = 0.0;
		    inner_subtree2.data[i][0][3] = 0.0;
		    inner_subtree2.data[i][0][dna2int[(int) c]] = 1.0;
		}
	    }
	}
    }

    // calc tree lengths
    int queue[tree->nnodes];
    int top = 0;

    double maintreelen = 0.0;
    // process in preorder maintree nodes
    queue[top++] = maintree_root;
    while (top > 0) {
        int node = queue[--top];
        if (node != maintree_root)
            maintreelen += max(tree->get_dist(node, model->times),
                               model->get_mintime(tree->nodes[node].age));
        if (!tree->nodes[node].is_leaf()) {
            queue[top++] = tree->nodes[node].child[0];
            queue[top++] = tree->nodes[node].child[1];
        }
    }

    double subtreelen = 0.0;
    if (internal) {
        // process in preorder subtree nodes
        queue[top++] = subtree_root;
        while (top > 0) {
            int node = queue[--top];
            if (node != subtree_root)
                subtreelen += max(tree->get_dist(node, model->times),
                                  model->get_mintime(tree->nodes[node].age));
            if (!tree->nodes[node].is_leaf()) {
                queue[top++] = tree->nodes[node].child[0];
                queue[top++] = tree->nodes[node].child[1];
            }
        }
    }


    // populate emission table
    for (int j=0; j<nstates; j++) {
        State state = states[j];

        // get nodes
        int node1 = internal ? subtree_root : 0;
        int node2 = state.node;
        int parent = tree->nodes[node2].parent;

        // get times
        double time1 = internal ? model->times[tree->nodes[node1].age] : 0.0;
        double time2 = model->times[tree->nodes[node2].age];
        double parent_time = (parent != -1) ?
            model->times[min(tree->nodes[parent].age, model->ntimes-1)] : 0.0;
        double coal_time = model->times[state.time];

        // get distances
	double dist[3], mut[3], nomut[3];
        double curr_mintime = model->get_mintime(state.time);
        dist[0] = max(coal_time - time1, curr_mintime);
        dist[1] = max(coal_time - time2, curr_mintime);
        dist[2] = max(parent_time - coal_time, curr_mintime);

        // get mutation probabilities
        mut[0] = prob_branch(dist[0], model->mu, true);
        mut[1] = prob_branch(dist[1], model->mu, true);
        mut[2] = prob_branch(dist[2], model->mu, true);
        nomut[0] = prob_branch(dist[0], model->mu, false);
        nomut[1] = prob_branch(dist[1], model->mu, false);
        nomut[2] = prob_branch(dist[2], model->mu, false);

        // get tree length
        double treelen;
        if (node2 == maintree_root)
            treelen = maintreelen + subtreelen
                + max(coal_time - time1, curr_mintime)
                + max(coal_time - model->times[tree->nodes[maintree_root].age],
                      curr_mintime);
        else
            treelen = maintreelen + subtreelen
                + max(coal_time - time1, curr_mintime);

        // calculate invariant_lk
        double invariant_lk = .25 * exp(- model->mu * treelen);

        // fill in row of emission table
        for (int i=0; i<seqlen; i++) {
            if (masked[i]) {
                // masked site
                emit[i][j] = 1.0;
            } else if (!variant[i]) {
                // invariant site
                emit[i][j] = invariant_lk;
            } else {
		emit[i][j] = calc_emit(inner.data[i], outer.data[i],
				       internal ? inner.data[i] : inner_subtree.data[i],
				       i, node1, node2, maintree_root,
				       nomut, mut);
                assert(!isnan(emit[i][j]));
		if (het != NULL && het[i]) {
		    double emit2 = calc_emit(inner2.data[i], outer2.data[i],
					    internal ? inner2.data[i] : inner_subtree2.data[i],
					    i, node1, node2, maintree_root,
					    nomut, mut);
		    phase_pr->add(i, j, emit[i][j]/(emit[i][j] + emit2), nstates);
		    emit[i][j] += emit2;
		    emit[i][j] *= 0.5;
                    assert(!isnan(emit[i][j]));
		}
            }
        }
    }

    // optionally enforce infinite sites model
    if (model->infsites_penalty < 1.0) {

        bool **valid_states = new_matrix<bool>(seqlen, nstates);
        get_infinite_sites_states(states, tree, seqs, nseqs, seqlen,
                                  variant, internal, valid_states,
                                  model->unphased ? phase_pr : NULL);
        for (int i=0; i<seqlen; i++) {
            if (variant[i]) {
                for (int j=0; j<nstates; j++)
                    if (!valid_states[i][j])
                        emit[i][j] *= model->infsites_penalty;
            }
        }
        delete_matrix<bool>(valid_states, seqlen);
    }


    // clean up
    delete [] variant;
    delete [] masked;
    if (het != NULL) delete [] het;
}

// calculate emissions for external branch resampling
void calc_emissions_external(const States &states, const LocalTree *tree,
                             const char *const *seqs,
                             const vector<vector<BaseProbs> > &base_probs,
                             int nseqs, int seqlen,
                             const ArgModel *model, double **emit,
			     PhaseProbs *phase_pr)
{
    calc_emissions(states, tree, seqs, base_probs, nseqs, seqlen, model, false,
                   emit, phase_pr);
}

// calculate emissions for internal branch resampling
void calc_emissions_internal(const States &states, const LocalTree *tree,
                             const char *const *seqs,
                             const vector<vector<BaseProbs> > &base_probs,
                             int nseqs, int seqlen,
                             const ArgModel *model, double **emit,
                             PhaseProbs *phase_pr)
{
    calc_emissions(states, tree, seqs, base_probs, nseqs, seqlen, model, true,
		   emit, phase_pr);
}




//=============================================================================
// counting non-compatiable sites


void parsimony_ancestral_set(const LocalTree *tree, const char * const *seqs,
                             int pos, int *postorder, int npostorder,
                             char *ancestral)
{
    const int nnodes = tree->nnodes;
    const LocalNode *nodes = tree->nodes;
    char sets[nnodes];

    // clear sets
    for (int node=0; node<nnodes; node++)
        sets[node] = 0;

    // do unweighted parsimony by postorder traversal
    int postorder2[nnodes];
    if (!postorder) {
        tree->get_postorder(postorder2);
        postorder = postorder2;
        npostorder = nnodes;
    }
    for (int i=0; i<npostorder; i++) {
        int node = postorder[i];
        if (nodes[node].is_leaf()) {
            char c = seqs[node][pos];
            if (c == 'N')
                sets[node] = 0;
            else
                sets[node] = 1 << dna2int[(int) c];
        } else {
            char lset = sets[nodes[node].child[0]];
            char rset = sets[nodes[node].child[1]];
            char intersect = lset & rset;
            if (intersect > 0)
                sets[node] = intersect;
            else
                sets[node] = lset | rset;
        }
    }

    // traceback
    int root = postorder[npostorder-1];
    ancestral[root] = sets[root];

    // traceback with preorder traversal
    for (int i=npostorder-2; i>=0; i--) {
        int node = postorder[i];
        char s = sets[node];
        char pset = ancestral[nodes[node].parent] & s;
        if (pset) {
            // use parent char if possible
            ancestral[node] = pset;
        } else {
            // otherwise do not refine set
            ancestral[node] = s;
        }
    }
}


void parsimony_ancestral_seq(const LocalTree *tree, const char * const *seqs,
                             int nseqs, int pos, char *ancestral,
                             int *postorder)
{
    const int nnodes = tree->nnodes;
    const LocalNode *nodes = tree->nodes;
    char sets[nnodes];
    int pchar;

    // clear sets
    for (int node=0; node<nnodes; node++)
        sets[node] = 0;

    // do unweighted parsimony by postorder traversal
    int postorder2[nnodes];
    if (!postorder) {
        tree->get_postorder(postorder2);
        postorder = postorder2;
    }
    for (int i=0; i<nnodes; i++) {
        int node = postorder[i];
        if (nodes[node].is_leaf()) {
            char c = seqs[node][pos];
            if (c == 'N')
                sets[node] = 1 + 2 + 4 + 8;
            else
                sets[node] = 1 << dna2int[(int) c];
        } else {
            char lset = sets[nodes[node].child[0]];
            char rset = sets[nodes[node].child[1]];
            char intersect = lset & rset;
            if (intersect > 0)
                sets[node] = intersect;
            else
                sets[node] = lset | rset;
        }
    }

    // traceback
    // arbitrary choose root base from set
    int root = postorder[nnodes-1];
    char rootset = sets[root];
    ancestral[root] = (rootset & 1) ? int2dna[0] :
        (rootset & 2) ? int2dna[1] :
        (rootset & 4) ? int2dna[2] : int2dna[3];

    // traceback with preorder traversal
    for (int i=nnodes-2; i>=0; i--) {
        int node = postorder[i];
        char s = sets[node];

        switch (s) {
        case 1: // just A
            ancestral[node] = int2dna[0];
            break;
        case 2: // just C
            ancestral[node] = int2dna[1];
            break;
        case 4: // just G
            ancestral[node] = int2dna[2];
            break;
        case 8: // just T
            ancestral[node] = int2dna[3];
            break;
        default:
            pchar = ancestral[nodes[node].parent];
            if (dna2int[pchar] & s) {
                // use parent char if possible
                ancestral[node] = pchar;
            } else {
                // use arbitrary char otherwise
                ancestral[node] = (s & 1) ? int2dna[0] :
                    (s & 2) ? int2dna[1] :
                    (s & 4) ? int2dna[2] : int2dna[3];
            }
        }
    }
}


int parsimony_cost_seq(const LocalTree *tree, const char * const *seqs,
                        int nseqs, int pos, int *postorder)
{
    const int nnodes = tree->nnodes;
    const LocalNode *nodes = tree->nodes;
    const int maxcost = 100000;
    int costs[nnodes][4];

    // do unweighted parsimony by postorder traversal
    int postorder2[nnodes];
    if (!postorder) {
        tree->get_postorder(postorder2);
        postorder = postorder2;
    }
    for (int i=0; i<nnodes; i++) {
        int node = postorder[i];
        if (tree->nodes[node].is_leaf()) {
            char c = seqs[node][pos];
            if (c == 'N') {
                for (int a=0; a<4; a++)
                    costs[node][a] = 0;
            } else {
                for (int a=0; a<4; a++)
                    costs[node][a] = maxcost;
                costs[node][dna2int[(int)c]] = 0;
            }
        } else {
            int *left_costs = costs[nodes[node].child[0]];
            int *right_costs = costs[nodes[node].child[1]];

            for (int a=0; a<4; a++) {
                int left_min = maxcost;
                int right_min = maxcost;
                for (int b=0; b<4; b++) {
                    left_min = min(left_min, int(a != b) + left_costs[b]);
                    right_min = min(right_min, int(a != b) + right_costs[b]);
                }
                costs[node][a] = left_min + right_min;
            }
        }
    }

    int root_min = maxcost;
    for (int a=0; a<4; a++)
        root_min = min(root_min, costs[tree->root][a]);

    return root_min;
}


int count_noncompat(const LocalTree *tree, const char * const *seqs,
                    int nseqs, int block_start, int block_len, int *postorder)
{
    // get postorder
    int postorder2[tree->nnodes];
    if (!postorder) {
        tree->get_postorder(postorder2);
        postorder = postorder2;
    }

    int noncompat = 0;
    for (int i=block_start; i<block_len; i++)
        if (!is_invariant_site(seqs, nseqs, i)) {
            int a = count_alleles(seqs, nseqs, i);
            int c = parsimony_cost_seq(tree, seqs, nseqs, i, postorder);
            noncompat += int(c > a - 1 );
        }

    return noncompat;
}




int count_noncompat(const LocalTrees *trees, const char * const *seqs,
                    int nseqs, int seqlen,
                    int start_coord, int end_coord)
{
    int noncompat = 0;
    if (start_coord == -1) start_coord = trees->start_coord;
    if (end_coord == -1) end_coord = trees->end_coord;

    int end = trees->start_coord;
    for (LocalTrees::const_iterator it=trees->begin();
         it != trees->end(); ++it)
    {
        int start = end;
        end += it->blocklen;

        //for my debugging purposes
        printLog(LOG_LOW, "start: %d, end: %d\n", start, end);

        if (end <= start_coord) continue;
        if (start >= end_coord) break;
        LocalTree *tree = it->tree;
        int blocklen = it->blocklen;
        int block_start = 0;
        if (start_coord > start) block_start = start_coord - start;
        int block_end = blocklen;
        if (end_coord < end) block_end = end_coord - start;

        // get subsequence block
        char const* subseqs[nseqs];
        for (int i=0; i<nseqs; i++)
            subseqs[i] = &seqs[i][start];

        noncompat += count_noncompat(tree, subseqs, nseqs, block_start,
                                     block_end, NULL);

    }

    return noncompat;
}


int count_noncompat(const LocalTrees *trees, const Sequences *sequences,
                        int start_coord, int end_coord) {
    int nseqs = trees->get_num_leaves();
    char *seqs[nseqs];
    for (int i=0; i < nseqs; i++)
        seqs[i] = sequences->seqs[trees->seqids[i]];
    return count_noncompat(trees, seqs, nseqs, sequences->length(),
                           start_coord, end_coord);
}




//=============================================================================
// slow literal emission calculation
// useful for testing against


void calc_emissions_external_slow(
    const States &states, const LocalTree *tree,
    const char *const *seqs, const vector<vector<BaseProbs> > &base_probs,
    int nseqs, int seqlen,
    const ArgModel *model, double **emit)
{
    const int nstates = states.size();
    const int newleaf = tree->get_num_leaves();
    bool *variant = new bool [seqlen];
    LikelihoodTable table(seqlen, tree->nnodes+2);

    // create local tree we can edit
    LocalTree tree2(tree->nnodes, tree->nnodes + 2);
    tree2.copy(*tree);

    // find invariant sites
    find_variant_sites(seqs, nseqs, seqlen, variant, base_probs);

    for (int j=0; j<nstates; j++) {
        State state = states[j];
        add_tree_branch(&tree2, state.node, state.time);

        likelihood_sites(&tree2, model, seqs, base_probs,
                         seqlen, j, variant, emit,
                         table.data, -1, -1);

        remove_tree_branch(&tree2, newleaf, model, NULL);
    }

    delete [] variant;
}


void calc_emissions_internal_slow(
    const States &states, const LocalTree *tree,
    const char *const *seqs, const vector<vector<BaseProbs> > &base_probs,
    int nseqs, int seqlen, const ArgModel *model, double **emit)
{
    const int nstates = states.size();
    const int subtree_root = tree->nodes[tree->root].child[0];
    const int subtree_root_age = tree->nodes[subtree_root].age;
    const int maxtime = model->ntimes + 1;

    // ignore fully specified local tree
    if (nstates == 0) {
        for (int i=0; i<seqlen; i++)
            emit[i][0] = 1.0;
        return;
    }

    bool *variant = new bool [seqlen];
    LikelihoodTable table(seqlen, tree->nnodes+2);

    // create local tree we can edit
    LocalTree tree2(tree->nnodes, tree->nnodes + 2);
    tree2.copy(*tree);

    // find variant sites
    find_variant_sites(seqs, nseqs, seqlen, variant, base_probs);

    for (int j=0; j<nstates; j++) {
        State state = states[j];
        assert(subtree_root != tree2.root);

        Spr add_spr(subtree_root, subtree_root_age, state.node, state.time);
        //note: pop_tree unnecessary for computing emissions
        apply_spr(&tree2, add_spr);

        likelihood_sites(&tree2, model, seqs, base_probs, seqlen, j, variant,
                         emit, table.data, -1, -1);

        Spr remove_spr(subtree_root, subtree_root_age,
                       tree2.root, maxtime);
        apply_spr(&tree2, remove_spr);
    }

    // clean up
    delete [] variant;
}



//=============================================================================
// assert emissions

bool assert_emissions(const States &states, const LocalTree *tree,
                      const char *const *seqs,
                      const vector<vector<BaseProbs> > &base_probs,
                      int nseqs, int seqlen,
                      const ArgModel *model)
{
    const int nstates = states.size();

    double **emit = new_matrix<double>(seqlen, nstates);
    double **emit2 = new_matrix<double>(seqlen, nstates);

    calc_emissions_external(states, tree, seqs, base_probs, nseqs, seqlen,
                            model, emit, NULL);
    calc_emissions_external_slow(states, tree, seqs, base_probs, nseqs, seqlen,
                                 model, emit2);

    // compare emission tables
    for (int j=0; j<nstates; j++) {
        for (int i=0; i<seqlen; i++) {
            bool invar = is_invariant_site(seqs, nseqs, i, base_probs);

            printLog(LOG_MEDIUM, ">> %d,%d: %e %e   invar=%d\n",
                     i, j, emit[i][j], emit2[i][j], invar);
            if (!fequal(emit[i][j], emit2[i][j], .0001, 1e-12))
                return false;
        }
    }

    delete_matrix<double>(emit, seqlen);
    delete_matrix<double>(emit2, seqlen);

    return true;
}


bool assert_emissions_internal(const States &states, const LocalTree *tree,
                               const char *const *seqs,
                               const vector<vector<BaseProbs> > &base_probs,
                               int nseqs, int seqlen,
                               const ArgModel *model)
{
    const int nstates = states.size();

    double **emit = new_matrix<double>(seqlen, nstates);
    double **emit2 = new_matrix<double>(seqlen, nstates);

    calc_emissions_internal(states, tree, seqs, base_probs, nseqs, seqlen,
                              model, emit);
    calc_emissions_internal_slow(states, tree, seqs, base_probs, nseqs, seqlen,
                                   model, emit2);

    // compare emission tables
    for (int j=0; j<nstates; j++) {
        for (int i=0; i<seqlen; i++) {
            printLog(LOG_MEDIUM, ">> %d,%d: %e %e\n",
                     i, j, emit[i][j], emit2[i][j]);
            if (!fequal(emit[i][j], emit2[i][j], .0001, 1e-12))
                return false;
        }
    }

    delete_matrix<double>(emit, seqlen);
    delete_matrix<double>(emit2, seqlen);

    return true;
}


//=============================================================================
// C interface
extern "C" {

double **new_emissions(intstate *istates, int nstates,
                       int *ptree, int nnodes, int *ages_index,
                       char **seqs, int nseqs, int seqlen,
                       double *times, int ntimes,
                       double mu)
{
    States states;
    make_states(istates, nstates, states);
    LocalTree tree(ptree, nnodes, ages_index);
    ArgModel model(ntimes, times, NULL, 0.0, mu);

    double **emit = new_matrix<double>(seqlen, nstates);
    vector<vector<BaseProbs> > base_probs;
    base_probs.clear();
    calc_emissions_external(states, &tree, seqs, base_probs,
                            nseqs, seqlen, &model, emit, NULL);

    return emit;
}


void delete_emissions(double **emit, int seqlen)
{
    delete_matrix<double>(emit, seqlen);
}


bool arghmm_assert_emit(
    LocalTrees *trees, int ntimes, double *times, double mu,
    char **seqs, int nseqs, int seqlen)
{
    ArgModel model(ntimes, times, NULL, 0.0, mu);
    States states;
    vector<vector<BaseProbs> > base_probs;
    base_probs.clear();

    int end = trees->start_coord;
    for (LocalTrees::const_iterator it=trees->begin(); it!=trees->end(); ++it) {
        int start = end;
        int blocklen = it->blocklen;
        end = start + blocklen;
        LocalTree *tree = it->tree;

        get_coal_states(tree, model.ntimes, states);

        // get subsequence
        char *seqs2[nseqs];
        for (int j=0; j<nseqs; j++)
            seqs2[j] = &seqs[j][start];

        if (!assert_emissions(states, tree, seqs2, base_probs,
                              nseqs, blocklen, &model))
            return false;
    }

    return true;
}

bool arghmm_assert_emit_internal(
    LocalTrees *trees, int ntimes, double *times, double mu,
    char **seqs, int nseqs, int seqlen)
{
    ArgModel model(ntimes, times, NULL, 0.0, mu);
    const int maxtime = model.ntimes + 1;
    States states;
    int *removal_path = new int [trees->get_num_trees()];
    vector<vector<BaseProbs> > base_probs;
    base_probs.clear();

    // randomly choose branch to remove
    LocalTrees trees2(*trees);

    // ramdomly choose a removal path
    sample_arg_removal_path_uniform(&trees2, removal_path);
    remove_arg_thread_path(&trees2, removal_path, maxtime);

    delete [] removal_path;

    int end = trees->start_coord;
    for (LocalTrees::const_iterator it=trees2.begin(); it!=trees2.end(); ++it)
    {
        int start = end;
        int blocklen = it->blocklen;
        end = start + blocklen;
        LocalTree *tree = it->tree;

        get_coal_states_internal(tree, model.ntimes, states);

        // get subsequence
        char *seqs2[nseqs];
        for (int j=0; j<nseqs; j++)
            seqs2[j] = &seqs[j][start];

        if (!assert_emissions_internal(states, tree, seqs2, base_probs,
                                       nseqs, blocklen, &model))
            return false;
    }

    return true;
}



} // extern "C"



} // namespace argweaver


