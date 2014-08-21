
#include <algorithm>
#include <vector>
#include <map>

#include "est_popsize.h"
#include "total_prob.h"
#include "logging.h"
#include "model.h"

#ifdef ARGWEAVER_MPI
#include "mpi.h"
#endif

namespace argweaver {

using namespace std;


void resample_popsizes(ArgModel *model, const LocalTrees *trees,
                       double heat) {

#ifdef ARGWEAVER_MPI
    printLog(0, "resample_popsizes %i\t%i\n", MPI::COMM_WORLD.Get_rank(),
             model->mc3.group_comm.Get_rank()); fflush(stdout);
    //    MPI::COMM_WORLD.Barrier();
    MPI::Intracomm comm = model->mc3.group_comm;
    int rank = comm.Get_rank();
    printLog(0, "rank=%i\n", rank);
    if (rank == 0) {
#endif
        int num_accept=0, total=0;
        list<PopsizeConfigParam> &l = model->popsize_config.params;
        double *num_coal, *num_nocoal;
        vector<double> lrs(2*model->ntimes), trans(2*model->ntimes),
            prior(2*model->ntimes), oldn(2*model->ntimes),
            newn(2*model->ntimes), praccept(2*model->ntimes);
        vector<int> accepted(2*model->ntimes);
        num_coal = (double*)malloc(2*model->ntimes*sizeof(double));
        num_nocoal = (double*)malloc(2*model->ntimes * sizeof(double));
        double curr_like = calc_arg_prior(model, trees, num_coal, num_nocoal);
#ifdef ARGWEAVER_MPI
        comm.Reduce(MPI_IN_PLACE, &curr_like, 1, MPI::DOUBLE, MPI_SUM, 0);
        comm.Reduce(MPI_IN_PLACE, num_coal, 2*model->ntimes-1, MPI::DOUBLE,
                    MPI_SUM, 0);
        comm.Reduce(MPI_IN_PLACE, num_nocoal, 2*model->ntimes-1, MPI::DOUBLE,
                    MPI_SUM, 0);
#endif
        for (int rep=0; rep < model->popsize_config.numsample; rep++) {
            for (list<PopsizeConfigParam>::iterator it=l.begin();
                          it != l.end(); it++) {
                if (it->sample == false) continue;
                int maxpop=-1;
                for (set<int>::iterator it2=it->pops.begin();
                     it2 != it->pops.end(); it2++)
                    if ((*it2) > maxpop) maxpop = *it2;
                double old_popsize = model->popsizes[maxpop];
                double s = min(500.0, old_popsize/2.0);
                s *= s;  //variance of gamma proposal from
                         // old_popsize to new_popsize

                double new_popsize = rand_gamma(old_popsize*old_popsize/s,
                                                s/old_popsize);
#ifdef ARGWEAVER_MPI
                comm.Bcast(&new_popsize, 1, MPI::DOUBLE, 0);
#endif
                double sp = min(500.0, new_popsize/2.0);
                sp *= sp;  //variance of proposal from new_popsize
                           // to old_popsize
                double logn = log(old_popsize); //log N
                double lognp = log(new_popsize); //log N'
                double nsquare = old_popsize*old_popsize;
                double npsquare = new_popsize*new_popsize;
                double trans_ratio = (npsquare/sp - nsquare/s - 1.0)*logn
                    + (1.0 - nsquare/s + npsquare/sp)*lognp
                    - old_popsize*new_popsize/sp + old_popsize*new_popsize/s
                    - npsquare/sp*log(sp) + nsquare/s*log(s)
                    - lgamma(npsquare/sp) + lgamma(nsquare/s);

                // using an uninformative gamma prior which slowly goes to zero
                // as you move out to infinity, still allows N to be at least a
                // million (maybe this is too big?)
                // has mean of 200000 (k*theta) and sd of 200000, and is pretty
                // flat from 0 to 400k or so
                double prior_ratio;
                double prior_theta = 200000;
                // double prior_k=1.0;  //this is the value but it is
                // commented-out since never used
                if (( ! model->popsize_config.neighbor_prior) ||
                    maxpop >= 2*model->ntimes-2)
                    prior_ratio = (old_popsize-new_popsize)/prior_theta;
                else {
                    double prev_popsize = model->popsizes[maxpop+1];
                    double pneighbor=0.99999;
                    static double neighbor_sigma = 50.0;
                    static double neighbor_sigma22 = 2.0*50.0*50.0;
                    static double neighbor_scale = 1.0/(neighbor_sigma *
                                                        sqrt(2.0*3.141593));
                    double newprior = (1.0 - pneighbor)*
                        (exp(-new_popsize/prior_theta)/prior_theta) +
                        pneighbor * neighbor_scale*
                        exp(-(new_popsize - prev_popsize)*
                            (new_popsize-prev_popsize)/neighbor_sigma22);
                    double oldprior = (1.0 - pneighbor)*
                        (exp(-old_popsize/prior_theta)/prior_theta) +
                        pneighbor * neighbor_scale *
                        exp(-(old_popsize - prev_popsize)*
                            (old_popsize - prev_popsize)/neighbor_sigma22);
                    prior_ratio = log(newprior/oldprior);
                }


                for (set<int>::iterator it2 = it->pops.begin();
                     it2 != it->pops.end(); it2++)
                    model->popsizes[*it2] = new_popsize;
                double new_like = calc_arg_prior(model, trees);
#ifdef ARGWEAVER_MPI
                comm.Reduce(MPI_IN_PLACE, &new_like, 1, MPI::DOUBLE, MPI_SUM, 0);
#endif
                double lr = new_like - curr_like;
                double ln_accept = trans_ratio + prior_ratio + lr;
                ln_accept *= heat;
                double pr_accept = (ln_accept > 0 ? 1.0 : exp(ln_accept));
                bool accept = (ln_accept > 0 || frand() < pr_accept);
#ifdef ARGWEAVER_MPI
                comm.Bcast(&accept, 1, MPI::BOOL, 0);
#endif
                for (set<int>::iterator it2=it->pops.begin();
                     it2 != it->pops.end(); it2++) {
                    lrs[*it2] = new_like - curr_like;
                    trans[*it2] = trans_ratio;
                    prior[*it2] = prior_ratio;
                    oldn[*it2] = old_popsize;
                    newn[*it2] = new_popsize;
                    accepted[*it2] = (accept == true);
                    praccept[*it2] = pr_accept;
                }
                if (accept) {
                    num_accept++;
                    curr_like = new_like;
                } else {
                    for (set<int>::iterator it2 = it->pops.begin();
                         it2 != it->pops.end(); it2++) {
                        model->popsizes[*it2] = old_popsize;
                    }
                }
                total++;
            }
        }
        printLog(LOG_LOW, "done resample_popsizes num_accept=%i/%i\n",
                 num_accept, total);
        for (int i=0; i < 2*model->ntimes-1; i++) {
            int found=0;
            for (list<PopsizeConfigParam>::iterator it=l.begin();
                 it != l.end(); it++) {
                if (it->pops.find(i) != it->pops.end()) {
                    found=1;
                    if (it->sample) {
                        printLog(LOG_LOW,
                                 "%i\t%.1f\t%.1f\t%f\t%f\t%f\t%f\t%f\t%s\n",
                                 i, num_coal[i], num_nocoal[i], oldn[i], newn[i],
                                 lrs[i], trans[i], prior[i],
                                 accepted[i]==1 ? "accept" : "reject");
                    } else {
                        printLog(LOG_LOW,
                                 "%i\t%.1f\t%.1f\t%f\tnot_sampled\n",
                                 i, num_coal[i], num_nocoal[i],
                                 model->popsizes[i]);
                    }
                }
            }
        }
        fflush(stdout);
        free(num_coal);
        free(num_nocoal);
#ifdef ARGWEAVER_MPI
    } else {
        list<PopsizeConfigParam> l = model->popsize_config.params;
        double *num_coal = (double*)malloc(2*model->ntimes*sizeof(double));
        double *num_nocoal = (double*)malloc(2*model->ntimes * sizeof(double));
        double curr_like = calc_arg_prior(model, trees, num_coal, num_nocoal);
        comm.Reduce(&curr_like, &curr_like, 1, MPI::DOUBLE, MPI_SUM, 0);
        comm.Reduce(num_coal, num_coal, 2*model->ntimes-1, MPI::DOUBLE,
                    MPI_SUM, 0);
        comm.Reduce(num_nocoal, num_nocoal, 2*model->ntimes-1, MPI::DOUBLE,
                    MPI_SUM, 0);

        for (int rep=0; rep < model->popsize_config.numsample; rep++) {
            for (list<PopsizeConfigParam>::iterator it=l.begin();
                 it != l.end(); it++) {
                if (it->sample == false) continue;
                double old_popsize = model->popsizes[*(it->pops.begin())];
                double new_popsize;
                comm.Bcast(&new_popsize, 1, MPI::DOUBLE, 0);
                for (set<int>::iterator it2 = it->pops.begin();
                     it2 != it->pops.end(); it2++)
                    model->popsizes[*it2] = new_popsize;
                double new_like = calc_arg_prior(model, trees);
                comm.Reduce(&new_like, &new_like, 1, MPI::DOUBLE, MPI_SUM, 0);
                bool accept;
                comm.Bcast(&accept, 1, MPI::BOOL, 0);
                if (!accept) {
                    for (set<int>::iterator it2 = it->pops.begin();
                         it2 != it->pops.end(); it2++) {
                        model->popsizes[*it2] = old_popsize;
                    }
                }
            }
        }
        free(num_coal);
        free(num_nocoal);
    }
    printLog(0, "done resample popsizes\n");
#endif
}

void est_popsize_arg(const ArgModel *model, const LocalTrees *trees,
                     double *popsizes)
{

}


void est_popsize_trees2(const ArgModel *model, const LocalTree *const *trees,
                        int ntrees, double *popsizes)
{
    assert(ntrees > 0);

    const int ntimes = model->ntimes;
    const int nleaves = trees[0]->get_num_leaves();
    LineageCounts lineages(ntimes);

    int total_ncoals[ntimes];
    int total_pairs[ntimes];
    int total_ncoals_pairs[ntimes];

    fill(total_ncoals, total_ncoals + ntimes, 0);
    fill(total_pairs, total_pairs + ntimes, 0);
    fill(total_ncoals_pairs, total_ncoals_pairs + ntimes, 0);

    printf("ntrees %d\n", ntrees);

    // count lineages
    for (int i=0; i<ntrees; i++) {
        lineages.count(trees[i]);

        for (int j=0; j<ntimes-1; j++) {
            int start = (j == 0) ? nleaves : lineages.nbranches[j-1];
            int end = lineages.nbranches[j];
            int ncoals = start - end;
            int pairs = start * (start-1) / 2;

            total_ncoals[j] += ncoals;
            total_pairs[j] += pairs;
            total_ncoals_pairs[j] += ncoals * pairs;
        }
    }

    for (int j=0; j<ntimes-1; j++) {
        if (total_ncoals[j] == 0)
            popsizes[j] = 0.0;
        else
            popsizes[j] = .5 * model->time_steps[j] *
                (total_ncoals_pairs[j] + total_pairs[j] - total_ncoals[j]) /
                double(total_ncoals[j]);
        printf("> %d %d %d\n", total_ncoals_pairs[j], total_pairs[j],
               total_ncoals[j]);
        printf("popsize %f\n", popsizes[j]);
    }
}



void est_popsize_trees(const ArgModel *model, const LocalTree *const *trees,
                       int ntrees, double *popsizes)
{
    assert(ntrees > 0);

    const int ntimes = model->ntimes;
    const int nleaves = trees[0]->get_num_leaves();
    LineageCounts lineages(ntimes);

    int total_ncoals[ntimes];
    int total_pairs[ntimes];

    fill(total_ncoals, total_ncoals + ntimes, 0);
    fill(total_pairs, total_pairs + ntimes, 0);

    printf("ntrees %d\n", ntrees);

    // count lineages
    for (int i=0; i<ntrees; i++) {
        lineages.count(trees[i]);

        for (int j=0; j<ntimes-1; j++) {
            int start = (j == 0) ? nleaves : lineages.nbranches[j-1];
            int end = lineages.nbranches[j];
            int ncoals = start - end;
            int pairs = start * (start-1) / 2;

            total_ncoals[j] += ncoals;
            total_pairs[j] += pairs;
        }
    }

    for (int j=0; j<ntimes-1; j++) {
        if (total_ncoals[j] == 0)
            popsizes[j] = 0.0;
        else
            popsizes[j] = .5 * model->time_steps[j] * total_pairs[j] /
                double(total_ncoals[j]);
        printf("> %d %d\n", total_pairs[j], total_ncoals[j]);
        printf("popsize %f\n", popsizes[j]);
    }
}


void est_popsize_trees(const ArgModel *model, const LocalTrees *trees,
                       int step, double *popsizes)
{
    vector<LocalTree*> indep_trees;

    int end = trees->start_coord;
    int pos = end;
    for (LocalTrees::const_iterator it=trees->begin();
         it != trees->end(); ++it)
        {
            int start = end;
            end += it->blocklen;

            while (start <= pos && pos < end) {
                // record tree
                indep_trees.push_back(it->tree);
                pos += step;
            }
        }

    est_popsize_trees(model, &indep_trees[0], indep_trees.size(), popsizes);
}


//=============================================================================
// C-interface

extern "C" {

// estimate population sizes
void arghmm_est_popsizes_trees(LocalTrees *trees, double *times, int ntimes,
                               int step, double *popsizes)
{
    // setup model, local trees, sequences
    ArgModel model(ntimes, times, NULL, 0.0, 0.0);
    est_popsize_trees(&model, trees, step, popsizes);
}


} // extern C

} // namespace argweaver

