/*=============================================================================

  Matt Rasmussen
  Copyright 2007-2012

  Common sequence functions

=============================================================================*/


#ifndef ARGWEAVER_SEQUENCES_H
#define ARGWEAVER_SEQUENCES_H

// c++ includes
#include <string>
#include <vector>
#include <map>
#include <set>

// arghmm includes
#include "track.h"
#include "common.h"
#include "tabix.h"
#include "seq.h"

namespace argweaver {

using namespace std;
 class LocalTrees;
 class ArgModel;

class SitesMapping;

class BaseProbs
{
 public:
    BaseProbs() {}
    BaseProbs(char c) {
        int pos = dna2int[(int)c];
        if (pos < 0 || pos >= 4) {
            for (int i=0; i < 4; i++) prob[i]=1;
        } else {
            for (int i=0; i < 4; i++) prob[i]=0;
            prob[pos] = 1;
        }
    }
    BaseProbs(double prob0[4]) {
        for (int i=0; i < 4; i++) prob[i] = prob0[i];
    }

    // pl is normalized phred-like score for genotypes REF/REF,
    // REF/ALT, ALT/ALT. hap_id is 0 or 1
    BaseProbs(const char refAllele, const char altAllele,
              const string &pl, const string &pp, const int hap_id) {
        if (pl != "" && pl != ".") {

        } else if (pp != "" && pp != ".") {

        } else if (pl == "." || pp == ".") {
            for (int i=0; i < 4; i++) prob[i]=1;
            return;
        }
        // if pl and pp are not specified do nothing; probs should

    }
    BaseProbs(const BaseProbs &other) {
        for (int i=0; i < 4; i++) prob[i] = other.prob[i];
    }

    // this works for PL or GL genotype probabilies; only difference
    // in VCF 4.2 specification is that PL is integers
    void set_by_pl(const char refAllele, const char altAllele,
                   const string &pl, int hap_id) {
        assert(hap_id == 0 || hap_id == 1);
        for (int i=0; i < 4; i++) prob[i]=0.0;
        vector<string> pl_score_strings;
        double pl_scores[3];
        split(pl.c_str(), ",", pl_score_strings);
        if (pl_score_strings.size() != 3) {
            printError("Error parsing PL string %s\n", pl.c_str());
            assert(0);
        }
        double sum=0.0;
        for (int i=0; i < 3; i++) {
            pl_scores[i] = atof(pl_score_strings[i].c_str());
            pl_scores[i] = pow(10, -pl_scores[i]/10.0);
            sum += pl_scores[i];
        }
        int refNum = dna2int[(int)refAllele];
        assert(refNum >= 0);
        prob[refNum] += pl_scores[0] / sum;
        if (hap_id == 0)
            prob[refNum] += pl_scores[1] / sum;

        int altNum = dna2int[(int)altAllele];
        if (altNum >= 0) {
            prob[altNum] += pl_scores[2] / sum;
            if (hap_id == 1)
                prob[altNum] += pl_scores[1] / sum;
        } else {
            for (int i=0; i < 4; i++)
                prob[i] += pl_scores[2] / sum / 4.0;
            if (hap_id == 1) {
                for (int i=0; i < 4; i++)
                    prob[i] += pl_scores[1] / sum / 4.0;
            }
        }
    }

    void set_by_pp(const string &pp, int hap_id) {
        assert(hap_id == 0 || hap_id == 1);
        for (int i=0; i < 4; i++) prob[i]=0.0;
        vector<string> pp_score_strings;
        double pp_scores[10];
        split(pp.c_str(), ",", pp_score_strings);
        if (pp_score_strings.size() != 10) {
            printError("Error parsing PP string %s\n", pp.c_str());
            assert(0);
        }
        double sum=0.0;
        for (int i=0; i < 10; i++) {
            pp_scores[i] = atof(pp_score_strings[i].c_str());
            pp_scores[i] = pow(10, -pp_scores[i]/10.0);
            sum += pp_scores[i];
        }
        // pp_score is given for genotypes:
        // AA, CC, GG, TT, AC, AG, AT, CG, CT, GT, in this order
        // so take the first allele if hap_id = 0, or else the second
        const char *alleles = ( hap_id == 0 ? "ACGTAAACCG" : "ACGTCGTGTT" );
        for (int i=0; i < 10; i++)
            prob[dna2int[(int)alleles[i]]] += pp_scores[i] / sum;
    }

    double maxProb() {
        double rv = prob[0];
        for (int i=1; i < 4; i++)
            if (prob[i] > rv) rv = prob[i];
        return rv;
    }
    void set_mask() {
        for (int i=0; i < 4; i++)
            prob[i] = 1.0;
    }
    void set_certain(char c) {
        if (c=='N') {
            set_mask();
            return;
        }
        for (int i=0; i < 4; i++) prob[i] = 0.0;
        prob[dna2int[(int)c]] = 1.0;
    }
    bool is_masked() {
        for (int i=0; i < 4; i++)
            if (prob[i] < 0.99) return false;
        return true;
    }
    bool is_certain(double tol=1e-8) const {
        int numPossible=0;
        for (int i=0; i < 4; i++)
            if (prob[i] > tol) {
                numPossible++;
                if (numPossible > 1) return false;
            }
        if (numPossible == 0) {
            printError("is_certain: got no allele possible\n");
            assert(0);
        }
        return true;
    }
    bool is_equal(const BaseProbs &other, double tol=1e-8) const {
        for (int i=0; i < 4; i++)
            if (fabs(prob[i] - other.prob[i]) > tol)
                return false;
        return true;
    }
    double prob[4];
};


 // The alignment of sequences
class Sequences
{
public:
    explicit Sequences(int seqlen=0) :
        seqlen(seqlen), owned(false)
    {}

    Sequences(char **_seqs, int nseqs, int seqlen) :
        seqlen(seqlen), owned(false)
    {
        extend(_seqs, nseqs);
    }

    // initialize from a subset of another Sequences alignment
    Sequences(const Sequences *sequences, int nseqs=-1, int _seqlen=-1,
              int offset=0) :
        seqlen(_seqlen), owned(false)
    {
        // use same nseqs and/or seqlen by default
        if (nseqs == -1)
            nseqs = sequences->get_num_seqs();
        if (seqlen == -1)
            seqlen = sequences->length();

        for (int i=0; i<nseqs; i++)
            append(sequences->names[i], &sequences->seqs[i][offset],
                   vector<BaseProbs>());

        ref = vector<char>(sequences->ref);

	if (sequences->pairs.size() > 0) {
	  pairs = vector<int>(nseqs);
	  for (int i = 0; i < nseqs; i++)
	    pairs[i] = sequences->pairs[i];
	}
        if (sequences->pops.size() > 0) {
            pops = vector<int>(nseqs);
            for (int i=0; i < nseqs; i++)
                pops[i] = sequences->pops[i];
        }
	ages = sequences->ages;
        if (sequences->base_probs.size() > 0) {
            for (int i=0; i < (int)base_probs.size(); i++) {
                base_probs.push_back(sequences->base_probs[i]);
            }
        }
    }

    ~Sequences()
    {
        clear();
    }

    inline int get_num_seqs() const
    {
        return seqs.size();
    }

    inline int length() const
    {
        return seqlen;
    }

    inline void set_length(int _seqlen)
    {
        seqlen = _seqlen;
    }

    inline char **get_seqs()
    {
        return &seqs[0];
    }

    inline const char * const *get_seqs() const
    {
        return &seqs[0];
    }


    void set_owned(bool _owned)
    {
        owned = _owned;
    }

    void extend(char **_seqs, int nseqs)
    {
        for (int i=0; i<nseqs; i++) {
            seqs.push_back(_seqs[i]);
            names.push_back("");
            pops.push_back(0);
        }
    }

    void extend(char **_seqs, char **_names, int nseqs)
    {
        for (int i=0; i<nseqs; i++) {
            seqs.push_back(_seqs[i]);
            names.push_back(_names[i]);
            pops.push_back(0);
        }
    }

    bool append(string name, char *seq, vector<BaseProbs> bp, int pop=0)
    {
        // check sequence length
        int new_seqlen = strlen(seq);
        if (new_seqlen > 0) {
            if (seqs.size() > 0) {
                if (new_seqlen != seqlen)
                    return false;
            } else
                seqlen = new_seqlen;
        }
        if (seqs.size() > 0) {
            if (base_probs.size() > 0) assert(bp.size() > 0);
            else assert(bp.size() == 0);
        }
        seqs.push_back(seq);
        if (bp.size() > 0) base_probs.push_back(bp);
        names.push_back(name);
        pops.push_back(pop);
	if (pairs.size() > 0) pairs.push_back(-1);
        return true;
    }

    void clear()
    {
        if (owned) {
            const int nseqs = get_num_seqs();
            for (int i=0; i<nseqs; i++)
                delete [] seqs[i];
        }
        seqs.clear();
        ref.clear();
        alt.clear();
        names.clear();
        pops.clear();
        pairs.clear();
        non_singleton_snp.clear();
        base_probs.clear();
    }


    //set pairs vector assuming that diploids are named XXXX_1 and XXXX_2
    bool set_pairs_by_name();
    void set_pairs_from_file(string fn);
    void set_pairs(const ArgModel *mod);
    void set_pops(vector<char*> seqname, vector<int> pop) {
        for (unsigned int i=0; i < seqname.size(); i++) {
            bool found=false;
            for (unsigned int j=0; j < names.size(); j++) {
                if (names[j].compare(seqname[i]) == 0) {
                    pops[j] = pop[i];
                    found=true;
                    break;
                }
            }
            if (!found) {
                printError("assign_pops: Did not find sequence %s in"
                           " sequences; skipping\n", seqname[i]);
            }
        }
    }

    void set_pops_from_file(string fn);
    int get_pop(int seq) const {
        if (pops.size() == 0) return 0;
        if ((int)pops.size() <= seq || seq < 0) {
            printError("get_pop: bad sequence argument\n");
            assert(0);
        }
        return pops[seq];
    }

    int get_pop(char *seqname) const {
        if (pops.size() == 0) return 0;
        for (unsigned int i=0; i < names.size(); i++) {
            if (names[i].compare(seqname) == 0) {
                return pops[i];
            }
        }
        printError("get_pop: Did not find sequence %s\n", seqname);
        assert(0);
        return -1;
    }


    void set_age(int sample, double time, int ntimes, const double *times);
    void set_age(string agefile, int ntimes, const double *times);
    void set_age();

    int get_pair(int i) {
	if ((int)pairs.size() < i) return -1;
	return pairs[i];
    }

    void switch_alleles(int coord, int seq1, int seq2) {
      char tmp = seqs[seq1][coord];
      seqs[seq1][coord] = seqs[seq2][coord];
      seqs[seq2][coord] = tmp;
      if (base_probs.size() > 0) {
          BaseProbs tmp = base_probs[seq1][coord];
          base_probs[seq1][coord] = base_probs[seq2][coord];
          base_probs[seq2][coord] = tmp;
      }
    }

    void randomize_phase(double frac);
    void add_switch_errors(double rate);

    // set vector to true for each snp which has frequency > 1
    bool get_non_singleton_snp(vector<bool> &nonsing);

    TrackNullValue get_masked_regions(string chrom,
                                      const SitesMapping *sites_mapping = NULL)
        const;

    vector <char*> seqs;
    vector<char> ref;
    vector<char> alt;
    vector <string> names;
    vector <int> pops;
    vector <int> pairs; // index of diploid pair partner
    vector <bool> non_singleton_snp;  //true if snp w frequency > 1
    vector <int> ages; // set to non-zero for ancient samples
    vector <double> real_ages;
    vector<vector<BaseProbs> > base_probs;

protected:
    int seqlen;
    bool owned;
};


class PhaseProbs
{
 public:
 PhaseProbs(int hap1, int _treemap1, Sequences *_seqs,
	    const LocalTrees *trees, const ArgModel *model);
  ~PhaseProbs() {}

  void add(int coord, int state, double pr, int nstate) {
      coord += offset;
      map<int, vector<double> >::iterator it = probs.find(coord);
      if (it == probs.end()) {
	  vector<double> tmp(nstate);
	  tmp[state] = pr;
	  probs[coord] = tmp;
      } else {
	  it->second[state] = pr;
      }
   }

   unsigned int size() {
     return probs.size();
   }

   void updateTreeMap2(const LocalTrees *tree);

   void sample_phase(int *thread_path);

   //map from coordinate to a vector of doubles giving the probability
   // of current phasing for each state
  map<int,vector<double> > probs;
  int hap1, hap2;
  int treemap1, treemap2;
  int offset;
  Sequences *seqs;
  vector<bool> non_singleton_snp;
};


// sites are represented internally as 0-index and end-exclusive
// file-format represents sites as 1-index and end-inclusive
class Sites
{
public:
    Sites(string chrom="", int start_coord=0, int end_coord=0) :
        chrom(chrom),
        start_coord(start_coord),
        end_coord(end_coord)
    {}
    ~Sites()
    {
        clear();
    }

    void append(int position, char* col, bool copy=false) {
        positions.push_back(position);
        if (copy) {
            unsigned int len = strlen(col);
            assert(len == names.size());
            char *col2 = new char [len + 1];
            strcpy(col2, col);
            cols.push_back(col2);
        } else {
            cols.push_back(col);
        }
    }

    void append_masked(int position, bool have_base_probs) {
        positions.push_back(position);
        int n = get_num_seqs();
        char *col = new char[n+1];
        col[n] = '\0';
        for (int i=0; i < n; i++) col[i] = 'N';
        cols.push_back(col);
        if (have_base_probs) {
            vector<BaseProbs> bp;
            for (int i=0; i < n; i++)
                bp.push_back(BaseProbs('N'));
            base_probs.push_back(bp);
        }
    }

    void clear()
    {
        for (unsigned int i=0; i<cols.size(); i++)
            delete [] cols[i];
        names.clear();
        positions.clear();
        cols.clear();
        ref.clear();
        alt.clear();
        pops.clear();
        base_probs.clear();
    }

    inline int length() const
    {
        return end_coord - start_coord;
    }

    inline int get_num_sites() const
    {
        return positions.size();
    }

    inline int get_num_seqs() const
    {
        return names.size();
    }

    // note: does not consider base_probs so just returns most likely bases
    bool is_snp(int i) const {
        if (i < 0 || i >= (int)cols.size()) return false;
        int len = get_num_seqs();
        int j=0;
        for ( ; j < len; j++)
            if (cols[i][j] != 'N') break;
        if (j == len-1) return false;
        for (int k=j+1; k < len; k++)
            if (cols[i][k] != 'N' && cols[i][k] != cols[i][j])
                return true;
        return false;
    }

    bool rename(string rename_file);
    int remove_invariant();
    int subset(vector<int> keep);
    int subset(set<string> names_to_keep);

    template<class T>
    int remove_overlapping(const Track<T> &track);
    TrackNullValue get_masked_regions() const;
    TrackNullValue remove_masked();

    // merge sites in other into this object.
    // Assumes that both original objects contain all positions where there
    // is a deviation from reference (including a non-reference allele that is
    // fixed in sample). Both sites must also have some chrom, start_coord, and
    // end_coord. Is appropriate to use to combine sites read from different
    // VCF files.
    bool merge(const Sites &other);

    string chrom;
    int start_coord;
    int end_coord;
    vector<string> names;
    vector<int> pops;
    vector<int> positions;
    vector<char*> cols;
    vector<char> ref;
    vector<char> alt;
    vector<vector<BaseProbs> > base_probs;
};


// mapping of sites between a compressed and uncompressed alignment
class SitesMapping
{
public:
    SitesMapping() {}
    ~SitesMapping() {}

    void init(const Sites *sites)
    {
        old_start = sites->start_coord;
        old_end = sites->end_coord;
        nsites = sites->get_num_sites();
        seqlen = sites->length();
    }

    //round_dir < 0: round down (returns lower bound on coordinate)
    //round_dir > 0: round up (returns upper bound on coordinate)
    int compress(int pos, int round_dir=0, int start=0) const {
        const int n = all_sites.size();
        if (start < 0) start=0;
        for (int pos2 = start; pos2<n; pos2++) {
            if (all_sites_start[pos2] <= pos &&
                all_sites_end[pos2] >= pos) {
                if (round_dir == 0) return pos2;
                if (round_dir < 0) {
                    if (all_sites_end[pos2] == pos) return pos2;
                    return pos2-1;
                } else {
                    if (all_sites_start[pos2] == pos) return pos2;
                    return pos2+1;
                }
            }
        }
        assert(0);
        return n;
    }

    int uncompress(int pos) const {
        return all_sites[pos];
    }

    int uncompress_start(int pos) const {
        return all_sites_start[pos];
    }
    int uncompress_end(int pos) const {
        return all_sites_end[pos];
    }

    vector<int> uncompress(const vector<int> pos,
                           vector<int> &newpos) const {
        newpos.clear();
        for (unsigned int i=0; i < pos.size(); i++)
            newpos.push_back(all_sites[pos[i]]);
        return newpos;
    }

    // compress a series of block lengths
    void compress_blocks(const vector<int> &blocks, vector<int> &blocks2) const
    {
        int cur = new_start;
        int new_seqlen = new_end - new_start;

        int end = old_start;
        for (vector<int>::const_iterator it=blocks.begin();
             it != blocks.end(); ++it)
        {
            end += *it;

            if (end < old_end) {
                int cur2 = cur;
                for (; cur2 < new_seqlen && all_sites[cur2] < end; cur2++) {}

                blocks2.push_back(cur2 - cur);
                cur = cur2;
            } else {
                // last block case
                blocks2.push_back(new_end - cur);
            }
        }
    }

    // uncompress a series of block lengths
    const void uncompress_blocks(const vector<int> &blocks,
                           vector<int> &blocks2) const
    {
        int cur = old_start;
        int end = new_start;
        for (vector<int>::const_iterator it=blocks.begin();
             it != blocks.end(); ++it)
        {
            end += *it;

            if (end < new_end) {
                // use median for placing block ends
                int cur2 = (all_sites[end-1] + 1 + all_sites[end]) / 2;
                blocks2.push_back(cur2 - cur);
                assert(cur2 > cur);
                cur = cur2;
            } else {
                // last block case
                blocks2.push_back(old_end - cur);
            }
        }
    }


    int old_start;
    int old_end;
    int new_start;
    int new_end;
    int nsites;
    int seqlen;

    vector<int> old_sites; // the original position of each variant site
    vector<int> new_sites; // the new position of each variant site

    // all_sites has length = the length of compressed sequence
    // all_sites[i] gives original position of compressed sequence i
    // if i is variant, then this should be exact
    // if i is invariant it will be mid-point of compressed region
    // The original coordinates of compressed site i will be
    // from all_sites_start[i] to all_sites_end[i]. These coordinates
    // are 0-based but end-inclusive. Therefore if compression is 1
    // then all_sites_start=all_sites_end=0..(numsite-1)
    vector<int> all_sites; // the original position of each site
    vector<int> all_sites_start;
    vector<int> all_sites_end;
};


// sequences functions
bool read_fasta(FILE *infile, Sequences *seqs);
bool read_fasta(const char *filename, Sequences *seqs);
bool write_fasta(const char *filename, Sequences *seqs);
void write_fasta(FILE *stream, Sequences *seqs);

bool check_sequences(Sequences *seqs);
bool check_sequences(int nseqs, int seqlen, char **seqs);
bool check_seq_names(Sequences *seqs);
bool check_seq_name(const char *name);
void resample_align(Sequences *aln, Sequences *aln2);

// sites functions
void write_sites(FILE *stream, Sites *sites, bool write_masked=false);
bool read_sites(FILE *infile, Sites *sites,
                int subregion_start=-1, int subregion_end=-1, bool quiet=false);
bool read_sites(const char *filename, Sites *sites,
                int subregion_start=-1, int subregion_end=-1, bool quiet=false);

bool read_vcf(FILE *infile, Sites *sites, double min_qual,
              const char *genotype_filter,
              bool parse_genotype_probs, double min_base_prob,
              bool add_ref=false, const set<string> keep_inds=set<string>());
bool read_vcf(const char *filename, Sites *sites, const char *region,
              double min_qual, const char *genotype_filter,
              bool parse_genotype_probs, double min_base_prob, bool add_ref=false,
              const char *tabix_dir=NULL, const set<string> keep_inds=set<string>());
bool read_vcf(const string filename, Sites *sites, const string region,
              double min_qual, const string genotype_filter,
              bool parse_genotype_probs, double min_base_prob, bool add_ref=false,
              const string tabix_dir="", set<string> keep_inds=set<string>());
bool read_vcfs(const vector<string> filenames, Sites* sites, const string region,
               double min_qual, const string genotype_filter,
               bool parse_genotype_probs, double min_base_prob,
               const string tabixdir, set<string> keep_inds=set<string>());
void make_sequences_from_sites(const Sites *sites, Sequences *sequencess,
                               char default_char='A');
void make_sites_from_sequences(const Sequences *sequences, Sites *sites);

void apply_mask_sequences(Sequences *sequences, const TrackNullValue &maskmap,
                          const char *ind=NULL);
void apply_mask_sites(Sites *sites, const TrackNullValue &mask,
                      const char *ind, vector<TrackNullValue> *ind_maskmap);

void print_masked_regions(const Sequences &sequences,
                          const SitesMapping *sites_mapping,
                          string chrom,
                          string filename);
void print_masked_sites_regions(const Sites &sites, string filename);

// sequence compression
bool find_compress_cols(const Sites *sites, int compress,
                        SitesMapping *sites_mapping);
void compress_sites(Sites *sites, const SitesMapping *sites_mapping);
void uncompress_sites(Sites *sites, const SitesMapping *sites_mapping);

// return track containing regions where sites have >= numN ns
TrackNullValue get_n_regions(const Sites &sites, int numN);
TrackNullValue get_snp_clusters(const Sites &sites, int numsnp, int window);

 void unmask_inds(Sites *sites, vector<TrackNullValue> *masked_regions);


 template<class T>
void compress_track(Track<T> &track, const SitesMapping *sites_mapping,
                    double compress_seq, bool is_rate)
{
    Track<T> track2;

    if (sites_mapping) {
        // get block lengths
        vector<int> blocks;
        for (unsigned int i=0; i<track.size(); i++)
            blocks.push_back(track[i].length());

        // compress block lengths
        vector<int> blocks2;
        sites_mapping->compress_blocks(blocks, blocks2);

        // build compress track
        int start = sites_mapping->new_start;
        for (unsigned int i=0; i<track.size(); i++) {
            int end = start + blocks2[i];
            if (end > start)
                track2.append(track[i].chrom, start, end, track[i].value);
            start = end;
        }

        // replace track
        track.clear();
        track.insert(track.begin(), track2.begin(), track2.end());
    }


    // compress rate
    if (is_rate) {
        for (unsigned int i=0; i<track.size(); i++)
            track[i].value *= compress_seq;
    }
}

template<class T>
void uncompress_track(Track<T> &track, const SitesMapping *sites_mapping,
                    double compress_seq, bool is_rate)
{
    Track<T> track2;

    if (sites_mapping) {
        // get block lengths
        vector<int> blocks;
        for (unsigned int i=0; i<track.size(); i++)
            blocks.push_back(track[i].length());

        // compress block lengths
        vector<int> blocks2;
        sites_mapping->uncompress_blocks(blocks, blocks2);

        // build compress track
        int start = sites_mapping->old_start;
        for (unsigned int i=0; i<track.size(); i++) {
            int end = start + blocks2[i];
            if (end > start)
                track2.append(track[i].chrom, start, end, track[i].value);
            start = end;
        }

        // replace track
        track.clear();
        track.insert(track.begin(), track2.begin(), track2.end());
    }


    // compress rate
    if (is_rate) {
        for (unsigned int i=0; i<track.size(); i++)
            track[i].value /= compress_seq;
    }
}


template<class T>
void compress_mask(Track<T> &track, SitesMapping *sites_mapping,
                       bool expand_mask=false)
{
    track.merge();
    if (sites_mapping) {
        int prev_start_orig = 0, prev_start_new = 0;
        int round_dir1 = expand_mask ? -1 : 1;
        int round_dir2 = expand_mask ? 1 : -1;
        for (unsigned int i=0; i<track.size(); i++) {
            if (track[i].start < prev_start_orig)
                prev_start_new = 0;
            prev_start_orig =track[i].start-1;
            track[i].start = sites_mapping->compress(track[i].start,
                                                     round_dir1, prev_start_new);
            prev_start_new = track[i].start-1;
            track[i].end = sites_mapping->compress(track[i].end-1,
                                                    round_dir2, prev_start_new)+1;
        }
    }
}



} // namespace argweaver

#endif // ARGWEAVER_SEQUENCES_H
