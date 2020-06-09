#include "common.h"
#include "logging.h"
#include "parsing.h"
#include "seq.h"
#include "sequences.h"
#include "local_tree.h"
#include "model.h"

// TODO: add sites validation
//       - positions should be sorted and unique
//       - bases should be acceptable characters


namespace argweaver {


//=============================================================================
// input/output: FASTA

bool read_fasta(FILE *infile, Sequences *seqs)
{
    // store lines until they are ready to discard
    class Discard : public vector<char*> {
    public:
        ~Discard() { clean(); }
        void clean() {
            for (unsigned int i=0; i<size(); i++)
                delete [] at(i);
            clear();
        }
    };


    // init sequences
    seqs->clear();
    seqs->set_owned(true);

    char *line;
    string key;
    vector<char*> seq;
    Discard discard;

    while ((line = fgetline(infile))) {
        chomp(line);

        if (line[0] == '>') {
            // parse key line

            if (seq.size() > 0) {
                // add new sequence
                char *full_seq = concat_strs(&seq[0], seq.size());
                if (!seqs->append(key, full_seq, vector<BaseProbs>())) {
                    printError("sequences are not the same length: %d != %d",
                               seqs->length(), strlen(full_seq));
                    delete [] full_seq;
                    return false;
                }
                seq.clear();
                discard.clean();
            }

            // new key found
            key = string(trim(&line[1]));
            delete [] line;
        } else {
            // parse sequence line

            seq.push_back(trim(line));
            discard.push_back(line);
        }
    }

    // add last sequence
    if (seq.size() > 0) {
        char *full_seq = concat_strs(&seq[0], seq.size());
        if (!seqs->append(key, full_seq, vector<BaseProbs>())) {
            printError("sequences are not the same length: %d != %d",
                       seqs->length(), strlen(full_seq));
            delete [] full_seq;
            return false;
        }
        discard.clean();
    }

    return true;
}


bool read_fasta(const char *filename, Sequences *seqs)
{
    FILE *infile = NULL;
    if ((infile = fopen(filename, "r")) == NULL) {
        printError("cannot read file '%s'", filename);
        return false;
    }

    bool result = read_fasta(infile, seqs);
    fclose(infile);

    return result;
}



bool write_fasta(const char *filename, Sequences *seqs)
{
    FILE *stream = NULL;

    if ((stream = fopen(filename, "w")) == NULL) {
        fprintf(stderr, "cannot open '%s'\n", filename);
        return false;
    }

    write_fasta(stream, seqs);

    fclose(stream);
    return true;
}

void write_fasta(FILE *stream, Sequences *seqs)
{
    for (int i=0; i<seqs->get_num_seqs(); i++) {
        fprintf(stream, ">%s\n", seqs->names[i].c_str());
        fprintf(stream, "%s\n", seqs->seqs[i]);
    }
}



//=============================================================================
// input/output: sites file format

void write_sites(FILE *stream, Sites *sites, bool write_masked) {
    bool have_base_probs = sites->base_probs.size() > 0;
    int nseq = (int)sites->names.size();
    fprintf(stream, "#NAMES");
    for (int i=0; i < nseq; i++)
        fprintf(stream, "\t%s", sites->names[i].c_str());
    fprintf(stream, "\n");
    fprintf(stream, "#REGION\t%s\t%i\t%i\n",
            sites->chrom.c_str(), sites->start_coord + 1, sites->end_coord);
    if (sites->positions.size() != sites->cols.size()) {
        fprintf(stderr, "Error in write_sites: positions.size()=%i cols.size=%i\n",
                (int)sites->positions.size(), (int)sites->cols.size());
        exit(-1);
    }
    for (unsigned int i=0; i < sites->positions.size(); i++) {
        int j;
        for (j=1; j < nseq; j++)
            if (sites->cols[i][j] != sites->cols[i][0]) break;
        if (j == nseq) {
            // this site is invariant though it may be masked
            if (!write_masked) continue;
            if (sites->cols[i][0] != 'N') continue;
        }
        fprintf(stream, "%i\t", sites->positions[i]+1);
        for (unsigned int k=0; k < sites->names.size(); k++)
            fprintf(stream, "%c", sites->cols[i][k]);
        if (have_base_probs)
            for (int k=0; k < (int)sites->names.size(); k++)
                for (int l=0; l<4; l++)
                    fprintf(stream, "\t%f", sites->base_probs[i][k].prob[l]);
        fprintf(stream, "%c", '\n');
    }
}

bool validate_site_column(char *col, int nseqs)
{
    for (int i=0; i<nseqs; i++) {
        col[i] = toupper(col[i]);
        if (col[i] != 'N' && dna2int[(int) col[i]] == -1)
            return false;
    }
    return true;
}


// Read a Sites stream
bool read_sites(FILE *infile, Sites *sites,
                int subregion_start, int subregion_end, bool quiet)
{
    const char *delim = "\t";
    char *line;
    int nseqs = 0;

    sites->clear();
    bool error = false;
    bool have_base_probs=false;

    // parse lines
    int lineno = 0;
    bool isHeader=false;
    char *line0;

    while (!error && (line = fgetline(infile))) {
        chomp(line);
        lineno++;
        if (line[0] == '#') {
            isHeader=true;
            unsigned int i=1;
            for (; i < strlen(line); i++)
                if (!(line[i] == '#' || isspace(line[i]))) break;
            line0 = &line[i];
        } else {
            isHeader=false;
            line0=line;
        }

        if (strncmp(line0, "NAMES\t", 6) == 0) {
            // parse NAMES line
            split(&line0[6], delim, sites->names);
            nseqs = sites->names.size();

            // assert every name is non-zero in length
            for (int i=0; i<nseqs; i++) {
                if (sites->names[i].length() == 0) {
                    if (!quiet)
                        printError(
                           "name for sequence %d is zero length (line %d)",
                           i + 1, lineno);
                    delete [] line;
                    return false;
                }
            }

        } else if (strncmp(line0, "REGION\t", 7) == 0) {
            // parse RANGE line
            char chrom[51];
            if (sscanf(line0, "REGION\t%50s\t%d\t%d",
                       chrom,
                       &sites->start_coord, &sites->end_coord) != 3) {
                if (!quiet) printError("bad REGION format");
                delete [] line;
                return false;
            }
            sites->chrom = chrom;
            sites->start_coord--;  // convert to 0-index

            // set region by subregion if specified
            if (subregion_start != -1)
                sites->start_coord = subregion_start;
            if (subregion_end != -1)
                sites->end_coord = subregion_end;


        } else if (strncmp(line0, "RANGE\t", 6) == 0) {
            // parse RANGE line
            if (!quiet)
                printError("deprecated RANGE line detected (use REGION instead)");
            delete [] line;
            return false;
        } else if (strncmp(line0, "POPS\t", 5) == 0) {
            if (nseqs == 0) {
                if (!quiet)
                    printError("NAMES line should come before POP line");
                delete [] line;
                return false;
            }
            vector<string> popstr;
            split(&line0[5], delim, popstr);
            if ((int)popstr.size() != nseqs) {
                if (!quiet)
                    printError("number of entries in POPS line should match entries in NAMES line");
                delete [] line;
                return false;
            }
            for (int i=0; i < nseqs; i++)
                sites->pops.push_back(atoi(popstr[i].c_str()));

        } else if (isHeader) {
            // no known tag; treat as comment
        } else {
            // parse a site line
            vector<string> fields;
            split(line, "\t", fields);
            assert(fields.size() >= 2);

            // parse site
            int position;
            if (sscanf(fields[0].c_str(), "%d", &position) != 1) {
                if (!quiet)
                    printError("first column is not an integer (line %d)", lineno);
                delete [] line;
                return false;
            }

            // skip site if not in region
            position--; //convert to 0-index
            if (position < sites->start_coord || position >= sites->end_coord) {
                delete [] line;
                continue;
            }

            // parse bases
            char col[fields[1].length()+1];
            strcpy(col, fields[1].c_str());
            unsigned int len = strlen(col);
            if (len != (unsigned int) nseqs) {
                if (!quiet)
                    printError(
                      "the number bases given, %d, does not match the "
                      "number of sequences %d (line %d)",
                      len, nseqs, lineno);
                delete [] line;
                return false;
            }
            if (!validate_site_column(col, nseqs)) {
                if (!quiet) {
                    printError("invalid sequence characters (line %d)", lineno);
                    printError("%s\n", line);
                }
                delete [] line;
                return false;
            }

            // validate site locations are unique and sorted.
            int npos = sites->get_num_sites();
            if (npos > 0 && sites->positions[npos-1] >= position) {
                if (!quiet) {
                    printError("invalid site location %d >= %d (line %d)",
                               sites->positions[npos-1], position, lineno);
                    printError("sites must be sorted and unique.");
                }
                delete [] line;
                return false;
            }

            // record site.
            sites->append(position, col, true);


            if (fields.size() == 2) {
                if (npos == 0) {
                    have_base_probs = false;
                } else if (have_base_probs) {
                    if (!quiet) {
                        printError("Error parsing line %d of sites file\n",
                                   lineno);
                    }
                    delete [] line;
                    return false;
                }
            } else {
                if (npos == 0) {
                    have_base_probs = true;
                } else if (!have_base_probs) {
                    if (!quiet) {
                        printError("Error parsing line %d of sites file\n",
                                   lineno);
                    }
                    delete [] line;
                    return false;
                }
                if ((int)fields.size() != 4*nseqs + 2) {
                    if (!quiet) {
                        printError("Error parsing base probs on line %i of sites file\n",
                                   lineno);
                    }
                    delete [] line;
                    return false;
                }
                vector<BaseProbs> bp_vec;
                bp_vec.clear();
                int pos=2;
                for (int i=0; i < nseqs; i++) {
                    BaseProbs bp;
                    for (int j=0; j < 4; j++)
                        sscanf(fields[pos++].c_str(), "%lf", &bp.prob[j]);
                    bp_vec.push_back(bp);
                }
                sites->base_probs.push_back(bp_vec);
            }

        }

        delete [] line;
    }

    return true;
}


// Read a Sites alignment file
bool read_sites(const char *filename, Sites *sites,
                int subregion_start, int subregion_end, bool quiet)
{
    CompressStream stream(filename);
    if (stream.stream == NULL) {
        if (!quiet)
            printError("cannot read file '%s'", filename);
        return false;
    }

    return read_sites(stream.stream, sites, subregion_start, subregion_end, quiet);
}


class GenoFilter {
public:
    string code;
    int cutoff;
    bool is_min;
    int index;

    GenoFilter(const char *str) {
        vector <string>tmp;
        split(str, ">", tmp);
        if (tmp.size() == 2) {
            is_min = false;
        } else {
            tmp.clear();
            split(str, "<", tmp);
            if (tmp.size() != 2) {
                printError("Could not parse genotype filter %s\n", str);
                assert(0);
            }
            is_min = true;
        }
        code = tmp[0];
        cutoff = atoi(tmp[1].c_str());
        index = -1;
    }
};


bool read_vcf(FILE *infile, Sites *sites, double min_qual,
              const char *genotype_filter, bool parse_genotype_probs,
              double min_base_prob, bool add_ref, const set<string> keep_inds) {
    const char *delim = "\t";
    char *line=NULL;
    int nseqs = 0, nsample=0;
    bool error = false;
    const char *headerStart = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t";
    string chrname = "";
    vector<GenoFilter> gf;
    int num_masked=0, total=0;
    int numIndel=0;
    static bool warnRefLen=false;
    static bool warnProbs=false;
    vector<bool> keep_ind;
    vector<string> sample_names;
    vector<int> ploidy;
    ploidy.clear();

    if (genotype_filter != NULL && strlen(genotype_filter) > 0) {
        vector<string> tmp;
        split(genotype_filter, ";", tmp);
        for (int i=0; i < (int)tmp.size(); i++) {
            gf.push_back(GenoFilter(tmp[i].c_str()));
        }
    }

    // note that this does not affect chrom, start_coord, end_coord
    sites->clear();

    int lineno = 1;
    while (!error) {
        if (line != NULL) delete [] line;
        line = fgetline(infile);
        if (line == NULL) break;
        chomp(line);
        lineno++;
        if (strncmp(line, "##", 2) == 0) {
            continue;
        }
        if (strncmp(line, headerStart, strlen(headerStart)) == 0) {
            split(&line[strlen(headerStart)], delim, sample_names);
            nsample = (int)sample_names.size();
            continue;
        }

        vector<string> fields;
        split(line, delim, fields);
        if ((int)fields.size() != 9 + nsample) {
            printError("Not enough fields in line %i of VCF file", lineno);
            return false;
        }
        if (chrname == "")
            chrname = fields[0];
        else if (chrname != fields[0]) {
            printError("VCF file contains multiple chromosomes. Must supply region str (chr:start-end)");
            return false;
        }
        int position;
        if (1 != sscanf(fields[1].c_str(), "%i", &position)) {
            printError("Error parsing position field in VCF\n");
            return false;
        }
        position--;  //convert to 0-index
        double qual = atof(fields[5].c_str());
        char alleles[5];  // alleles can only be A,C,G,T,N
        int num_alleles=1;
        if (fields[3].length() != 1) {
            if (!warnRefLen) {
                printWarning("Reference allele is not length one on line %i of VCF... skipping this and future similar lines",
                             lineno);
                warnRefLen=true;
            }
            numIndel++;
            continue;
        }
        alleles[0] = fields[3].c_str()[0];
        vector<string> alt;
        split(fields[4].c_str(), ",", alt);
        static bool badAlleleWarn=false;
        if (alt.size() > 4) {
            if (!badAlleleWarn) {
                printError("length of ALT allele should not be more than 4 on line %i of VCF\n",
                           lineno);
                badAlleleWarn=true;
            }
            continue;
        }
        bool badAllele=false;
        for (int i=0; i < (int)alt.size(); i++) {
            if (alt[i].length() != 1) {
                if (!badAlleleWarn) {
                    printWarning("ReadVCF can only handle alleles A,C,G,T,N currently;"
                                 " got allele %s on line %i; skipping this line and"
                                 " other similar ones",
                               alt[i].c_str(), lineno);
                    badAlleleWarn=true;
                }
                badAllele=true;
                numIndel++;
                break;
            }
            alleles[num_alleles++] = alt[i].c_str()[0];
        }
        if (badAllele) continue;
        // next: parse FORMAT in fields[8] and figure out where to find
        // GT
        vector<string> format;
        split(fields[8].c_str(), ":", format);
        int gt_idx=-1;
        int pl_idx=-1;
        int gl_idx=-1;
        int pp_idx=-1;
        for (int i=0; i < (int)format.size(); i++) {
            if (strcmp(format[i].c_str(), "GT")==0) {
                gt_idx=i;
            }
            if (parse_genotype_probs) {
                if (strcmp(format[i].c_str(), "PL")==0) {
                    pl_idx = i;
                }
                if (strcmp(format[i].c_str(), "GL")==0) {
                    gl_idx = i;
                }
                if (strcmp(format[i].c_str(), "PP")==0) {
                    pp_idx = i;
                }
            }
        }
        if (gt_idx == -1) {
            printError("Did not find GT in format field in VCF file line %i",
                       lineno);
            return false;
        }

        vector<BaseProbs> base_probs;
        // get positions for genotype filter(s)
        for (int i=0; i < (int)gf.size(); i++) {
            gf[i].index = -1;
            for (int j=0; j < (int)format.size(); j++) {
                if (format[j] == gf[i].code) {
                    gf[i].index = j;
                    break;
                }
            }
        }


        vector<string> seqfields;
        string gtstr;
        // on first input line, process sample names and figure out ploidy
        // (only ploidy 1 or two supported)
        if (ploidy.size() == 0) {
            nseqs = 0;
            for (int i=0; i < nsample; i++) {
                if (keep_inds.size() > 0 && keep_inds.find(sample_names[i]) == keep_inds.end()) {
                    keep_ind.push_back(false);
                }
                keep_ind.push_back(true);
                split(fields[9+i].c_str(), ":", seqfields);
                gtstr = seqfields[gt_idx];
                if (gtstr.length() == 1) {
                    ploidy.push_back(1);
                    if (keep_ind[i]) {
                        nseqs++;
                        sites->names.push_back(sample_names[i]);
                    }
                } else if (gtstr.length() == 3) {
                    ploidy.push_back(2);
                    if (keep_ind[i]) {
                        nseqs += 2;
                        for (int j=0; j < 2; j++) {
                            // append _1 and _2 suffix to distinguish two haplotypes of a diploid individual
                            char tmp[sample_names[i].length()+3];
                            sprintf(tmp, "%s_%i", sample_names[i].c_str(), j+1);
                            sites->names.push_back(string(tmp));
                        }
                    }
                } else {
                    printError("Bad genotype on line %i of VCF", lineno);
                    return false;
                }
            }
            if (add_ref) {
                sites->names.push_back("REF");
                nseqs++;
            }
            printf("nseqs = %i\n", nseqs - add_ref);
            // otherwise this line contains a variant
        }
        if (nseqs - add_ref  <= 0) {
            printError("Did not find sequences to keep in VCF file\n");
            return false;
        }

        char col[nseqs+1];
        col[nseqs] = '\0';
        int idx=0;
        for (int i=0; i < nsample; i++) {
            if (!keep_ind[i]) continue;
            bool masked = ( num_alleles > 2 || qual < min_qual );
            split(fields[9+i].c_str(), ":", seqfields);
            if (seqfields.size() != format.size()) {
                if (gt_idx < (int)seqfields.size() && seqfields[gt_idx] == "./.")
                    masked=true;
                else {
                    printError("Field %i does not match format string on line %i of VCF file\n",
                               9+i+1, lineno);
                    return false;
                }
            } else {
                for (int j=0; j < (int)gf.size(); j++) {
                    if (gf[j].index >= 0) {
                        int val = atoi(seqfields[gf[j].index].c_str());
                        if ((  gf[j].is_min  && val < gf[j].cutoff) ||
                            ((!gf[j].is_min) && val > gf[j].cutoff)) {
                            masked=true;
                            break;
                        }
                    }
                }
            }

            if (parse_genotype_probs) {
                BaseProbs bp = BaseProbs('N');
                base_probs.push_back(bp);
                if (ploidy[i] == 2) base_probs.push_back(bp);
            }
            total++;
            if (masked) {
                col[idx] = 'N';
                if (ploidy[i] == 2) col[idx+1] = 'N';
                idx += ploidy[i];
                num_masked+= ploidy[i];
                continue;
            }
            if (parse_genotype_probs && pl_idx == -1 &&
                gl_idx == -1 && pp_idx == -1) {
                if (!warnProbs) {
                    printWarning("Did not find PL, GL, or PP in format field in VCF file line %i",
                                 lineno);
                    warnProbs=true;
                }
            }
            gtstr = seqfields[gt_idx];
            if (ploidy[i]==2 && gtstr.length() != 3) {
                printError("genotype not length three on line %i of VCF",
                           lineno);
                return false;
            }
            if (ploidy[i]==1 && gtstr.length() != 1) {
                printError("genotype not length one on line %i of VCF for haploid sample",
                           lineno);
                return false;
            }
            if (ploidy[i] == 2) {
                if (gtstr.c_str()[1] != '|' &&
                    gtstr.c_str()[1] != '/') {
                    printError("genotype middle character not '|' or '/' on line %i",
                               lineno);
                    return false;
                }
            }
            for (int j=0; j < ploidy[i]; j++) {
                char allele = gtstr.c_str()[j*2];
                if (allele == '.') {
                    col[idx] = 'N';
                    if (parse_genotype_probs)
                        base_probs[idx].set_mask();
                } else {
                    int ia = allele - '0';
                    if (ia < 0 || ia >= num_alleles) {
                        printError("Bad GT in field %i,line %i of VCF",
                                   i+9+1, lineno);
                        return false;
                    }
                    col[idx] = alleles[ia];
                    if (parse_genotype_probs) {
                        // PL and GL are the same except GL is float;
                        // set_by_pl treats input as float anyway
                        if (gl_idx >= 0) pl_idx = gl_idx;
                        if (pl_idx >= 0)
                            base_probs[idx].set_by_pl(alleles[0], alleles[1],
                                                      seqfields[pl_idx], j);
                        else if (pp_idx >= 0)
                            base_probs[idx].set_by_pp(seqfields[pp_idx], j);
                        else base_probs[idx].set_certain(alleles[ia]);
                        if (base_probs[idx].maxProb() < min_base_prob) {
                            col[idx] = 'N';
                            base_probs[idx].set_mask();
                            num_masked++;
                        }
                    }
                }
                idx++;
            }
        }
        if (add_ref) {
            assert(idx == nseqs-1);
            col[idx] = alleles[0];
            if (parse_genotype_probs)
                base_probs.push_back(BaseProbs(alleles[0]));
        }
        sites->append(position, col, true);
        if (parse_genotype_probs)
            sites->base_probs.push_back(base_probs);
    }

    printLog(LOG_LOW, "Read %i sites from %i lines of VCF file (num skipped indels=%i)\n",
             sites->get_num_sites(), lineno, numIndel);
    if (gf.size() > 0) printLog(LOG_LOW, "Masked %.1f out of %i genotypes\n",
                                (double)num_masked/2, total);
    return true;
}


bool read_vcf(const char *filename, Sites *sites, const char *region,
              double min_qual, const char *genotype_filter,
              bool parse_genotype_probs, double min_base_prob, bool add_ref,
              const char *tabixdir, const set<string> keep_inds) {
    TabixStream ts(filename, region, tabixdir);
    char chr[10000];
    int start_coord, end_coord;
    if (region == NULL || strlen(region) == 0) {
        printError("read_vcf requires --region string\n");
        return false;
    }
    vector<string> tmp;
    split(region, ":", tmp);
    if (tmp.size() != 2) {
        printError("Error parsing region string %s. Must be in format chr:start-end\n", region);
        return false;
    }
    strcpy(chr, tmp[0].c_str());
    if (2 != (sscanf(tmp[1].c_str(), "%i-%i", &start_coord, &end_coord))) {
        printError("Error parsing region string %s. Must be in format chr:start-end\n", tmp[1].c_str());
        return false;
    }
    if (ts.stream == NULL) {
        return false;
    }
    printLog(LOG_LOW, "Reading %s %s\n", filename, region);
    sites->start_coord = start_coord - 1;
    sites->end_coord = end_coord;
    sites->chrom = string(chr);
    if ( ! read_vcf(ts.stream, sites, min_qual, genotype_filter,
                    parse_genotype_probs, min_base_prob, add_ref, keep_inds))
        return false;
    return true;
}

bool read_vcf(const string filename, Sites *sites, const string region,
              double min_qual, const string genotype_filter,
              bool parse_genotype_probs, double min_base_prob, bool add_ref,
              const string tabixdir, const set<string> keep_inds) {
    return read_vcf(filename.c_str(), sites, region.c_str(),
                    min_qual, genotype_filter.c_str(), parse_genotype_probs,
                    min_base_prob, add_ref, tabixdir.c_str(), keep_inds);
}

bool read_vcfs(const vector<string> filenames, Sites* sites, const string region,
               double min_qual, const string genotype_filter,
               bool parse_genotype_probs, double min_base_prob,
               const string tabixdir, const set<string> keep_inds) {
    if (filenames.size() == 0) {
        fprintf(stderr, "Read_vcfs expects at least one filename\n");
        return false;
    }
    if (!read_vcf(filenames[0], sites, region, min_qual,
                  genotype_filter, parse_genotype_probs, min_base_prob,
                  true, tabixdir, keep_inds))
        return false;
    for (int i=1; i < (int)filenames.size(); i++) {
        Sites s;
        if (!read_vcf(filenames[i], &s, region, min_qual,
                      genotype_filter, parse_genotype_probs, min_base_prob,
                      true, tabixdir, keep_inds))
            return false;
        if (!sites->merge(s)) return false;
    }
    // need to remove REF
    vector<int> keep;
    for (int i=0; i < sites->get_num_seqs(); i++) {
        if (sites->names[i] == "REF") continue;
        bool duplicated=false;
        for (int j=0; j < i; j++)
            if (sites->names[j] == sites->names[i]) {
                fprintf(stderr, "read_vcfs: found multiple sequences named %s; removing one\n",
                        sites->names[j].c_str());
                duplicated=true;
                break;
            }
        if (!duplicated)
            keep.push_back(i);
    }
    sites->subset(keep);
    return true;
}



// Converts a Sites alignment to a Sequences alignment
void make_sequences_from_sites(const Sites *sites, Sequences *sequences,
                               char default_char)
{
    int nseqs = sites->names.size();
    int seqlen = sites->length();
    int start = sites->start_coord;
    int nsites = sites->get_num_sites();
    bool have_pops = ( sites->pops.size() > 0);
    bool have_base_probs = ( sites->base_probs.size() > 0 );

    sequences->clear();
    sequences->set_owned(true);

    for (int i=0; i<nseqs; i++) {
        char *seq = new char [seqlen+1];
        seq[seqlen] = '\0';
        vector<BaseProbs> base_probs;
        int col = 0;
        for (int j=0; j<seqlen; j++) {
            if (col < nsites && start+j == sites->positions[col]) {
                // variant site
                seq[j] = sites->cols[col][i];
                if (have_base_probs)
                    base_probs.push_back(sites->base_probs[col][i]);
                col++;
            } else {
                seq[j] = default_char;
                if (have_base_probs)
                    base_probs.push_back(BaseProbs(default_char));
            }
        }

        sequences->append(sites->names[i], seq, base_probs,
                          have_pops ? sites->pops[i] : 0);
    }

    sequences->set_length(seqlen);
}


int Sites::remove_invariant() {
    int nsite=0;
    int orig_nsite = (int)positions.size();
    int nseq = get_num_seqs();
    bool have_base_probs = ( base_probs.size() > 0 );
    for (int i=0; i < orig_nsite; i++) {
        bool variant=false;
        char allele1='N';
        for (int j=0; j < nseq; j++) {
            if (cols[i][j] == 'N') {
                variant=true;
            } else if (have_base_probs && !base_probs[i][j].is_certain()) {
                variant=true;
            } else if (allele1=='N') {
                allele1 = cols[i][j];
            } else if (allele1 != cols[i][j]) {
                variant=true;
            }
            if (variant) break;
        }
        if (variant) {
            if (i != nsite) {
                positions[nsite] = positions[i];
                strncpy(cols[nsite], cols[i], nseq);
                if (have_base_probs)
                    base_probs[nsite] = base_probs[i];
            }
            nsite++;
        }
    }
    if (nsite < orig_nsite) {
        positions.resize(nsite);
        cols.resize(nsite);
        if (have_base_probs)
            base_probs.resize(nsite);
    }
    return orig_nsite - nsite;
}


bool Sites::rename(string rename_file) {
    FILE *infile = fopen(rename_file.c_str(), "r");
    if (infile == NULL) {
        printError("Error opening %s", rename_file.c_str(), "\n");
        assert(0);
    }
    char oldname[1000], newname[1000];
    int nseq = get_num_seqs();
    while (EOF != fscanf(infile, "%s %s", oldname, newname)) {
        // first check for exact matches
        bool found=false;
        for (int i=0; i < nseq; i++) {
            if (strcmp(names[i].c_str(), oldname) == 0) {
                names[i] = string(newname);
                found=true;
                break;
            }
        }
        // now check for oldname_1, oldname_2
        if (!found) {
            char tmpstr[2][3]={"_1", "_2"};
            char tmptarget[2][1002];
            char newname2[1002];
            for (int j=0; j < 2; j++)
                sprintf(tmptarget[j], "%s%s", oldname, tmpstr[j]);
            for (int i=0; i < nseq; i++) {
                for (int j=0; j < 2; j++) {
                    if (strcmp(tmptarget[j], names[i].c_str())==0) {
                        sprintf(newname2, "%s%s", newname, tmpstr[j]);
                        names[i] = string(newname2);
                    }
                }
            }
        }
    }
    fclose(infile);
    return true;
}

int Sites::subset(vector<int> keep) {
    bool have_base_probs = ( base_probs.size() > 0 );
    std::sort(keep.begin(), keep.end());
    vector<string> new_names;
    for (unsigned int i=0; i < keep.size(); i++)
        new_names.push_back(names[keep[i]]);
    names = new_names;
    if (pops.size() != 0) {
        vector<int> new_pops;
        for (unsigned int i=0; i < keep.size(); i++)
            new_pops.push_back(pops[keep[i]]);
        pops = new_pops;
    }
    vector<int> new_positions;
    vector<char*> new_cols;
    vector<vector<BaseProbs> > new_base_probs;
    for (unsigned int i=0; i < positions.size(); i++) {
        char *tmp = new char[keep.size()+1];
        bool variant=false;
        vector<BaseProbs> bp;
        for (unsigned int j=0; j < keep.size(); j++) {
            tmp[j] = cols[i][keep[j]];
            if (have_base_probs) {
                if (!base_probs[i][keep[j]].is_certain())
                    variant=true;
                bp.push_back(base_probs[i][keep[j]]);
            }
            if (tmp[j]=='N' || tmp[j] != tmp[0]) variant=true;
        }
        tmp[keep.size()] = '\0';
        delete [] cols[i];
        if (variant) {
            new_cols.push_back(tmp);
            new_positions.push_back(positions[i]);
            if (have_base_probs) new_base_probs.push_back(bp);
        } else delete [] tmp;
    }
    cols = new_cols;
    positions = new_positions;
    if (have_base_probs) base_probs = new_base_probs;
    printLog(LOG_LOW, "subset sites (nseqs=%i, nsites=%i)\n",
             names.size(), positions.size());
    return 0;
}

int Sites::subset(set<string> names_to_keep) {
    vector<int> keep;
    for (unsigned int i=0; i < names.size(); i++) {
        if (names_to_keep.find(names[i]) != names_to_keep.end())
            keep.push_back(i);
    }
    if (keep.size() > 0) {
        if (keep.size() != names_to_keep.size()) {
            fprintf(stderr, "Error in subset: not all names found in sites\n");
            return 1;
        }
        return subset(keep);
    }

    // if nothing found, these may be individual names rather than haploid
    // names
    set<string> hapkeep;
    for (set<string>::iterator it=names_to_keep.begin();
         it != names_to_keep.end(); it++) {
        hapkeep.insert((*it) + string("_1"));
        hapkeep.insert((*it) + string("_2"));
    }
    for (unsigned int i=0; i < names.size(); i++) {
        if (hapkeep.find(names[i]) != hapkeep.end())
            keep.push_back(i);
    }
    if (keep.size() != hapkeep.size()) {
        fprintf(stderr, "Error in subset: not all names found in sites\n");
        return 1;
    }
    return subset(keep);
}

template<>
int Sites::remove_overlapping(const TrackNullValue &track) {
    bool have_base_probs = ( base_probs.size() > 0 );
    int numhap = get_num_seqs();
    int idx=0, mask_idx=0;
    assert(track.is_sorted());
    for (int i=0; i < (int)positions.size(); i++) {
        bool overlapping = track.find(positions[i], &mask_idx, true);
        if (overlapping) continue;
        if (i != idx) {
            positions[idx] = positions[i];
            strncpy(cols[idx], cols[i], numhap);
            if (have_base_probs)
                base_probs[idx] = base_probs[i];
        }
        idx++;
    }
    for (unsigned int i=idx; i < cols.size(); i++)
        delete [] cols[i];
    positions.resize(idx);
    cols.resize(idx);
    if (have_base_probs)
        base_probs.resize(idx);
    return 0;
}


TrackNullValue Sites::get_masked_regions() const {
    TrackNullValue masked_regions;
    int numhap = get_num_seqs();
    for (int i=0; i < (int)positions.size(); i++) {
        bool masked=true;
        for (int j=0; j < numhap; j++) {
            if (cols[i][j] != 'N') {
                masked=false;
                break;
            }
        }
        if (masked) {
            masked_regions.push_back(RegionNullValue(chrom, positions[i],
                                                     positions[i]+1, ' '));
        }
    }
    masked_regions.merge();
    return masked_regions;
}


TrackNullValue Sequences::get_masked_regions(string chrom,
                                             const SitesMapping *sites_mapping)
    const {
    TrackNullValue masked_regions;
    int numhap = get_num_seqs();
    for (int i=0; i < length(); i++) {
        bool masked=true;
        for (int j=0; j < numhap; j++) {
            if (seqs[j][i] != 'N') {
                masked=false;
                break;
            }
        }
        if (masked) {
            int start = ( sites_mapping ?
                          sites_mapping->uncompress_start(i) :
                          i );
            int end = ( sites_mapping ?
                        sites_mapping->uncompress_end(i)+1 :
                        ( i+1 ));
            masked_regions.push_back(RegionNullValue(chrom, start, end, ' '));
        }
    }
    masked_regions.merge();
    return masked_regions;
}


bool Sites::merge(const Sites &other) {
    if (chrom != other.chrom) {
        fprintf(stderr, "Error: Cannot merge sites from different chromosomes\n");
        return false;
    }
    if (start_coord != other.start_coord ||
        end_coord != other.end_coord) {
        fprintf(stderr, "Error: regions do not match in sites.merge()\n");
        return false;
    }
    int old_nseq =  get_num_seqs();
    int other_nseq = other.get_num_seqs();
    int old_ref = -1;
    int other_ref = -1;
    for (int i=0; i < old_nseq; i++)
        if (names[i] == "REF") {
            old_ref = i;
            break;
        }
    for (int i=0; i < other_nseq; i++)
        if (other.names[i] == "REF") {
            other_ref = i;
            break;
        }
    if (old_ref == -1 || other_ref == -1) {
        fprintf(stderr, "Sites.merge() requires both sites to have an sequence named REF");
        return false;
    }
    if (pops.size() != 0 || other.pops.size() != 0) {
        if ((int)pops.size() != old_nseq || (int)other.pops.size() != other_nseq) {
            fprintf(stderr, "Warning in sites.merge(): both sites do not "
                    "have population information; dropping");
            pops.clear();
        }
    }

    int old_num_sites = get_num_sites();
    int other_num_sites = other.get_num_sites();
    for (int i=0; i < other_nseq; i++) {
        if (i != other_ref) {
            names.push_back(other.names[i]);
            if (pops.size() != 0)
                pops.push_back(other.pops[i]);
        }
    }

    vector<int> old_positions = positions;
    vector<char*> old_cols = cols;
    vector<vector<BaseProbs> > old_base_probs = base_probs;
    bool have_base_probs=false;
    if (old_base_probs.size() > 0 || other.base_probs.size() > 0)  {
        if ((int)old_base_probs.size() != old_num_sites ||
            (int)other.base_probs.size() != other_num_sites) {
            fprintf(stderr, "Error in Sites.merge(); both sites should have"
                    " base probabilities, or neither");
            return false;
        }
        have_base_probs=true;
    }
    positions.clear();
    cols.clear();
    base_probs.clear();
    int i1=0, i2=0;
    int num_seq = names.size();
    char col[num_seq + 1];
    col[num_seq] = '\0';
    vector<BaseProbs>bp;
    while (i1 < old_num_sites || i2 < other_num_sites) {
        int pos=-1;
        bp.clear();
        if (i1 == old_num_sites ||
            (i2 < other_num_sites && other.positions[i2] < old_positions[i1])) {
            // next site is from "other"
            int i=0;
            char ref = other.cols[i2][other_ref];
            char fixedAllele = ref;
            if (fixedAllele == 'N') fixedAllele='A';
            for (; i < old_nseq; i++) {
                col[i] = ( i == old_ref ? ref : fixedAllele );
                if (have_base_probs)
                    bp.push_back(BaseProbs(col[i]));
            }
            for (int j=0; j < other_nseq; j++) {
                if (j != other_ref) {
                    col[i++] = other.cols[i2][j];
                    if (have_base_probs)
                        bp.push_back(other.base_probs[i2][j]);
                }
            }
            assert(i == num_seq);
            pos = other.positions[i2];
            i2++;
        } else if (i2 == other_num_sites || old_positions[i1] < other.positions[i2]) {
            // next site is from orig object
            int i=0;
            char fixedAllele = old_cols[i1][old_ref];
            if (fixedAllele == 'N') fixedAllele = 'A';
            for ( ; i < old_nseq; i++) {
                col[i] = old_cols[i1][i];
                if (have_base_probs)
                    bp.push_back(old_base_probs[i1][i]);
            }
            while (i < num_seq) {
                col[i++] = fixedAllele;
                if (have_base_probs)
                    bp.push_back(BaseProbs(fixedAllele));
            }
            pos = old_positions[i1];
            i1++;
        } else {
            // both objects have this site; merge
            assert(old_positions[i1] == other.positions[i2]);
            int i=0;
            for (; i < old_nseq; i++) {
                col[i] = old_cols[i1][i];
                if (have_base_probs)
                    bp.push_back(old_base_probs[i1][i]);
            }
            for (int j=0; j < other_nseq; j++) {
                if (j != other_ref) {
                    col[i++] = other.cols[i2][j];
                    if (have_base_probs)
                        bp.push_back(other.base_probs[i2][j]);
                }
            }
            assert(i == num_seq);
            pos = old_positions[i1];
            i1++;
            i2++;
        }
        append(pos, col, true);
        if (have_base_probs) {
            assert((int)bp.size() == num_seq);
            base_probs.push_back(bp);
        }
    }
    for (int i=0; i < (int)old_cols.size(); i++)
        delete [] old_cols[i];
    return true;
}


TrackNullValue Sites::remove_masked() {
    TrackNullValue masked_regions;
    int numhap = get_num_seqs();
    bool have_base_probs = ( base_probs.size() > 0 );
    int idx=0;
    for (int i=0; i < (int)positions.size(); i++) {
        bool masked=true;
        for (int j=0; j < numhap; j++) {
            if (cols[i][j] != 'N') {
                masked=false;
                break;
            }
        }
        if (masked) {
            masked_regions.push_back(RegionNullValue(chrom, positions[i],
                                                     positions[i]+1, ' '));
        } else {
            if (idx != i) {
                positions[idx] = positions[i];
                strncpy(cols[idx], cols[i], numhap);
                if (have_base_probs)
                    base_probs[idx] = base_probs[i];
            }
            idx++;
        }
    }
    for (unsigned int i=idx; i < cols.size(); i++)
        delete [] cols[i];
    positions.resize(idx);
    cols.resize(idx);
    if (have_base_probs)
        base_probs.resize(idx);
    masked_regions.merge();
    return masked_regions;
}

void apply_mask_sequences(Sequences *sequences,
                          const TrackNullValue &maskmap,
                          const char *ind)
{
    const char maskchar = 'N';
    int maskind[sequences->get_num_seqs()];
    int num_mask=0;
    bool have_base_probs = ( sequences->base_probs.size() > 0 );

    if (ind == NULL) {
        num_mask = (int)sequences->get_num_seqs();
        for (int i=0; i < num_mask; i++)
            maskind[i]=i;
    } else {
        // first see if there is an exact match
        for (int i=0; i < (int)sequences->get_num_seqs(); i++) {
            if (strcmp(sequences->names[i].c_str(), ind)==0) {
                num_mask = 1;
                maskind[0] = i;
                printLog(LOG_HIGH, "Applying individual mask to %s\n", ind);
                break;
            }
        }
        // else check to see if it is a diploid match
        if (num_mask == 0) {
            for (int j=0; j <= 2; j++) {
                char hapname[strlen(ind)+3];
                sprintf(hapname, "%s_%i", ind, j);
                for (int i=0; i < (int)sequences->get_num_seqs(); i++) {
                    if (strcmp(sequences->names[i].c_str(), hapname)==0) {
                        maskind[num_mask++] = i;
                        printLog(LOG_LOW, "Applying individual mask to %s\n",
                                 sequences->names[i].c_str());
                        break;
                    }
                }
            }
            if (num_mask != 2) {
                printLog(LOG_LOW, "Warning: No sequences matched ind mask name %s; mask not applied\n", ind);
                return;
            }
        }
    }

    for (unsigned int k=0; k<maskmap.size(); k++) {
        for (int i=maskmap[k].start; i<maskmap[k].end; i++) {
            for (int j=0; j < num_mask; j++) {
                sequences->seqs[maskind[j]][i] = maskchar;
                if (have_base_probs)
                    sequences->base_probs[maskind[j]][i].set_mask();
            }
        }
    }
}




/* Apply track maskmap to samples indicated by "ind" (either on haplotype
   or two diploid lineages named ind_1 and ind_2). If a masked region falls
   within a variant position of the sites object, then the allele is changed to
   'N' for this position. However, invariant sites are not added to the object.
   Instead, the maskmap is added to ind_masks for the samples indicated by ind
 */
void apply_mask_sites(Sites *sites,
                      const TrackNullValue &mask,
                      const char *ind,
                      vector<TrackNullValue> *ind_masks)
{
    const char maskchar = 'N';
    int maskind[sites->get_num_seqs()];
    int num_mask=0;
    bool have_base_probs = ( sites->base_probs.size() > 0 );

    // first see if there is an exact match to ind
    for (int i=0; i < (int)sites->get_num_seqs(); i++) {
        if (strcmp(sites->names[i].c_str(), ind)==0) {
            num_mask = 1;
            maskind[0] = i;
            printLog(LOG_HIGH, "Applying individual mask to %s\n", ind);
            break;
        }
    }
    // else check to see if it is a diploid match
    if (num_mask == 0) {
        for (int j=0; j <= 2; j++) {
            char hapname[strlen(ind)+3];
            sprintf(hapname, "%s_%i", ind, j);
            for (int i=0; i < (int)sites->get_num_seqs(); i++) {
                if (strcmp(sites->names[i].c_str(), hapname)==0) {
                    maskind[num_mask++] = i;
                    printLog(LOG_LOW, "Applying individual mask to %s\n",
                             sites->names[i].c_str());
                    break;
                }
            }
        }
        if (num_mask != 2) {
            printLog(LOG_LOW, "Warning: No sequences matched ind mask name %s; mask not applied\n", ind);
            return;
        }
    }
    for (int j=0; j < num_mask; j++)
        (*ind_masks)[maskind[j]].merge_tracks(mask);

    unsigned int last_pos=0;
    int last_mask_pos=0;
    for (unsigned int k=0; k<mask.size(); k++) {
        if (mask[k].start < last_mask_pos) last_pos=0;
        for (int i=mask[k].start; i<mask[k].end; i++) {
            while (last_pos < sites->positions.size() &&
                   sites->positions[last_pos] < i)
                last_pos++;
            if (last_pos != sites->positions.size() &&
                sites->positions[last_pos] == i) {
                for (int j=0; j < num_mask; j++) {
                    sites->cols[last_pos][maskind[j]] = maskchar;
                    if (have_base_probs)
                        sites->base_probs[last_pos][maskind[j]].set_mask();
                }
            }
            last_mask_pos = i;
        }
    }
}



// Converts a Sequences alignment to a Sites alignment
void make_sites_from_sequences(const Sequences *sequences, Sites *sites)
{
    int nseqs = sequences->get_num_seqs();
    int seqlen = sequences->length();
    const char * const *seqs = sequences->get_seqs();
    bool have_base_probs = ( sequences->base_probs.size() > 0 );

    sites->clear();
    sites->start_coord = 0;
    sites->end_coord = seqlen;
    sites->names.insert(sites->names.begin(),
                        sequences->names.begin(), sequences->names.end());

    for (int i=0; i<seqlen; i++) {
        bool isSite=false;  //need to make site if position is N or varaint or
                            // has baseprobs and is not certain
        for (int j=0; j < nseqs; j++) {
            if (seqs[j][i] == 'N' ||
                seqs[j][i] != seqs[0][i] ||
                (have_base_probs && !sequences->base_probs[j][i].is_certain())) {
                isSite=true;
                break;
            }
        }
        if (!isSite) continue;

        char *col = new char [nseqs+1];
        col[nseqs]='\0';
        for (int j=0; j<nseqs; j++)
            col[j] = seqs[j][i];
        sites->append(i, col);
        if (have_base_probs) {
            vector<BaseProbs> bp;
            for (int j=0; j < nseqs; j++)
                bp.push_back(sequences->base_probs[j][i]);
            sites->base_probs.push_back(bp);
        }
    }
}


bool Sequences::get_non_singleton_snp(vector<bool> &nonsing) {
    if (seqs.size() == 0) return false;
    for (int i=0; i < seqlen; i++) {
        char allele[4]={'N','N','N','N'};
        int count[4]={0,0,0,0};
        for (unsigned int j=0; j < seqs.size(); j++) {
            char c = seqs[j][i];
            if (c == 'N') continue;
            for (int k=0; k < 4; k++) {
                if (allele[k] == 'N')
                    allele[k] = c;
                if (allele[k] == c) {
                    count[k]++;
                    break;
                }
            }
        }
        int total = 0;
        for (int k=0; k < 4; k++)
            total += (count[k] > 1);
        nonsing[i] = (total >= 2);
    }
    return true;
}

bool Sequences::set_pairs_by_name() {
    int numfound=0;
    pairs.resize(names.size());
    for (unsigned int i=0; i < names.size(); i++) pairs[i] = -1;
    for (unsigned int i=0; i < names.size(); i++) {
        if (pairs[i] != -1) continue;
        int len = names[i].length();
        if (len < 3) continue;
        string basename = names[i].substr(0, len-2);
        string ext = names[i].substr(len-2, 2);
        string target_ext;
        if (ext == "_1") {
            target_ext = "_2";
        } else if (ext == "_2") {
            target_ext = "_1";
        } else {
            continue;
        }
        for (unsigned int j=i+1; j < names.size(); j++) {
            if (names[i].substr(0, len-2) == basename &&
                names[j].substr(len-2, 2) == target_ext) {
                pairs[i] = j;
                pairs[j] = i;
                numfound += 2;
                break;
            }
        }
    }
    return (numfound > 0);
}

void Sequences::set_pairs_from_file(string fn) {
  FILE *infile = fopen(fn.c_str(), "r");
  char p1[1000], p2[1000];
  if (infile == NULL) {
    fprintf(stderr, "Error opening %s\n", fn.c_str());
    exit(-1);
  }
  pairs = vector<int>(names.size());
  for (unsigned int i=0; i < names.size(); i++)
    pairs[i] = -1;
  int x;
  while (EOF != (x=fscanf(infile, "%s %s", p1, p2))) {
    if (x != 2) {
      fprintf(stderr, "Error: bad format for file %s\n", fn.c_str());
      exit(-1);
    }
    int x1=-1, x2=-1;
    for (unsigned int i=0; i < names.size(); i++) {
      if (names[i] == p1) {
	x1=i;
	if (x2 >= 0) break;
      } else if (names[i] == p2) {
	x2=i;
	if (x1 >= 0) break;
      }
    }
    if (x1 >=0 && x2 >= 0) {
      pairs[x1] = x2;
      pairs[x2] = x1;
    }
  }
  fclose(infile);
}


void Sequences::set_pops_from_file(string fn) {
    FILE *infile = fopen(fn.c_str(), "r");
    char seqname[1000];
    int pop;
    if (infile == NULL) {
        fprintf(stderr, "Error opening %s\n", fn.c_str());
        exit(-1);
    }
    pops = vector<int>(names.size());
    for (unsigned int i=0; i < names.size(); i++) pops[i]=-1;
    while (2 == fscanf(infile, "%s %i", seqname, &pop)) {
        unsigned int i=0;
        for (i=0; i < names.size(); i++) {
            if (names[i] == seqname) {
                pops[i] = pop;
                break;
            }
        }
        if (i == names.size()) {
            int found=0;
            for (int hap=1; hap <=2; hap++) {
                char tmpname[1003];
                sprintf(tmpname, "%s_%i", seqname, hap);
                for (int i=0; i < (int)names.size(); i++)
                    if (names[i] ==tmpname) {
                        pops[i] = pop;
                        found++;
                    }
            }
            if (found != 2) {
                printError("set_pops_from_file: do not see sequence %s\n", seqname);
            }
        }
    }
    for (unsigned int i=0; i < names.size(); i++) {
        if (pops[i] == -1) {
            printError("set_pops_from_file: sequence %s does not have assignment\n",
                       names[i].c_str());
            assert(0);
        }
    }
}



void Sequences::set_pairs(const ArgModel *mod) {
    if (mod->unphased_file != "") {
        set_pairs_from_file(mod->unphased_file);
        return;
    }
    if (set_pairs_by_name()) return;

    printWarning("set_pairs: --unphased-file not given and haploid names do not"
                 " appear to follow convention <ind>_1, <ind>_2. Assuming that"
                 " the two sequences for each individual are adjacent in the"
                 " sequence file");
    pairs = vector<int>(names.size());
    for (unsigned int i=0; i < names.size(); i++) {
        if (i%2==0) {
            if (i+1 < names.size()) pairs[i] = i+1;
        } else pairs[i] = i-1;
    }
}

void PhaseProbs::sample_phase(int *thread_path) {
    int sing_tot=0, sing_switch=0, non_sing_tot=0, non_sing_switch=0;
    if (probs.size() == 0)
        return;
    bool have_base_probs = (seqs->base_probs.size() > 0);
    for (map<int,vector<double> >::iterator it=probs.begin(); it != probs.end();
       it++) {
        int coord = it->first;
        vector<double> prob = it->second;
        if (seqs->seqs[hap1][coord] != seqs->seqs[hap2][coord] ||
            (have_base_probs &&
             !seqs->base_probs[hap1][coord].is_equal(seqs->base_probs[hap2][coord]))) {
            if (frand() > prob[thread_path[coord]]) {
                seqs->switch_alleles(coord, hap1, hap2);
                if (non_singleton_snp[coord])
                    non_sing_switch++;
                else sing_switch++;
            }
            if (non_singleton_snp[coord])
                non_sing_tot++;
            else sing_tot++;
        }
    }
    printLog(LOG_LOW, "sample_phase %s %s size=%i %i frac_switch=%f %f\n",
             seqs->names[hap1].c_str(), seqs->names[hap2].c_str(),
             sing_tot, non_sing_tot,
             (double)sing_switch/sing_tot,
             (double)non_sing_switch/non_sing_tot);
}


void Sequences::randomize_phase(double frac) {
    printLog(LOG_LOW, "randomizing phase (frac=%f)\n", frac);
    int count=0, total=0;
    for (int i=0; i < seqlen; i++) {
        for (int j=0; j < (int)seqs.size(); j++) {
            if (pairs[j] < j) continue;
            total++;
            if (frand() < frac) {
                if (frand() < 0.5) {
                    switch_alleles(i, j, pairs[j]);
                    count++;
                }
            }
        }
    }
    printLog(LOG_LOW, "switched %i of %i (%f)\n", count, total,
             (double)count/(double)total);
}

void Sequences::add_switch_errors(double rate) {
    printLog(LOG_LOW, "Adding phase switch errors at rate %.2g\n", rate);
    int num_switches = 0, num_pairs=0;
    for (int i1=0; i1 < (int)seqs.size(); i1++) {
        int i2=pairs[i1];
        if (i2 < i1) continue;
        bool switch_state = false;
        num_pairs++;
        for (int j=0; j < seqlen; j++) {
            if (seqs[i1][j] == seqs[i2][j])
                continue;
            if (seqs[i1][j] != 'N' && seqs[i2][j] != 'N') {
                if (frand() < rate) {
                    switch_state = !switch_state;
                    num_switches++;
                }
            }
            if (switch_state)
                switch_alleles(j, i1, i2);
        }
    }
    printLog(LOG_LOW, "introduced %i switches to %i haploid pairs\n", num_switches, num_pairs);
}

void Sequences::set_age() {
    int nseqs = names.size();
    ages.resize(nseqs, 0);
    real_ages.resize(nseqs, 0);
}


void Sequences::set_age(int sample, double time, int ntimes, const double *times) {
    real_ages[sample] = time;
    double mindif = fabs(times[0] - time);
    int whichmin = 0;
    for (int j=1; j < ntimes - 1; j++) {
        double tempdif = fabs(times[j] - time);
        if (tempdif < mindif) {
            whichmin = j;
            mindif = tempdif;
        }
        if (times[j] > time) break;
    }
    ages[sample] = whichmin;
    printLog(LOG_LOW, "rounded age for sample %s to %f\n",
             names[sample].c_str(), times[whichmin]);
}

void Sequences::set_age(string agefile, int ntimes, const double *times) {
    char currseq[1000];
    FILE *infile = fopen(agefile.c_str(), "r");
    double time;
    int nseqs=names.size();
    if (infile == NULL) {
	fprintf(stderr, "Error opening %s\n", agefile.c_str());
	exit(-1);
    }
    ages.resize(nseqs, 0);
    real_ages.resize(nseqs, 0);
    while (EOF != fscanf(infile, "%s %lf", currseq, &time)) {
        int found = 0;
	for (int i=0; i < nseqs; i++) {
	    if (strcmp(names[i].c_str(), currseq)==0) {
                set_age(i, time, ntimes, times);
                found++;
            }
	}
        if (found > 1) {
            fprintf(stderr, "Found too many matches for %s from agefile\n", currseq);
            exit(0);
        }
	if (found==0) {
            char tmpstr[strlen(currseq)+3];
            for (int hap=1; hap <= 2; hap++) {
                found=0;
                sprintf(tmpstr, "%s_%i", currseq, hap);
                for (int i=0; i < nseqs; i++) {
                    if (strcmp(names[i].c_str(), tmpstr) == 0) {
                        set_age(i, time, ntimes, times);
                        found++;
                    }
                }
                if (found != 1) {
                    fprintf(stderr, "WARNING: could not find sequence %s (from %s) in sequences\n",
		    currseq, agefile.c_str());
                    exit(0);
                }
            }
	}
    }
    fclose(infile);
}

// Compress the sites by a factor of 'compress'.
//
// Return true if compression is successful.
//
// The compression maintains the following properties:
//
// - The coorindate system (start_coord, end_coord) will be adjusted to
//   roughly (0, seqlen / compress).
// - Every variant column is kept.
//
bool find_compress_cols(const Sites *sites, int compress,
                        SitesMapping *sites_mapping)
{
    const int ncols = sites->get_num_sites();

    int blocki = 0;
    int next_block = sites->start_coord + compress;
    int half_block = compress / 2;

    // record old coords
    sites_mapping->init(sites);

    // special case
    if (compress == 1) {
        for (int i=sites->start_coord; i<=sites->end_coord; i++) {
            sites_mapping->all_sites.push_back(i);
            sites_mapping->all_sites_start.push_back(i);
            sites_mapping->all_sites_end.push_back(i);
        }

        for (int i=0; i<ncols; i++) {
            int col = sites->positions[i];
            sites_mapping->old_sites.push_back(col);
            sites_mapping->new_sites.push_back(col - sites->start_coord);
        }

        // record new coords
        sites_mapping->new_start = 0;
        sites_mapping->new_end = sites->length();
        return true;
    }

    // iterate through variant sites
    for (int i=0; i<ncols; i++) {
        int col = sites->positions[i];

        // find next block with variant site
        while (col >= next_block) {
            sites_mapping->all_sites.push_back(next_block - half_block);
            next_block += compress;
            blocki++;
        }

        // record variant site.
        sites_mapping->old_sites.push_back(col);
        sites_mapping->new_sites.push_back(blocki);
        sites_mapping->all_sites.push_back(col);
        next_block += compress;
        blocki++;

        // each original site should be unique
        const int n = sites_mapping->all_sites.size();
        if (n > 1)
            assert(sites_mapping->all_sites[n-1] !=
                   sites_mapping->all_sites[n-2]);

        // Check whether compression is not possible
        if (next_block - compress > sites->end_coord && i != (ncols-1))
            return false;
    }

    // record non-variants at end of alignment
    while (sites->end_coord >= next_block) {
        sites_mapping->all_sites.push_back(next_block - half_block);
        next_block += compress;
        blocki++;
    }

    // record new coords
    sites_mapping->new_start = 0;
    int new_end = sites->length() / compress;
    if (ncols > 0)
        sites_mapping->new_end = max(sites_mapping->new_sites[ncols-1]+1,
                                     new_end);
    else
        sites_mapping->new_end = new_end;

    sites_mapping->all_sites_start.clear();
    sites_mapping->all_sites_end.clear();
    sites_mapping->all_sites_start.push_back(sites_mapping->old_start);
    for (unsigned int i=0; i < sites_mapping->all_sites.size()-1; i++) {
        int pos = ( sites_mapping->all_sites[i] +
                    sites_mapping->all_sites[i+1] ) / 2;
        sites_mapping->all_sites_end.push_back(pos);
        sites_mapping->all_sites_start.push_back(pos+1);
    }
    sites_mapping->all_sites_end.push_back(sites_mapping->old_end);

    return true;
}


// Apply compression using sites_mapping.
void compress_sites(Sites *sites, const SitesMapping *sites_mapping)
{
    const int ncols = sites->cols.size();
    sites->start_coord = sites_mapping->new_start;
    sites->end_coord = sites_mapping->new_end;

    for (int i=0; i<ncols; i++)
        sites->positions[i] = sites_mapping->new_sites[i];
}


// Uncompress sites using sites_mapping.
void uncompress_sites(Sites *sites, const SitesMapping *sites_mapping)
{
    const int ncols = sites->cols.size();
    sites->start_coord = sites_mapping->old_start;
    sites->end_coord = sites_mapping->old_end;

    for (int i=0; i<ncols; i++)
        sites->positions[i] = sites_mapping->uncompress(sites->positions[i]);
}


TrackNullValue get_n_regions(const Sites &sites, int numN) {
    TrackNullValue track;
    int numhap = sites.get_num_seqs();
    for (int i=0; i < sites.get_num_sites(); i++) {
        const char *cols = sites.cols[i];
        int count=0;
        for (int j=0; j < numhap; j++) {
            if (cols[j] == 'N') {
                count++;
                if (count >= numN) break;
            }
        }
        if (count >= numN) {
            track.push_back(RegionNullValue(sites.chrom, sites.positions[i],
                                            sites.positions[i]+1, ' '));
        }
    }
    track.merge();
    return track;
}

TrackNullValue get_snp_clusters(const Sites &sites, int numsnp, int window) {
    TrackNullValue track;
    if (numsnp < 0 || numsnp > window) return track;
    if (numsnp == 0) {  // mask everything
        track.push_back(RegionNullValue(sites.chrom, sites.start_coord,
                                        sites.end_coord, ' '));
        return track;
    }
    int numsite = sites.get_num_sites();
    bool isSnp[numsite];

    for (int i=0; i < numsite; i++)
        isSnp[i] = sites.is_snp(i);

    for (int i=0; i < numsite; i++) {
        if (!isSnp[i]) continue;
        int count=1;
        int start_pos = sites.positions[i];
        int j=i+1;
        for ( ; j < numsite; j++) {
            if (sites.positions[j] - start_pos + 1 > window) break;
            if (isSnp[j]) count++;
        }
        j--;
        if (count >= numsnp) {
            int end_pos = sites.positions[j];
            assert(start_pos <= end_pos);
            end_pos++;  // sites positions are zero-based; want end non-inclusive
            int diff = window - (end_pos - start_pos);
            track.push_back(RegionNullValue(sites.chrom, start_pos - diff,
                                            end_pos + diff, ' '));
        }
    }
    track.merge();
    return track;
}


// This is used to help with compression. We do not want sites without any
// known variants to inhibit ability to compress. So first use this function
// to remove N's from individuals but create a map of which regions were
// Ns in each. Then after compression can reapply this mask to each
// individual. The alleles should be changed to the allele of some other
// non-N individual at the same site so that the site is only variant
// if there is an actual known variant.
TrackNullValue unmask_ind(Sites *sites, int ind) {
    TrackNullValue masked;
    int numind = sites->get_num_seqs();
    for (int i=0; i < sites->get_num_sites(); i++) {
        if (sites->cols[i][ind] == 'N') {
            bool variant=false;
            char c='N';
            for (int j=0; j < numind; j++) {
                if (sites->cols[i][j] != 'N') {
                    if (c == 'N')
                        c = sites->cols[i][j];
                    else if (c != sites->cols[i][j]) {
                        variant=true;
                    }
                }
            }
            if (!variant) {
                sites->cols[i][ind] = c;
                masked.push_back(RegionNullValue(sites->chrom,
                                                 sites->positions[i],
                                                 sites->positions[i]+1, ' '));
            }
        }
    }
    masked.merge();
    return masked;
}

void unmask_inds(Sites *sites, vector<TrackNullValue> *masked_regions) {
    if (masked_regions->size() == 0) {
        for (int i=0; i < sites->get_num_seqs(); i++)
            masked_regions->push_back(TrackNullValue());
    }
    for (int i=0; i < sites->get_num_seqs(); i++) {
        TrackNullValue curr_mask = unmask_ind(sites, i);
        (*masked_regions)[i].merge_tracks(curr_mask);
    }
}


void print_masked_regions(const Sequences &sequences,
                          const SitesMapping *sites_mapping,
                          string chrom,
                          string filename) {
    TrackNullValue masked_regions = sequences.get_masked_regions(chrom,
                                                                 sites_mapping);
    masked_regions.write_track_regions(filename);
}


void print_masked_sites_regions(const Sites &sites, string filename) {
    TrackNullValue masked_regions = sites.get_masked_regions();
    masked_regions.write_track_regions(filename);
}


//=============================================================================
// assert functions


// Returns True if sequences pass all sanity checks.
bool check_sequences(Sequences *seqs)
{
    return check_sequences(
        seqs->get_num_seqs(), seqs->length(), seqs->get_seqs()) &&
        check_seq_names(seqs);
}


// Ensures that all characters in the alignment are sensible.
bool check_sequences(int nseqs, int seqlen, char **seqs)
{
    // check seqs
    for (int i=0; i<nseqs; i++) {
        for (int j=0; j<seqlen; j++) {
            char x = seqs[i][j];
            if (strchr("NnRrYyWwSsKkMmBbDdHhVv", x))
                // treat Ns as gaps
                x = '-';
            if (x != '-' && dna2int[(int) (unsigned char) x] == -1) {
                // an unknown character is in the alignment
                printError("unknown char '%c' (char code %d)\n", x, x);
                return false;
            }
        }
    }

    return true;
}


// Return True if all gene names are valid.
bool check_seq_names(Sequences *seqs)
{
    for (unsigned int i=0; i<seqs->names.size(); i++) {
        if (!check_seq_name(seqs->names[i].c_str())) {
            printError("sequence name has illegal characters '%s'",
                       seqs->names[i].c_str());
            return false;
        }
    }

    return true;
}


//
// A valid gene name and species name follows these rules:
//
// 1. the first and last characters of the ID are a-z A-Z 0-9 _ - .
// 2. the middle characters can be a-z A-Z 0-9 _ - . or the space character ' '.
// 3. the ID should not be purely numerical characters 0-9
// 4. the ID should be unique within a gene tree or within a species tree
//

bool check_seq_name(const char *name)
{
    int len = strlen(name);

    if (len == 0) {
        printError("name is zero length");
        return false;
    }

    // check rule 1
    if (name[0] == ' ' || name[len-1] == ' ') {
        printError("name starts or ends with a space '%c'");
        return false;
    }

    // check rule 2
    for (int i=0; i<len; i++) {
        char c = name[i];
        if (!((c >= 'a' && c <= 'z') ||
              (c >= 'A' && c <= 'Z') ||
              (c >= '0' && c <= '9') ||
              c == '_' || c == '-' || c == '.' || c == ' ')) {
            printError("name contains illegal character '%c'", c);
            return false;
        }
    }

    // check rule 3
    int containsAlpha = false;
    for (int i=0; i<len; i++) {
        if (name[i] < '0' || name[i] > '9') {
            containsAlpha = true;
            break;
        }
    }
    if (!containsAlpha) {
        printError("name is purely numeric '%s'", name);
        return false;
    }

    return true;
}


//=============================================================================
// Misc

void resample_align(Sequences *aln, Sequences *aln2)
{
    assert(aln->get_num_seqs() == aln2->get_num_seqs());
    char **seqs = aln->get_seqs();
    char **seqs2 = aln2->get_seqs();
    if (aln->base_probs.size() > 0)
        printError("Warning: resample_align not implemented for base_probs\n");

    for (int j=0; j<aln2->length(); j++) {
        // randomly choose a column (with replacement)
        int col = irand(aln->length());

        // copy column
        for (int i=0; i<aln2->get_num_seqs(); i++) {
            seqs2[i][j] = seqs[i][col];
        }
    }
}

PhaseProbs::PhaseProbs(int _hap1, int _treemap1, Sequences *_seqs,
			    const LocalTrees *trees, const ArgModel *model) {
  if (!model->unphased) return;
  hap1 = _hap1;
  treemap1 = _treemap1;
  seqs = _seqs;
  if ((int)seqs->pairs.size() != seqs->get_num_seqs()) {
    seqs->set_pairs(model);
  }
  hap2 = seqs->pairs[hap1];
  if (hap2 != -1) {
    updateTreeMap2(trees);
  } else treemap2 = -1;
  if (seqs->get_num_seqs() > 0) {
      non_singleton_snp = vector<bool>(seqs->length());
      seqs->get_non_singleton_snp(non_singleton_snp);
  }
}



void PhaseProbs::updateTreeMap2(const LocalTrees *tree) {
    for (unsigned int i=0; i < tree->seqids.size(); i++) {
	if (tree->seqids[i] == hap2) {
	    treemap2 = i;
	    return;
	}
    }
    treemap2 = -1;
    return;
}


//=============================================================================
// C interface
extern "C" {


Sites *arghmm_read_sites(const char *filename,
                         int subregion_start=-1, int subregion_end=-1)
{
    Sites *sites = new Sites();
    bool results = read_sites(filename, sites, subregion_start, subregion_end);
    if (!results) {
        delete sites;
        return NULL;
    }
    return sites;
}


void arghmm_delete_sites(Sites *sites)
{
    delete sites;
}


} // extern "C"
} // namespace argweaver
