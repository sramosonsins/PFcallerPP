#ifndef PSPARSE_H
#define	PSPARSE_H

#include "common.h"
#include "main.h"
#include "ran1.h"
#include "PSstats.h"
#include "PScaller.h"
#include "getdelim.h"
#include "zutil.h"
#include "zindex.h"
#include "progressbar.h"
#include "cdflib.h"

#if noMacOS==1
#include <malloc.h>
#endif

#ifdef	__cplusplus
extern "C" {
#endif
    
    void printfasta(char *file,char **seqs, char *chr_name, unsigned long pmc,unsigned long lentotal);
    void printtfasta(FILE *fasta_file, SGZip *fasta_file_gz, char *chr_name,unsigned long first, unsigned long lentotal, unsigned long pmc, struct typepos *typep, struct featpos *featp, unsigned long min_qual, unsigned long iter, unsigned long cnt2, double theta, long int *seed, unsigned long poolsize, struct list_chr_pileup *chr_data, unsigned long num_scaffolds_pileup, int argc, char **argv);
    void printlikelihoods(char *file,struct typepos *typep, struct featpos *featp,char *chr_name,unsigned long first, unsigned long iter, unsigned long cnt, unsigned long cnt2, unsigned long pmc, int biallelic);
    void printgVCF(char *file, char *chr_name, struct typepos *typep, struct featpos *featp,unsigned long first,long int lentotal,unsigned long min_qual, unsigned long cnt2, double theta, unsigned long pmc, unsigned long poolsize, unsigned long iter, long int *seed, struct list_chr_pileup *chr_data, unsigned long num_scaffolds_pileup);
    void Freqs2fasta(struct typepos *typep, struct featpos *featp, char **seqs, char *prfasta,char *chr_name, unsigned long min_qual, unsigned long iter, unsigned long lentotal, unsigned long cnt2, double theta, long int *seed, unsigned long pmc, unsigned long poolsize, struct list_chr_pileup *chr_data, unsigned long num_scaffolds_pileup);
    void Rawdata2fasta(struct typepos *typep, struct featpos *featp, char **seqs, char *prfasta,char *chr_name, unsigned long min_qual, unsigned long iter, unsigned long lentotal, unsigned long cnt2, double theta, long int *seed, unsigned long pmc, unsigned long poolsize, struct list_chr_pileup *chr_data, unsigned long num_scaffolds_pileup,unsigned long max_cov);
    int  read_line_pileup(FILE * bam_file, SGZip * bam_file_gz, unsigned long * cov, unsigned long * pos_base, unsigned long *freq, char cc, float *error_mean,char *cchrom, unsigned long max_cov, unsigned long min_cov, unsigned long min_qual, double min_error, double *noerror_prod, unsigned long nscaffolds, char **chr_name_array, long int *number_scaffold, unsigned char *crefbase, int teststrands);
    int test_strands(int *frqb, double thr);
    int base_to_num(int base);
    float Bquality_error(float perr);
    int base_to_num_char(int base,char *alt);
    void sortObs (struct Obs *tmp, struct Obs *sorted, char nucl[5], int outg);
    void ReduceInput(char *bampath, struct typepos **typep,struct featpos **featp, unsigned long max_cov, unsigned long min_cov, unsigned long min_qual, long int *ntotalreads, long int *ntotalreads_novariant, long int *npositions_novariant, long int *npositions_variant, double *n_nr, double min_error, double *sum_errors, double *sum_errors_variant, double *npositions_noerror_exp, unsigned long iter, unsigned long *cnt, unsigned long *cnt2, int biallelic, float qerr_range, struct list_chr_pileup **chr_data, unsigned long *num_scaffolds_pileup, unsigned long nscaffolds, char **chr_name_array, char **chr_length_array, double *sum_per, long int *npositions_per, double *thetant, unsigned long *thata_pos, unsigned long p,int outg, double *divergence, unsigned long *div_pos,int teststrands);
    void usage(void);
    int check_param(int argc, char **argv, unsigned long *max_cov, unsigned long *min_cov, unsigned long *min_qual, char *chr_length_all, char *chr_name_all, int *biallelic, char *file_in, char *file_out, int *priortype, long int *samplemhmcmc,long int *burnin_par, unsigned long *iter, unsigned long *nitermc, char *prfasta, unsigned long *lentotal, double *theta, long int *seed, unsigned long *poolsize, float *qerr_range, int *outg, int *teststrands);
    void function_gvcf_results_null(unsigned long jj1, FILE *gvcf_file, SGZip *gvcf_file_gz, char *seqs, char *chr_name, unsigned long pmc);
    void function_gvcf_results(unsigned long jj1, struct list_chr_pileup *chr_data, struct typepos *typep, struct featpos *featp, FILE *gvcf_file, SGZip *gvcf_file_gz, char *seqs, char *chr_name, unsigned long cnt2, unsigned long pmc, unsigned long iter, long int *seed, long int lentotal);
    int read_index_file(char *chr_name_all, unsigned long *nscaffolds,char ***chr_name_array,char ***chr_length_array);
#ifdef	__cplusplus
}
#endif

#endif	/* PSPARSE_H */
