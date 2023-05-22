#ifndef PSSTATS_H
#define	PSSTATS_H

#include "common.h"
#include "main.h"
#include "ran1.h"
#include "PScaller.h"
#include "PSparse.h"

#if noMacOS==1
#include <malloc.h>
#endif


#ifdef	__cplusplus
extern "C" {
#endif
    
    double ln(double x);
    void lnfactorial(double *lnfact, unsigned long pmc2);
    double lnfactorial1(unsigned long j);
    void stirling2(long double **s2, unsigned long max_cov, unsigned long poolsize);
    double mind(double x,double y);
    double maxd(double x,double y);
    unsigned long minl(unsigned long x,unsigned long y);
    double minint(unsigned long int x,unsigned long int y);
    double maxint(int x,int y);
    void prob_nrb_st(double **combin, double **lncombin, double *lnfact,long double **s2,unsigned long pR, unsigned long max_cov);
    void prob_nrb_mc(double **combin, double **lncombin, unsigned long nitermc,unsigned long pR, unsigned long max_cov, long int *seed, unsigned long poolsize, unsigned long pmc);
    double lnmultinomialfr(unsigned long realnR, unsigned long realnA, unsigned long realnE, double pR, double pe1, double pe2, double pee, double *lnfact);
    double lnbinomial(unsigned long realnR, unsigned long realnA, double pR, double *lnfact);
    double lnmultinomial4(unsigned long n1,unsigned long n2, unsigned long n3,unsigned long n4,double p1,double p2,double p3,double *lnfact);
    double lnmultinomial3(unsigned long n1,unsigned long n2, unsigned long n3,double p1,double p2,double *lnfact);
    double lnbinomial_seqbial(unsigned long n1, unsigned long n2, unsigned long n3, unsigned long obs_p1, unsigned long obs_p2, unsigned long obs_p3, double pe1, double pe2, double pee, double *lnfact);
    void binomialner_seqbial_sample(unsigned long obs_p1, unsigned long obs_p2, unsigned long obs_p3, double pe1, double pe2, double pee,  double *lnfact, struct MCMCsample *next, double random0, double random1, double random2, unsigned long pmc, double theta, double *lnprob);
    unsigned long binomial_seqbial_sample(unsigned long nc, double freq, unsigned long min, double random4, double *lnfact, double *lnprob);
    void multinomialfr_seqbial_sample(unsigned long obs_p1, unsigned long obs_p2, unsigned long obs_p3, unsigned long obs_p4 ,unsigned char obs_outg,double pe1, double pe2,double pe3, double pe4, double pee,double *lnprior,double *lnfact, struct MCMCsample *mcmc_par,unsigned long pmc, int outg,double *lnpval,double *sumpval);
    unsigned long combin_seq_sample(unsigned long obs_p1, unsigned long exp_frp1, unsigned long exp_ner1, unsigned long exp_ner2, unsigned long exp_ner3,  double ***combin, double random1, unsigned long poolsize, int diffnt, double *lnprob);
    double error0_calc(double nr, double e, double ***combin, double theta, unsigned long pmc);
    double maxval(double x,double y);
    int compunsign(const void *i, const void *j);
    void freq_prior(double *lnprior, struct MCMCsample *next, double random6);
    
#ifdef	__cplusplus
}
#endif

#endif	/* PSSTATS_H */
