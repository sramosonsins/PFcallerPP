#ifndef PSCALLER_H
#define	PSCALLER_H

#include "common.h"
#include "main.h"
#include "ran1.h"
#include "PSstats.h"
#include "PSparse.h"

#if noMacOS==1
#include <malloc.h>
#endif


#ifdef	__cplusplus
extern "C" {
#endif
    
    void acceptance(struct MCMCsample *current, struct MCMCsample *next);
    void doMCMCiterbial(struct MCMCsample *current, double *lnfact, double ***lncombin,double ***combin, double *lnprior, struct featpos *featp, unsigned long max_cov, double theta, long int *seed, unsigned long pmc, unsigned long poolsize, double pe1, double pe2, double pe3, double pe4, double pee, int priortype,int outg, unsigned long it, unsigned long iter);
    int PScallerMCMCbial(double *lnfact, double ***lncombin,double ***combin, double *lnprior, struct featpos *featp, unsigned long iter, unsigned long max_cov, long int samplemhmcmc,long int burnin_par, double theta, long int *seed, unsigned long pmc, unsigned long poolsize/*,double omega_obs*/, double mean_ml_err, int priortype, int outg);
    void anp(unsigned long max, double *an_m, int priortype);
    void ComputePrior(double *lnprior,unsigned long max, double *an_m, int priortype, double theta, int outg, double divergence);
    double ComputePriorFreq(unsigned long j, int priortype, double theta, unsigned long poolsize, int outg);
    void find_common_bases(struct MCMCsample *post_together,double threshlik,unsigned long maxthin,int *cp1, int *cp2, int *cp3, int *cp4,int *cf1,int *cf2,int *cer1,int *cer2,int *cer3,int *cer4);
    int compare_(const void *i,const void *j);
    void newParamProposal_bial(struct MCMCsample *current, struct MCMCsample *next,double ***combin, struct Obs observed, double *lnfact,double *lnprior, double random0, double random1, double random2, double random3, double random4, double random5, /*double random6,*/ long max_cov, double theta, unsigned long pmc, unsigned long poolsize, double pe1, double pe2, double pe3, double pe4, double pee, int priortype, int outg, unsigned long it, unsigned long iter, double *lnprob);
    double calculate_prob_bial(unsigned long obs_p1, unsigned long obs_p2, unsigned long obs_p3, unsigned char obs_outg, unsigned long exp_p1, unsigned long exp_p2,unsigned long exp_frp1, unsigned long exp_frp2, unsigned long exp_popfr, unsigned long exp_ner1,unsigned long exp_ner2,unsigned long exp_ner3,unsigned long exp_ner4,double ***lncombin,double ***combin,double *lnprior,double *lnfact, unsigned long max_cov, double theta, unsigned long pmc, unsigned long poolsize, double pe1, double pe2, double pee, int priortype, int outg);
#ifdef	__cplusplus
}
#endif

#endif	/* PSCALLER_H */
