
#include "PSstats.h"

/*///////////////////////////////////////*/
double ln(double x){
  if (x==0){
    return (-1e+300);
  }else{
      if(x<0) {
          printf("\nERROR: Negative log.\n");
          return (-1e+300);
      }
  }
  return (log(x));
}
/*/////////////////////////////////////*/
double lnmultinomialfr(unsigned long realnR, unsigned long realnA, unsigned long realnE, double pR, double pe1, double pe2, double pee, double *lnfact)
{
    double lnbin=0.;
    lnbin = lnfact[realnA+realnR+realnE] - lnfact[realnR] - lnfact[realnA] - lnfact[realnE] + realnR*ln(pR*(1.0-pe2-pee)+(1.0-pR)*pe1) + realnA*ln((1.0-pR)*(1.0-pe1-pee)+pR*pe2);
    if(realnE > 0) {
        lnbin += realnE * ln(pee);
    }
    return (lnbin);
}
/*/////////////////////////////////////*/
double lnbinomial(unsigned long realnR, unsigned long realnA, double pR, double *lnfact)
{
    double lnbin=0.;
    lnbin = lnfact[realnA+realnR] - lnfact[realnR] - lnfact[realnA] + realnR*ln(pR) + realnA*ln(1.0-pR);
    return (lnbin);
}
/*///////////////////////////////////////*/
double lnbinomial_seqbial(unsigned long n1, unsigned long n2, unsigned long n3, unsigned long obs_p1, unsigned long obs_p2, unsigned long obs_p3, double peR, double peA, double peE, double *lnfact) {
    
    double lnprob1,lnprob2,lnprob3,lnprob;
    double lnpow1,lnpow2,lnpow3,lnpow4,lnpow5,lnpow6;
    unsigned long nobs;
    
    nobs = obs_p1+obs_p2+obs_p3;
    if(nobs==0)
        return -1E+300;
    
    lnpow1  = (double)n1 * (double)ln(peR);
    lnpow2  = (double)(obs_p1-n1) * (double)ln(1.0-peR);
    lnpow3  = (double)n2 * (double)ln(peA);
    lnpow4  = (double)(obs_p2-n2) * (double)ln(1.0-peA);
    lnpow5=lnpow6=0;
    if(obs_p3) {
        lnpow5  = (double)n3 * (double)ln(peE);
        lnpow6  = (double)(obs_p3-n3) * (double)ln(1.0-peE);
    }

    lnprob1 = lnfact[obs_p1] - lnfact[n1] - lnfact[obs_p1-n1] + lnpow1 + lnpow2;
    lnprob2 = lnfact[obs_p2] - lnfact[n2] - lnfact[obs_p2-n2] + lnpow3 + lnpow4;
    lnprob3 = lnfact[obs_p3] - lnfact[n3] - lnfact[obs_p3-n3] + lnpow5 + lnpow6;

    lnprob =  lnprob1 + lnprob2 + lnprob3;
    
    return(lnprob);
}
/*///////////////////////////////////////*/
double lnmultinomial4(unsigned long n1,unsigned long n2, unsigned long n3, unsigned long n4, double p1,double p2,double p3,double *lnfact)
{
    double lnprob;
    double lnpow1,lnpow2,lnpow3,lnpow4;
    
    if(n1+n2+n3+n4==0)
        return -1E+300;
    
    lnpow1 = lnpow2 = lnpow3 = lnpow4 = 0.0;
    
    if(n1 != 0) lnpow1  = (double)n1 * (double)ln(p1);
    if(n2 != 0) lnpow2  = (double)n2 * (double)ln(p2);
    if(n3 != 0) lnpow3  = (double)n3 * (double)ln(p3);
    if(n4 != 0) lnpow4  = (double)n4 * (double)ln(1.0-(p1+p2+p3));
    
    
    lnprob =  lnfact[n1+n2+n3+n4] -
    lnfact[n1] - lnfact[n2] - lnfact[n3] - lnfact[n4] +
    lnpow1 + lnpow2 + lnpow3 + lnpow4;
    
    return(lnprob);
}
/*///////////////////////////////////////*/
double lnmultinomial3(unsigned long n1,unsigned long n2, unsigned long n3,double p1,double p2,double *lnfact)
{
    double lnprob;
    double lnpow1,lnpow2,lnpow3;
    
    if(n1+n2+n3==0)
        return -1E+300;
    
    lnpow1  = (double)n1 * (double)ln(p1);
    lnpow2  = (double)n2 * (double)ln(p2);
    if(n3 > 0) lnpow3  = (double)n3 * (double)ln(1.0-(p1+p2));
    else lnpow3 = 0.;
    
    
    lnprob =  lnfact[n1+n2+n3] -
    lnfact[n1] - lnfact[n2] - lnfact[n3] +
    lnpow1 + lnpow2 + lnpow3;
    
    return(lnprob);
}
/*/////////////////////////////////////*/    /*...20230905...*/
void multinomialfr_seqbial_folded(unsigned long obs_p1, unsigned long obs_p2, unsigned long obs_p3,
                                  unsigned long obs_p4, unsigned char obs_outg,double pe1, double pe2, double pe3, double pe4,
                                  double pee, double *lnprior,double *lnfact, struct MCMCsample *mcmc_par,
                                  unsigned long pmc, int outg, double *lnpval,double *sumpval, double ****pbfold) {
    
    unsigned long i,k,j;
    double lnpar,lnp,popf;
    double sumtot/*,pR,pA,pE*/;
    
    *sumpval = 0.;
    for(i=0;i<pmc+1;i++) {
        lnpar = pbfold[minl(pmc-i,i)][obs_p2][obs_p3][obs_p4] ;
        lnpval[i] = 0.;
        for(k=0;k<NPOPINTV+1;k++) { /*pop frequency (popfr) used as a non-inferred parameter*/ /*CHANGE AND INCLUDE IN VCF?, then add new dimension to lnpval?*/
            lnp = lnbinomial(i,pmc-i,(double)k/(double)NPOPINTV,lnfact);
            if(outg){
                if(obs_outg==0) j=k;
                else j=NPOPINTV-k;
            }
            else j = minint(k,NPOPINTV-k);
            popf = lnprior[j];
            sumtot = exp(lnpar + lnp + popf);
            
            *sumpval += (sumtot);
            lnpval[i] += (sumtot);
        }
        lnpval[i] = log(lnpval[i]);
        
        mcmc_par[i].frp1 = i;
        mcmc_par[i].frp2 = pmc - mcmc_par[i].frp1;
        mcmc_par[i].popfr = 0;
    }

    return;
}
/*/////////////////////////////////////*/    /*...20230905...*/
void multinomialfr_seqbial_sample(unsigned long obs_p1, unsigned long obs_p2, unsigned long obs_p3,
                                  unsigned long obs_p4, unsigned char obs_outg,double pe1, double pe2,double pe3, double pe4,
                                  double pee, double *lnprior,double *lnfact, struct MCMCsample *mcmc_par,
                                  unsigned long pmc, int outg, double *lnpval,double *sumpval) {
    
    unsigned long i,k,j;
    double lnpar,lnp,popf;
    double sumtot,pR,pA/*,pE*/;
    
    *sumpval = 0.;
    for(i=0;i<pmc+1;i++) {
        pR = (double)i/pmc * (1.0-pe2-pee) + (1.0-(double)i/pmc)*pe1;
        pA = (1.0-(double)i/pmc) * (1.0-pe1-pee) + ((double)i/pmc)*pe2;
        //pE = pee;
        lnpar = lnmultinomial3(obs_p1,obs_p2,obs_p3,pR,pA,lnfact); /*to modify, including sorting in a value precalculated*/
        lnpval[i] = 0.;
        for(k=0;k<NPOPINTV+1;k++) { /*pop frequency (popfr) used as a non-inferred parameter*/ /*CHANGE AND INCLUDE IN VCF?, then add new dimension to lnpval?*/
            lnp = lnbinomial(i,pmc-i,(double)k/(double)NPOPINTV,lnfact);
            if(outg){
                if(obs_outg==0) j=k;
                else j=NPOPINTV-k;
            }
            else j = minint(k,NPOPINTV-k);
            popf = lnprior[j];
            sumtot = exp(lnpar + lnp + popf);
            
            *sumpval += (sumtot);
            lnpval[i] += (sumtot);
        }
        lnpval[i] = log(lnpval[i]);
        
        mcmc_par[i].frp1 = i;
        mcmc_par[i].frp2 = pmc - mcmc_par[i].frp1;
        mcmc_par[i].popfr = 0;
    }

    return;
}
/*///////////////////////////////////////*/
void binomialner_seqbial_sample(unsigned long obs_p1, unsigned long obs_p2, unsigned long obs_p3, double pe1, double pe2, double pee, double *lnfact, struct MCMCsample *next, double random0, double random1, double random2, unsigned long pmc, double theta, double *lnprob){
    
    unsigned long i,j,k;
    double pval0,pval1,pval2;
    double prob1,prob2;
    double pA;
    unsigned long beg,end,inc;
    
    pA = 1.0 - (double)next->frp1/(double)pmc;
    if(pmc > 1 && theta > 0.0) {
        pval0 = 0.;
        for(k=0;k<obs_p3;k++) {
            pval0 += exp(lnbinomial(k,obs_p3-k,1.0-pA,lnfact));
            if(pval0 >= random0) break;
        }
        pval0 = 0.;
        pval1 = ((1.0-pA)*pe2)/(pA*(1-pe1-pee)+(1.0-pA)*pe2);
        beg = 0; end = obs_p2; inc = 1;
        if(pval1 > 0.5) {beg = obs_p2; end = 0; inc = -1;}
        for(j=beg;j!=end;j+=inc) {
            pval0 += exp(lnbinomial(j,obs_p2-j,pval1,lnfact));
            if(pval0 >= random1) break;
        }
        pval0 = 0.;
        pval2 = (pA*pe1)/((1.0-pA)*(1-pe2-pee)+pA*pe1);
        beg = 0; end = obs_p1; inc = 1;
        if(pval2 > 0.5) {beg = obs_p1; end = 0; inc = -1;}
        for(i=beg;i!=end;i+=inc) {
            pval0 += exp(lnbinomial(i,obs_p1-i,pval2,lnfact));
            if(pval0 >= random2) break;
        }
        *lnprob = lnbinomial(k,obs_p3-k,1.0-pA,lnfact) +
                  lnbinomial(j,obs_p2-j,pval1,lnfact) +
                  lnbinomial(i,obs_p1-i,pval2,lnfact);
    }
    else {
        //if all variants come from pR, the error will be:
        pval1 = 0.; pval2 = 0.;
        prob1 = (1.0-pA)*pe2 / ((1.0-pA)*pe2 + pA*(1.0-pe1-pee));
        if(obs_p2+obs_p3) pval1 = exp(lnbinomial(obs_p2+obs_p3,0,prob1,lnfact));
        //if all variants come from pA, the error will be:
        pval0 = 0.;
        prob2 = pA*pe1 / (pA*pe1 + (1.0-pA)*(1.0-pe2-pee));
        if(obs_p1) pval0 = exp(lnbinomial(obs_p1,0,prob2,lnfact));
        //decide
        if(pval0/(pval0+pval1) >= random1) {
            //sequence error of pA -> frp1...
            k = 0;
            j = 0;
            i = obs_p1;
            *lnprob = lnbinomial(obs_p1,0,prob2,lnfact);
        }
        else {
            //sequence error of pR -> frp2
            k = obs_p3;
            j = obs_p2;
            i = 0;
            *lnprob = lnbinomial(obs_p2+obs_p3,0,prob1,lnfact);
        }
    }

    next->ner1 = i;
    next->ner2 = j;
    next->ner3 = k;
    next->ner4 = obs_p3-k;

    return;
}
/*/////////////////////////////////////*/
unsigned long binomial_seqbial_sample(unsigned long nc, double freq, unsigned long min, double random4, double *lnfact, double *lnprob)
{
    unsigned long ncA;
    double pval0,pT;
    
    if(nc > 0) {
        if(min==0) pT = 1.;
        else { //min=1;
            pT = 1. - exp(lnbinomial(0,nc-0,freq,lnfact));
        }
        pval0 = 0.;
        for(ncA=min;ncA<nc;ncA++) {
            pval0 += exp(lnbinomial(ncA,nc-ncA,freq,lnfact));
            if(pval0 >= random4) break;
        }
        *lnprob = lnbinomial(ncA,nc-ncA,freq,lnfact)-log(pT);
    }
    else {
        ncA = 0;
        *lnprob = 0;
    }
   
    return ncA;
}

unsigned long combin_seq_sample(unsigned long obs_p1, unsigned long exp_frp1, unsigned long exp_ner1, unsigned long exp_ner2, unsigned long exp_ner3,  double ***combin, double random1, unsigned long poolsize, int difnt, double *lnprob) {
    
    double pval0;
    unsigned long nc;
    double sumtot;
    unsigned long minrange;
    
    if(poolsize > LIMIT_COMB) {
        nc = (unsigned long) ceil(exp_frp1 * (1.0 - pow(1.0 - 1.0/exp_frp1,obs_p1-exp_ner1+exp_ner2+exp_ner3)));
        return(nc);
    }
    
    minrange = minl(exp_frp1,obs_p1-exp_ner1+exp_ner2+exp_ner3);
    sumtot = 0.;
    for(nc=difnt;nc<=minrange;nc++) {
        sumtot += combin[exp_frp1][nc][obs_p1-exp_ner1+exp_ner2+exp_ner3];
    }

    pval0 = 0.;
    for(nc=difnt;nc<minrange;nc++) {
        pval0 += combin[exp_frp1][nc][obs_p1-exp_ner1+exp_ner2+exp_ner3]/sumtot;
        if(pval0 >= random1) break;
    }
    if(nc>minrange) nc=minrange; /*in the strange case random1=1?*/
    
    *lnprob = ln(combin[exp_frp1][nc][obs_p1-exp_ner1+exp_ner2+exp_ner3]);
    
    return(nc);
}
/*///////////////////////////////////*/
void stirling2(long double **s2, unsigned long max_cov, unsigned long poolsize){
  unsigned long  nr,nc;
  for (nr=0;nr<=max_cov;nr++){
    for (nc=0;nc<=poolsize;nc++){
      /*//fill first row/col, diagonal, and upper matrix*/
      if ((nr == 0 && nc == 0) || (nr>0 && nr == nc)){s2[nr][nc]=1;}
      else if (nc == 0 || nr < nc){s2[nr][nc]=0;}
      else{ /*//fill lower matrix*/
        long double temp1=s2[nr-1][nc]*nc;
        long double temp2=s2[nr-1][nc-1];
        s2[nr][nc]=temp1+temp2;
      }
    }
  }
}
/*///////////////////////////////////////*/
void lnfactorial(double *lnfact, unsigned long pmc2){
  unsigned long j;
  double lf;
  
  for (j = 0; j <= pmc2; j++) {
    if (j==0){
      lf=0;
    }else if (j<80){
      lf=lnfact[j-1]+ln(j);
    }else{
      lf=j*ln(j)-j;/*Stirling approx for numerical stability and to speed up computation*/
    }
    lnfact[j] = lf;
  }
}
/*///////////////////////////////////////*/
double lnfactorial1(unsigned long j){
    return(j*ln(j)-j);
}
/*///////////////////////////////////////*/
double mind(double x,double y) {
    if(x<y)
        return x;
    else
        return y;
}
/*///////////////////////////////////////*/
double maxd(double x,double y) {
    if(x>y)
        return x;
    else
        return y;
}
/*///////////////////////////////////////*/
unsigned long minl(unsigned long x,unsigned long y) {
    if(x<y)
        return x;
    else
        return y;
}
/*/////////////////////////////////////*/
double minint(unsigned long int x,unsigned long int y) {
    if(x<y)
        return x;
    else
        return y;
}
/*/////////////////////////////////////*/
double maxint(int x,int y) {
    if(x>y)
        return x;
    else
        return y;
}
/*///////////////////////////////////////*/
double deltak(unsigned long i, unsigned long j) {
    if(i==j) return(1);
    return(0);
}
/*///////////////////////////////////////*/
void prob_nrb_st(double **combin, double **lncombin, double *lnfact,long double **s2,unsigned long pR, unsigned long max_cov){
    
    unsigned long nr,nc0;
    double ll,lll;
    
    for (nr=0;nr<=max_cov;nr++){
        lll=-1e+300;
        for (nc0=0;nc0<=pR;nc0++){
            if(pR == 0 && nc0==0 && nr==0) ll = 0.0;
            else if (nr == 0) {
                if (nc0 == 0) ll = 0.0;
                else ll = -1e+300;
            }
            else if (nr < nc0) {
                ll = -1e+300;
            }
            else { /*//fill lower matrix*/
                if(pR == 0 || nc0 == 0) ll = lll;
                else ll = lnfact[pR]-lnfact[pR-nc0]+ln(s2[nr][nc0])-(double)nr*ln(pR);
            }
            if(ll < -1e+100) ll = lll;
            if(ll > 0. && ll < 0.01) ll = 0.;
            lncombin[nc0][nr]= ll;
            combin[nc0][nr]=exp(ll);
            //printf("%.3e\t",combin[nc0][nr]);
        }
        //printf("\n");
    }
    return;
}
/*///////////////////////////////////////*/
void prob_nrb_mc2(double **combin, double **lncombin,unsigned long nitermc,unsigned long pR, unsigned long max_cov, long int *seed, unsigned long poolsize, unsigned long pmc){
    unsigned long nr,nr2,nc0,it,it1;
    double ll,sum;
    unsigned long *numberc;
    unsigned long *nunique;
    
    for (nr=0;nr<=max_cov;nr++){
        nunique =  (unsigned long *)calloc(nitermc, sizeof(unsigned long));
        for(it=0;it<nitermc;it++) {
            numberc = (unsigned long *)calloc(poolsize+1, sizeof(unsigned long));
            for(nr2=0;nr2<nr;nr2++) {
                numberc[it1=(unsigned long)(ran1(seed)*poolsize)] += 1;
                if(numberc[it1]>1) {
                    nunique[it] += 1;
                }
            }
            free(numberc);
        }
        
        sum = 0;
        for (nc0=0;nc0<=pmc;nc0++){
            if(pR == 0 && nc0==0 && nr==0) ll = 1.0;
            else if (nr == 0){
                if (nc0 == 0) ll = 1.0;
                else ll = 1E-300;
            }
            else if (nr < nc0){
                ll = 1E-300;
            }
            else if (nc0 == 0) {
                ll = 0.0;
            }
            else{ /*//fill lower matrix*/
                ll = 0;
                for(it=0;it<nitermc;it++) {
                    if(nunique[it] == nr - nc0)
                        ll += 1.0;
                }
                ll = ll/(double)nitermc;
                sum += ll;
            }
            if(ll > 0.) {
                combin[nc0][nr]= ll;
                lncombin[nc0][nr]=log(ll);
                if(sum >= 1.0)
                    break;
            }
            else {
                combin[nc0][nr]= 0.;
                lncombin[nc0][nr]=-1e+300;
            }
            //printf("%.3e\t",combin[nc0][nr]);
        }
        //printf("\n");
        free(nunique);
    }
    return;
}
/*///////////////////////////////////////*/
void prob_nrb_mc(double **combin, double **lncombin,unsigned long nitermc,unsigned long pR, unsigned long max_cov, long int *seed, unsigned long poolsize, unsigned long pmc){
    unsigned long nr,nc0,it,it1;
    double ll,sum;
    unsigned long *numberc;
    unsigned long *nunique;
    
    unsigned long ppsize = (unsigned long) ((double)pR * (double)poolsize / (double) pmc);
    
    for(nr=0;nr<=max_cov;nr++){
        nunique = (unsigned long *)calloc(nitermc, sizeof(unsigned long));
        numberc = (unsigned long *)calloc(nr+1, sizeof(unsigned long));
        if(pR > 0) {
            for(it=0;it<nitermc;it++) {
                for(it1=0;it1<nr;it1++)
                    numberc[it1] = (unsigned long)floor(ran1(seed)*ppsize);
                qsort(numberc,nr,sizeof(unsigned long),compunsign);
                if(nr>0) nunique[it] = 1;
                for(it1=1;it1<nr;it1++) {
                    if(numberc[it1] != numberc[it1-1])
                        nunique[it] = nunique[it] + 1;
                }
            }
        }
        sum = 0;
        for (nc0=0;nc0<=pmc;nc0++){
            if(pR == 0 && nc0==0 && nr==0) ll = 1.0;
            else if (nr == 0){
                if (nc0 == 0) ll = 1.0;
                else ll = 1E-300;
            }
            else if (nr < nc0){
                ll = 1E-300;
            }
            else if (nc0 == 0) {
                ll = 0.0;
            }
            else{ /*//fill lower matrix*/
                ll = 0;
                for(it=0;it<nitermc;it++) {
                    if(nunique[it] == nc0)
                        ll += 1.0;
                }
                ll = ll/(double)nitermc;
                sum += ll;
            }
            /*calculate combin*/
            if(ll > 0.) {
                combin[nc0][nr]= ll;
                lncombin[nc0][nr]=log(ll);
                if(sum >= 1.0)
                    break;
            }
            else {
                combin[nc0][nr]= 0.;
                lncombin[nc0][nr]=-1e+300;
            }
            //printf("%.3e\t",combin[nc0][nr]);
        }
        //printf("\n");
        free(numberc);
        free(nunique);
    }
    return;
}
/*/////////////////////////////////////*/
double maxval(double x,double y) {
    return(x>y?x:y);
}
/*/////////////////////////////////////*/
int compunsign(const void *i, const void *j)
{
    if(*(unsigned long *)i < *(unsigned long *)j) return(-1);
    if(*(unsigned long *)i > *(unsigned long *)j) return(+1);
    return 0;
}
/*/////////////////////////////////////*/


