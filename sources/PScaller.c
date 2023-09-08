
#include "PScaller.h"
#include "main.h"

/*///////////////////////////////////////*/
void acceptance (struct MCMCsample *current, struct MCMCsample *next){
    current->nc=next->nc;
    current->p1=next->p1;
    current->p2=next->p2;
    current->p3=next->p3;
    current->p4=next->p4;
    current->frp1=next->frp1;
    current->frp2=next->frp2;
    current->popfr=next->popfr;
    current->ner1=next->ner1;
    current->ner2=next->ner2;
    current->ner3=next->ner3;
    current->ner4=next->ner4;
    current->lnprob=next->lnprob;
}
/*///////////////////////////////////////*/
void doMCMCiterbial(struct MCMCsample *current, double *lnfact, double ***lncombin,double ***combin, double *lnprior, struct featpos *featp, unsigned long max_cov, double theta, long int *seed, unsigned long pmc, unsigned long poolsize, double pe1, double pe2, double pe3, double pe4, double pee, int priortype, int outg, unsigned long it, unsigned long iter, double ****pbfold){
    
    struct MCMCsample next;
    double random0,random1,random2,random3,random4,random5/*,random6*/;
    double nq,dq;
    double ratio, random,ratio_p;
    double lnprob;
    
    random0=ran1(seed);
    random1=ran1(seed);
    random2=ran1(seed);
    random3=ran1(seed);
    random4=ran1(seed);
    random5=ran1(seed);
    //random6=ran1(seed);
    
    newParamProposal_bial(current,&next,combin,featp->observed,lnfact,lnprior,
                          random0,random1,random2,random3,random4,random5/*,random6*/,
                          max_cov, theta,pmc,poolsize, pe1, pe2, pe3, pe4, pee, priortype,
                          outg, it, iter, &lnprob, pbfold);
    next.lnprob = lnprob;
    
    /*
     next.lnprob = calculate_prob_bial(featp->observed.freq[0],featp->observed.freq[1],
     featp->observed.freq[2]+featp->observed.freq[3],
     featp->observed.outg,
     next.p1, next.p2,
     next.frp1,next.frp2,next.popfr,
     next.ner1,next.ner2,next.ner3,next.ner4,
     lncombin,combin,lnprior,lnfact, max_cov,
     theta,pmc,poolsize, pe1, pe2, pee, priortype, outg);
    */
    /**/
    if(next.lnprob > -1e+299) {
        dq = nq = 0; /*metropolis algorithm. Less efficient than mhmcmc but simple*/
        ratio = next.lnprob - current->lnprob  + nq - dq ;
        
        if (ratio>=0){
            acceptance(current,&next);
        }else{
            random = ran1(seed);
            ratio_p = exp(ratio);
            if (ratio_p > random)
                acceptance(current,&next);
        }
    }
    /**/
    //acceptance(current,&next);
}
/*///////////////////////////////////////*/
int PScallerMCMCbial(double *lnfact, double ***lncombin, double ***combin,double *lnprior, struct featpos *featp, unsigned long iter, unsigned long max_cov, long int samplemhmcmc, long int burnin_par, double theta,long int *seed, unsigned long pmc, unsigned long poolsize, double mean_ml_error, int priortype, int outg){
    
    /*MCMC behaviour, I just do a few iterations to approx the variance*/
    unsigned long thinin=samplemhmcmc; /*arbitrary, just to have an independent sampled distribution, works well?*/
    unsigned long MCMCiter= iter*thinin + burnin_par*thinin; /*1.0/4.0 * (iter*thinin); *//*667;*//*500*/
    unsigned long burnin= burnin_par*thinin; /*(unsigned long) 1.0/4.0 * (iter*thinin)*/ /*MCMCiter/4.0*/;
    double cnt_samples=0; /*count the number of samples from posterior*/
    int done=0;
    
    /*posteriors*/
    unsigned long i,j,k;
    unsigned long *post_nc;
    unsigned long *post_p1, *post_p2;
    unsigned long *post_ner1, *post_ner2, *post_ner3, *post_ner4;
    unsigned long *post_frp1, *post_frp2;
    float *post_popf;
    post_nc = (unsigned long *)calloc(pmc+1,sizeof(unsigned long));
    post_p1 = (unsigned long *)calloc(pmc+1,sizeof(unsigned long));
    post_p2 = (unsigned long *)calloc(pmc+1,sizeof(unsigned long));
    post_ner1 = (unsigned long *)calloc(max_cov+1,sizeof(unsigned long));
    post_ner2 = (unsigned long *)calloc(max_cov+1,sizeof(unsigned long));
    post_ner3 = (unsigned long *)calloc(max_cov+1,sizeof(unsigned long));
    post_ner4 = (unsigned long *)calloc(max_cov+1,sizeof(unsigned long));
    post_frp1 = (unsigned long *)calloc(pmc+1,sizeof(unsigned long));
    post_frp2 = (unsigned long *)calloc(pmc+1,sizeof(unsigned long));
    post_popf = (float *)calloc(NPOPINTV+1,sizeof(float));

    struct MCMCsample *post_together;
    unsigned long maxthin = floor((iter*thinin)/(double)thinin);
    post_together = (struct MCMCsample *)calloc(maxthin,sizeof(struct MCMCsample));
    
    /*MCMC initialization*/
    struct MCMCsample current;
    double pe1,pe2,pee,pe3,pe4;
    double lnprob;
    unsigned long nr;

    nr  = featp->observed.freq[0] + featp->observed.freq[1] + featp->observed.freq[2] + featp->observed.freq[3];
    
    if(mean_ml_error == 0.0) {
        pe2 = 1e-300;
        pe3 = 1e-300;
        pe4 = 1e-300;
        pee = 1e-300;
        pe1 = featp->observed.error_mean[0]; /*P(observed1 | real n2) = phred score*/
        if(featp->observed.freq[1]>0) pe2 = featp->observed.error_mean[1];/*P(observed2 | real n1)*/
        if(featp->observed.freq[2]>0) pee = pe3 = featp->observed.error_mean[2];
        if(featp->observed.freq[3]>0) {pe4 = featp->observed.error_mean[3]; pee += pe4; }
        if(featp->observed.freq[2]>0 && featp->observed.freq[3]>0) pee /= 2.0;
        /*Raineri approach*/
        pe1 = pe2 * (1.-pe1) / (1.-pe1-pe2); /*P(real n1 | observed2) modified 28052023 (Raineri et al. 2012 approach)*/
        pe2 = pe1 * (1.-pe2) / (1.-pe1-pe2); /*P(real n2 | observed1) modified 28052023 (Raineri et al. 2012 approach)*/
    }
    else {
        pe1 = pe2 = pe3 = pe4 = pee = mean_ml_error;
    }
    
    /*///////////////////////////////////*/
    /*...20230905...*/
    /* Include here the function to calculate all pool + population frequencies before sorting reads by size: include restrictions */
    /* The array can be calculated here (before entering in next functions), because the parameters will change in any "reduced" case */
    
    double phiA, phiR,prob_bin_p;
    long int minorE1, minorE2,prE1,prE2;
    long int pr;
    unsigned long nrm[4];

    double ****pbfold;
    /*The size of dimensions is reduced to reduce memory and increase speed: ok?*/
    minorE1 = minl(featp->observed.cov*(pe1 + pe2) + ADDMHITS1, featp->observed.cov/2); /*approx*//* +ADDMHITS1 is an arbitrary (large) number for additional mhits over mean*/
    minorE2 = minl(featp->observed.cov*(pe1*pe1 + pe2*pe2) + ADDMHITS2, featp->observed.cov/2); /*approx*//* +ADDMHITS2 is an arbitrary (large) number  for additional mhits over mean*/
    pbfold = (double ****)calloc((unsigned int)pmc/2+1,sizeof(double ***));
    for(i=0;i<pmc/2+1;i++) {
        pbfold[i] = (double ***)calloc((unsigned int)nr/2+1,sizeof(double **));
        for(j=0;j<nr/2+1;j++) {
            pbfold[i][j] = (double **)calloc((unsigned int)minorE1+1,sizeof(double *));
            for(k=0;k<minorE1+1;k++) {
                pbfold[i][j][k] = (double *)calloc((unsigned int)minorE2+1,sizeof(double));
            }
        }
    }
    
    for(i=0;i<pmc+1;i++) {
        //pool
        phiA = (double)i/pmc * (1.0-pe1) + (1.0-(double)i/pmc)*pe1;
        phiR = (1.0-(double)i/pmc) * (1.0-pe2) + ((double)i/pmc)*pe2;
        phiA = phiA/(phiA+phiR+pe1+pe2);
        phiR = phiR/(phiA+phiR+pe1+pe2);
        pe1  = pe1/(phiA+phiR+pe1+pe2);
        pe2  = pe2/(phiA+phiR+pe1+pe2);
        for(prE1=0; prE1<minorE1+1; prE1++) {
            for(prE2=0; prE2<minorE2+1; prE2++) {
                for(pr=nr-minorE1-minorE2; pr>=0; pr--) {
                    //unfolded but order known
                    nrm[0] = nr-pr-prE1-prE2; nrm[1] = pr; nrm[2] = prE1; nrm[3] = prE2;
                    prob_bin_p = lnmultinomial4(nrm[0],nrm[1],nrm[2],nrm[3],phiA,phiR,pe1,lnfact);
                    qsort(nrm,4,sizeof(unsigned long),compunsign);
                    //folded and order unknown
                    pbfold[minl(pmc-i,i)][nrm[2]][nrm[1]][nrm[0]] =  pbfold[minl(pmc-i,i)][nrm[2]][nrm[1]][nrm[0]] + exp(prob_bin_p) / (2.0-deltak(pmc-i,i));
                }
            }
        }
    }
    /*...20230905...*/
    /*///////////////////////////////////*/
    newParamProposal_bial(&current,&current,combin,featp->observed,lnfact,lnprior,0.5, 0.5,0.5,0.5,0.5,0.5/*,0.5*/, max_cov,theta, pmc,poolsize, pe1, pe2, pe3, pe4, pee, priortype, outg, 0, iter, &lnprob, pbfold);
    
    current.lnprob = lnprob;
    /*
    current.lnprob = calculate_prob_bial(featp->observed.freq[0],featp->observed.freq[1],
                                         featp->observed.freq[2]+featp->observed.freq[3],
                                         featp->observed.outg,
                                         current.p1, current.p2,
                                         current.frp1,current.frp2,current.popfr,
                                         current.ner1,current.ner2,current.ner3,current.ner4,
                                         lncombin,combin,lnprior,lnfact, max_cov,
                                         theta,pmc,poolsize,pe1,pe2, pee, priortype, outg);
    */
    i=1;
    do{ /*MCMC*/
        doMCMCiterbial(&current,lnfact,lncombin,combin,lnprior,featp,max_cov,theta,seed,pmc,poolsize, pe1, pe2, pe3, pe4, pee, priortype, outg, i, iter, pbfold);
        /*add to posterior, each i replicates and after discarding burnin*/
        if (i>burnin){ /*burn-in*/
            if (i % thinin == 0 && current.lnprob > -10000){
                post_nc[current.nc]++;
                post_p1[current.p1]++;
                post_p2[current.p2]++;
                post_ner1[current.ner1]++;
                post_ner2[current.ner2]++;
                post_ner3[current.ner3]++;
                post_ner4[current.ner4]++;
                post_frp1[(unsigned long)(current.frp1)]++;
                post_frp2[(unsigned long)(current.frp2)]++;
                post_popf[(unsigned long)(current.popfr)]++;

                post_together[(unsigned long)cnt_samples].nc = current.nc;
                post_together[(unsigned long)cnt_samples].p1 = current.p1;
                post_together[(unsigned long)cnt_samples].p2 = current.p2;
                post_together[(unsigned long)cnt_samples].p3 = current.p3;
                post_together[(unsigned long)cnt_samples].ner1 = current.ner1;
                post_together[(unsigned long)cnt_samples].ner2 = current.ner2;
                post_together[(unsigned long)cnt_samples].ner3 = current.ner3;
                post_together[(unsigned long)cnt_samples].ner4 = current.ner4;
                post_together[(unsigned long)cnt_samples].frp1 = current.frp1;
                post_together[(unsigned long)cnt_samples].frp2 = current.frp2;
                post_together[(unsigned long)cnt_samples].popfr = current.popfr;
                post_together[(unsigned long)cnt_samples].lnprob = current.lnprob;
                
                featp[0].dist_estimated.nr[(unsigned long)cnt_samples] = current.nc;
                featp[0].dist_estimated.nrFreq[0][(unsigned long)cnt_samples] = current.p1;
                featp[0].dist_estimated.nrFreq[1][(unsigned long)cnt_samples] = current.p2;
                //featp[0].dist_estimated.nrFreq[2][(unsigned long)cnt_samples] = current.p3;
                //featp[0].dist_estimated.nrFreq[3][(unsigned long)cnt_samples] = current.p4;
                featp[0].dist_estimated.nerr[0][(unsigned long)cnt_samples] = current.ner1;
                featp[0].dist_estimated.nerr[1][(unsigned long)cnt_samples] = current.ner2;
                featp[0].dist_estimated.nerr[2][(unsigned long)cnt_samples] = current.ner3;
                featp[0].dist_estimated.nerr[3][(unsigned long)cnt_samples] = current.ner4;
                featp[0].dist_estimated.freqpool[0][(unsigned long)cnt_samples] = current.frp1;
                featp[0].dist_estimated.freqpool[1][(unsigned long)cnt_samples] = current.frp2;
                featp[0].dist_estimated.freqpop[(unsigned long)cnt_samples] = current.popfr;
                featp[0].dist_estimated.lnprob[(unsigned long)cnt_samples] = current.lnprob;

                cnt_samples++;
            }
        }
        if (i>=MCMCiter){done=1;}
        i++;
    }while(done==0);
    
    if(cnt_samples == 0) {
        free(post_ner1);
        free(post_ner2);
        free(post_ner3);
        free(post_ner4);
        free(post_frp1);
        free(post_frp2);
        free(post_popf);
        free(post_p2);
        free(post_p1);
        free(post_nc);
        free(post_together);
        return 1;
    }
    
    /*free 4-d matrix*/
    for(i=0;i<pmc/2+1;i++) {
        for(j=0;j<nr/2+1;j++) {
            for(k=0;k<minorE1;k++) {
                free(pbfold[i][j][k]);
            }
            free(pbfold[i][j]);
        }
        free(pbfold[i]);
    }
    free(pbfold);

    /*integrate posterior*/
    double integ_nc, integ_p1,integ_p2;
    double integ_fr1, integ_fr2;
    double integ_ner1,integ_ner2;
    double integ_ner3,integ_ner4;
    double integ_popf;
    integ_nc = integ_p1 = integ_p2 = 0;
    integ_fr1 = integ_fr2 = 0;
    integ_ner1 = integ_ner2 = 0;
    integ_ner3 = integ_ner4 = 0;
    integ_popf = 0;
    
    /*all distributions already kept*/
    for (i=1;i<=pmc;i++){
        integ_nc+=i*post_nc[i];
        integ_p1+=i*post_p1[i];
        integ_p2+=i*post_p2[i];
        integ_fr1+=i*post_frp1[i];
        integ_fr2+=i*post_frp2[i];
    }
    for (i=1;i<=NPOPINTV;i++){/*say 100 divisions, for example, */
        integ_popf+=i*post_popf[i];
    }
    for (i=1;i<=max_cov;i++){
        integ_ner1+=i*post_ner1[i];
        integ_ner2+=i*post_ner2[i];
        integ_ner3+=i*post_ner3[i];
        integ_ner4+=i*post_ner4[i];
    }
    /*return results*/
    featp->estimated.nr=(integ_nc/cnt_samples);
    featp->estimated.nrFreq[0]=(integ_p1/cnt_samples);
    featp->estimated.nrFreq[1]=(integ_p2/cnt_samples);
    featp->estimated.nerr[0]=(integ_ner1/cnt_samples);
    featp->estimated.nerr[1]=(integ_ner2/cnt_samples);
    featp->estimated.nerr[2]=(integ_ner3/cnt_samples);
    featp->estimated.nerr[3]=(integ_ner4/cnt_samples);
    featp->estimated.freqpool[0]=(integ_fr1/cnt_samples)/(double)pmc;/*<- here becomes frequency!!! [0,1]*/
    featp->estimated.freqpool[1]=(integ_fr2/cnt_samples)/(double)pmc;/*<- here becomes frequency!!! [0,1]*/
    featp->estimated.freqpop=(integ_popf/cnt_samples);/*this is already frequency*/

    free(post_ner1);
    free(post_ner2);
    free(post_ner3);
    free(post_ner4);
    free(post_frp1);
    free(post_frp2);
    free(post_popf);
    free(post_p2);
    free(post_p1);
    free(post_nc);
    free(post_together);
    /**/
    
    return 1;
}
/*///////////////////////////////////////*/
void anp(unsigned long pmc, double *an_m, int priortype)
{
    unsigned long j;
    double tmp,tmp2;
    
    if(priortype == 1){
        for (j = 1; j < pmc; j++) {/*assuming a neutral SFS*/
            tmp=1.0/(double)j;
            an_m[j] = an_m[j-1] + tmp;
        }
    }
    if(priortype == 2){
        tmp2=(1.0/(double)pmc);//assuming a uniform SFS
        for (j = 1; j < pmc; j++) {
            tmp = tmp2;
            an_m[j] = an_m[j-1] + tmp;
        }
    }
    if(priortype == 3){
        for (j = 1; j < pmc; j++) {/*assuming an exponential SFS*/
            tmp=1.0/((double)j*(double)(j+1));
            an_m[j] = an_m[j-1] + tmp;
        }
    }
    if(priortype == 4){
        tmp2=(1.0/(double)(pmc+1));//uniform, also for freq[0]
        an_m[0] = tmp2;
        for (j = 1; j < pmc+1; j++) {
            tmp = tmp2;
            an_m[j] = an_m[j-1] + tmp;
        }
    }
}
/*///////////////////////////////////////*/
void ComputePrior(double *lnprior,unsigned long pmc, double *an_m, int priortype, double theta, int outg, double divergence){
    unsigned long j;
    double tmp,tmp2;
    double fN;
    
    if(priortype == 1){ /*SNM*/
        /*for large values an_m(x) = 0.578 + ln(x-1)*/
        for (j = 1; j < pmc; j++) {/*assuming a neutral SFS*/
            tmp=1.0/((double)j);
            if(outg==0) {
                //if((double)j/(double)pmc != 0.5) {
                    tmp+=1.0/((double)(pmc-j));
                    tmp/=2.0;
                //}
            }
            lnprior[j] = ln(theta*tmp);
        }
        /*invariant*/
        if(pmc > LIMIT_FACT) {
            if(theta*(0.578+ln((double)pmc-1.0)) < 1.0){
                if(outg==0) {lnprior[0]=ln((1.0-theta*(0.578+ln((double)pmc-1.0)))/2.0);}
                else lnprior[0]=ln((1.0-theta*(0.578+ln((double)pmc-1.0))-divergence));
            }else lnprior[0]=-1e+300;
        } else {
            if(theta*an_m[pmc-1] < 1.0) {
                if(outg==0) {lnprior[0]=ln((1.0-theta*an_m[pmc-1])/2.0);}
                else lnprior[0]=ln((1.0-theta*an_m[pmc-1]-divergence));
            }else lnprior[0]=-1e+300;
        }
        if(outg==0) fN = lnprior[0]; else fN = ln(divergence);
        lnprior[pmc]=fN;
    }
    if(priortype == 2){
        tmp2=(1.0/(double)(pmc-1.0));/*assuming a uniform SFS*/
        for (j = 1; j < pmc; j++) {
            tmp = tmp2;
            if(outg==0) {
                //if((double)j/(double)pmc != 0.5) {
                    tmp+=tmp2;
                    tmp/=2.0;
                //}
            }
            lnprior[j] = ln(theta*tmp);
        }
        /*invariant**/
        if(theta*an_m[pmc-1] < 1.0) {
            if(outg==0) {lnprior[0]=ln((1.0-theta*an_m[pmc-1])/2.0);}
            else lnprior[0]=ln(1.0-theta*an_m[pmc-1]-divergence);
        } else lnprior[0]=-1e+300;
        if(outg==0) fN = lnprior[0]; else fN = ln(divergence);
        lnprior[pmc]=fN;
    }
    if(priortype == 3){
        /* include exponential (as in Ohtsuki and Innan 2017) 1/((i*(i+1)): tmp=1.0/((double)j*(double)(j+1));*/
        /*for large values an_m(x) = 1.0 + ln(x) - ln(x+1)*/
        for (j = 1; j < pmc; j++) {/*assuming a neutral SFS*/
            tmp=1.0/((double)j*(double)(j+1));
            if(outg==0) {
                //if((double)j/(double)pmc != 0.5) {
                    tmp+=1.0/((double)(pmc-j)*(double)(pmc-(j+1)));
                    tmp/=2.0;
                //}
            }
            lnprior[j] = ln(theta*tmp);
        }
        /*invariant**/
        if(pmc > LIMIT_FACT) {
            if(theta*(1.0+ln((double)pmc-1)-ln(pmc-1+1)) < 1.0) {
                if(outg==0) lnprior[0]=ln((1.0-theta*(1.0+ln((double)pmc-1)-ln(pmc-1+1)))/2.0);
                else lnprior[0]=ln(1.0-theta*(1.0+ln((double)pmc-1)-ln(pmc-1+1))-divergence);
            } else lnprior[0]=-1e+300;
        } else {
            if(theta*an_m[pmc-1] < 1.0) {
                if(outg==0)  lnprior[0]=ln((1-theta*an_m[pmc-1])/2.0);
                else lnprior[0]=ln(1-theta*an_m[pmc-1]-divergence);
            } else lnprior[0]=-1e+300;
        }
        if(outg==0) fN = lnprior[0]; else fN = ln(divergence);
        lnprior[pmc]=fN;
    }
    if(priortype == 4){
        for (j = 0; j < pmc; j++) {/*assuming all NULL*/
            lnprior[j]=0.;//log(1.0/(double)pmc);
        }
    }
}
/*///////////////////////////////////////*//*NOT USED*/
double ComputePriorFreq(unsigned long j, int priortype, double theta, unsigned long poolsize,int outg){
    double tmp=0.0;
    double tmp2=0.0;
    double lnprior=-1e+300;

    if(priortype == 1){ /*SNM*/
        /*for large values an_m(x) = 0.578 + ln(x-1)*/
        if(j>0 && j<poolsize) {
            tmp=1.0/(double)j;
            if(outg==0) {
                //if((double)j/(double)poolsize != 0.5) {
                    tmp+=1.0/(double)(poolsize-j);
                //}
            }
            lnprior = (ln(theta*tmp*0.5));
        }
        else { /*invariant: do not consider outgroup**/
            if(theta*(0.578+ln((double)poolsize-1.0)) < 1.0)
                lnprior = ln(1.0-theta*(0.578+ln((double)poolsize-1.0)));
            else
                lnprior = -1e+300;
        }
    }
    if(priortype == 2){ /*assuming a uniform SFS*/
        if(j>0 && j<poolsize) {
            tmp2=(1.0/(double)(poolsize-1.0));
            if(outg==0) {
                //if((double)j/(double)poolsize != 0.5) {
                    tmp+=tmp2;
                //}
            }
            lnprior = (ln(theta*tmp));
        }
        else { /*invariant: do not consider outgroup**/
            if(theta*tmp2*(poolsize-1) < 1.0)
                lnprior=ln(1.0-theta*tmp2*(poolsize-1));
            else
                lnprior=-1e+300;
        }
    }
    if(priortype == 3){
        /* include exponential (as in Ohtsuki and Innan 2017) 1/((i*(i+1)): tmp=1.0/((double)j*(double)(j+1));*/
        /*for large values an_m(x) = 1.0 + ln(x) - ln(x+1)*/
        if(j>0 && j<poolsize) {
            tmp=1.0/((double)j*(double)(j+1));
            if(outg==0) {
                //if((double)j/(double)poolsize != 0.5) {
                    tmp+=1.0/((double)(poolsize-j)*(double)(poolsize-(j+1)));
                //}
            }
            lnprior = (ln(theta*tmp));
        }
        else { /*invariant: do not consider outgroup**/
            if(theta*(1.0+ln((double)poolsize-1)-ln(poolsize)) < 1.0)
                lnprior=ln(1.0-theta*(1.0+ln((double)poolsize-1)-ln(poolsize)));
            else
                lnprior=-1e+300;
        }
    }
    if(priortype == 4){
        lnprior = 0;
    }
    return(lnprior);
}
/*///////////////////////////////////////*/
void find_common_bases(struct MCMCsample *post_together,double threshlik,unsigned long maxthin, int *cp1, int *cp2, int *cp3, int *cp4,int *cf1,int *cf2,int *cer1,int *cer2,int *cer3,int *cer4) {
    int i=0;
    //int ccp4;
    *cp1 = (int)post_together[i].p1;
    *cp2 = (int)post_together[i].p2;
    *cp3 = (int)post_together[i].p3;
    *cp4 = (int)post_together[i].p4;

    *cf1 = (int)post_together[i].frp1;
    *cf2 = (int)post_together[i].frp2;
    *cer1 = (int)post_together[i].ner1;
    *cer2 = (int)post_together[i].ner2;
    *cer3 = (int)post_together[i].ner3;
    *cer4 = (int)post_together[i].ner4;

    for(i=1;i< maxthin;i++) {
        if(post_together[i].lnprob >= threshlik) {
            if(*cp1 > (int)post_together[i].p1) {*cp1=(int)post_together[i].p1;}
            if(*cp2 > (int)post_together[i].p2) {*cp2=(int)post_together[i].p2;}
            if(*cp3 > (int)post_together[i].p3) {*cp3=(int)post_together[i].p3;}
            if(*cp4 > (int)post_together[i].p4) {*cp4=(int)post_together[i].p4;}
            if(*cf1 > (int)post_together[i].frp1) {*cf1=(int)post_together[i].frp1;}
            if(*cf2 > (int)post_together[i].frp2) {*cf2=(int)post_together[i].frp2;}
            if(*cer1 > (int)post_together[i].ner1) {*cer1=(int)post_together[i].ner1;}
            if(*cer2 > (int)post_together[i].ner2) {*cer2=(int)post_together[i].ner2;}
            if(*cer3 > (int)post_together[i].ner3) {*cer3=(int)post_together[i].ner3;}
            if(*cer4 > (int)post_together[i].ner4) {*cer4=(int)post_together[i].ner4;}
            i++;
        }
    }
    return;
}
/*///////////////////////////////////////*/
/*compare two double numbers in a long int list*/
int compare_(const void *i,const void *j)
{
    if(*(double *)i < *(double *)j) return +1;
    if(*(double *)i > *(double *)j) return -1;
    return 0;
}

/*///////////////////////////////////////*/
void newParamProposal_bial(struct MCMCsample *current, struct MCMCsample *next,double ***combin, struct Obs observed, double *lnfact,double *lnprior, double random0, double random1, double random2, double random3, double random4, double random5, /*double random6,*/ long max_cov, double theta, unsigned long pmc, unsigned long poolsize, double pe1, double pe2, double pe3, double pe4, double pee, int priortype, int outg, unsigned long it, unsigned long iter, double *lnprob, double ****pbfold) {

    double pval0;
    unsigned long i/*,k,j*/;
    double Sprob;
    
    /*static*/ struct MCMCsample *mcmc_par = 0;
    /*static*/ double *lnpval=0;
    static double sumpval=0;
    
    next->nc=0;
    next->p1=next->p2=0;
    next->p3=next->p4=0;/*not used*/
    next->frp1=next->frp2=next->lnprob=0.0;
    next->ner1=next->ner2=0;
    next->ner3=next->ner4=0;
    //next->popfr=0;
    
    /*NO USING CURRENT (previous iteration) parameters...(finally NOT MHMCMMC) */

    /*WE MUST USE THE GENERATION OF PARAMETERS FROM A 7-DIMENSION MATRIX PROBABILITY*/
    /*AN APPROACH IS TO COMBINE ONLY THE FIRST 2-D WITH THE REST */
    /*sampling ... considering the first 2 pars and the rest independently*/
    /* calculate the fr1 and the popfr parameters using observations*/
    
    if(it == 0) {
        sumpval = 0.;
    }
    
    mcmc_par = (struct MCMCsample *)calloc((pmc+1),sizeof(struct MCMCsample));
    lnpval =(double *)calloc((pmc+1),sizeof(double));
    /*
    double **lnpval;
     lnpval = (double **)calloc((unsigned int)pmc+1,sizeof(double ***));
    for(k=0;k<NPOPINTV+1;k++) {
        //frequency
        lnpval[k] = (double *)calloc((unsigned int)k+1,sizeof(double *));
    }
    */

    /*calculate the 2-D grid-matrix probabilities for population freq and for pool freq but not sample until having the total prob.*/
    //multinomialfr_seqbial_sample(observed.freq[0], observed.freq[1], observed.freq[2],observed.freq[3],observed.outg,pe1,pe2,pe3,pe4,pee,lnprior,lnfact, mcmc_par, pmc, outg,lnpval,&sumpval); /*to include sorting to calculate lnpval and sumpval*/
    multinomialfr_seqbial_folded(observed.freq[0], observed.freq[1], observed.freq[2],observed.freq[3],observed.outg, pe1,pe2,pe3,pe4,pee,lnprior,lnfact, mcmc_par, pmc, outg,lnpval,&sumpval,pbfold); /*modify to include sorting to calculate lnpval and sumpval*/
    //}
    /*estimate the probabilities for each parameter (in three sections): */
    Sprob = 0.;
    /*chose random parameters according the probabilities in three sections (heuristic): seq_errors, combinatorics and total(given grid 2D freqa)*/
    for(i=0;i<pmc+1;i++) {
        /* 0. include the probability of having (freq k in the pop and) freq i in the pool*/
        mcmc_par[i].lnprob = lnpval[i] - ln(sumpval);
        
        /* 1. sampling the ner1 and ner2 values from the multibinomial distribution using observations*/
        binomialner_seqbial_sample(observed.freq[0], observed.freq[1], observed.freq[2]+observed.freq[3], pe1, pe2, pee,lnfact, &mcmc_par[i], random0, random1, random2, pmc, theta, lnprob);
        mcmc_par[i].lnprob += (*lnprob);
        
        /* 2. p1: sampling from the P(nc1|nr1,pmc) using observations. Combinatorial. Reference nt*/
        mcmc_par[i].p1 = combin_seq_sample(observed.freq[0], mcmc_par[i].frp1, mcmc_par[i].ner1, mcmc_par[i].ner2, mcmc_par[i].ner3, combin,random3, poolsize, (observed.freq[0]>0), lnprob);
        /* total nc: sampling from the P(nc|nr,pmc) the nc = nc1 + nc2 values using observations. Combinatorial. Reference nt*/
        //mcmc_par[i].nc = combin_seq_sample(observed.freq[0]+observed.freq[1],mcmc_par[i].frp1+mcmc_par[i].frp2, 0, 0, 0, combin,random3, poolsize, (observed.freq[0]>0)+(observed.freq[1]>0), lnprob);
        mcmc_par[i].lnprob += (*lnprob);
        
        /* 3. p2: sampling from the P(nc2|nr1,pmc) using observations. Combinatorial. Reference nt*/
        mcmc_par[i].p2 = combin_seq_sample(observed.freq[1], pmc - mcmc_par[i].frp1, mcmc_par[i].ner2, mcmc_par[i].ner1, mcmc_par[i].ner4, combin,random4, poolsize, (observed.freq[1]>0), lnprob);
        /*sampling from the P(nc2|nc,fr2,pmc), sample Binomial, the nc2 value (and therefore the nc1) using observations. Alternative nt*/
        //mcmc_par[i].p2 = binomial_seqbial_sample(mcmc_par[i].nc, mcmc_par[i].frp2/pmc,(observed.freq[1]>0),random4,lnfact,lnprob);
        mcmc_par[i].lnprob += (*lnprob);
        
        /*calculate nc*/
        mcmc_par[i].nc = mcmc_par[i].p1 + mcmc_par[i].p2;
        /*calculate p1*/
        //mcmc_par[i].p1 = mcmc_par[i].nc + mcmc_par[i].p2;
        
        /*sums of probs*/
        Sprob += exp(mcmc_par[i].lnprob);
    }
    if(Sprob == 0.0){
        next->frp1 = ((double)observed.freq[0]/(double)(observed.freq[0]+observed.freq[1]+observed.freq[2]+observed.freq[3]) * (1.0-pe2-pee) +
                      (1.0-(double)observed.freq[0]/(double)(observed.freq[0]+observed.freq[1]+observed.freq[2]+observed.freq[3]))*pe1) * pmc;
        next->frp2 = pmc - next->frp1;
        next->popfr = 0;
        next->ner1 = 0;
        next->ner2 = 0;
        next->ner3 = 0;
        next->ner4 = 0;
        next->p1 = next->frp1;
        next->p2 = next->frp2;
        next->nc = pmc;
    }
    else {
        //Choose the SNP frequency in the pool, given the frequency
        pval0 = 0.;
        for(i=0;i<pmc;i++) {
            pval0 += exp(mcmc_par[i].lnprob)/Sprob;
            if(pval0 >= random5) break;
        }

        if(i > pmc) {
            printf("\nError: multinomial assigns an impossible value.\n");
            exit(1);
        }
    }
    
    /*Keep the mcmc_par (heuristic approach):*/
    next->popfr = 0;
    next->frp1 = i;
    next->frp2 = pmc - next->frp1;
    next->ner1 = mcmc_par[i].ner1;
    next->ner2 = mcmc_par[i].ner2;
    next->ner3 = mcmc_par[i].ner3;
    next->ner4 = mcmc_par[i].ner4;
    next->p1 = mcmc_par[i].p1;
    next->p2 = mcmc_par[i].p2;
    next->nc = mcmc_par[i].nc;
    next->lnprob = mcmc_par[i].lnprob;
    *lnprob = next->lnprob;
    
    //if(it==iter) {
        free(mcmc_par);
        free(lnpval);
    //}

    if((next->p1>0) + (next->p2>0) < (observed.freq[0] - next->ner1 + next->ner2 + next->ner3 > 0) + (observed.freq[1] - next->ner2 + next->ner1 + next->ner4 >0)) {
        printf("\nError: More different variants observed than estimated (without error)");
        exit(1);
    }
    if((next->p1)  + (next->p2) > pmc) {
        printf("\nError: next->p1 + next->p2 > pmc \n %ld + %ld > %ld \n",next->p1,next->p2,pmc);
        exit(1);
    }
    if((next->frp1 == pmc && next->p2>0) || (next->frp2 == pmc && next->p1>0)) {
        printf("\nError: (next->frp1 or next->frp2) == pmc but p1 or p2 exist: \n (frp1=%ld, p2=%ld) or (frp2=%ld, p1=%ld), pmc=%ld \n",next->frp1,next->p2,next->frp2,next->p1,pmc);
        exit(1);
    }

    return;
}
/*///////////////////////////////////////*/
double calculate_prob_bial(unsigned long obs_p1, unsigned long obs_p2, unsigned long obs_p3, unsigned char obs_outg, unsigned long exp_p1, unsigned long exp_p2, unsigned long exp_frp1, unsigned long exp_frp2, unsigned long exp_popfr, unsigned long exp_ner1,unsigned long exp_ner2,unsigned long exp_ner3,unsigned long exp_ner4,double ***lncombin,double ***combin,double *lnprior,double *lnfact, unsigned long max_cov, double theta, unsigned long pmc, unsigned long poolsize, double pe1, double pe2, double pee, int priortype, int outg)
{
    double lnprob;
    double cc1,cc2,ff1n,ff1d,popf,ee;
    unsigned long i,j,k;
    double lnp;
    double peR,peA,peE,pR;
    
    //cc1 = lnbinomial(exp_p1,exp_p2,(double)exp_frp1/(double)poolsize,lnfact);/*binomial*/
    /*combinatorial*/
    if(poolsize > LIMIT_COMB) {
        cc1 = cc2 = 0.0;
        //cc2 = 0.0;
    }
    else {
        cc1 = lncombin[exp_frp1][exp_p1][obs_p1-exp_ner1+exp_ner2+exp_ner3];
        cc2 = lncombin[exp_frp2][exp_p2][obs_p2-exp_ner2+exp_ner1+exp_ner4];
        //cc2 = lncombin[exp_frp1+exp_frp2][exp_p1+exp_p2][obs_p1+obs_p2];
    }
    
    /*sequence errors*/
    pR = (double)exp_frp1/((double)exp_frp1+(double)exp_frp2);
    peR = (1.0-pR)*pe1 / ((    pR)*(1.0-pe2-pee)+(1.0-pR)*pe1);
    peA = (    pR)*pe2 / ((1.0-pR)*(1.0-pe1-pee)+(    pR)*pe2);
    peE = pR;
    ee  = lnbinomial_seqbial(exp_ner1, exp_ner2, exp_ner3, obs_p1, obs_p2, obs_p3, peR, peA, peE, lnfact);

    /*Pool and population frequency*/
    /*num*/
    lnp = lnbinomial(exp_frp1,pmc-exp_frp1,(double)exp_popfr/(double)NPOPINTV,lnfact);
    if(outg) {if(obs_outg==0) j=exp_popfr; else j=(int)NPOPINTV-(unsigned long int)exp_popfr;}
    else j = minint((unsigned long int)exp_popfr,(int)NPOPINTV-(unsigned long int)exp_popfr);
    popf = lnprior[j];
    ff1n = lnmultinomialfr(obs_p1, obs_p2, obs_p3,pR,pe1,pe2,pee,lnfact) + lnp + popf;
    /*den*/
    ff1d = 0.;
    for(i=0;i<pmc+1;i++) {
        for(k=0;k<NPOPINTV+1;k++) {
            lnp = lnbinomial(i,pmc-i,(double)k/(double)NPOPINTV,lnfact);
            if(outg) {if(obs_outg==0) j = k; else j=NPOPINTV-k;}
            else j = minint(k,NPOPINTV-k);
            popf = lnprior[j];
            ff1d += exp(lnmultinomialfr(obs_p1, obs_p2, obs_p3, (double)i/pmc, pe1, pe2, pee,lnfact) + lnp + popf);
        }
    }
    ff1d = ln(ff1d);

    /*total*/
    lnprob = cc1 + cc2 + ee + ff1n - ff1d;
    return(lnprob);
}
/*///////////////////////////////////////*/

