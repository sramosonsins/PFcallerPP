/* 
 * File:   main.c
 * Author: Pablo Librado
 *
 * Created on 5 de noviembre de 2014, 14:03
 */

#include "main.h"
#include "progressbar.h"

/*////////////////////////////////////////////////////////////////////*/
int main(int argc, char** argv) {
    unsigned long i;
    long double **s2 = 0; /*matrix of stirling numbers second kind S(nr,nc)*/
    double *lnfact=0;       /*logn of the factorial number*/
    double ***combin=0;        /*P(nc|nr,p)*/
    double ***lncombin=0;      /*ln(P(nc|nr,p))*/
    double *lnprior=0;         /*lnP(pa|p)*/
    char file_fasta[LEN_FILENAME];
    //char file_fasta2[LEN_FILENAME];
    char **seqs;
    char file_liklh[LEN_FILENAME];
    
    unsigned long nscaffolds;
    char **chr_name_array;
    char **chr_length_array;
    char chrname[MSP_MAX_NAME];
    unsigned long first = 0;
    unsigned long pp;
    
    unsigned long cnt,cnt2;
    long int seed,output_seed;

    unsigned long max_cov, min_cov;
    unsigned long min_qual;
    int biallelic;
    int outg;
    double divergence=0;
    unsigned long div_pos=0;
    double min_error;
    double sum_errors = 0;
    double sum_errors_variant = 0;
    int priortype;
    long int samplemhmcmc;
    long int burnin_par;
    unsigned long iter;
    unsigned long nitermc;
    char *prfasta;
    unsigned long lentotal;
    
    double *an_m=0;
    double theta=0.0;
    double thetant=0.0;
    unsigned long theta_pos = 0;
    unsigned long pmc/*,pmc2*/;
    unsigned long pmm;
    unsigned long poolsize;
    float qerr_range;
    int teststrands;

    char file_in[LEN_FILENAME];
    char file_out[LEN_FILENAME];
    char *chr_name_all;
    char *chr_length_all;
    
    long int ntotalreads = 0;
    long int ntotalreads_novariant = 0;
    long int npositions_novariant = 0;
    long int npositions_variant = 0;
    double   npositions_noerror_exp = 0;
    
    double sum_per = 0.0;
    long int npositions_per = 0;
    
    double *n_nr;
    struct list_chr_pileup *chr_data=0;
    unsigned long num_scaffolds_pileup=0;
    
    double mean_ml_err;
 
    char file_outw[LEN_FILENAME];
    FILE *fasta_file = 0;
    SGZip fasta_file_gz;
    struct SGZIndex fasta_file_gz_index;          /* This is the index for the output gz file. */
    init_gzindex_structure(&fasta_file_gz_index); /* IMPORTANT TO INITIALIZE!*/
    
    chr_name_all=(char *)calloc(LEN_FILENAME,sizeof(char));
    chr_length_all=(char *)calloc(LEN_FILENAME,sizeof(char));
    prfasta=(char *)calloc(5,sizeof(char));

    file_in[0] =  '\0';
    file_out[0] =  '\0';
    file_fasta[0] =  '\0';
    file_liklh[0] =  '\0';
    
    check_param(argc,argv, &max_cov, &min_cov, &min_qual, chr_length_all, chr_name_all, &biallelic,file_in,file_out, &priortype, &samplemhmcmc, &burnin_par, &iter, &nitermc, prfasta, &lentotal, &theta, &seed, &poolsize,&qerr_range, &outg, &teststrands);
    
    printf(PSCALLER);
    int arg = 1;
    printf("PFcaller ");
    if( argc > 1 ) {
        while(arg < argc) {
            printf("%s ",argv[arg]);
            arg++;
        }
    }
    printf("\n");
    
    min_error=1.0-(pow(10,-(double)(min_qual)/10.0));
    if(file_out[0] != '\0') printf("OK\n");
    
    /*////////////////////////////////////////////////////////*/
    /*read pileup file, and save all types of different positions*/
    if(file_out[0] != '\0') printf ("Reading pileup file...\n");
    struct typepos *typep = (struct typepos *) calloc(1,sizeof(struct typepos));
    struct featpos *featp = (struct featpos *) calloc(1,sizeof(struct featpos));
    /*separate all values of the list chr_name_all in chr_name_array: */
    /* Only do the list if input and output is tfa*/
    if(read_index_file(chr_name_all,&nscaffolds,&chr_name_array,&chr_length_array)) {
        printf("\nError reading the index file %s",chr_name_all);
        exit(1);
    }
    
    n_nr = (double *)calloc((unsigned int)max_cov+1,sizeof(double));/*number of positions with nr reads*/
    
    //printf ("\nEach dot means 10Gb positions readed.\n");
    ReduceInput(file_in,&typep,&featp, max_cov, min_cov, min_qual, &ntotalreads, &ntotalreads_novariant, &npositions_novariant, &npositions_variant, n_nr,min_error,&sum_errors, &sum_errors_variant, &npositions_noerror_exp, iter, &cnt, &cnt2, biallelic, qerr_range, &chr_data, &num_scaffolds_pileup,nscaffolds,chr_name_array,chr_length_array, &sum_per, &npositions_per, &thetant, &theta_pos,poolsize,outg,&divergence,&div_pos,teststrands);
    
    /*fraction_error = (double)ntotalreads / (double)(ntotalreads - ntotalreads_novariant);*/

    if(file_out[0] != '\0') printf ("Number of different patterns retained: %lu (from %lu annotated positions)\n",cnt,cnt2);
    if(file_out[0] != '\0') printf ("OK\n");
    
    if(theta==-1) {
        if(theta_pos>0)
            thetant = thetant / (double)theta_pos;
        else
            thetant = 0.0;
    }
    else {
        if(poolsize > 1) thetant = theta;
        else thetant = 0.;
    }
    printf("\n\nthetant value is set to %.4f\n",thetant);
    if(outg) {
        divergence = divergence/(double)div_pos;
        printf("Divergence is set to %.4f\n",divergence);
    }
    /*////////////////////////////////////////////////////////*/
    /*pre-compute combinatorics*/
    
    pmm = (max_cov > (poolsize < LIMIT_FACT ? poolsize : LIMIT_FACT) ? max_cov : (poolsize < LIMIT_FACT ? poolsize : LIMIT_FACT));
    lnfact = (double *)calloc(pmm+1, sizeof(double));
    lnfactorial(lnfact,pmm);
    
    /*Combinatorics for frequency and read depth moderate, less than LIMIT_COMB*/
    /*if more, we calculate until maximum and use MonteCarlo for the rest*/
    if(poolsize <= LIMIT_COMB) {
        pmc = poolsize;
        if(file_out[0] != '\0') printf("Computing combinatorics...");
        combin = (double ***)calloc((unsigned int)pmc+1,sizeof(double **));
        lncombin = (double ***)calloc((unsigned int)pmc+1,sizeof(double **));
        
        /*stirling*/
        s2 = (long double **)calloc((unsigned int)max_cov+1,sizeof(long double *));
        for(i=0;i<max_cov+1;i++) {
            s2[i] = (long double *)calloc((unsigned int)pmc/*poolsize*/+1,sizeof(long double));
        }
        stirling2(s2,max_cov,pmc /*poolsize*/);
        
        for(pp=0;pp<pmc+1;pp++) {
            /*combinatorics*/
            lncombin[pp] = (double **)calloc((unsigned int)pmc+1,sizeof(double *));
            for(i=0;i<pmc+1;i++) {
                lncombin[pp][i] = (double *)calloc((unsigned int)max_cov+1,sizeof(double));
            }
            combin[pp] = (double **)calloc((unsigned int)pmc+1,sizeof(double *));
            for(i=0;i<pmc+1;i++) {
                combin[pp][i] = (double *)calloc((unsigned int)max_cov+1,sizeof(double));
            }
            prob_nrb_st(combin[pp],lncombin[pp],lnfact,s2,pp,max_cov);
        }
    }
    else {
        pmc = (max_cov < poolsize ? max_cov : poolsize);//LIMIT_COMB;//
    }
    
    if(file_out[0] != '\0') printf ("OK\n");
    
    /*////////////////////////////////////////////////////////*/
    /*pre-compute prior*/
    if (priortype == 1){
        if(file_out[0] != '\0') printf("Computing prior for the std. neutral model...");
    }
    if (priortype == 2){
        if(file_out[0] != '\0') printf("Computing uniform prior...");
    }
    if (priortype == 3){
        if(file_out[0] != '\0') printf("Computing exponential prior...");
    }
    if (priortype == 4){
        if(file_out[0] != '\0') printf("Computing uniform prior not considering variability...");
    }

    //if(poolsize <= LIMIT_COMB) {
        lnprior = (double *)calloc(NPOPINTV+1,sizeof(double));
        an_m = (double *)calloc(NPOPINTV+1, sizeof(double));
        anp(NPOPINTV,an_m,priortype);
        ComputePrior(lnprior,NPOPINTV,an_m,priortype,thetant,outg,divergence);
    //}

    if(npositions_per > THR_ERR) {
        mean_ml_err = sum_per/(double)npositions_per;
        printf("\nML Sequence misread estimate obtained from sequence data. per = %.3e\n",mean_ml_err);
    }
    else {
        mean_ml_err = 0.0;
        printf("\nNot enough information to estimate the sequence misread using ML. \b Misread sequence assumed from Phred scores.\n");
    }
        
    /* ///////////////////////////// */
    
    if(file_out[0] != '\0') printf ("OK\n");
    
    /* /\*\/////////////////////////////////////////////////////////\*\/ */
    /* /\*\/////////////////////////////////////////////////////////\*\/ */
    /* /\*\///////////////////   BEGIN CORE    /////////////////////\*\/ */
    /* /\*\/////////////////////////////////////////////////////////\*\/ */
    /* /\*\/////////////////////////////////////////////////////////\*\/ */
    /*do calling */
    if(file_out[0] != '\0') printf ("Estimating posterior probabilities of the allele frequencies...\n");
    for (i=0;i<cnt;i++){
        if (featp[i].observed.cov>0){
                PScallerMCMCbial(lnfact,lncombin,combin,lnprior,&featp[i],iter,max_cov,samplemhmcmc,burnin_par, thetant,&seed,pmc,poolsize,mean_ml_err, priortype,outg);
        }
        if(file_out[0] != '\0') {if(cnt>0) ShowProgressBar(i, cnt, 30); }
    }

    printf("\n");

    if(file_out[0] != '\0') printf ("OK\n");
    /* /\*\/////////////////////////////////////////////////////////\*\/ */
    /* /\*\/////////////////////////////////////////////////////////\*\/ */
    /* /\*\/////////////////////   END CORE    /////////////////////\*\/ */
    /* /\*\/////////////////////////////////////////////////////////\*\/ */
    /* /\*\/////////////////////////////////////////////////////////\*\/ */
        
 
    output_seed = seed; /*to have the same results for each output*/
    /* print fasta */
    if(strchr(prfasta,'f')!=NULL) {
        char prf[1];
        prf[0] = 'f';
        printf ("\nWriting %ld samples ...",pmc);
        printf ("\n");
        for(first=0;first<nscaffolds;first++) {
            lentotal = atol(chr_length_array[first]);
            memset(chrname, '\0', MSP_MAX_NAME);
            strcpy(chrname,chr_name_array[first]);
            
            if(file_out[0] != '\0') {
                sprintf(file_fasta, "%s",file_out);
                printf ("\nWriting fasta for scaffold %s ...",chrname);
                printf ("\n");
            }
            seqs=(char **)calloc(pmc,sizeof(char *));
            for(i=0;i<pmc;i++)
                seqs[i]=(char *)calloc((lentotal+1),sizeof(char));
            Freqs2fasta(typep,featp,seqs,prf,chrname,min_qual,iter,lentotal,cnt2, thetant, &seed,pmc, poolsize,chr_data,num_scaffolds_pileup);
            printfasta(file_fasta, seqs,chrname,pmc,lentotal);
            
            /*TESTING*//*
            char **seqs2;
            seqs2=(char **)calloc(max_cov,sizeof(char *));
            for(i=0;i<max_cov;i++)
                seqs2[i]=(char *)calloc((lentotal+1),sizeof(char));
            if(file_out[0] != '\0') {
                sprintf(file_fasta2, "%s_RAW",file_out);
                printf ("\nWriting RAW fasta for scaffold %s ...",chrname);
                printf ("\n");
            }
            Rawdata2fasta(typep,featp,seqs2,prf,chrname,min_qual,iter,lentotal,cnt2, thetant, &seed,pmc, poolsize,chr_data,num_scaffolds_pileup,max_cov);
            printfasta(file_fasta2, seqs2,chrname,max_cov,lentotal);
            for(i=0;i<max_cov;i++)
                free(seqs2[i]);
            free(seqs2);
            *//*END TESTING*/
            
            for(i=0;i<pmc;i++)
                free(seqs[i]);
            free(seqs);
            if(file_out[0] != '\0') printf("\nOK");
        }
        if(file_out[0] != '\0') printf("\nOK all fasta prints\n");
    }

    seed =  output_seed; /*to have the same results for each output*/
    /* print TFASTA */
    if(strchr(prfasta,'t')!=NULL) {
        printf ("\nWriting %ld samples ...",pmc);
        printf ("\n");
        
        if(*file_out=='\0')
            fasta_file = stdout;
        else {
            sprintf(file_outw, "%s.tfa.gz",file_out);
            fasta_file=fzopen(file_outw,"w",&fasta_file_gz);
            if (fasta_file == NULL){
                printf ("\nUnable to write tfasta file: %s\n",file_out);
                exit(1);
            }
            fasta_file_gz.index = &fasta_file_gz_index;
        }
        
        for(first=0;first<nscaffolds;first++) {
            lentotal = atol(chr_length_array[first]);
            memset(chrname, '\0', MSP_MAX_NAME);
            strcpy(chrname,chr_name_array[first]);
            
            if(file_out[0] != '\0') {
                printf ("\nWriting tfasta for scaffold %s in %s...",chrname,file_outw);
                printf ("\n");
            }
            printtfasta(fasta_file,&fasta_file_gz,chrname,first,lentotal,pmc,typep, featp, min_qual,iter,cnt2, thetant,&seed, poolsize,chr_data, num_scaffolds_pileup,argc,argv);
            if(file_out[0] != '\0') {printf("\nOK");}
        }
        if(file_out[0] != '\0') {
            fzclose(fasta_file, &fasta_file_gz);
        }
        if(file_out[0] != '\0') printf("\nOK all tfasta prints\n");
    }
    
    seed =  output_seed; /*to have the same results for each output*/
    /* print gVCF */
    if(strchr(prfasta,'g')!=NULL) {
        printf ("\nWriting %ld samples ...",pmc);
        printf ("\n");
        for(first=0;first<nscaffolds;first++) {
            lentotal = atol(chr_length_array[first]);
            memset(chrname, '\0', MSP_MAX_NAME);
            strcpy(chrname,chr_name_array[first]);
            
            if(file_out[0] != '\0') {
                sprintf(file_fasta, "%s.gvf.gz",file_out);
                printf ("\nWriting gVCF for scaffold %s in %s...",chrname,file_fasta);
                printf ("\n");
            }
            printgVCF(file_fasta,chrname,typep,featp,first,lentotal,min_qual, cnt2,theta,pmc,poolsize,iter,&seed,chr_data,num_scaffolds_pileup);
            if(file_out[0] != '\0') printf("\nOK");
        }
        if(file_out[0] != '\0') printf("\nOK all gVCF prints\n");
    }

    /* print NUMERICAL */
    if(strchr(prfasta,'0')!=NULL) {
        for(first=0;first<nscaffolds;first++) {
            lentotal = atol(chr_length_array[first]);
            memset(chrname, '\0', MSP_MAX_NAME);
            strcpy(chrname,chr_name_array[first]);
            
            if(file_out[0] != '\0') {
                sprintf(file_liklh, "%s",file_out);
                printf ("\nWriting posterior estimates for scaffold %s in %s...",chrname,file_liklh);
                printf ("\n");
            }
            printlikelihoods(file_liklh,typep,featp,chrname,first,iter,cnt,cnt2,pmc,biallelic);
            if(file_out[0] != '\0') printf("\nOK");
       }
        if(file_out[0] != '\0') printf("\nOK all tabulated prints\n");
    }
    
    /*/////////////////////////////////////////////////////////*/
    /*//free*/
    free(typep);
    /*calloc the doubles with iter values of nr and nrFreq[4] de DistEst inside featpos*/
    for(i=0;i<cnt;i++) {
        free(featp[i].dist_estimated.nr);
        free(featp[i].dist_estimated.nrFreq[0]);
        free(featp[i].dist_estimated.nrFreq[1]);
        free(featp[i].dist_estimated.nerr[0]);
        free(featp[i].dist_estimated.nerr[1]);
        free(featp[i].dist_estimated.nerr[2]);
        free(featp[i].dist_estimated.nerr[3]);
        free(featp[i].dist_estimated.freqpool[0]);
        free(featp[i].dist_estimated.freqpool[1]);
        free(featp[i].dist_estimated.freqpop);
        free(featp[i].dist_estimated.lnprob);
    }
    free(featp);
    if(file_out[0] != '\0') printf ("\nFinished\n");
  
    if(poolsize <= LIMIT_COMB) {
         for(i=0;i<max_cov+1;i++)
            free(s2[i]);
        free(s2);
        for(pp=0;pp<pmc+1;pp++){
            for(i=0;i<pmc+1;i++) {
                free(lncombin[pp][i]);
                free(combin[pp][i]);
            }
            free(lncombin[pp]);
            free(combin[pp]);
        }
        free(lncombin);
        free(combin);
    }
    free(lnprior);
    free(an_m);
    free(lnfact);
    free(n_nr);
    
    for(i=0;i<nscaffolds;i++) free(chr_name_array[i]);
    for(i=0;i<nscaffolds;i++) free(chr_length_array[i]);
    
    free(chr_name_array);
    free(chr_length_array);
    free(chr_name_all);
    free(chr_length_all);
    free(prfasta);
    free(chr_data);

    return (EXIT_SUCCESS);
}

