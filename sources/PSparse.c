
#include "PSparse.h"

/*///////////////////////////////////////////////////////////////////*/
void printfasta(char *file, char **seqs,char *chr_name, unsigned long pmc,unsigned long lentotal){
    unsigned long i,j;
    FILE *fasta_file;
    SGZip fasta_file_gz;
    char filechr[LEN_FILENAME];
    
    if(*file=='\0')
        fasta_file = stdout;
    else {
        filechr[0]=0;
        strncat(filechr,file,LEN_FILENAME*sizeof(char));
        strncat(filechr,"_",LEN_FILENAME*sizeof(char));
        strncat(filechr,chr_name,LEN_FILENAME*sizeof(char));
        strncat(filechr,".fa.gz",LEN_FILENAME*sizeof(char));
        fasta_file=fzopen(filechr,"w",&fasta_file_gz);
        if (fasta_file == NULL){
            fprintf (stderr,"Unable to write fasta file: %s\n",file);
            exit(1);
        }
    }
    for (i=0;i<pmc;i++){
        fzprintf(fasta_file,&fasta_file_gz,">%lu\n",i);
        //fzprintf(fasta_file,&fasta_file_gz/*,"%s"*/,seqs[i]);/*the print is unfinished with this option! include a buffer.*/
        for(j=0;j<lentotal;j++) {
            fzprintf(fasta_file,&fasta_file_gz,"%c",seqs[i][j]);
        }
        fzprintf(fasta_file,&fasta_file_gz,"\n");
        if(pmc>10) ShowProgressBar(i,pmc, 30);
    }
    if(*file!='\0')
        fzclose(fasta_file, &fasta_file_gz);
}
/*///////////////////////////////////////////////////////////////////*/
void printtfasta(FILE *fasta_file, SGZip *fasta_file_gz, char *chr_name,unsigned long first, unsigned long lentotal, unsigned long pmc, struct typepos *typep, struct featpos *featp, unsigned long min_qual, unsigned long iter, unsigned long cnt2, double theta, long int *seed, unsigned long poolsize, struct list_chr_pileup *chr_data, unsigned long num_scaffolds_pileup, int argc, char **argv){
    
    unsigned long i,j;
    unsigned long jj1=0;
    unsigned long pos;
    long int val;
    struct MCMCsample next;
    unsigned long f01;
    char *seqs;
    int flag=0;
    char tmp[1];
    tmp[0]='N';
    
    if(first == 0) {
        int arg = 1;
        fzprintf(fasta_file,fasta_file_gz,"#PFcaller ");
        if( argc > 1 ) {
            while(arg < argc) {
                fzprintf(fasta_file,fasta_file_gz,"%s ",argv[arg]);
                arg++;
            }
        }
        fzprintf(fasta_file,fasta_file_gz,"\n");
        fzprintf(fasta_file,fasta_file_gz,"#NAMES: ");
        for (i=0;i<pmc;i++){
            fzprintf(fasta_file,fasta_file_gz,">%lu\t",i);
        }
        fzprintf(fasta_file,fasta_file_gz,"\n#CHR:POSITION\tGENOTYPES\n");
    }
    
    pos = 0;
    seqs = (char *)calloc(pmc+1, sizeof(char));
    for(j=0;j<num_scaffolds_pileup;j++) {
        if(strcmp(chr_data[j].cchrom,chr_name)==0) {
            flag=1;
            if(typep[chr_data[j].init].pos-1 >0) {
                memset(seqs,tmp[0],pmc); seqs[pmc] = '\0';
                for(jj1=0;jj1<typep[chr_data[j].init].pos-1;jj1++) {
                    pos=jj1+1;
                    fzprintf(fasta_file,fasta_file_gz,"%s:%ld\t%s\n",chr_name,pos,seqs);
                    if(jj1%BLOCKBASES == 0) {
                        if(chr_data[j].end-chr_data[j].init > 10)
                            ShowProgressBar(jj1,chr_data[j].end-chr_data[j].init, 30);
                    }
                }
            }
            for(jj1=chr_data[j].init;jj1<=chr_data[j].end;jj1++) {
                memset(seqs,tmp[0],pmc); seqs[pmc] = '\0';
                while(typep[jj1].pos != pos+1) {
                    pos += 1;
                    fzprintf(fasta_file,fasta_file_gz,"%s:%ld\t%s\n",chr_name,pos,seqs);
                }
                pos=typep[jj1].pos;
                if(featp[typep[jj1].t].estimated.nr > 0){
                    next.p1 = next.p2 = next.p3 = 0; // added Initialization
                    
                    val = (long int)floor(ran1(seed)*(double)(iter-1));
                    next.nc = (long int)(double)floor(featp[typep[jj1].t].dist_estimated.nr[val]);
                    next.p1 = (long int)(double)floor(featp[typep[jj1].t].dist_estimated.nrFreq[0][val]);
                    next.p2 = (long int)(double)floor(featp[typep[jj1].t].dist_estimated.nrFreq[1][val]);
                    
                    f01=(next.p1+next.p2);
                    for (i=0;i<next.p1;i++)  {seqs[i] = typep[jj1].nucl[0];}
                    for (i=next.p1;i<f01;i++){seqs[i] = typep[jj1].nucl[1];}
                }
                fzprintf(fasta_file,fasta_file_gz,"%s:%ld\t%s\n",chr_name,pos,seqs);
                if(jj1%BLOCKBASES == 0) {
                    if(chr_data[j].end-chr_data[j].init > 10)
                        ShowProgressBar(jj1,chr_data[j].end-chr_data[j].init, 30);
                }
                fflush(stdout);
            }
            if(typep[chr_data[j].end].pos-1 < lentotal) {
                memset(seqs,tmp[0],pmc); seqs[pmc] = '\0';
                for(jj1=typep[chr_data[j].end].pos;jj1<lentotal;jj1++) {
                    pos=jj1+1;
                    fzprintf(fasta_file,fasta_file_gz,"%s:%ld\t%s\n",chr_name,pos,seqs);
                }
                if(jj1%BLOCKBASES == 0) {
                    if(chr_data[j].end-chr_data[j].init > 10)
                        ShowProgressBar(jj1,chr_data[j].end-chr_data[j].init, 30);
                }
            }
        }
    }
    if(flag==0) {/*print required chromosomes although they are not in the dataset*/
        memset(seqs,tmp[0],pmc); seqs[pmc] = '\0';
        for(jj1=0;jj1<lentotal;jj1++) {
            pos=jj1+1;
            fzprintf(fasta_file,fasta_file_gz,"%s:%ld\t%s\n",chr_name,pos,seqs);
        }
    }
    free(seqs);
}
/*///////////////////////////////////////////////////////////////////*/
void printgVCF(char *file, char *chr_name, struct typepos *typep, struct featpos *featp,unsigned long first,long int lentotal,unsigned long min_qual, unsigned long cnt2, double theta, unsigned long pmc, unsigned long poolsize, unsigned long iter, long int *seed, struct list_chr_pileup *chr_data, unsigned long num_scaffolds_pileup) {
    
    unsigned long jjc;
    FILE *gvcf_file;
    SGZip gvcf_file_gz;
    unsigned long jj1=0;
    char alt[4];
    char alt2[4];
    char *seqs;
    int flag=0;
    char tmp[1];
    tmp[0]='N';
    
    alt[0]=alt[1]=alt[2]=alt[3]=0;
    alt2[0]=alt2[1]=alt2[2]=alt2[3]=0;
    
    if(*file=='\0')
        gvcf_file = stdout;
    else {
        if(first==0) {
            gvcf_file=fzopen(file,"wb",&gvcf_file_gz);
            if (gvcf_file == NULL){
                fprintf (stderr,"Unable to write gvcf file: %s\n",file);
                exit(1);
            }
        }
        else {
            gvcf_file=fzopen(file,"ab",&gvcf_file_gz);
            if (gvcf_file == NULL){
                fprintf (stderr,"Unable to write gvcf file: %s\n",file);
                exit(1);
            }
        }
    }
    
    if(first==0) {
        fzprintf(gvcf_file,&gvcf_file_gz,"##fileformat=VCFv4.1\n");
        fzprintf(gvcf_file,&gvcf_file_gz,"##FILTER=<ID=LowRDQ,Description=\"Low read depth/quality\">\n");
        fzprintf(gvcf_file,&gvcf_file_gz,"##INFO=<ID=TP,Number=1,Type=Integer,Description=\"Number of Total Samples\">\n");
        fzprintf(gvcf_file,&gvcf_file_gz,"##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">\n");
        fzprintf(gvcf_file,&gvcf_file_gz,"##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n");
        fzprintf(gvcf_file,&gvcf_file_gz,"##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">\n");
        fzprintf(gvcf_file,&gvcf_file_gz,"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
        fzprintf(gvcf_file,&gvcf_file_gz,"##ALT=<ID=NON_REF,Description=\"Represents any possible alternative allele at this location\">\n");
        fzprintf(gvcf_file,&gvcf_file_gz,"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tPOOL00000\n");
    }
    
    seqs = (char *)calloc(pmc+1, sizeof(char));
    for(jjc=0;jjc<num_scaffolds_pileup;jjc++) {
        if (strcmp(chr_data[jjc].cchrom,chr_name)==0) {
            flag=1;
            if(typep[chr_data[jjc].init].pos-1 >0) {
                memset(seqs,tmp[0],pmc); seqs[pmc] = '\0';
                for(jj1=0;jj1<typep[chr_data[jjc].init].pos-1;jj1++) {
                    function_gvcf_results_null(jj1, gvcf_file, &gvcf_file_gz, seqs, chr_name, pmc);
                }
            }
            for(jj1=chr_data[jjc].init;jj1<=chr_data[jjc].end;jj1++) {
                function_gvcf_results(jj1, chr_data, typep, featp, gvcf_file, &gvcf_file_gz, seqs, chr_name, cnt2, pmc,iter,seed,lentotal);
                if(jj1%BLOCKBASES == 0) {
                    if(chr_data[jjc].end-chr_data[jjc].init > 10)
                        ShowProgressBar(jj1,chr_data[jjc].end-chr_data[jjc].init, 30);
                }
                fflush(stdout);
            }
            if(typep[chr_data[jjc].end].pos-1 < lentotal) {
                memset(seqs,tmp[0],pmc); seqs[pmc] = '\0';
                for(jj1=typep[chr_data[jjc].end].pos;jj1<lentotal;jj1++) {
                    function_gvcf_results_null(jj1, gvcf_file, &gvcf_file_gz, seqs, chr_name, pmc);
                }
            }
        }
    }
    if(flag==0) {/*print required chromosomes although they are not in the dataset*/
        memset(seqs,tmp[0],pmc); seqs[pmc] = '\0';
        for(jj1=0;jj1<lentotal;jj1++) {
            function_gvcf_results_null(jj1, gvcf_file, &gvcf_file_gz, seqs, chr_name, pmc);        }
    }
    free(seqs);
    
    if(*file != '\0')
        fzclose(gvcf_file, &gvcf_file_gz);
    
}
/*/////////////////////////////////////////*/
void function_gvcf_results_null(unsigned long jj1, FILE *gvcf_file, SGZip *gvcf_file_gz, char *seqs, char *chr_name, unsigned long pmc){

    int j,q;
    unsigned long pos;
    char tmp[1];
    tmp[0]='N';

    pos=jj1+1;;
    memset(seqs,tmp[0],pmc); seqs[pmc] = '\0';
    fzprintf(gvcf_file,gvcf_file_gz,"%s\t",chr_name);/*CHROM*/
    fzprintf(gvcf_file,gvcf_file_gz,"%ld\t",pos);/*POS*/
    fzprintf(gvcf_file,gvcf_file_gz,".\t");/*ID*/
    fzprintf(gvcf_file,gvcf_file_gz,"%c\t",seqs[0]);/*REF NULL*/
    fzprintf(gvcf_file,gvcf_file_gz,".\t");
    fzprintf(gvcf_file,gvcf_file_gz,".\t");/*QUAL:TO COMPLETE LATER*/
    /*FILTER*/
    fzprintf(gvcf_file,gvcf_file_gz,"LowRDQ\t");
    /*INFO*/
    fzprintf(gvcf_file,gvcf_file_gz,"TP=%ld,NS=%ld,DP=0,AF=0.000",pmc,0);
    fzprintf(gvcf_file,gvcf_file_gz,"\t");
    /*FORMAT*/
    fzprintf(gvcf_file,gvcf_file_gz,"GT\t");
    /*POOL00000*/
    q=0;
    for(j=0;j<pmc;j++) {
        if(q>0) fzprintf(gvcf_file,gvcf_file_gz,"/");
        fzprintf(gvcf_file,gvcf_file_gz,".");
        q++;
    }
    fzprintf(gvcf_file,gvcf_file_gz,"\n");
    
    return;
}
/*/////////////////////////////////////////*/
void function_gvcf_results(unsigned long jj1, struct list_chr_pileup *chr_data, struct typepos *typep, struct featpos *featp, FILE *gvcf_file, SGZip *gvcf_file_gz, char *seqs, char *chr_name, unsigned long cnt2, unsigned long pmc, unsigned long iter, long int *seed, long int lentotal){
    
    unsigned long i;
    int j,k,q;
    char alt[4];
    char alt2[4];
    unsigned long pos;
    long int val=0;
    struct MCMCsample next;
    unsigned long f01=0;
     char tmp[1];
    tmp[0]='N';
    
    next.nc = next.p1 = next.p2 = next.p3 = 0; // added Initialization
    pos=typep[jj1].pos;
    memset(seqs,tmp[0],pmc); seqs[pmc] = '\0';
     if(featp[typep[jj1].t].estimated.nr > 0){
        val = (long int)floor(ran1(seed)*(double)(iter-1));
        next.nc = (long int)(double)floor(featp[typep[jj1].t].dist_estimated.nr[val]);
        next.p1 = (long int)(double)floor(featp[typep[jj1].t].dist_estimated.nrFreq[0][val]);
        next.p2 = (long int)(double)floor(featp[typep[jj1].t].dist_estimated.nrFreq[1][val]);
        f01=(next.p1+next.p2);
        if(pos > lentotal) {
            printf("\nERROR: The given length of the chromosome is shorter than included in the mpileup data.\n");
            exit(1);
        }
        for (i=0;i<next.p1;i++)  {seqs[i] = typep[jj1].nucl[0];}
        for (i=next.p1;i<f01;i++){seqs[i] = typep[jj1].nucl[1];}
    }
    fzprintf(gvcf_file,gvcf_file_gz,"%s\t",chr_name);/*CHROM*/
    fzprintf(gvcf_file,gvcf_file_gz,"%ld\t",pos);/*POS*/
    fzprintf(gvcf_file,gvcf_file_gz,".\t");/*ID*/
    fzprintf(gvcf_file,gvcf_file_gz,"%c\t",seqs[0]);/*REF IN SEQ[0]!!!*/
    /*ALT*/
    alt[0]='A';alt[1]='T';alt[2]='C';alt[3]='G';
    alt2[0]=alt2[1]=alt2[2]=alt2[3]=0;
    q=0;
    for(j=0;j<pmc;j++) {
        for(k=0;k<4;k++) {
            if(seqs[j] == alt[k]) {
                if(j>0) {/*exclude the first sample because is in REF*/
                    if(q>0) fzprintf(gvcf_file,gvcf_file_gz,",");
                    fzprintf(gvcf_file,gvcf_file_gz,"%c",seqs[j]);
                    q+=1;
                }
                alt2[q] = alt[k];
                alt[k]='0';
                break;
            }
        }
    }
    if(q==0) fzprintf(gvcf_file,gvcf_file_gz,".");
    fzprintf(gvcf_file,gvcf_file_gz,"\t");
    fzprintf(gvcf_file,gvcf_file_gz,".\t");/*QUAL:TO COMPLETE LATER*/
    /*FILTER*/
    if(seqs[0]=='N') fzprintf(gvcf_file,gvcf_file_gz,"LowRDQ\t");
    else fzprintf(gvcf_file,gvcf_file_gz,"PASS\t");
    /*INFO*/
    if(featp[typep[jj1].t].estimated.nr > 0) {
        fzprintf(gvcf_file,gvcf_file_gz,"TP=%ld,NS=%ld,DP=%ld,",pmc,next.nc,featp[typep[jj1].t].observed.cov);
        fzprintf(gvcf_file,gvcf_file_gz,"AF=%.3f",featp[typep[jj1].t].estimated.nrFreq[0]); /*CHECK*/
        if(next.p2) {
            fzprintf(gvcf_file,gvcf_file_gz,",%.3f",featp[typep[jj1].t].estimated.nrFreq[1]);
        }
    }
    else {
        fzprintf(gvcf_file,gvcf_file_gz,"TP=%ld,NS=%ld,DP=0,AF=0.000",pmc,0);
    }
    fzprintf(gvcf_file,gvcf_file_gz,"\t");
    /*FORMAT*/
    fzprintf(gvcf_file,gvcf_file_gz,"GT\t");
    /*POOL00000*/
    q=0;
    for(j=0;j<pmc;j++) {
        if(q>0) fzprintf(gvcf_file,gvcf_file_gz,"/");
        int k=0;
        while(k < 4 && seqs[j] != alt2[k]) k++;
        if(k < 4 && seqs[j] == alt2[k]) fzprintf(gvcf_file,gvcf_file_gz,"%d",k);
        else fzprintf(gvcf_file,gvcf_file_gz,".");
        q++;
    }
    fzprintf(gvcf_file,gvcf_file_gz,"\n");
    
    return;

}
/*///////////////////////////////////////////////////////////////////*/
void printlikelihoods(char *file, struct typepos *typep, struct featpos *featp,char *chr_name, unsigned long first, unsigned long iter, unsigned long cnt, unsigned long cnt2, unsigned long pmc, int biallelic){
    unsigned long i/*,j*/;
    FILE *lik_file;
    SGZip lik_file_gz;
    char file_complete_table[LEN_FILENAME];
    char filelik[1024];
    long int jj1 = 0;
    long int jj2 = 0;
    long int jj3 = 0;
    int flag=0;
    
    filelik[0] = '\0';
    if(*file=='\0')
        lik_file = stdout;
    else {
        strcat(filelik,file);
        filelik[0]=0;
        filelik[0] = '\0';
        strncat(filelik,file,LEN_FILENAME*sizeof(char));

        strncat(filelik,"_lik.txt.gz",LEN_FILENAME*sizeof(char));
        if(first==0) {
            lik_file=fzopen(filelik,"wb",&lik_file_gz);
            if (lik_file == NULL){
                fprintf (stderr,"Unable to write lik file: %s\n",file);
                exit(1);
            }
        }
        else {
            lik_file=fzopen(filelik,"ab",&lik_file_gz);
            if (lik_file == NULL){
                fprintf (stderr,"Unable to write lik file: %s\n",file);
                exit(1);
            }
        }
    }
    
    if(first==0) {
        fzprintf(lik_file,&lik_file_gz, "#CHR:POS\tnEff\tAs\tCs\tGs\tTs\tfrA\tfrC\tfrG\tfrT\tfrPop\n"); /*\tErrA\tErrC\tErrG\tErrT\n");*/
    }
    
    while(jj1<cnt2) {
        if (featp[typep[jj1].t].estimated.nr>0 && strcmp(typep[jj1].cchrom,chr_name)==0){
            flag=1;
            fzprintf(lik_file,&lik_file_gz,"%s:%lu\t%.3lf\t",chr_name,typep[jj1].pos,featp[typep[jj1].t].estimated.nr);
            double *probs;
            probs = (double *)calloc(9,sizeof(double));
            for (i=0;i<4;i++){
                int index=base_to_num(typep[jj1].nucl[i]);
                if(index != 999)
                    probs[index]=featp[typep[jj1].t].estimated.nrFreq[i];
            }
            for (i=4;i<8;i++){
                int index=base_to_num(typep[jj1].nucl[i-4]);
                if(index != 999)
                     probs[index+4]=featp[typep[jj1].t].estimated.freqpool[i-4];
            }
            probs[8]=featp[typep[jj1].t].estimated.freqpop/(double)NPOPINTV;
            double pt=probs[4]+probs[5]+probs[6]+probs[7];
            if( pt < 1.0) { /*some errors of unobserved nt positions */
                double p4,p5,p6,p7;
                p4 = probs[4]/pt;
                p5 = probs[5]/pt;
                p6 = probs[6]/pt;
                p7 = probs[7]/pt;
                probs[4] = p4;
                probs[5] = p5;
                probs[6] = p6;
                probs[7] = p7;
            }
            fzprintf(lik_file,&lik_file_gz,"%.3lf\t%.3lf\t%.3lf\t%.3lf\t%.3lf\t%.3lf\t%.3lf\t%.3lf\t%.3lf\n",
                     probs[0],probs[1],probs[2],probs[3],
                     probs[4],probs[5],probs[6],probs[7],
                     probs[8]); //t%.3lf\n",
            free(probs);
        }
        jj1++;
        if(jj1%BLOCKSECTS == 0) if(cnt>10) ShowProgressBar(jj1,cnt2, 30);
    }
    if(cnt2==0) {
        fzprintf(lik_file, &lik_file_gz, "NO DATA\n");
    }
    if(*file != '\0')
        fzclose(lik_file, &lik_file_gz);
    
    if(first == 0) {
        flag = 0;
        if(*file=='\0')
            lik_file = stdout;
        else {
            file_complete_table[0]=0;
            file_complete_table[0] = '\0';
            strncat(file_complete_table,file,LEN_FILENAME*sizeof(char));
            strncat(file_complete_table,"_complete_table.txt.gz",LEN_FILENAME*sizeof(char));
            
            lik_file=fzopen(file_complete_table,"w",&lik_file_gz);
            if (lik_file == NULL){
                fprintf (stderr,"Unable to write lik file: %s\n",file_complete_table);
                exit(1);
            }
        }
        while(jj2 < cnt) {
            if(featp[jj2].observed.cov > 0) {
                flag = 1;
                /*show the nucleotide positions for each different case*/
                fzprintf(lik_file, &lik_file_gz, "npos: ");
                jj3=0;
                while(jj3 < cnt2) {
                    if(typep[jj3].t == jj2)
                        fzprintf(lik_file, &lik_file_gz, "%s:%ld ",typep[jj3].cchrom,typep[jj3].pos);
                    jj3++;
                }
                fzprintf(lik_file, &lik_file_gz, "\n");
                if(biallelic) {
                    //show observed values for each different case
                    fzprintf(lik_file, &lik_file_gz, "RD\tNR1\tNR2\tNR3\tNR4\tpERR1\tpERR2\tpERR3\tpERR4\n");
                    float *err;
                    err = (float *)calloc(4,sizeof(float));
                    for(i=0;i<4;i++) {
                        if(featp[jj2].observed.freq[i]>0) err[i]=featp[jj2].observed.error_mean[i];
                    }
                    fzprintf(lik_file, &lik_file_gz, "%ld\t%ld\t%ld\t%ld\t%ld\t%.5f\t%.5f\t%.5f\t%.5f\n", featp[jj2].observed.cov, featp[jj2].observed.freq[0], featp[jj2].observed.freq[1], featp[jj2].observed.freq[2], featp[jj2].observed.freq[3], err[0], err[1], err[2], err[3]);
                    free(err);
                    //show results of the MHMCMC for each different case
                    if(featp[jj2].observed.cov) {
                        fzprintf(lik_file, &lik_file_gz, "ITER\tNC\tNC1\tNC2\tFRPOOL\tFRPOP\tNERR1\tNERR2\tNERR1->3\tNERR2->3\tlnPROB\n");
                        for (i=0;i<iter;i++){
                            fzprintf(lik_file, &lik_file_gz, "%lu\t%lu\t%lu\t%lu\t%.5f\t%.5f\t%lu\t%lu\t%lu\t%lu\t%.5f\n", i, featp[jj2].dist_estimated.nr[i], featp[jj2].dist_estimated.nrFreq[0][i], featp[jj2].dist_estimated.nrFreq[1][i], featp[jj2].dist_estimated.freqpool[0][i]/(double)pmc,featp[jj2].dist_estimated.freqpop[i]/(double)NPOPINTV, featp[jj2].dist_estimated.nerr[0][i], featp[jj2].dist_estimated.nerr[1][i], featp[jj2].dist_estimated.nerr[2][i], featp[jj2].dist_estimated.nerr[3][i], featp[jj2].dist_estimated.lnprob[i]);
                        }
                    }
                }
            }
            jj2++;
            if(jj2%BLOCKSECTS == 0) if(cnt>10) ShowProgressBar(jj2,cnt, 30);

            fzprintf(lik_file, &lik_file_gz, "##############################################################");
            fzprintf(lik_file, &lik_file_gz, "##############################################################");
            fzprintf(lik_file, &lik_file_gz, "##############################################################\n");
            fzprintf(lik_file, &lik_file_gz, "##############################################################");
            fzprintf(lik_file, &lik_file_gz, "##############################################################");
            fzprintf(lik_file, &lik_file_gz, "##############################################################\n");
        }
        if(flag==0) {
            fzprintf(lik_file, &lik_file_gz, "NO DATA\n");
            fzprintf(lik_file, &lik_file_gz, "##############################################################");
            fzprintf(lik_file, &lik_file_gz, "##############################################################");
            fzprintf(lik_file, &lik_file_gz, "##############################################################\n");
            fzprintf(lik_file, &lik_file_gz, "##############################################################");
            fzprintf(lik_file, &lik_file_gz, "##############################################################");
            fzprintf(lik_file, &lik_file_gz, "##############################################################\n");
        }
        if(*file != '\0')
            fzclose(lik_file, &lik_file_gz);
    }
}
/*/////////////////////////////////////////*/
void Freqs2fasta(struct typepos *typep, struct featpos *featp, char **seqs, char *prfasta,char *chr_name, unsigned long min_qual, unsigned long iter, unsigned long lentotal, unsigned long cnt2, double theta, long int *seed, unsigned long pmc, unsigned long poolsize, struct list_chr_pileup *chr_data, unsigned long num_scaffolds_pileup){
    unsigned long i,j;
    char tmp[1];
    unsigned long jj1=0;
    unsigned long pos;
    long int val;
    struct MCMCsample next;
    unsigned long f01;
    *tmp='N';
    
    if(prfasta[0]=='f') {
        for (i=0;i<pmc;i++){
            memset(seqs[i],*tmp,lentotal+1);
            memset(&seqs[i][lentotal],'\0',1);
        }
    }
    if(prfasta[0]=='t' || prfasta[0]=='g') {
        for (i=0;i<lentotal+1;i++){
            memset(seqs[i],*tmp,pmc);
            memset(&seqs[i][pmc],'\0',1);
        }
    }
    
    for(j=0;j<num_scaffolds_pileup;j++) {
        if (strcmp(chr_data[j].cchrom,chr_name)==0) {
            jj1 = chr_data[j].init;
            while(jj1<=chr_data[j].end) {
                if (strcmp(typep[jj1].cchrom,chr_name)==0){
                    if(featp[typep[jj1].t].estimated.nr > 0){
                        next.p1 = next.p2 = next.p3 = 0; // added Initialization
                        pos=typep[jj1].pos-1;
                        
                        val = (long int)floor(ran1(seed)*(double)(iter-1));
                        next.nc = (long int)(double)floor(featp[typep[jj1].t].dist_estimated.nr[val]);
                        next.p1 = (long int)(double)floor(featp[typep[jj1].t].dist_estimated.nrFreq[0][val]);
                        next.p2 = (long int)(double)floor(featp[typep[jj1].t].dist_estimated.nrFreq[1][val]);
                        
                        f01=(next.p1+next.p2);
                        if(prfasta[0]=='f') {
                            if(pos > lentotal) {
                                printf("\nERROR: The given length of the chromosome is shorter than included in the mpileup data.\n");
                                exit(1);
                            }
                            for (i=0;i<next.p1;i++)  {seqs[i][pos] = typep[jj1].nucl[0];}
                            for (i=next.p1;i<f01;i++){seqs[i][pos] = typep[jj1].nucl[1];}
                        }
                        if(prfasta[0]=='t' || prfasta[0]=='g') {
                            if(pos > lentotal) {
                                printf("\nERROR: The given length of the chromosome is shorter than included in the mpileup data.\n");
                                exit(1);
                            }
                            //BE CAREFUL IF THE USER INCLUDE A SHORTER LENGTH OF CHR THAN EXPECTED
                            for (i=0;i<next.p1;i++)  {seqs[pos][i] = typep[jj1].nucl[0];}
                            for (i=next.p1;i<f01;i++){seqs[pos][i] = typep[jj1].nucl[1];}
                        }
                    }
                }
                jj1++;
            }
        }
    }
}
/*/////////////////////////////////////////*/
void Rawdata2fasta(struct typepos *typep, struct featpos *featp, char **seqs, char *prfasta,char *chr_name, unsigned long min_qual, unsigned long iter, unsigned long lentotal, unsigned long cnt2, double theta, long int *seed, unsigned long pmc, unsigned long poolsize, struct list_chr_pileup *chr_data, unsigned long num_scaffolds_pileup,unsigned long max_cov){
    unsigned long i,j;
    char tmp[1];
    unsigned long jj1=0;
    unsigned long pos;
    struct MCMCsample next;
    unsigned long f01;
    *tmp='N';
    
    if(prfasta[0]=='f') {
        for (i=0;i<max_cov;i++){
            memset(seqs[i],*tmp,lentotal+1);
            memset(&seqs[i][lentotal],'\0',1);
        }
    }
    if(prfasta[0]=='t' || prfasta[0]=='g') {
        for (i=0;i<lentotal+1;i++){
            memset(seqs[i],*tmp,max_cov);
            memset(&seqs[i][max_cov],'\0',1);
        }
    }
    
    for(j=0;j<num_scaffolds_pileup;j++) {
        if (strcmp(chr_data[j].cchrom,chr_name)==0) {
            jj1 = chr_data[j].init;
            while(jj1<=chr_data[j].end) {
                if (strcmp(typep[jj1].cchrom,chr_name)==0){
                    if(featp[typep[jj1].t].observed.cov > 0){
                        next.p1 = next.p2 = next.p3 = 0; // added Initialization
                        pos=typep[jj1].pos-1;
                        
                        next.nc = (long int)(double)floor(featp[typep[jj1].t].observed.cov);
                        next.p1 = (long int)(double)floor(featp[typep[jj1].t].observed.freq[0]);
                        next.p2 = (long int)(double)floor(featp[typep[jj1].t].observed.freq[1]);
                        
                        f01=(next.p1+next.p2);
                        if(prfasta[0]=='f') {
                            if(pos > lentotal) {
                                printf("\nERROR: The given length of the chromosome is shorter than included in the mpileup data.\n");
                                exit(1);
                            }
                            for (i=0;i<next.p1;i++)  {seqs[i][pos] = typep[jj1].nucl[0];}
                            for (i=next.p1;i<f01;i++){seqs[i][pos] = typep[jj1].nucl[1];}
                        }
                        if(prfasta[0]=='t' || prfasta[0]=='g') {
                            if(pos > lentotal) {
                                printf("\nERROR: The given length of the chromosome is shorter than included in the mpileup data.\n");
                                exit(1);
                            }
                            //BE CAREFUL IF THE USER INCLUDE A SHORTER LENGTH OF CHR THAN EXPECTED
                            for (i=0;i<next.p1;i++)  {seqs[pos][i] = typep[jj1].nucl[0];}
                            for (i=next.p1;i<f01;i++){seqs[pos][i] = typep[jj1].nucl[1];}
                        }
                    }
                }
                jj1++;
            }
        }
    }
}
/*/////////////////////////////////////////*/
int  read_line_pileup(FILE * bam_file, SGZip * bam_file_gz, unsigned long * cov, unsigned long * pos_base, unsigned long *freq, char cc, float *error_mean, char *cchrom, unsigned long max_cov, unsigned long min_cov, unsigned long min_qual, double min_error, double *noerror_prod, unsigned long nscaffolds, char **chr_name_array, long int *number_scaffold, unsigned char  *crefbase, int teststrands){
    
    int count_i,flag_scaffold;
    unsigned char *cline,*cpileup, *cqual;
    unsigned long count_j, n_ins, nseq;
    size_t nline;
    unsigned long nsc;
    int frqb[8];
    
    frqb[0]=frqb[1]=frqb[2]=frqb[3]=frqb[4]=frqb[5]=frqb[6]=frqb[7]=0;
    noerror_prod[0]=noerror_prod[1]=noerror_prod[2]=noerror_prod[3]=0.0;
    (*cov)=0;
    
    cpileup=(unsigned char *) malloc(40);
    cqual=(unsigned char *) malloc(20);
    cline=(unsigned char *) malloc(1);
    
    nseq=-1; 
    *pos_base=-1;
    *number_scaffold = -1;

    memset( cchrom, '\0', MSP_MAX_NAME );
    memset( cqual, '\0', 20 );
    memset( cpileup, '\0', 40 );
    
    nline=1;
    
    getdelim2(&cline,&nline,9,bam_file, bam_file_gz);
    cchrom[0]=cc;
    sscanf((char *)cline,"%s\t",cchrom+1);
    cchrom[MSP_MAX_NAME-1] = '\0';
    getdelim2(&cline,&nline,9,bam_file, bam_file_gz);
    sscanf((char *)cline,"%lu\t",pos_base);
    getdelim2(&cline,&nline,9,bam_file, bam_file_gz);
    sscanf((char *)cline,"%c\t",crefbase);
    getdelim2(&cline,&nline,9,bam_file, bam_file_gz);
    sscanf((char *)cline,"%lu\t",&nseq);
    getdelim2(&cline,&nline,9,bam_file, bam_file_gz);
    if (strlen((char *)cline)>=40) cpileup=(unsigned char *)realloc(cpileup,strlen((char *)cline)+1);
    sscanf((char *)cline,"%s\t",cpileup);
    getdelim2(&cline,&nline,10,bam_file, bam_file_gz);
    if (strlen((char *)cline)>=20) cqual=(unsigned char *)realloc(cqual,strlen((char *)cline)+1);
    sscanf((char *)cline,"%s\n",cqual);

    flag_scaffold = 0;
    for(nsc=0;nsc<nscaffolds;nsc++) {
        if(strcmp(cchrom, chr_name_array[nsc])==0) {
            flag_scaffold = 1;
            *number_scaffold = nsc;
            break;
        }
    }
    if(flag_scaffold == 0) {
        (*cov)=0;
        free(cpileup);
        free(cqual);
        free(cline);
        return 0;
    }

    if (nseq == 0){
        (*cov)=0;
        free(cpileup);
        free(cqual);
        free(cline);
        return 0;
    }

    if (nseq  == -1 || cchrom[0] == '\0' || cqual[0] == '\0' || cpileup[0] == '\0' || (*crefbase != 'A' && *crefbase != 'a' && *crefbase != 'C' && *crefbase != 'c' && *crefbase != 'G' && *crefbase != 'g' && *crefbase != 'T' && *crefbase != 't' && *crefbase != 'N' && *crefbase != 'n' && *crefbase != 'N' && *crefbase != 'W' && *crefbase != 'S' && *crefbase != 'Y' && *crefbase != 'R' && *crefbase != 'M' && *crefbase != 'K' && *crefbase != 'B' && *crefbase != 'D' && *crefbase != 'H' && *crefbase != 'V') || *pos_base == -1){
        fprintf(stderr,"\t\tWarning: corrupted line at scaffold %s and position %lu (%c) in mpileup file (skipping it)\n",cchrom,*pos_base, *crefbase);
        (*cov)=0;
        free(cpileup);
        free(cqual);
        free(cline);
        return 0;
    }
    
    count_j=0;
    if((nseq>=min_cov)&&(nseq<=max_cov)) {
        for(count_i=0;count_i<strlen((char *)cpileup);count_i++) {
            float error=pow(10.0,-(double)(cqual[count_j]-33)/10.0);
            //printf("%f\n",error);
            switch (cpileup[count_i]) {
                case '^': count_i++; break;
                case '$': break;
                case '*': count_j++; break;
                case '+': {
                    count_i++;
                    for(n_ins=0;isdigit(cpileup[count_i])!=0;count_i++){
                        n_ins=n_ins*10+(cpileup[count_i]-48);
                    };
                    for(n_ins=n_ins-1;n_ins>0;n_ins--) {
                        count_i++;
                    };
                }; break;
                case '-': {
                    count_i++;
                    for(n_ins=0;isdigit(cpileup[count_i])!=0;count_i++){
                        n_ins=n_ins*10+(cpileup[count_i]-48);
                    };
                    for(n_ins=n_ins-1;n_ins>0;n_ins=n_ins-1) {
                        count_i++;
                    };
                }; break;
                case 'N': count_j++; break;
                case 'n': count_j++; break;
                case '.': {
                    if(cqual[count_j]>=min_qual+33){
                        if (*crefbase == 'A' || *crefbase == 'a'){
                            error_mean[0]+=(error);
                            noerror_prod[0] += log(1.0-error);
                            freq[0]++;
                            frqb[0]++;
                        }else if (*crefbase == 'C' || *crefbase == 'c'){
                            error_mean[1]+=(error);
                            noerror_prod[1] += log(1.0-error);
                            freq[1]++;
                            frqb[1]++;
                        }else if (*crefbase == 'G'|| *crefbase == 'g'){
                            error_mean[2]+=(error);
                            noerror_prod[2] += log(1.0-error);
                            freq[2]++;
                            frqb[2]++;
                        }else if (*crefbase == 'T' || *crefbase == 't'){
                            error_mean[3]+=(error);
                            noerror_prod[3] += log(1.0-error);
                            freq[3]++;
                            frqb[3]++;
                        }
                    }
                    count_j++;
                }; break;
                case ',': {
                    if(cqual[count_j]>=min_qual+33){
                        if (*crefbase == 'A' || *crefbase == 'a'){
                            error_mean[0]+=(error);
                            noerror_prod[0] += log(1.0-error);
                            freq[0]++;
                            frqb[4]++;
                        }else if (*crefbase == 'C' || *crefbase == 'c'){
                            error_mean[1]+=(error);
                            noerror_prod[1] += log(1.0-error);
                            freq[1]++;
                            frqb[5]++;
                        }else if (*crefbase == 'G'|| *crefbase == 'g'){
                            error_mean[2]+=(error);
                            noerror_prod[2] += log(1.0-error);
                            freq[2]++;
                            frqb[6]++;
                        }else if (*crefbase == 'T' || *crefbase == 't'){
                            error_mean[3]+=(error);
                            noerror_prod[3] += log(1.0-error);
                            freq[3]++;
                            frqb[7]++;
                        }
                    }
                    count_j++;
                }; break;
                case 'A': {
                    if(cqual[count_j]>=min_qual+33){
                        error_mean[0]+=(error);
                        noerror_prod[0] += log(1.0-error);
                        freq[0]++;
                        frqb[0]++;
                    }
                    count_j++;
                }; break;
                case 'C': {
                    if(cqual[count_j]>=min_qual+33){
                        error_mean[1]+=(error);
                        noerror_prod[1] += log(1.0-error);
                        freq[1]++;
                        frqb[1]++;
                    }
                    count_j++;
                }; break;
                case 'G': {
                    if(cqual[count_j]>=min_qual+33){
                        error_mean[2]+=(error);
                        noerror_prod[2] += log(1.0-error);
                        freq[2]++;
                        frqb[2]++;
                    }
                    count_j++;
                }; break;
                case 'T': {
                    if(cqual[count_j]>=min_qual+33) {
                        error_mean[3]+=(error);
                        noerror_prod[3] += log(1.0-error);
                        freq[3]++;
                        frqb[3]++;
                    }
                    count_j++;
                }; break;
                case 'a': {
                    if(cqual[count_j]>=min_qual+33){
                        error_mean[0]+=(error);
                        noerror_prod[0] += log(1.0-error);
                        freq[0]++;
                        frqb[4]++;
                    }
                    count_j++;
                }; break;
                case 'c': {
                    if(cqual[count_j]>=min_qual+33){
                        error_mean[1]+=(error);
                        noerror_prod[1] += log(1.0-error);
                        freq[1]++;
                        frqb[5]++;
                    }
                    count_j++;
                }; break;
                case 'g': {
                    if(cqual[count_j]>=min_qual+33){
                        error_mean[2]+=(error);
                        noerror_prod[2] += log(1.0-error);
                        freq[2]++;
                        frqb[6]++;
                    }
                    count_j++;
                }; break;
                case 't': {
                    if(cqual[count_j]>=min_qual+33) {
                        error_mean[3]+=(error);
                        noerror_prod[3] += log(1.0-error);
                        freq[3]++;
                        frqb[7]++;
                    }
                    count_j++;
                }; break;
            }
        }
        /*Sign Test included as indicated in Meacham et al. BMC Bioinformatics 2011, 12:451*/
        //Introduce Sign Test function for both strands and for all nucleotides:
        // if one nt no pass then eliminate that nt.
        if(teststrands==1) {
            if(test_strands(frqb,0.005)) {
                if(frqb[0] + frqb[4] == 0) {
                    error_mean[0] = 0.; noerror_prod[0] = 0.;
                    freq[0] = 0;
                }
                if(frqb[1] + frqb[5] == 0) {
                    error_mean[1] = 0.; noerror_prod[1] = 0.;
                    freq[1] = 0;
                }
                if(frqb[2] + frqb[6] == 0) {
                    error_mean[2] = 0.; noerror_prod[2] = 0.;
                    freq[2] = 0;
                }
                if(frqb[3] + frqb[7] == 0) {
                    error_mean[3] = 0.; noerror_prod[3] = 0.;
                    freq[3] = 0;
                }
            }
        }
        
        if((freq[0]+freq[1]+freq[2]+freq[3]>=min_cov)&&(freq[0]+freq[1]+freq[2]+freq[3]<=max_cov)) {
            int var=0;
            for (count_i=0;count_i<4;count_i++){
                if (freq[count_i]>0){
                    if (freq[count_i]>=1){var++;}
                    error_mean[count_i] = mind((float)(error_mean[count_i]/((double)freq[count_i])), min_error);
                }else{
                    error_mean[count_i] = 1e-300/*min_error*/;
                }
                (*cov)+=freq[count_i];
            }
            free(cpileup);
            free(cqual);
            free(cline);

            if ((*cov)>0){
                return 1;
            }else{
                return 0;
            }
        }
        else {
            (*cov)=0;
            free(cpileup);
            free(cqual);
            free(cline);
        }
    }
    else {
        (*cov)=0;
        free(cpileup);
        free(cqual);
        free(cline);

        return 0;
    }
    return 0;
}
/*////////////////////////////////////////////////////////////////////*/
int test_strands(int *frqb,double thr) {
    /*cdfbin: Milton Abramowitz and Irene Stegun, Handbook of Mathematical Functions 1966, Formula 26.5.24.*/
    int which=1;
    double pr = 0.5;
    double ompr = 0.5;
    double p,q;
    double s,xn;
    int status;
    double bound;
    int j;
    for(j=0;j<4;j++) {
        s = frqb[j];
        xn= frqb[j]+frqb[j+4];
        if(xn>5)  { /*smaller total will never be significant*/
            cdfbin(&which,&p,&q,&s,&xn,&pr,&ompr,&status,&bound);
            if(p<thr || q<thr) {
                frqb[j] = frqb[j+4] = 0;
            }
        }
    }
    return 1;
}
/*////////////////////////////////////////////////////////////////////*/
int base_to_num(int base){
    if (base == 'A'){
        return 0;
    }else if (base == 'C'){
        return 1;
    }else if (base == 'G'){
        return 2;
    }else if (base == 'T'){
        return 3;
    }
    return 999;
}
/*///////////////////////////////////////////////////////////////////*/
int base_to_num_char(int base,char *alt){
    if (base == alt[0]){
        return 0;
    }else if (base == alt[1]){
        return 1;
    }else if (base == alt[2]){
        return 2;
    }else if (base == alt[3]){
        return 3;
    }
    return 999;
}
/*///////////////////////////////////////////////////////////////////*/
void sortObs (struct Obs *tmp, struct Obs *sorted, char nucl[5], int outg){
    unsigned long  i,j;
    for (i=0;i<3;i++){
        for (j=(i+1);j<4;j++){
            if (sorted->freq[i]<sorted->freq[j]){
                unsigned long tmpfrq=sorted->freq[i];
                float tmperror=sorted->error_mean[i];
                char tmpnucl=*(nucl+i);
                sorted->freq[i]= sorted->freq[j];
                sorted->error_mean[i]= sorted->error_mean[j];
                sorted->freq[j]= tmpfrq;
                sorted->error_mean[j]= tmperror;
                memset(nucl+i,*(nucl+j),1);
                memset(nucl+j,tmpnucl,1);
            }
        }
    }
    for(j=1;j<4;j++) {
        if(sorted->freq[j]==0) {
            memset(nucl+4,*(nucl+j),1);
            memset(nucl+j,'N',1);
            break;
        }
    }
    if(outg) {
        if(sorted->outg == nucl[0]) sorted->outg = 0;
        else if(sorted->outg == nucl[1] || nucl[1]=='N') sorted->outg = 1;
        else sorted->outg = 2;
    }
    else
        sorted->outg = 0;
    return;
}
/*///////////////////////////////////////////////////////////////////*/
float Bquality_error(float perr) {
    float BQerr;
    BQerr = -10.0*log10(1.0-perr);
    return BQerr;
}
/*///////////////////////////////////////////////////////////////////*/
void ReduceInput(char *bampath, struct typepos **typep, struct featpos **featp, unsigned long max_cov, unsigned long min_cov, unsigned long min_qual, long int *ntotalreads, long int *ntotalreads_novariant, long int *npositions_novariant, long int *npositions_variant, double *n_nr, double min_error, double *sum_errors, double *sum_errors_variant, double *npositions_noerror_exp, unsigned long iter, unsigned long *cnt, unsigned long *cnt2, int biallelic, float qerr_range, struct list_chr_pileup **chr_data, unsigned long *num_scaffolds_pileup, unsigned long nscaffolds, char **chr_name_array, char **chr_length_array, double *sum_per, long int *npositions_per, double *thetant, unsigned long *theta_pos, unsigned long p,int outg, double *divergence, unsigned long *div_pos, int teststrands){
    
    FILE *bam_file;
    SGZip bam_file_gz;
    double noerror_prod[4];
    char *cchrom_;
    char *cchrom;
    unsigned long i,j;
    long int remain_scaffolds_to_read;
    long int number_scaffold = -1;
    cchrom =(char *) malloc(MSP_MAX_NAME);
    cchrom_=(char *) malloc(MSP_MAX_NAME);
    memset(cchrom_, '\0', MSP_MAX_NAME);
    
    char cc;
    *cnt = *cnt2 = 0;
    if( bampath[0] == '\0' ) {
        bam_file = stdin;
    }
    else {
        bam_file=fzopen(bampath,"r", &bam_file_gz);
    }
    if (bam_file == NULL){
        fprintf (stderr,"Unable to open the %s mpileup file\n",bampath);
        exit(0);
    }
    remain_scaffolds_to_read = nscaffolds;
    cc=fzgetc(bam_file, &bam_file_gz);
    
    for(i=1;(cc!=EOF && cc!=GZ_EOF); i++){
        if (cc != '\n'){
            struct Obs observed_tmp;
            struct Obs observed_sorted;
            char *nucl/*[5]="ACGT"*/;
            unsigned long pos;
            double theta_i;
            memset(cchrom, '\0', MSP_MAX_NAME);
            
            nucl = (char *)calloc(5,sizeof(char));
            //nucl[0]='N';nucl[1]='A';nucl[2]='C';nucl[3]='G';nucl[4]='T';
            nucl[0]='A';nucl[1]='C';nucl[2]='G';nucl[3]='T';//nucl[4]='T';

            *typep=(struct typepos *)realloc(*typep,sizeof(struct typepos) * (*cnt2+1));
            int repeat;
            
            observed_tmp.freq[0]=0;  observed_tmp.freq[1]=0; observed_tmp.freq[2]=0; observed_tmp.freq[3]=0;
            observed_tmp.error_mean[0]=0.0; observed_tmp.error_mean[1]=0.0; observed_tmp.error_mean[2]=0.0; observed_tmp.error_mean[3]=0.0;
            observed_tmp.outg=0;
            observed_tmp.cov=0;
            
            observed_sorted.freq[0]=0;  observed_sorted.freq[1]=0; observed_sorted.freq[2]=0; observed_sorted.freq[3]=0;
            observed_sorted.error_mean[0]=0.0; observed_sorted.error_mean[1]=0.0; observed_sorted.error_mean[2]=0.0; observed_sorted.error_mean[3]=0.0;
            observed_sorted.outg=0;
            observed_sorted.cov=0;

            (*typep)[*cnt2].pos=0; (*typep)[*cnt2].t=0; (*typep)[*cnt2].nucl[0]=0; (*typep)[*cnt2].nucl[1]=0; (*typep)[*cnt2].nucl[2]=0; (*typep)[*cnt2].nucl[3]=0; (*typep)[*cnt2].nucl[4]=0;
            (*typep)[*cnt2].cchrom = 0;//memset((*typep)[*cnt2].cchrom, 0, MSP_MAX_NAME);
            
            if(read_line_pileup(bam_file, &bam_file_gz, &observed_tmp.cov, &pos, observed_tmp.freq, cc, observed_tmp.error_mean,cchrom, max_cov, min_cov, min_qual, min_error,noerror_prod, nscaffolds, chr_name_array, &number_scaffold,&observed_tmp.outg, teststrands)){
                repeat=0;
                if (observed_tmp.cov == 0){
                    observed_tmp.error_mean[0]=1e-300/*min_error*/;
                    observed_tmp.error_mean[1]=1e-300/*min_error*/;
                    observed_tmp.error_mean[2]=1e-300/*min_error*/;
                    observed_tmp.error_mean[3]=1e-300/*min_error*/;
                    observed_tmp.outg=0;
                }
                
                for (j=0;j<4;j++){
                    observed_sorted.freq[j]=observed_tmp.freq[j];
                    observed_sorted.error_mean[j]=observed_tmp.error_mean[j];
                }
                observed_sorted.cov=observed_tmp.cov;
                observed_sorted.outg=observed_tmp.outg;

                sortObs (&observed_tmp, &observed_sorted, nucl, (int)outg);/*also modify the REF/outg to 0/1/2*/
                
                n_nr[(long int)observed_sorted.cov]++;
                
                if (observed_tmp.cov >= min_cov) {
                    /*introduce lynch approach: */
                    double per,per2,pest2,pest/*,llp,llm,lr*/;
                    double div_raw=0;
                    /*
                    unsigned long obs0,obs1;
                    per = 3.0/2.0 * (double)(observed_sorted.freq[1] + observed_sorted.freq[2] + observed_sorted.freq[3])/(double)observed_sorted.cov;
                    pest = 1.0;
                    */
                    per2 = 3.0/2.0 * (double)(observed_sorted.freq[2] + observed_sorted.freq[3])/(double)observed_sorted.cov;
                    pest2 = (((double)observed_sorted.freq[0]/(double)(observed_sorted.freq[0]+observed_sorted.freq[1])) * (1.0 - 2.0/3.0 * per2) - per2/3.0) / (1.0 - 4.0/3.0*per2);
                    /*estimation of divergence versus reference*/
                    if((int)observed_sorted.outg==0)
                        div_raw = (double)(observed_sorted.cov - observed_sorted.freq[0]) / (double)observed_sorted.cov;
                    if((int)observed_sorted.outg==1)
                        div_raw = (double)(observed_sorted.cov - observed_sorted.freq[1]) / (double)observed_sorted.cov;
                    if((int)observed_sorted.outg==2)
                        div_raw = (double)(observed_sorted.freq[0] + observed_sorted.freq[1]) / (double)observed_sorted.cov;
                    *divergence += div_raw;
                    /*
                    obs0 = observed_sorted.freq[0] + (double)observed_sorted.freq[0]/(observed_sorted.freq[0]+observed_sorted.freq[1]) *(observed_sorted.freq[2]+observed_sorted.freq[3]);
                    obs1 = observed_sorted.freq[1] + (double)observed_sorted.freq[1]/(observed_sorted.freq[0]+observed_sorted.freq[1]) *(observed_sorted.freq[2]+observed_sorted.freq[3]);
                    
                    if(pest2<1.0) {
                        llp = observed_sorted.freq[0] * log(pest2 * (1.0 - 4.0/3.0*per2)+per2/3.0) +
                              observed_sorted.freq[1] * log(pest2 * (4.0/3.0 * per2 - 1.0) + (1.0-per2));
                        if(per2) llp += (observed_sorted.freq[2] + observed_sorted.freq[3]) * log(per2/3.0);
                        llm = observed_sorted.freq[0] * log(1.0-per) + (observed_sorted.cov - observed_sorted.freq[0]) * log(per/3.0);
                        lr = 2.0 * (llp - llm);
                        //introduce a funcion that gives the threshold value depending on the minimum seq error
                        if(lr >= 6.635) { //error2 0.01
                            per = per2; pest = pest2;
                        }
                        else {
                            obs0 = observed_sorted.cov;
                            obs1 = 0;
                        }
                    }
                    */
                    /*estimation of theta from Tajima*/
                    per = per2; pest = pest2;
                    double nC;
                    nC = round(p * (1.0-pow(1.0-1.0/(double)p,(double)observed_sorted.cov))); //from ferretti: last eq. in combinatorics
                    theta_i = 2.0 * pest * (1.0-pest) * ((double)nC-1.0)/(double)nC;
                    /*
                    double  ncR,ncA;
                    unsigned long pRp,pAp;
                    pRp = (unsigned long)ceil(p*pest);
                    pAp = p - pRp;
                    ncR = (pRp * (1.0-pow(1.0-1.0/((double)pRp),(double)obs0))); //from ferretti: last eq. in combinatorics
                    ncA = (pAp * (1.0-pow(1.0-1.0/((double)pAp),(double)obs1))); //from ferretti: last eq. in combinatorics
                    theta_i = 2.0 * pest * (1.0-pest) * ((double)(ncR + ncA)-1.0)/(double)(ncR + ncA);
                    if(theta_i < 0.0) theta_i = 0.;
                    nC = ncR+ncA;
                    */
                    if(nC >= 2.0) {
                        *thetant += theta_i;
                        *theta_pos += 1;
                    }
                    if(nC >= 1.0)
                        *div_pos += 1;
                    
                    if(pest > 1.0/nC && pest < 1.0 - 1.0/nC) {
                        *sum_per += per;
                        *npositions_per += 1;
                        //printf("\n%.3f\t%.3f\t%ld",per,pest,observed_tmp.cov);
                    }
                    
                    *ntotalreads += (long int)observed_sorted.cov;
                    *sum_errors += observed_sorted.freq[0] * (observed_sorted.error_mean[0]) + observed_sorted.freq[1] * (observed_sorted.error_mean[1]) + observed_sorted.freq[2] * (observed_sorted.error_mean[2]) + observed_sorted.freq[3] * (observed_sorted.error_mean[3]);
                    *npositions_noerror_exp += exp(noerror_prod[0] + noerror_prod[1] + noerror_prod[2] + noerror_prod[3]);//expected invariant positions
                    if(observed_sorted.freq[0] == observed_sorted.cov) {
                        *ntotalreads_novariant += (long int)observed_sorted.cov;
                        if(observed_sorted.cov > 1) *npositions_novariant += (long int)1;
                    } else {
                        *npositions_variant += (long int)1;
                        *sum_errors_variant += observed_sorted.freq[0] * (observed_sorted.error_mean[0]) + observed_sorted.freq[1] * (observed_sorted.error_mean[1]) + observed_sorted.freq[2] * (observed_sorted.error_mean[2]) + observed_sorted.freq[3] * (observed_sorted.error_mean[3]);
                    }
                     /**/
                }
                else {
                    observed_sorted.cov=0;
                    observed_sorted.freq[0]=0;
                    observed_sorted.freq[1]=0;
                    observed_sorted.freq[2]=0;
                    observed_sorted.freq[3]=0;
                }
                
                /*construct a list with all required scaffolds and their positions initial and end*/
                if(*num_scaffolds_pileup == 0) {
                    remain_scaffolds_to_read -=1;
                    num_scaffolds_pileup[0]=num_scaffolds_pileup[0]+1;
                    *chr_data = (struct list_chr_pileup *)calloc(*num_scaffolds_pileup,sizeof(struct list_chr_pileup));
                    (*chr_data)[*num_scaffolds_pileup-1].cchrom = chr_name_array[number_scaffold];//memcpy((*chr_data)[*num_scaffolds_pileup-1].cchrom, cchrom, MSP_MAX_NAME);
                    (*chr_data)[*num_scaffolds_pileup-1].init = *cnt2;
                    (*chr_data)[*num_scaffolds_pileup-1].end = *cnt2;
                }
                else {
                    if (strcmp((*chr_data)[*num_scaffolds_pileup-1].cchrom,cchrom)==0) {
                        (*chr_data)[*num_scaffolds_pileup-1].end = *cnt2;
                    }
                    else {
                        remain_scaffolds_to_read -=1;
                        if(remain_scaffolds_to_read>=0) {
                            num_scaffolds_pileup[0]=num_scaffolds_pileup[0]+1;
                            *chr_data = (struct list_chr_pileup *)realloc(*chr_data,*num_scaffolds_pileup * sizeof(struct list_chr_pileup));
                            (*chr_data)[*num_scaffolds_pileup-1].cchrom = chr_name_array[number_scaffold];//memcpy((*chr_data)[*num_scaffolds_pileup-1].cchrom, cchrom, MSP_MAX_NAME);
                            (*chr_data)[*num_scaffolds_pileup-1].init = *cnt2;
                            (*chr_data)[*num_scaffolds_pileup-1].end = *cnt2;
                        }
                    }
                }
                
                for (j=0;j<*cnt;j++){
                    if ((*featp)[j].observed.cov==observed_sorted.cov){
                        if ((*featp)[j].observed.freq[0]==observed_sorted.freq[0] &&
                            (*featp)[j].observed.freq[1]==observed_sorted.freq[1] &&
                            (*featp)[j].observed.freq[2]==observed_sorted.freq[2] &&
                            (*featp)[j].observed.freq[3]==observed_sorted.freq[3] &&
                            (*featp)[j].observed.outg==observed_sorted.outg &&
                            /**/
                            Bquality_error((*featp)[j].observed.error_mean[0]) <= Bquality_error(observed_sorted.error_mean[0]) + qerr_range &&
                            Bquality_error((*featp)[j].observed.error_mean[0]) >= Bquality_error(observed_sorted.error_mean[0]) - qerr_range &&
                            Bquality_error((*featp)[j].observed.error_mean[1]) <= Bquality_error(observed_sorted.error_mean[1]) + qerr_range &&
                            Bquality_error((*featp)[j].observed.error_mean[1]) >= Bquality_error(observed_sorted.error_mean[1]) - qerr_range &&
                            Bquality_error((*featp)[j].observed.error_mean[2]) <= Bquality_error(observed_sorted.error_mean[2]) + qerr_range &&
                            Bquality_error((*featp)[j].observed.error_mean[2]) >= Bquality_error(observed_sorted.error_mean[2]) - qerr_range &&
                            Bquality_error((*featp)[j].observed.error_mean[3]) <= Bquality_error(observed_sorted.error_mean[3]) + qerr_range &&
                            Bquality_error((*featp)[j].observed.error_mean[3]) >= Bquality_error(observed_sorted.error_mean[3]) - qerr_range ){
                            (*typep)[*cnt2].pos=pos;
                            (*typep)[*cnt2].t=j;
                            (*typep)[*cnt2].nucl[0]=nucl[0];
                            (*typep)[*cnt2].nucl[1]=nucl[1];
                            (*typep)[*cnt2].nucl[2]=nucl[2];
                            (*typep)[*cnt2].nucl[3]=nucl[3];
                            (*typep)[*cnt2].cchrom = chr_name_array[number_scaffold];//memcpy((*typep)[*cnt2].cchrom, cchrom, MSP_MAX_NAME);
                            repeat=1;
                            /*printf("observed_error_mean[0]=%f\tobserved_error_mean[1]=%f\tobserved_error_mean[2]=%f\tobserved_error_mean[3]=%f\n",
                             (*featp)[j].observed.error_mean[0],(*featp)[j].observed.error_mean[1],(*featp)[j].observed.error_mean[2],(*featp)[j].observed.error_mean[3]);*/
                            break;
                        }
                    }
                }
                if (repeat == 0){
                    *featp=(struct featpos *)realloc(*featp,sizeof(struct featpos) * (*cnt+1));
                    /*calloc the doubles with iter values of nr and nrFreq[4] de DistEst inside featpos*/
                    featp[0][*cnt].dist_estimated.nr = (unsigned long *) calloc(iter,sizeof(unsigned long));
                    featp[0][*cnt].dist_estimated.nrFreq[0] =(unsigned long *) calloc(iter,sizeof(unsigned long));
                    featp[0][*cnt].dist_estimated.nrFreq[1] =(unsigned long *) calloc(iter,sizeof(unsigned long));
                    featp[0][*cnt].dist_estimated.nerr[0] =(unsigned long *) calloc(iter,sizeof(unsigned long));
                    featp[0][*cnt].dist_estimated.nerr[1] =(unsigned long *) calloc(iter,sizeof(unsigned long));
                    featp[0][*cnt].dist_estimated.nerr[2] =(unsigned long *) calloc(iter,sizeof(unsigned long));
                    featp[0][*cnt].dist_estimated.nerr[3] =(unsigned long *) calloc(iter,sizeof(unsigned long));
                    featp[0][*cnt].dist_estimated.freqpool[0] =(unsigned long *) calloc(iter,sizeof(unsigned long));
                    featp[0][*cnt].dist_estimated.freqpool[1] =(unsigned long *) calloc(iter,sizeof(unsigned long));
                    featp[0][*cnt].dist_estimated.freqpop =(unsigned long  *) calloc(iter,sizeof(unsigned long ));
                    featp[0][*cnt].dist_estimated.lnprob = (double *) calloc(iter,sizeof(double));
                    
                    (*featp)[*cnt].observed.freq[0]=observed_sorted.freq[0];
                    (*featp)[*cnt].observed.freq[1]=observed_sorted.freq[1];
                    (*featp)[*cnt].observed.freq[2]=observed_sorted.freq[2];
                    (*featp)[*cnt].observed.freq[3]=observed_sorted.freq[3];
                    (*featp)[*cnt].observed.error_mean[0]=observed_sorted.error_mean[0];
                    (*featp)[*cnt].observed.error_mean[1]=observed_sorted.error_mean[1];
                    (*featp)[*cnt].observed.error_mean[2]=observed_sorted.error_mean[2];
                    (*featp)[*cnt].observed.error_mean[3]=observed_sorted.error_mean[3];
                    (*featp)[*cnt].observed.cov=observed_sorted.cov;
                    (*featp)[*cnt].observed.outg=observed_sorted.outg;
                    /**/
                    (*featp)[*cnt].estimated.nr=0;
                    (*featp)[*cnt].estimated.nrFreq[0]=0.;
                    (*featp)[*cnt].estimated.nrFreq[1]=0.;
                    (*featp)[*cnt].estimated.nrFreq[2]=0.;
                    (*featp)[*cnt].estimated.nrFreq[3]=0.;
                    (*featp)[*cnt].estimated.nerr[0]=0.;
                    (*featp)[*cnt].estimated.nerr[1]=0.;
                    (*featp)[*cnt].estimated.nerr[2]=0.;
                    (*featp)[*cnt].estimated.nerr[3]=0.;
                    (*featp)[*cnt].estimated.freqpool[0]=0.;
                    (*featp)[*cnt].estimated.freqpool[1]=0.;
                    (*featp)[*cnt].estimated.freqpool[2]=0.;
                    (*featp)[*cnt].estimated.freqpool[3]=0.;
                    (*featp)[*cnt].estimated.freqpop=0.;
                    /**/
                    (*typep)[*cnt2].pos=pos;
                    (*typep)[*cnt2].t=*cnt;
                    (*typep)[*cnt2].nucl[0]=nucl[0];
                    (*typep)[*cnt2].nucl[1]=nucl[1];
                    (*typep)[*cnt2].nucl[2]=nucl[2];
                    (*typep)[*cnt2].nucl[3]=nucl[3];
                    (*typep)[*cnt2].cchrom = chr_name_array[number_scaffold];//memcpy((*typep)[*cnt2].cchrom, cchrom, MSP_MAX_NAME);
                    
                    *cnt += 1;
                }
                free(nucl);
                *cnt2 += 1;
            }else{
                if(*num_scaffolds_pileup == 0) {
                    remain_scaffolds_to_read -=1;
                    num_scaffolds_pileup[0]=num_scaffolds_pileup[0]+1;
                    *chr_data = (struct list_chr_pileup *)calloc(*num_scaffolds_pileup,sizeof(struct list_chr_pileup));
                    (*chr_data)[*num_scaffolds_pileup-1].cchrom = chr_name_array[number_scaffold];
                    (*chr_data)[*num_scaffolds_pileup-1].init = *cnt2;
                    (*chr_data)[*num_scaffolds_pileup-1].end = *cnt2;
                }
                else {
                    if(strcmp((*chr_data)[*num_scaffolds_pileup-1].cchrom,cchrom)!=0) {
                        remain_scaffolds_to_read -=1;
                        if(remain_scaffolds_to_read>=0) {
                            num_scaffolds_pileup[0]=num_scaffolds_pileup[0]+1;
                            *chr_data = (struct list_chr_pileup *)realloc(*chr_data,*num_scaffolds_pileup * sizeof(struct list_chr_pileup));
                            (*chr_data)[*num_scaffolds_pileup-1].cchrom = chr_name_array[number_scaffold];
                            (*chr_data)[*num_scaffolds_pileup-1].init = *cnt2;
                            (*chr_data)[*num_scaffolds_pileup-1].end = *cnt2;
                        }
                    }
                }
                if(number_scaffold >= 0) {
                    repeat=0;
                    for (j=0;j<*cnt;j++){
                        if ((*featp)[j].observed.cov==0){
                            (*typep)[*cnt2].pos=pos;
                            (*typep)[*cnt2].t=j;
                            (*typep)[*cnt2].nucl[0]=nucl[0];
                            (*typep)[*cnt2].nucl[1]=nucl[1];
                            (*typep)[*cnt2].nucl[2]=nucl[2];
                            (*typep)[*cnt2].nucl[3]=nucl[3];
                            (*typep)[*cnt2].cchrom = chr_name_array[number_scaffold];//memcpy((*typep)[*cnt2].cchrom, cchrom, MSP_MAX_NAME);
                            
                            repeat=1;
                            break;
                        }
                    }
                    if (repeat == 0){
                        *featp=(struct featpos *)realloc(*featp,sizeof(struct featpos) * (*cnt+1));
                         //calloc the doubles with iter values of nr and nrFreq[4] de DistEst inside featpos
                        featp[0][*cnt].dist_estimated.nr = (unsigned long *) calloc(iter,sizeof(unsigned long));
                        featp[0][*cnt].dist_estimated.nrFreq[0] =(unsigned long *) calloc(iter,sizeof(unsigned long));
                        featp[0][*cnt].dist_estimated.nrFreq[1] =(unsigned long *) calloc(iter,sizeof(unsigned long));
                        featp[0][*cnt].dist_estimated.nerr[0] =(unsigned long *) calloc(iter,sizeof(unsigned long));
                        featp[0][*cnt].dist_estimated.nerr[1] =(unsigned long *) calloc(iter,sizeof(unsigned long));
                        featp[0][*cnt].dist_estimated.nerr[2] =(unsigned long *) calloc(iter,sizeof(unsigned long));
                        featp[0][*cnt].dist_estimated.nerr[3] =(unsigned long *) calloc(iter,sizeof(unsigned long));
                        featp[0][*cnt].dist_estimated.freqpool[0] =(unsigned long *) calloc(iter,sizeof(unsigned long));
                        featp[0][*cnt].dist_estimated.freqpool[1] =(unsigned long *) calloc(iter,sizeof(unsigned long));
                        featp[0][*cnt].dist_estimated.freqpop =(unsigned long  *) calloc(iter,sizeof(unsigned long ));
                        featp[0][*cnt].dist_estimated.lnprob = (double *) calloc(iter,sizeof(double));

                        (*featp)[*cnt].observed.cov=0;
                        (*featp)[*cnt].observed.error_mean[0]=observed_sorted.error_mean[0];
                        (*featp)[*cnt].observed.error_mean[1]=observed_sorted.error_mean[1];
                        (*featp)[*cnt].observed.error_mean[2]=observed_sorted.error_mean[2];
                        (*featp)[*cnt].observed.error_mean[3]=observed_sorted.error_mean[3];
                        (*featp)[*cnt].observed.freq[0]=0;
                        (*featp)[*cnt].observed.freq[1]=0;
                        (*featp)[*cnt].observed.freq[2]=0;
                        (*featp)[*cnt].observed.freq[3]=0;
                        (*featp)[*cnt].observed.outg=0;
                        (*featp)[*cnt].estimated.nr=0;
                        (*featp)[*cnt].estimated.nrFreq[0]=0.;
                        (*featp)[*cnt].estimated.nrFreq[1]=0.;
                        (*featp)[*cnt].estimated.nrFreq[2]=0.;
                        (*featp)[*cnt].estimated.nrFreq[3]=0.;
                        (*featp)[*cnt].estimated.nerr[0]=0.;
                        (*featp)[*cnt].estimated.nerr[1]=0.;
                        (*featp)[*cnt].estimated.nerr[2]=0.;
                        (*featp)[*cnt].estimated.nerr[3]=0.;
                        (*featp)[*cnt].estimated.freqpool[0]=0.;
                        (*featp)[*cnt].estimated.freqpool[1]=0.;
                        (*featp)[*cnt].estimated.freqpool[2]=0.;
                        (*featp)[*cnt].estimated.freqpool[3]=0.;
                        (*featp)[*cnt].estimated.freqpop=0.;
                        (*typep)[*cnt2].pos=pos;
                        (*typep)[*cnt2].t=*cnt;
                        (*typep)[*cnt2].cchrom = chr_name_array[number_scaffold];//memcpy((*typep)[*cnt2].cchrom, cchrom, MSP_MAX_NAME);
                        
                        *cnt += 1;
                    }
                }
                free(nucl);
                *cnt2 += 1;
            }
        }
        do{
            cc=fzgetc(bam_file, &bam_file_gz);
        }while (cc == '\n');
        
        if(strcmp(cchrom, cchrom_)!=0) {
            if(remain_scaffolds_to_read < 0) {
                break; /*exit from reading pileup file when there are no more required scaffolds*/
            }
            if(cchrom_[0]!='\0' && *num_scaffolds_pileup/* && number_scaffold >= 0*/)
                printf("\nNumber of different patterns retained: %ld / %ld ",*cnt,*cnt2);
            if(number_scaffold >= 0) {
                printf(" Scaffold %s \n",cchrom);
            }
            else {
                printf("\n Searching for demanded scaffolds (now parsing through scaffold %s \n) ",cchrom);
            }
            strcpy(cchrom_, cchrom);
        }
        if(i%BLOCKBASES == 0) {
            if(atol(chr_length_array[nscaffolds-1]) > 10)
            ShowProgressBar(i,atol(chr_length_array[nscaffolds-1]), 30);
        }
        fflush(stdout);
    }
    free(cchrom);
    free(cchrom_);
    printf("\n");
    if(bampath[0] !='\0')
        fzclose(bam_file, &bam_file_gz);
}
/*///////////////////////////////////////////////////////////////////*/
void usage(void)
{
  printf(PSCALLER);
  printf("\nUsage of PFcaller:\n");

  printf("\nPFcaller -h [help and exit]\n");
  printf("#INPUT/OUTPUT:\n");
  printf("         -i [input pileup file (nothing is stdin)]\n");
  printf("         -o [output file name (nothing is stdout, if extension is '.gz' then zipped)]\n");
  printf("#MANDATORY PARAMETERS:\n");
  printf("         -m [minimum read depth]\n");
  printf("         -M [maximum read depth]\n");
  printf("         -q [minumum Base Quality. (Phred)]\n");
  printf("         -p [ploidy]\n");
  printf("         -s [seed]\n");
  printf("         -S [name of the file containing the name(s) of scaffold(s) and their length (separated by a tab), one per line (ex. fai file)]\n");
  printf("#OPTIONAL PARAMETERS:\n");
  printf("         -f [Output format: '0':numeric, 'f':fasta, 't':transposed fasta, 'g':gVCF. DEFAULT: 't']\n");
  printf("         -r [prior dist: '1':snm, '2':uniform, '3':exponential, '4':uniform-no-theta. DEFAULT: '1']\n");/*'0':uniform over theta, */
  printf("         -t [prior value: theta/nt (-1 means auto-inferred from whole data). DEFAULT: -1]\n");
  printf("         -n [# MC iterations. DEFAULT: 100]\n");
/*printf("         -b [#burnin iterations in MCMC. DEFAULT: 10]\n");*/
/*printf("         -v [#interval between iterations in MCMC. DEFAULT: 5]\n");*/
  printf("         -e [sections in which are divided the SNPs in relation to BaseQuality (higher means lower precision). DEFAULT: 1]\n");
/*printf("         -C [combinatorics prob: # iterations in case ploidy>100. DEFAULT: 200]\n");*/
  printf("         -g [outgroup defined at reference (1/0). DEFAULT: 0]\n");
  printf("         -d [testing for having reads at both strands (1/0). DEFAULT: 0]\n");
  printf("\n");
}
/*///////////////////////////////////////////////////////////////////*/
int check_param(int argc, char **argv, unsigned long *max_cov, unsigned long *min_cov, unsigned long *min_qual, char *chr_length_all, char *chr_name_all, int *biallelic, char *file_in, char *file_out, int *priortype, long int *samplemhmcmc, long int *burnin_par, unsigned long *iter, unsigned long *nitermc, char *prfasta, unsigned long *lentotal, double *theta, long int *seed, unsigned long *poolsize, float *qerr_range, int *outg, int *teststrands){
  FILE *bam_file;
  SGZip bam_file_gz;

  /*mandatory*/
  *min_cov = 0;
  *max_cov = 0;
  *poolsize = 0;
  *seed = 0;
  *lentotal = 0;
  *min_qual= 123456;

  /*optional default*/
  *theta = -1.0;
  *iter = 100;
  *priortype = 1;
  *samplemhmcmc = 1;
  *burnin_par = 0;
  prfasta[0] = 't';
  *nitermc = NITER_MC;
  *biallelic = 1;
  *qerr_range = 1;
  *outg = 0;
  *teststrands = 0;

  int arg;
  if( argc > 1 )
    {
      arg = 1;
      while(arg < argc)
	{
	  if( argv[arg][0] != '-' )
            {
	      if(argv[arg][0] == '>')
		break;
	      printf(" argument should be -%s ?\n", argv[arg]);
	      usage();
	      exit(1);
            }

	  switch (argv[arg][1])
            {
	    case 'i': /* i Input File, el path */
	      arg++;
	      strcpy( file_in, argv[arg] );
            if( file_in[0] == '\0' ) {
                bam_file = stdin;
            }
            else {
              bam_file=fzopen(file_in,"r", &bam_file_gz);
              if (bam_file == NULL){
                fprintf (stderr,"Unable to open mpileup file for reading: %s\n",file_in);
                return(0);
              }
              fzclose (bam_file, &bam_file_gz);
            }
	      break;

	    case 'o' : /* o Output type*/
	      arg++;
	      strcpy( file_out, argv[arg] );
	      break;

	    case 'm' : /* minimum read depth */
	      arg++;
	      *min_cov = (int)atoi(argv[arg]);
	      if(*min_cov < 1 || *max_cov > LIMIT_FACT) {
              printf("\n Error in -m argument: The value must be a positive integer not larger than %d.\n",LIMIT_FACT);
              usage();
              exit(1);
	      }
        break;

	    case 'M' : /* maximum read depth */
	      arg++;
	      *max_cov = (int)atoi(argv[arg]);
	      if(*max_cov < 1 || *max_cov > LIMIT_FACT) {
              printf("\n Error in -M argument: The value must be a positive integer not larger than %d.\n",LIMIT_FACT);
              usage();
              exit(1);
	      }
	      break;

	    case 'p' : /* p Ploidy */
	      arg++;
	      *poolsize = (int)atoi(argv[arg]);
	      if(*poolsize < 1) {
              printf("\n Error in -p argument: The value must be a positive integer.");
              usage();
              exit(1);
	      }
	      break;
        case 't' : /* t theta */
            arg++;
            *theta = atof(argv[arg]);
            if(*theta > 1.0 || (*theta < 0 && *theta !=-1)) {
                printf("\n Error in -t argument: The value must be positive [0,1] or -1 (undefined).");
                usage();
                exit(1);
            }
            break;
	    case 's' : /* s Seed, positive value */
	      arg++;
	      *seed = (long int)atol(argv[arg]);
	      if(*seed < 1) {
		printf("\n Error in -s argument: The value must be a positive long integer.");
		usage();
		exit(1);
	      }
	      break;

	    case 'q' : /* minimum Base Quality */
	      arg++;
	      *min_qual = (int)atoi(argv[arg]);
	      if(*min_qual < 1) {
		printf("\n Error in -q argument: The value must be a positive integer.");
		usage();
		exit(1);
	      }
	      break;

	    case 'n' : /* n number of iterations*/
	      arg++;
	      *iter = (long int)atol( argv[arg] );
	      if(*iter < 1) {
		printf("\n Error in -n argument: The value must be a positive integer.");
		usage();
		exit(1);
	      }
	      break;
	    case 'r' : /* prior method  */
	      arg++;
	      /* if format is not recognized, error and out */
	      if(/*strcmp( argv[arg], "0" ) !=0 &&*/
             strcmp( argv[arg], "1" ) !=0 &&
             strcmp( argv[arg], "2" ) !=0 &&
             strcmp( argv[arg], "3" ) !=0 &&
             strcmp( argv[arg], "4" ) !=0)
		{
		  printf("\n Error: the argument -r has only the choices '1', '2', '3' or '4'."); /*'0', '*/
		  usage();
		  exit(1);
		}
	      /*else if(strcmp(argv[arg],"0")==0) {
		*priortype = 0;
	      }*/
	      else if(strcmp(argv[arg],"1")==0) {
		      *priortype = 1;
	      }
          else if(strcmp(argv[arg],"2")==0) {
              *priortype = 2;
          }
          else if(strcmp(argv[arg],"3")==0) {
              *priortype = 3;
          }
          else if(strcmp(argv[arg],"4")==0) {
              *priortype = 4;
          }

	      break;

	    case 'v' : /* v interval between iterations*/
	      arg++;
	      *samplemhmcmc = (long int)atol( argv[arg] );
	      if(*samplemhmcmc < 1) {
		printf("\n Error in -v argument: The value must be a positive integer.");
		usage();
		exit(1);
	      }
	      break;
	    case 'b' : /* b burning*/
	      arg++;
	      *burnin_par = (long int)atol( argv[arg] );
	      if(*burnin_par < 0) {
		printf("\n Error in -b argument: The value must be a zero or a positive integer.");
		usage();
		exit(1);
	      }
	      break;

        case 'f' : /* print fasta or tfasta files  */
            arg++;
            *prfasta = 0;
            /* if format is not recognized, error and out */
            if(strchr( argv[arg], '0' ) == NULL &&
               strchr( argv[arg], 'f' ) == NULL &&
               strchr( argv[arg], 'g' ) == NULL &&
               strchr( argv[arg], 't' ) == NULL)
            {
                printf("\n Error: the argument -f has only the choices '0', 'f', 't' or 'g'.");
                usage();
                exit(1);
            }
            prfasta[0]='\0';
            strcat(prfasta,argv[arg]);
            break;
        case 'B':
            arg++;
            *biallelic = (int)atoi(argv[arg]);
            if(*biallelic < 0 && *biallelic != 1) {
                printf("\n Error in -B argument: The value must be 0 or 1.");
                usage();
                exit(1);
            }
            break;
       case 'S' : // name of the file containg the scaffold(s) and the lengths to analyze
            arg++;
            strcpy( chr_name_all, argv[arg] );
            break;
        case 'C':
            arg++;
            *nitermc = (long int)atol(argv[arg]);
            if(*nitermc < 100) {
                printf("\n Error in -C argument: The value must be larger or equal than 100.");
                usage();
                exit(1);
            }
            break;
        case 'e' : /* Base Quality error range for classify SNPs*/
            arg++;
            *qerr_range = (int)atoi(argv[arg] );
            if(*qerr_range < 0) {
                printf("\n Error in -e argument: The value must be positive");
                usage();
                exit(1);
            }
            break;
        case 'g' : /* Base Quality error range for classify SNPs*/
            arg++;
            *outg = (int)atoi(argv[arg] );
            if(*outg != 0 && *outg != 1) {
                printf("\n Error in -g argument: The value must be 1 or 0.");
                usage();
                exit(1);
            }
            break;
        case 'd' : /* Testing for having reads at both strands*/
            arg++;
            *teststrands = (int)atoi(argv[arg] );
            if(*teststrands != 0 && *teststrands != 1) {
                printf("\n Error in -d argument: The value must be 1 or 0.");
                usage();
                exit(1);
            }
            break;
        case 'h' : /* h HEEEEEEEEEEELPPPPPPPPPPPPPP!!! */
            usage();
            exit(0);
            break;
      }
	  arg++;
    }
  }
  else {
    usage();
    exit(1);
  }
  if (*burnin_par>*iter || *iter<*samplemhmcmc){
    printf("\n Error: Bad configuration of the MCMC parameters. \n\n");
    usage();
    exit(1);
  }
  if(strcmp(chr_name_all,"") == 0) {
        printf("\nError: the name of the scaffold(s) (option -S) must be defined\n");
        usage();
        exit(1);
  }
  if(/*file_in[0]=='\0' || file_out[0]=='\0' || */*min_cov == 0 || *max_cov == 0 || *poolsize == 0 || *seed == 0 || /*chr_length_all[0] == 0 ||*/ chr_name_all[0] == 0 || *min_qual == 123456) {
    printf("\n\nError: Mandatory PScaller parameters not defined. \n\n");
    usage();
    exit(/*(file_in[0]=='\0') +
	 (file_out[0]=='\0') +*/
	 (*min_cov == 0) +
	 (*max_cov == 0) +
	 (*poolsize == 0) +
	 (*seed == 0) +
	 (*lentotal == 0) +
	 (*min_qual == 123456));
  }
  return 1;
}

int read_index_file(char *chr_name_all, unsigned long *nscaffolds,char ***chr_name_array,char ***chr_length_array) {
    
    FILE *file_scaffolds;
    char *buf;
    int c;
    int k;
    
    *nscaffolds = 1;
    chr_name_array[0] = (char **)calloc(*nscaffolds,sizeof(char *));
    chr_name_array[0][0] = (char *)calloc(MSP_MAX_NAME,sizeof(char));
    chr_length_array[0] = (char **)calloc(*nscaffolds,sizeof(char *));
    chr_length_array[0][0] = (char *)calloc(MSP_MAX_NAME,sizeof(char));
    
    if (!(file_scaffolds = fopen(chr_name_all,"r"))) {
        printf("Error reading the input file %s\n",chr_name_all);
        return(1);
    }
    if(!(buf = (char *)malloc(BUFSIZ))) {
        puts("\nError: Not enough memory to read  the input file.\n");
        return(1);
    }
    setbuf(file_scaffolds,buf);
    c=fgetc(file_scaffolds);
    while(c != EOF) {
        k=0;
        chr_name_array[0][*nscaffolds-1][k] = c; k++;
        while((c=fgetc(file_scaffolds))!= 9 && c!= 10 && c!=13 && c!=-1 && c!=0 && k<MSP_MAX_NAME-1) {
            chr_name_array[0][*nscaffolds-1][k] = c; k++;
        }
        chr_name_array[0][*nscaffolds-1][k] = '\0';
        if(c!= 9 && c!= 32) {
            printf("Error reading the input file %s:\n scaffold (%s) without length information.\n",chr_name_all, chr_name_array[0][*nscaffolds-1]);
            return(1);
        }
        do {
            c=fgetc(file_scaffolds);
        }while(!(c!= 9 && c!= 32 && c!= 10 && c!=13 && c!=-1 && c!=EOF));
        if(c==EOF) {
            printf("Error reading the input file %s:\n scaffold (%s) without length information.\n",chr_name_all, chr_name_array[0][*nscaffolds-1]);
            return(1);
        }
        k=0;
        chr_length_array[0][*nscaffolds-1][k] = c; k++;
        while((c=fgetc(file_scaffolds)) != 9 && c != 32 && c!= 10 && c!=13 && c!=-1 && c!=0 && k<MSP_MAX_NAME-1) {
            chr_length_array[0][*nscaffolds-1][k] = c; k++;
        }
        chr_length_array[0][*nscaffolds-1][k] = '\0';
        /*check next line, if exist*/
        if(c==32 || c==9) {
            do {
                c=fgetc(file_scaffolds);
            }while(c!=10 && c!= 13 && c!=EOF);
        }
        while(c==10 || c==13) {
            c=fgetc(file_scaffolds);
        }
        if(c==EOF)
            break;
        /*if exist, prepare new row in arrays*/
        *nscaffolds += 1;
        chr_name_array[0] = (char **)realloc(chr_name_array[0],*nscaffolds*sizeof(char *));
        chr_name_array[0][*nscaffolds-1] = (char *)calloc(MSP_MAX_NAME,sizeof(char));
        chr_length_array[0] = (char **)realloc(chr_length_array[0],*nscaffolds*sizeof(char *));
        chr_length_array[0][*nscaffolds-1] = (char *)calloc(MSP_MAX_NAME,sizeof(char));
    }
    fclose(file_scaffolds);
    free(buf);
    
    return(0);
}
