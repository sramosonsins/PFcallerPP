#ifndef COMMON_H_
#define COMMON_H_

#ifdef	__cplusplus
extern "C" {
#endif
    
#define _CRT_SECURE_NO_DEPRECATE
    
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>

/* For compatibility for some old compilers */
#ifndef NULL
#define NULL	0
#endif
    
#define PSCALLER "PFcaller: A Frequency Caller for Polyploid and Pooled Sequences. version b20220518.\n"

#ifndef noMacOS
    #define noMacOS	0
#endif

#if noMacOS==1
    #include <malloc.h>
#endif
    
#define LIMIT_COMB 100
#define LIMIT_FACT 30000
#define NITER_MC 200
#define MSP_MAX_NAME 21
#define LEN_FILENAME 1024
#define BLOCKBASES 1000000
#define BLOCKSECTS 100000
#define NPOPINTV 20
    
#define THR_ERR 100
    
    /*declare struct*/
    struct typepos {
        char *cchrom/*[MSP_MAX_NAME]*/;
        unsigned long pos;
        unsigned long t;
        char nucl[5];
    };
    
    struct Obs {
        unsigned long freq[4];
        float error_mean[4];
        unsigned long cov;
        unsigned char outg;
    };
    
    struct Est {
        float nr;
        float nrFreq[4];
        float nerr[4];
        float freqpool[4];
        float freqpop;
    };
    
    struct DistEst {
        unsigned long *nr;
        unsigned long *nrFreq[2];
        unsigned long *nerr[4];
        unsigned long *freqpool[2];
        unsigned long *freqpop;
        double *lnprob;
    };
    
    struct featpos {
        struct Obs observed;
        struct Est estimated;
        struct DistEst dist_estimated;
    };
    
    struct fasta {
        char *name;
        char *seq;
    };
    
    struct MCMCsample {
        unsigned long nc;
        unsigned long p1;
        unsigned long p2;
        unsigned long p3;/*not used*/
        unsigned long p4;/*not used*/
        unsigned long frp1; /*times pmc*/
        unsigned long frp2;
        unsigned long popfr;/*not evaluated*/
        unsigned long ner1;
        unsigned long ner2;
        unsigned long ner3;
        unsigned long ner4;
        double lnprob;
    };
    
    struct list_chr_pileup {
        char *cchrom/*[MSP_MAX_NAME]*/;
        unsigned long init;
        unsigned long end;
    };
    
#ifdef	__cplusplus
}
#endif

#endif /* COMMON_H_ */
