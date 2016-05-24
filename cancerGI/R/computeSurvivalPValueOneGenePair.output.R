computeSurvivalPValueOneGenePair.output <- function (file.out, genes.info, data.mut, data.surv, colTime=2, colStatus=3, groups=c("All", "Two"), PRINT=FALSE) {
    #############################################################################################
    # For a pair of genes, estimate the survival curves for each group and 
    # compute the p value for the log rank test of same survival prospect for the groups.  There 
    # may be two or four groups.
    #
    # Input:
    # data.mut: integer matrix of cases by 2 genes; 0 as wildtype, 1 amplification and -1 deletion.
    # data.surv: matrix of data frame containing case ID, survival time and survival status.  
    #            Cases should be ordered in the same way as data.mut.
    # Output:
    # vector with elements: nSingleMut (no. of cases with single mutation), nDoubleMut (no. of cases with 
    # double mutation), obsSingleMut (no. of deceased cases with single mutation), obsDoubleMut 
    # (no. of deceased cases with double mutation), expSingleMut (expected no. of deceased cases
    # with single mutation), expDbouleMut (expected no. of deceased cases with double mutation), 
    # medianSingleMut (estimated median survival time for single mutation), medianDoubleMut (estimated 
    # median survival time for double mutation)
    #############################################################################################
    
    groups <- match.arg (groups)
    
    data.mut <- data.mut[,as.numeric(genes.info[5:6])]
    
    # remove NAs in survival data
    which.na <- which (is.na (data.surv[,colTime]) | is.na (data.surv[,colStatus]))
    if (length (which.na)) {
        data.mut <- data.mut[-which.na,]
        data.surv <- data.surv[-which.na,]        
    }
    
    if (PRINT) {
        cat ("computing p values...\n")
    }
    # convert to binary vectors based on type of mutation
    x <- as.numeric (data.mut[,1]==as.numeric(genes.info[2]))
    y <- as.numeric (data.mut[,2]==as.numeric(genes.info[4]))
    #print (table (x,y))
    if (groups=="Two") {
        z <- x+y
        # consider cases with at least one mutation
        z.inc <- which (z>0)   
        results <- rep (NA, 9)
    }
    if (groups=="All") {
        z <- 2*y+x
        z.inc <- 1:length (z)
        results <- rep (NA, 17)
    }
    if (PRINT) {
        cat ("summary of z\n")
        print (table (z[z.inc]))
    }
    if (length (unique (z[z.inc]))>1) {
        # estimate survival times
        z.surv.sum <- summary (survfit (Surv (data.surv[z.inc,colTime], data.surv[z.inc,colStatus]) ~ z[z.inc]))
        # perform logrank test for survival times between groups
        z.surv <- survdiff (Surv (data.surv[z.inc,colTime], data.surv[z.inc,colStatus]) ~ z[z.inc])
        if (PRINT) {
            print (z.surv.sum)
            print (z.surv)
            #print (z.surv$obs)
            #print (z.surv$exp)                        
        }
        # compute chi-square p value for log rank test, and extract total counts, observed and expected counts, and estimated median (or minimum
        # survival time
        results <- c(as.vector (z.surv$n), z.surv$obs, z.surv$exp, findSurvMedian (z.surv.sum), pchisq (z.surv$chisq, df=nlevels (as.factor(z[z.inc]))-1, lower.tail=FALSE))
        write.table (genes.info[1:4], file.out, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE, append=TRUE, eol="\t")
        write (results, file.out, ncolumns=length (results), append=TRUE, sep="\t")
    }
    #else {
    #    print (genes.info[1:4])
    #}
}

