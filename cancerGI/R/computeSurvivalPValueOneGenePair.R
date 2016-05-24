computeSurvivalPValueOneGenePair <- function (data.mut, data.surv, colTime=2, colStatus=3, type.gene1=(-1), type.gene2=(-1), groups=c("All", "Two"), compare=c("Both", "Gene1", "Gene2"), PLOT=FALSE, PRINT=FALSE, pvalue.text.x=10, pvalue.text.y=0.1, legend.x=150, legend.y=1.0) {
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
    if (groups=="Two") {
        compare <- match.arg (compare)        
    }
    
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
    x <- as.numeric (data.mut[,1]==type.gene1)
    y <- as.numeric (data.mut[,2]==type.gene2)
    #print (table (x,y))
    if (groups=="All") {
        z <- 2*y+x
        z.inc <- 1:length (z)
        #results <- rep (NA, 18)
    } else {
        z <- x+y
        results <- rep (NA, 10)
        if (compare=="Gene1"){
            z.x <- x
            z.x[which(z==2)] <- 2
            z <- z.x
        } 
        if (compare=="Gene2"){
            z.y <- y
            z.y[which(z==2)] <- 2
            z <- z.y
        }
        # consider cases with at least one mutation
        z.inc <- which (z>0)   
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
        results <- c(nlevels (as.factor(z[z.inc])), as.vector (z.surv$n), z.surv$obs, z.surv$exp, findSurvMedian (z.surv.sum), pchisq (z.surv$chisq, df=nlevels (as.factor(z[z.inc]))-1, lower.tail=FALSE))
        if (PLOT) {
            if (groups=="All") {
                plot (survfit (Surv (data.surv[z.inc,colTime], data.surv[z.inc,colStatus]) ~ z[z.inc]), xlim=c(0,250), lty=1, col=c("red", "blue", "magenta", "green"), xlab="Estimated survival time (months)", ylab="Survival probability", main=paste (colnames(data.mut)[1], "-", colnames(data.mut)[2], sep=""))
                legend (x=legend.x, y=legend.y, legend = c("WT", paste (colnames(data.mut)[1], ifelse (type.gene1==1, "_gain", "_loss"), sep=''), paste (colnames(data.mut)[2], ifelse (type.gene2==1, "_gain", "_loss"), sep=''), "Double mutation"), bty="n", lty=1, col=c("red", "blue", "magenta", "green"))
                text (pvalue.text.x, pvalue.text.y, paste ("p=", sprintf ("%1.1e", pchisq (z.surv$chisq, df=nlevels (as.factor(z[z.inc]))-1, lower.tail=FALSE)), sep=""))                
            }
            else {
                plot (survfit (Surv (data.surv[z.inc,colTime], data.surv[z.inc,colStatus]) ~ z[z.inc]), xlim=c(0,250), lty=1, col=c("magenta", "green"), xlab="Estimated survival time (months)", ylab="Survival probability", main=paste (colnames(data.mut)[1], "-", colnames(data.mut)[2], sep=""))
                legend (x=legend.x, y=legend.y, legend = c("Single mutation", "Double mutation"), bty="n", lty=1, col=c("magenta", "green"))
                text (pvalue.text.x, pvalue.text.y, paste ("p=", sprintf ("%1.1e", pchisq (z.surv$chisq, df=nlevels (as.factor(z[z.inc]))-1, lower.tail=FALSE)), sep=""))                
            }
        }
    }
    
    return (results)
}

