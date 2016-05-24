#'@title Do a permutation test for functional regression
#' @description Performs a permutation F test (Ramsay, Hooker, and Graves, 2009, p. 145)
#' for the significance of a functional covariate, and a permutation
#' likelihood ratio test.  The permutation test function currently doesn't allow models
#'  with multiple functional covariates, but 
#' subject-level covariates are allowed.
#' @param object An object of class funreg
#' @param num.permute The number of permutations to use.  Ramsay, Hooker and Graves (2009) 
#' recommended ``several hundred'' (p. 145), but for a quicker initial look
#' it might suffice to use 100.  
#' @param seed An optional random number seed.
#' @references
#' Onghena, P., & May, R. B. (1995). Pitfalls in computing and interpreting
#'          randomization test p values: A commentary on Chen and Dunlap.
#'          Behavior Research Methods, Instruments, & Computers, 27(3), 408-411.
#' 
#' Phipson, Belinda and Smyth, Gordon K. (2010) "Permutation P-values 
#'          Should Never Be Zero:  Calculating Exact P-values When Permutations
#'          Are Randomly Drawn," Statistical Applications
#'         in Genetics and Molecular Biology: Vol. 9: Iss. 1, Article 39.    
#' 
#' Ramsay, J. O., Hooker, G., & Graves, S. (2009). Functional
#'           data analysis with R and MATLAB. NY: Springer.
#' 
#' Sen, S. (2013) Permutation Tests. Available at
#'  \url{http://www.epibiostat.ucsf.edu/biostat/sen/statgen/permutation.html}
#' 
#' @return Returns a list with several components. First, 
#' \code{pvalue.F} is the p-value for the F test.  Second, \code{conf.int.for.pvalue.F} is
#' the confidence interval for estimating the p-value that would
#' be obtained from the dataset as \code{num.permute} approached infinity. 
#' The idea of a confidence interval for a p-value is explained further by 
#' Sen (2013), with a STATA example.  See "Permutation Tests" by Saunak Sen (2013) at
#' http://www.epibiostat.ucsf.edu/biostat/sen/statgen/permutation.html.
#' Third, \code{orig.F} is the F statistic calculated on the original
#' dataset. Last, \code{permuted.F} is the vector of F statistics calculated
#' on each of the random permuted datasets.  Also included are \code{pvalue.LR},
#' \code{conf.int.for.pvalue.LR},  \code{orig.LR}, \code{permuted.LR} for
#' the permutation test with a likelihood ratio statistic.
#' A more conservative alternative formula for the p-value is used in 
#' \code{pvalue.F.better} and \code{pvalue.LR.better}.  
#' It is not obvious whether to define the p-value as the 
#' proportion of permuted datasets with statistics less than or equal to
#' the original, or simply less than the original.  This should usually not 
#' matter, as a tie is not likely.  We made the arbitrary decision to use 
#' the former here because it was presented in this way in the 
#' Wikipedia article for permutation tests.                                     
#' The conservative alternative formula is 
#' the number of less extreme permuted datasets plus one,
#' over the total number of datasets plus one.  Adding one to the numerator
#' and denominator is suggested by some authors, partly  in order to prevent
#' a nonsensical zero p-value (Onghena & May, 1995; Phipson, Belinda & Smyth, 2010).
#'@export
funreg.permutation <- function(object,num.permute=500,seed=NULL) {
    stopifnot(class(object)=="funreg");
    p <- num.functional.covs.in.model(object);
    if (p>1) {
        stop(paste("funreg.permutation currently does not work for",
                   "multiple functional covariates."));
    }
    if (!is.null(seed)) {set.seed(seed);}
    result <- object;
    id <- result$data$id;
    y <- result$data$response;
    s <- result$data$other.covariates;
    x <- result$data$x;
    time <- result$data$time;
    short.id <- result$subject.info[,"id"];
    short.y <- result$subject.info[,"response"];
    short.s <- NULL;
    if (ncol(result$subject.info)>3) {
        short.s <- as.matrix(result$subject.info[,4:ncol(result$subject.info),drop=FALSE]);
        stopifnot(length(short.y)==nrow(short.s));
    }
    stopifnot(length(short.y)==length(short.id));
    stopifnot(identical(as.numeric(fitted(object)$id),as.numeric(short.id)));
    short.yhat <- fitted(object)$fitted;
    Fstat.orig <- var(short.yhat)/mean((short.y-short.yhat)^2);
    Fstat.permuted <- rep(NA,num.permute);
    loglikH1 <- logLik(result$intermediate.work.results$fit.model);
    loglikH1.permuted <- rep(NA,num.permute);
    for (this.permutation in 1:num.permute) {
        # Permute the covariates and responses once;
        shuffled.indices <- sample(1:length(short.y));
        shuffled.short.y <- short.y[shuffled.indices]; # do random permutation;
        shuffled.y <- NA*result$data$response;
        if (!is.null(short.s)) {
            shuffled.short.s <- short.s[shuffled.indices,,drop=FALSE];
            shuffled.s <- NA*s;
        } else {
            shuffled.short.s <- NULL;
            shuffled.s <- NULL;
        }
        for (i in short.id) {
            shuffled.y[which(id==i)] <- shuffled.short.y[which(short.id==i)];
            if (!is.null(short.s)) {
               for (j in 1:ncol(shuffled.short.s)) {
                    shuffled.s[which(id==i),j] <-
                             shuffled.short.s[which(short.id==i),j];
               }
            }
        }
        result.permuted <- redo.funreg(object=result,
                                 id=id,
                                 response=shuffled.y,
                                 time=time,
                                 other.covariates=shuffled.s,
                                 x=x);
        Fstat.permuted[this.permutation] <- var(fitted(result.permuted)$fitted)/
                     (mean((fitted(result.permuted)$response-fitted(result.permuted)$fitted)^2));
        loglikH1.permuted[this.permutation] <-
               logLik(result.permuted$intermediate.work.results$fit.model);
    }
    permutation.pvalue.F <- mean(Fstat.orig<=Fstat.permuted);
    permutation.pvalue.F.better <- (sum(Fstat.orig<=Fstat.permuted)+1)/num.permute;
    conf.int.for.pvalue.F <- binom.test(x=sum(Fstat.orig<=Fstat.permuted),
                                      n=length(Fstat.permuted))$conf.int;
    permutation.pvalue.LR <- mean(loglikH1<=loglikH1.permuted);                 
    permutation.pvalue.LR.better <- (sum(loglikH1<=loglikH1.permuted)+1)/num.permute;
        
    conf.int.for.pvalue.LR <- binom.test(x=sum(loglikH1<=loglikH1.permuted),
                                      n=length(loglikH1.permuted))$conf.int;
    cat(paste("Permutation F test two-sided p-value:  ",permutation.pvalue.F.better),"\n");
    cat(paste("Permutation LR test two-sided p-value:  ",permutation.pvalue.LR.better),"\n");
    return(list(pvalue.F=permutation.pvalue.F,
                conf.int.for.pvalue.F=conf.int.for.pvalue.F,
                pvalue.F.better=permutation.pvalue.F.better,
                orig.F=Fstat.orig,
                permuted.F=Fstat.permuted,      
                pvalue.LR=permutation.pvalue.LR,
                pvalue.LR.better=permutation.pvalue.LR.better,
                conf.int.for.pvalue.LR=conf.int.for.pvalue.LR,
                orig.loglikH1=loglikH1,
                permuted.loglikH1 = loglikH1.permuted));
}



