
pvals.cat <- function(pest,pval.test) 
{

# check whether the number of columns in 'pest' is correct (it must be a multiple of 3 plus 1):
if (((ncol(pest)-1)%%3)>0) {
  warning("Invalid number of columns in 'pest'! The vector with p-values was set to NULL.")
  pval <- NULL
} else {
  pval <- array(NA,nrow(pest)) # allocation of a vector with p-values
  st <- array(NA,(ncol(pest)-1)/3)
  cit <- array(NA,(ncol(pest)-1)/3)
  ivalid <- c(1:nrow(pest))[(rowSums(is.na(pest))==0)&(rowSums(pest[,seq(3,ncol(pest),3)]==0)==0)&(rowSums(pest[,seq(3,ncol(pest),3)]==1)==0)] # creates a vector of indexes of time points in which all stratified curves have defined values and survival estimates are > 0 and < 1
  for (i in ivalid) {
     st <- pest[i,seq(3,ncol(pest),3)] # survival function point estimates
     if (sum(pest[i,]<0.1)>0) {
        cit <- pest[i,seq(4,ncol(pest),3)] # upper bound of the confidence interval (upper bound is taken because lower bound could be biased if it is too small)
     } else {
        cit <- pest[i,seq(2,ncol(pest),3)] # lower bound of the confidence interval (lower bound is taken because curves are sometimes >1 which means that the estimated sigma would not be correct)
     }
     Vt <- ((st-cit)/1.96)^2
     ot.sq <- Vt / (st^2) # estimated variance

     # computation of the chi-square test statistic:
     chisq <- switch(pval.test, 
                     # naive test:
                     naive = chisq.naive(st,ot.sq), 
                     # test based on log transformation of the survival estimates:
                     log = chisq.log(st,ot.sq), 
                     # test based on log-log transformation of the survival estimates:
                     loglog = chisq.loglog(st,ot.sq) ) 

     # computation of p-value:
     if (!is.null(chisq)) {
        pval[i] <- 1-pchisq(chisq,df=1)
     }
  }
}

pval <- cbind(pest[,1],pval)
colnames(pval) <- c("Month","p-value")

pvals.cat <- pval

}




