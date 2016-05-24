
pvals.2cat <- function(pest,pval.test) 
{

pval <- array(NA,nrow(pest)) # allocation of a vector with p-values
ivalid <- c(1:nrow(pest))[(rowSums(is.na(pest))==0)&(rowSums(pest[,seq(3,ncol(pest),3)]==0)==0)&(rowSums(pest[,seq(3,ncol(pest),3)]==1)==0)] # creates a vector of indexes of time points in which all stratified curves have defined values and survival estimates are > 0 and < 1

# computation of p-values:
for (i in ivalid) {
   s1t <- as.vector(pest[i,3]) # survival estimate
   s2t <- as.vector(pest[i,6]) # survival estimate
   if (sum(pest[i,]<0.1)>0) {
      ci1t <- as.vector(pest[i,4]) # upper bound of the confidence interval (upper bound is taken because lower bound could be biased if it is too small))
      ci2t <- as.vector(pest[i,7]) # upper bound of the confidence interval (upper bound is taken because lower bound could be biased if it is too small)
   } else {
      ci1t <- as.vector(pest[i,2]) # lower bound of the confidence interval (lower bound is taken because curves are sometimes >1 which means that the estimated sigma would not be correct)
      ci2t <- as.vector(pest[i,5]) # lower bound of the confidence interval (lower bound is taken because curves are sometimes >1 which means that the estimated sigma would not be correct)
   }
   V1t <- ((s1t-ci1t)/1.96)^2 
   V2t <- ((s2t-ci2t)/1.96)^2 
   o1t.sq <- V1t / (s1t^2) # estimated variance
   o2t.sq <- V2t / (s2t^2) # estimated variance

   # computation of the chi-square test statistic:
   chisq <- switch(pval.test, 
                   # naive test:
                   naive = ((s1t-s2t)^2) / ((s1t^2)*o1t.sq + (s2t^2)*o2t.sq), 
                   # test based on log transformation of the survival estimates:
                   log = ((log(s1t)-log(s2t))^2)/(o1t.sq+o2t.sq), 
                   # test based on log-log transformation of the survival estimates:
                   loglog = ((log(-log(s1t)) - log(-log(s2t)))^2) / ((o1t.sq/((log(s1t))^2))+(o2t.sq/((log(s2t))^2))) ) 

   # computation of p-value:
   if (!is.null(chisq)) {
     pval[i] <- 1-pchisq(chisq,df=1)
   }
}

pval <- cbind(pest[,1],pval)
colnames(pval) <- c("Month","p-value")

pvals.2cat <- pval

}


