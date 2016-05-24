"cond.test"<-function(data, groups = c(1:nrow(data)), kappa = NULL, biasCorrect = TRUE){
# DATE WRITTEN: 4 Ago 2010          LAST REVISED:  26 Sep 2012
# AUTHOR: Sandra Barragan based on the SAS routines written by Miguel A. Fernandez.
# DESCRIPTION: 
# data are the unordered estimator
# replic=T means that data have replication => kappa estimated
# kappa is the value of kappa when is known
# REFERENCE: Fernandez, M. A., Rueda, C. and Shyamal, P. (2011). 
#             Isotropic order among core set of orthologs conserved between 
#             budding and fission yeasts. Preprint.
# SEE ALSO: mrl, CIRE.

if(is.vector(data)){data <- cbind(data)}
if(ncol(data) == 1 && missing(kappa)){stop("If data do not have replications, kappa must be known")}

q <- nrow(data)
n <- as.vector(table(groups))
posi <- prod(factorial(n))
resultCIRE <- CIRE(data, groups, circular=TRUE)
SCE <- resultCIRE$SCE
ordermeans <- unlist(resultCIRE$CIRE)
data1 <- sort(ordermeans)
data2 <- c(data1[2:q],data1[1])
difference <- ifelse((data1 - data2) == 0, 0, 1)
m <- sum(difference)

if(!is.null(kappa)& SCE!=0){
  kap <- kappa
  T <- 2*kap*SCE
  pvalue <- pchisq(T, q-m, lower.tail=FALSE)*(1 - ((posi)/(factorial(q - 1))))
  } # end kappa known

if(is.null(kappa)){
  Rs<-apply(data,1,mrl)
  kap<-A1inv(mean(Rs))
  if(biasCorrect == TRUE) {
    N<-ncol(data)*nrow(data)
    if(kap < 2){
      kap<-max(kap-(2/(N*kap)), 0)
    }else{
      kap<-(((N-1)^3)*kap)/((N^3)+N)
    }
  }
if(SCE!=0){
  T <- 2*kap*SCE
  pvalue <- pf(T, q-m, q-1, lower.tail=FALSE)*(1 - ((posi)/(factorial(q - 1))))
}
} # end kappa unknown and replications
if(SCE==0){
  pvalue<-1
  warning("SCE=0")
}
lexit <- resultCIRE
class(lexit)<-"isocir"
lexit$pvalue <- pvalue
lexit$kappa <- kap
attr(lexit,"estkappa")<-ifelse(is.null(kappa),"Kappa has been estimated","Kappa was known")

return(lexit)
} # end cond.test
