#' @importFrom stats ks.test
#' @importFrom stats chisq.test
#' @importFrom stats dpois
#' @importFrom stats pchisq
#' @importFrom stats rnorm
#' @importFrom stats runif
#' @importFrom parallel parLapply
#' @importFrom parallel clusterExport
#' @importFrom parallel setDefaultCluster
#' @importFrom parallel stopCluster
#' @importFrom parallel makePSOCKcluster
#' @importFrom Rmpfr mpfrArray
#' @importFrom Rmpfr mpfr
#' @importFrom gmp Stirling2
#' @importFrom LambertW W
#' @importFrom MissMech AndersonDarling
#' @importFrom kSamples ad.test
#' @importFrom tseries jarque.bera.test
#' @importFrom sfsmisc digitsBase
#' @importFrom methods new
#' @importFrom utils write.table

#' @S3method adaptive.chi.square default
#' @S3method birthday.spacings default
#' @S3method book.stack default
#' @S3method GCD.test default
#' @S3method random.walk.tests default
#' @S3method topological.binary default
#' @S3method adaptive.chi.square main
#' @S3method birthday.spacings main
#' @S3method book.stack main
#' @S3method GCD.test main
#' @S3method random.walk.tests main
#' @S3method topological.binary main
#' @S3method print CryptRndTest

#' @export adaptive.chi.square
#' @export birthday.spacings
#' @export book.stack
#' @export GCD.test
#' @export random.walk.tests
#' @export topological.binary
#' @export GCD
#' @export GCD.big
#' @export GCD.q
#' @export Strlng2
#' @export TBT.criticalValue
#' @export toBaseTen
#' @export toBaseTwo

KSADdga=function(e,alfa,n,m,lambda,num.class=10){
	z=0
	expected=0
  p=0
  p[1:(num.class+1)]=dpois(0:num.class,lambda=lambda) 
  p[num.class+1]=p[num.class+1]+sum(dpois((num.class+1):1000,lambda=lambda)) 
	N=length(e)
	expected=round(p*N)
	if (sum(expected)!=N){
		expected[which.max(expected)]=expected[which.max(expected)]-(sum(expected)-N)
	}
	z=rep(0:(length(expected)-1),expected) 
	test=kSamples::ad.test(e,z,method="simulated",dist=FALSE,Nsim=1000)
	ADtest=test$ad
	if (min(ADtest[,3])<alfa){ 
		sonucAD=0
	} else{ 
		sonucAD=1
	}
  
	KStest=0
	test2=ks.test(e,z)
	KStest[1]=test2$statistic
	KStest[2]=test2$p.value
	if (KStest[2]<alfa){ 
		sonucKS=0
	} else{ 
		sonucKS=1
	}
  
  k=length(expected)
  observed=array(0,k)
  observed[1:length(table(e))]=table(e)
  KKtest=0 
  observed[which(observed==0)]=10^-5
  expected[which(expected==0)]=10^-5
  KKtest[1]=sum(((observed-expected)^2)/expected,na.rm=TRUE) #gives test statistic
  KKtest[2]=pchisq(KKtest[1],(k-1),lower.tail=FALSE) #gives p.value of the chi-square test
	  if (KKtest[2]<alfa){ 
  		sonucKK=0
	  } else{ 
		  sonucKK=1
	  }
 	result=list(sonucAD=sonucAD,ADtest=ADtest,sonucKS=sonucKS,KStest=KStest,sonucKK=sonucKK,KKtest=KKtest)
	return(result)
}