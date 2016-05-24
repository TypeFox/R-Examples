#' Max test for detecting simultaneous signals
#'
#' Given two sequences of paired test statistics, tests whether any simultaneous signals exist
#'
#' @param T1,T2 paired vectors of test statistics, both must be the same length; must be stochastically larger under the alternative than under the null; must contain only positive values
#'
#' @return
#' \item{p}{p-value}
#' \item{M}{value of max statistic}
#'
#' @examples
#' ## generate paired test statistics
#' p <- 10^6; ## total number of pairs
#' X <- c(rep(0,p-30),rep(1,10),rep(2,10),rep(3,10));
#' ## X=0: no signal in either sequence of tests
#' ## X=1: signal in sequence 1 only
#' ## X=2: signal in sequence 2 only
#' ## X=3: simultaneous signal
#' set.seed(1);
#' Z1 <- rnorm(p,0,1); Z1[X==1|X==3] <- rnorm(20,3,1);
#' Z2 <- rnorm(p,0,1); Z2[X==2|X==3] <- rnorm(20,4,1);
#' maxtest(abs(Z1),abs(Z2));
#'
#' @import stats
#' @export

maxtest <- function(T1,T2){
  if(length(T1)!=length(T2)){
    stop("Test statistic vectors must be of same length");
  }
  if(sum(c(T1,T2)<0)>0){
    stop("Test statistics must be positive");
  }
  if(sum(is.na(T1)||is.na(T2))>0){
    stop("No missing data allowed");
  }
  
  T <- pmin(T1,T2);
  M <- max(T);
  n <- length(T);
  k <- sum(T1>=M);
  m <- sum(T2>=M);
  p <- 1-phyper(0,k,n-k,m);
  
  return(list(p=p,M=M));
}
