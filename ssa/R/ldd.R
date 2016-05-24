#' Latent dependency detection
#'
#' Given two sequences of paired test statistics, tests whether the latent indicators of significance are positively dependent.
#'
#' @param T1,T2 paired vectors of test statistics, both must be the same length; can be p-values or otherwise; if not p-values, must be stochastically larger under the null
#' @param m1,m2 search only up the m1th (m2th) most significant test statistic in T1 (T2); NULL to search through all statistics
#' @param p1,p2 TRUE if T1 (T2) is a vector of p-values
#' @param perm the indices of T1 will be randomly permuted \code{perm} times
#' @param jitter NULL if no jittering is desired to resolve ties, otherwise a jitter of \code{runif(0,jitter)} will be added to all entries of T1 and T2
#'
#' @return
#' \item{D}{value of the test statistic}
#' \item{p.perm}{permutation p-value}
#' \item{p.asymp}{asymptotic approximate p-value}
#'
#' @examples
#' ## generate paired test statistics
#' \donttest{
#' p <- 10^6; ## total number of pairs
#' X <- c(rep(0,p-30),rep(1,10),rep(2,10),rep(3,10));
#' ## X=0: no signal in either sequence of tests
#' ## X=1: signal in sequence 1 only
#' ## X=2: signal in sequence 2 only
#' ## X=3: simultaneous signal
#' set.seed(1);
#' Z1 <- rnorm(p,0,1); Z1[X==1|X==3] <- rnorm(20,3,1);
#' Z2 <- rnorm(p,0,1); Z2[X==2|X==3] <- rnorm(20,4,1);
#' ## convert to p-value
#' P1 <- 2*pnorm(-abs(Z1));
#' P2 <- 2*pnorm(-abs(Z2));
#' ## run different version of ldd()
#' out.pp <- ldd(P1,P2,perm=100);
#' out.zp <- ldd(abs(Z1),P2,p1=FALSE,perm=100);
#' out.pz <- ldd(P1,abs(Z2),p2=FALSE,perm=100);
#' out.zz <- ldd(abs(Z1),abs(Z2),p1=FALSE,p2=FALSE,perm=100);
#' }
#'
#' @import stats
#' @useDynLib ssa
#' @export

ldd <- function(T1,T2,m1=1000,m2=1000,perm=0,p1=TRUE,p2=TRUE,jitter=NULL)
{
  if(length(T1)!=length(T2)){
    stop("Test statistic vectors must be of same length");
  }
  if(sum(is.na(T1)||is.na(T2))>0){
    stop("No missing data allowed");
  }
  if(sum(c(T1,T2)<0)>0){
    warning("Test statistics contain negative values");
  }
  
  m <- length(T1);
  if(is.null(m1)){ m1 <- m; } else { m1 <- min(m,m1); }
  if(is.null(m2)){ m2 <- m; } else { m2 <- min(m,m2); }
  m1 <- as.integer(m1); m2 <- as.integer(m2);

  ## deal with ties
  if(sum(duplicated(T1))>0||sum(duplicated(T2))>0){
    if(!is.null(jitter)){
      if(p1){
        T1 <- 10^-(-log10(T1)+runif(m,0,jitter));
      } else {
        T1 <- T1+runif(m,0,jitter);
      }
      if(p2){
        T2 <- 10^-(-log10(T2)+runif(m,0,jitter));
      } else {
        T2 <- T2+runif(m,0,jitter);
      }
      ## if still tied after jittering
      if(sum(duplicated(T1))>0||sum(duplicated(T2))>0){
        warning("Need to increase jitter");
      }
    } else {
      warning("Some test statistics are tied");
    }
  }
  
  ## the C code is written in terms of empirical CDFs instead of
  ## empirical survival functions
  ## if not a p-value, take negatives to make non-null statistics
  ## stochastically smaller
  if(!p1){
    T1 <- -as.numeric(T1);
  }
  if(!p2){
    T2 <- -as.numeric(T2);
  }
  
  ## calculate test statistics statistics
  Ds <- rep(NA,(perm+1));
  for(i in 1:(perm+1)){
    if(i==1){
      UU <- T1; ## don't permute; this is the real statistic
    } else {
      UU <- T1[sample(1:m,m,replace=FALSE)]; ## permute the T1
    }
    ord <- sort.list(UU,method="quick",na.last=NA);
    UU <- UU[ord];
    VV <- T2[ord];
    Vord <- sort.list(VV,method="quick",na.last=NA)-1; ## subtract 1 bc C indices start at 0
    ret <- .Call("ldd",as.numeric(UU),as.integer(Vord),m1,m2);
    Ds[i] <- ret[1];
  }
  Ds <- Ds*sqrt(m);
  
  ## permutation p-value
  if(perm==0){
    p.perm <- NA;
  } else {
    p.perm <- mean(Ds[-1]>Ds[1]);
  }
  
  ## closed-form asymptotic p-value
  p.asymp <- 1-exp(-(Ds[1]/sqrt(log(m)))^(-2));
  
  return(list(D=Ds[1],p.perm=p.perm,p.asymp=p.asymp));
}
