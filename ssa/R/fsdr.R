#' False simultaneous discovery rate control
#'
#' Given two sequences of paired test statistics, returns the optimal rectangular rejection that identifies the largest number of simultaneous signals while controlling the false discovery rate.
#'
#' @param T1,T2 paired vectors of test statistics, both must be the same length; can be p-values or otherwise; if not p-values, must be stochastically larger under the null; must contain only positive values
#' @param alpha desired false simultaneous discovery rate
#' @param m1,m2 search only up the m1th (m2th) most significant test statistic in T1 (T2); NULL to search through all statistics
#' @param p1,p2 TRUE if T1 (T2) is a vector of p-values
#' @param jitter NULL if no jittering is desired to resolve ties, otherwise a jitter of \code{runif(0,jitter)} will be added to all entries of T1 and T2
#'
#' @return two-component vector; the first component is the optimal threshold for T1 and the second is for T2
#'
#' @examples
#' \donttest{
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
#' ## convert to p-value
#' P1 <- 2*pnorm(-abs(Z1));
#' P2 <- 2*pnorm(-abs(Z2));
#' ## run different version of fsdr()
#' out.pp <- fsdr(P1,P2,alpha=0.05);
#' out.zp <- fsdr(abs(Z1),P2,p1=FALSE,alpha=0.05);
#' out.pz <- fsdr(P1,abs(Z2),p2=FALSE,alpha=0.05);
#' out.zz <- fsdr(abs(Z1),abs(Z2),p1=FALSE,p2=FALSE,alpha=0.05);
#' ## discovered simultaneous features
#' R1 <- which(P1<=out.pp[1]&P2<=out.pp[2]);
#' R2 <- which(abs(Z1)>=out.zp[1]&P2<=out.zp[2]);
#' R3 <- which(P1<=out.pz[1]&abs(Z2)>=out.pz[2]);
#' R4 <- which(abs(Z1)>=out.zz[1]&abs(Z2)>=out.zz[2]);
#' }
#'
#' @import stats
#' @useDynLib ssa
#' @export

fsdr <- function(T1,T2,alpha,m1=10000,m2=10000,p1=TRUE,p2=TRUE,jitter=NULL){ 
  if(length(T1)!=length(T2)){
    stop("Test statistic vectors must be of same length");
  }
  if(sum(c(T1,T2)<0)>0){
    stop("Test statistics must be positive");
  }
  if(sum(is.na(T1)||is.na(T2))>0){
    stop("No missing data allowed");
  }
  if(alpha>1){
    stop("alpha must be less than 1");
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
  
  if(p1&&p2){
    ord <- sort.list(T1,method="quick",na.last=NA,decreasing=FALSE);
    TT1 <- as.numeric(T1[ord]);
    TT2 <- as.numeric(T2[ord]);
    ## subtract 1 bc C indices start at 0
    T2ord <- as.integer(sort.list(TT2,method="quick",na.last=NA,decreasing=FALSE)-1);
    
    t <- .Call("fsdr_p",TT1,TT2,T2ord,m1,m2,as.numeric(alpha));
  } else {
    ## transform if the T1 and T2 are p-values
    if(p1){ T1 <- -log10(T1); }
    if(p2){ T2 <- -log10(T2); }
    
    ord <- sort.list(T1,method="quick",na.last=NA,decreasing=TRUE);
    TT1 <- as.numeric(T1[ord]);
    TT2 <- as.numeric(T2[ord]);
    ## subtract 1 bc C indices start at 0
    T2ord <- as.integer(sort.list(TT2,method="quick",na.last=NA,decreasing=TRUE)-1);
    
    t <- .Call("fsdr",TT1,TT2,T2ord,as.numeric(max(T1)*max(T2)),
               m1,m2,as.numeric(alpha));
    
    if(p1){
      if(t[1]==-99){ t[1] <- 0; } else { t[1] <- 10^-t[1]; }
    } else {
      if(t[1]==-99){ t[1] <- max(T1)+1; }
    }
    if(p2){
      if(t[2]==-99){ t[2] <- 0; } else { t[2] <- 10^-t[2]; }
    } else {
      if(t[2]==-99){ t[2] <- max(T2)+1; }
    }
  }
  return(t);
}

#' False simultaneous discovery rate control -- report all thresholds
#'
#' Given two sequences of paired test statistics, returns all rectangular rejections that identify the largest number of simultaneous signals while also controlling the false discovery rate.
#'
#' @param T1,T2 paired vectors of test statistics, both must be the same length; must be stochastically larger under the alternative than under the null; must contain only positive values
#' @param alpha desired false simultaneous discovery rate
#' @param m1,m2 search only up the m1th (m2th) most significant test statistic in T1 (T2)
#'
#' @return k x 2 matrix, where k is the number of rectangular regions found; the first column is the threshold for T1 and the second column is for T2
#'
#' @examples
#' \donttest{
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
#' ## all rectangular rejection regions
#' out.zz <- fsdr_all(abs(Z1),abs(Z2),alpha=0.05,m1=1000,m2=1000);
#' ## 10 sets of identified simultaneous signals
#' R <- apply(out.zz[1:10,],1,function(x){
#'              which(abs(Z1)>=x[1]&abs(Z2)>=x[2]);
#'            });
#' }
#'
#' @import stats
#' @useDynLib ssa
#' @export

fsdr_all<- function(T1,T2,alpha,m1=5000,m2=5000){
  if(length(T1)!=length(T2)){
    stop("Test statistic vectors must be of same length");
  }
  if(sum(c(T1,T2)<0)>0){
    stop("Test statistics must be positive");
  }
  if(sum(is.na(T1)||is.na(T2))>0){
    stop("No missing data allowed");
  }
  if(alpha>1){
    stop("alpha must be less than 1");
  }
  if(sum(duplicated(T1))>0||sum(duplicated(T2))>0){
    warning("Some test statistics are tied");
  }
  
  m <- length(T1);
  if(is.null(m1)){ m1 <- m; } else { m1 <- min(m,m1); }
  if(is.null(m2)){ m2 <- m; } else { m2 <- min(m,m2); }
  
  ord <- sort.list(T1,method="quick",na.last=NA,decreasing=TRUE);
  U <- T1[ord];
  V <- T2[ord];
  Vord <- sort.list(V,method="quick",na.last=NA,decreasing=TRUE)-1; ## subtract 1 bc C indices start at 0
  
  t <- .Call("fsdr_all",as.numeric(U),as.integer(Vord),
             as.integer(m1),as.integer(m2),as.numeric(alpha));
  
  if(t[1]==-99){
    return(NA);
  } else {
    nt <- length(t)/2; ## number of thresholds that attain max fdr
    ret <- cbind(t[(1:nt)*2-1],V[t[(1:nt)*2]+1]);
    return(ret);
  }
}
