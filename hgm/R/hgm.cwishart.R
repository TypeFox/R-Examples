# $OpenXM: OpenXM/src/R/r-packages/hgm/R/hgm.cwishart.R,v 1.8 2014/04/18 05:04:40 takayama Exp $
"hgm.tk.pwishart" <-
function(m=3,n=5,beta=c(1,2,3),q0=0.2,approxdeg=-1,h=0.01,dp=20,q=10,
         mode=c(1,1,0),method="a-rk4",err=c(-1.0,-1.0),
         automatic=1,assigned_series_error=0.00001,verbose=0) { 
  x<-q; x0<-q0;
  nstrategy<-0;
  mm<-charmatch(method,c("rk4","a-rk4","a-rk4-m"));
  if (!is.na(mm)) nstrategy<- (mm-1);

  if ((m>=200) || (m<=0)) stop("Invalid size of m."); #200 is M_m_MAX in jack-n.c

  if (!is.loaded("hgm")) library.dynam("hgm",package="hgm",lib.loc=NULL);

  .C("Rmh_set_strategy",as.integer(nstrategy),as.double(err),retv=double(1),
     package="hgm")$retv ;

  if (m<1) m=1;
  rank <- 2^m;
  rsize <- rank+1;
  if (mode[3] > 0) rsize <- rsize+mode[3]; 
  if (approxdeg < 0) approxdeg <- 6;
##argchecks
  if (class(beta) != "numeric") stop("beta must be a vector.");
  if (length(beta) != m) stop("The length beta must be equal to m.");
  for (i in 1:m) {
    if (beta[i] <= 0) stop("beta[i] must be positive.");
  }
  if (m != 1) {
   for (i in 1:(m-1)) {
    for (j in  (i+1):m) {
       if (beta[i] == beta[j]) stop("beta's must be different.");
    }
   } 
  }
  if (class(mode) != "numeric") stop("mode must be a vector of length 3.");
  if (length(mode) != 3) stop("mode must be a vector of length 3.");
##end of argchecks

  .C("Rmh_cwishart_gen",as.integer(m),as.integer(n),as.double(beta),as.double(x0),
     as.integer(approxdeg),
     as.double(h),
     as.integer(dp),as.double(x),
     as.integer(mode),as.integer(rank),
     as.integer(automatic),as.double(assigned_series_error),as.integer(verbose),
     retv=double(rsize),PACKAGE="hgm")$retv
}
