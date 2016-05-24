errorTesting <- function(N,
                         nevents,
                         tmedian,
                         meantf,
                         p,
                         t,
                         SigmaTheta,
                         sigmae_2,
                         ordtraj,
                         beta, 
                         alpha,
                         tol){

  # Error checking for input parameters
  # ****************************************
  N <- as.integer(N)
  if( N < 20L ) stop("N - total sample size must be >= 20.")

  nevents <- as.integer(nevents)
  if( (nevents < 20L) || (nevents > N) ){
    stop("nevents - number of events must be [20, N].")
  }

  if( tmedian < tol ) stop("tmedian - median survival time must be positive.")

  if( meantf < tol || meantf > max(t) + tol ) {
    stop("meantf - mean follow-up time must be (0, max(t)].")
  }

  if( !is.numeric(p) ) {
    stop("p - subject proportions per measurements must be a numeric vector.")
  }

  p <- sort(p)
  n <- length(p)

  if( p[1] < -tol || p[n] > 1.0-tol) {
    stop("Subject proportions (p) must be [0.0, 1.0); 0.0 proportion allowed.")
  }

  if( !isTRUE(all.equal(sum(p),1.0,tolerance=tol)) ) {
    stop("Subject proportions (p) must sum to 1.")
  }

  if( !is.numeric(t) ) {
    stop("t - measurement times must be entered as a numeric vector.")
  }

  t <- as.array(t)
  if( length(t) <= 2) stop("length(t)>2, - at least three measurements.")

  if( n != length(t) ) stop("Lengths of p and t must agree.")

  tst <- apply(X=array(1:n),
               MARGIN=1,
               FUN=function(x){
                     tst <- apply(X=array(1:n),
                                  MARGIN=1,
                                  FUN=function(x,y){
                                        if(x == y) return(FALSE)
                                        isTRUE(all.equal(t[x],t[y],tolerance=tol))
                                      }, y=x)
                     return(any(tst))
                   })

  if( any(tst) ) stop("Measurement times (t) must be distinct.")

  if( any(t < tol) ) stop("All measurement times (t) must be positive.")


  m <- dim(SigmaTheta)
  if( m[1] != m[2] ) stop("SigmaTheta must be a square matrix.")

  if( !isSymmetric(SigmaTheta) ) stop("SigmaTheta must be a symmetric matrix.")

  eig <- eigen(SigmaTheta, symmetric=TRUE, only.values=TRUE)
  if( any(eig$values < tol) ) {
    stop("SigmaTheta must be a positive-definite matrix.")
  }

  if( sigmae_2 < tol ) stop("Measurement error, sigmae_2, must be positive.")

  if( ordtraj >= length(eig$values) || ordtraj < 0) {
    stop("Trajectory order, ordtraj, must be [0,ncol(SigmaTheta)-1].")
  }

  if( (alpha < tol) || (alpha > 1.0-tol) ) {
    stop("The significance level, alpha, must be (0.0, 1.0)." )
  }

  return(NULL)

} # End of error.powerLongSurv

