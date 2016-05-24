#    Model Implementation:
#       sse()         Return BK's 'stupid simple estimator': % ever found
#       recur()       Recursive calculation   needed in Rst() below
#       simp()        Simpson-rule integrator needed in Rst() below
#       Fs()          Integrand function      needed in Rst() below
##################################################################
#
# Scripts to IMPLEMENT efficient evaluation of terms comprising
# R* from "1eqn" document, and also the "expansion factor"
# kappa = Iij/R*. 

##################################################################
#
sse <- function(rd) {
  # Stupid Simple Estimator
  # Brown, Smallwood, & Karas (2013), section 2.3.2.2
  # "Overall detection rate" D
  ids <- rd$scav$Id;      # List of carcass Id strings
  sch <- rd$srch[,-2];    # Search results
  fnd <- logical(no <- length(ids));
  for(i in 1:no) {
    ok     <- sch$Id == ids[i];
    fnd[i] <- any(sch$Found[ok]);
  }
  return( sum(fnd)/no );
}  

##################################################################
#
recur <- function(kmax=3) {
  # Which terms (k,m,n) occur in expression of R* as sum of
  #             bt^k * (-1)^(m+1) * Q[k,m,n] ???
  
  kmn <- matrix(NA,nrow=(nr <- 2^kmax-1),ncol=3);
  dimnames(kmn)[[2]]<-c("k","m","n");
  kmn[1,] <- c(0,1,0);
  old <- 1;
  new <- 2;
# From (k,m,n) we get:  (k+1,m,n+1) and (k+1,m+1,n+k+1)
  while(new <= nr) {
    k <- kmn[old,1]+1;          # Note this is old (k+1)
    kmn[new,] <- kmn[old,] + c(1,0,1); new <- new+1;
    kmn[new,] <- kmn[old,] + c(1,1,k); new <- new+1;
    old <- old+1;
  }
  return(kmn);
}

##################################################################
#
simp <- function(f, npts=100, ...) {
  # Simpson's Rule for integral on [0,1]
  # Assumes f() can take (and return) vectors.
  
   npts <- 2 * (mpts <- ceiling(npts/2));
   Even <- seq(0,1,,mpts+1);        # = {0/n, 2/n,   ...,    n/n}
   Odd  <- seq(1/npts,1,1/mpts);    # = {  1/n,  3/n,...,(n-1)/n}
   rv   <- mean(f(Odd,...))*(2/3) + # = Odds*4/3 + Evens*2/3 - Ends/3
          (2*sum(f(Even,...))-sum(f(0:1,...)))/(3*npts);
   return(rv);
}

##################################################################
#
Fs <- function(x, kmn=c(k=0,m=1,n=0), pars=c(a=1,bI=0,alp=1,rI=1)) {
  # Weibull integrand for Qst calculation.  Args are 0<x<1 (vector okay),
  # kmn=c(k,m,n), and pars=c(a,bI, alp,rI, the) where bI=b*I, rI=r*I
  
  a <- pars[1]; bI <- pars[2]; alp <- pars[3]; rI <- pars[4];
  return(exp( -(rI*(kmn["k"]+x))^alp - kmn %*% rbind(0,a+bI*x,bI)));
}