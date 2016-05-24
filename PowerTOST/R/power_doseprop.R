#-----------------------------------------------------------------------------
# power for dose proportionality studies via power model
# ----------------------------------------------------------------------------
# degrees of freedom:
# crossover: total    = prds*nt-1
#            subjects = nt-1
#            periods  = prds-1
#            regr.    = 1
#            error    = prds*nt-1 - (nt-1) - (prds-1) -1 = prds*nt - nt - prds
# parallel:  total    = nt-1
#            regr     = 1
#            error    = nt-2
#library(PowerTOST)

power.dp <- function(alpha=0.05, CV, doses, n, beta0, theta1=0.8, theta2=1/theta1,
                     design=c("crossover", "parallel", "IBD"), dm=NULL, CVb)
{
  desi <- match.arg(design)
  
  grps <- length(doses) # dose groups and periods in case of crossover
  if (desi=="IBD"){
    # check
    if(!is.matrix(dm)) stop("Design matrix has to be given.")
    if(ncol(dm)>length(doses)) stop("Design matrix must have <=", length(doses),
                                    "columns (periods).")
    if(any(dm>grps) | any(dm<1)) stop("Wrong content of 'design' matrix.")
    grps <- nrow(dm) # sequence groups
    if (missing(CVb)) CVb <- mse2CV(0.8) # sethurahman et al.
  } else {
    if (missing(CVb)) CVb <- 0           # don't need CVb for the other designs
  }
  if (length(doses)<=1) stop("At least two doses have to be given.")
  
  if (CV<=0) stop("CV must be greater then zero.")
  s2   <- CV2mse(CV)

  if (theta1<=0 | theta2<=0) stop("theta1/theta2 must be greater then zero.")
  
  if (length(n)==1){
    # then we assume n=ntotal
    # for unbalanced designs we divide the n's by ourself
    # to have only small imbalance
    n  <- nvec(n=n, grps=grps)
    # give a message?
    if (n[1]!=n[length(n)]){
      message("Unbalanced design. n(i)=", paste(n, collapse="/"), " assumed.")
    } 
  } 
  # else n is the vector of subjects in (sequence) groups
  if (length(n)!=grps) stop("n as vector must have length ", grps,".")
  nt   <- sum(n)
  if(nt<2) stop("n(total) has to be >2.")
  # periods
  prds <- ifelse(desi=="parallel", 1, grps)
  prds <- ifelse(desi=="IBD", ncol(dm), prds)
  # degrees of freedom
  df   <- ifelse(desi=="parallel", nt-2, (nt*prds)-(nt+prds-1)-1)
  # range of doses as ratio 
  rd   <- max(doses)/min(doses)
  # acceptance range for beta
  bl   <- 1+log(theta1)/log(rd)
  bu   <- 1+log(theta2)/log(rd)
  if (missing(beta0)) beta0 <- 1+log(0.95)/log(rd)
  if (beta0<=0) stop("beta0 must be greater then zero.")
  # log doses corrected sum of squares
  Sdd   <- .css3(doses=doses, design=desi, dm=dm, n=n, s02=CV2mse(CV), 
                 omega2=CV2mse(CVb))
  # variance of slope
  vbeta <- s2/Sdd
  
  tval  <- qt(1-alpha, df)
  
  # non-centrality parms according to Patterson/Jones
#  nc1 <- (sqrt(nt))*((beta0-bl)/sqrt(vbeta))
#  nc2 <- (sqrt(nt))*((beta0-bu)/sqrt(vbeta))
  # question: where sqrt(nt) comes from?
  # only without sqrt(nt) the 'power' calculations of Hummel et. al 
  # (large sample approx. ?) will be obtained!
  # and only then the results of power.dp() and power.TOST() in case of two doses 
  # coincide if beta0=1 / theta0=1 are used
  nc1 <- (beta0-bl)/sqrt(vbeta)
  nc2 <- (beta0-bu)/sqrt(vbeta)
  
  # nct approximation
  pwr <- pmax(pt(-tval, df, nc2) - pt(tval, df, nc1), 0)

  return(pwr)
  
}

# --------------------------------------------------------------------------
# power function of Hummel et.al, large sample approx., parallel group
# seems alpha has to be set to 2*alpha
# --------------------------------------------------------------------------
power.dpLS <- function(alpha=0.05, CV, doses, n, beta0=1, theta1=0.8, 
                       theta2=1/theta1)
{
  s2   <- CV2mse(CV)
  
  rd   <- max(doses)/min(doses)
  bl   <- 1+log(theta1)/log(rd)
  bu   <- 1+log(theta2)/log(rd)
  
  grps <-length(doses)
  if (length(n)==1){
    # then we assume n=ntotal
    # for unbalanced designs we divide the ns by ourself
    # to have only small imbalance
    n  <- nvec(n=n, grps=grps)
  }  
  ld    <- log(doses)
  meand <- mean(ld)
  Sdd   <- sum(n*(ld-mean(ld))^2) 
  
  w   <- Sdd/s2
  
  u   <- qnorm(1-alpha) # original was alpha/2
  pwr <- pnorm(-u - (beta0-bu)*sqrt(w)) - pnorm(u - (beta0-bl)*sqrt(w))
  
  return(pmax(pwr,0))
  
} # end function
