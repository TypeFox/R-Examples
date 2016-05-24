# -------------------------------------------------------------------------
# helper functions for power.dp / sampleN.dp and others
# -------------------------------------------------------------------------
# distribute n into grps with 'best' balance between grps
# this function is also used in other power functions
nvec <- function(n, grps)
{
  ni <- trunc(n/grps)
  nv <- rep.int(ni, times=grps)
  rest <- n-grps*ni
  if(rest!=0){
    nv <- nv + c(rep.int(1,rest), rep.int(0,grps-rest))
  }
  nv
}
# --------------------------------------------------------------------------
# css without ni, assuming all ni equal
# --------------------------------------------------------------------------
.css <- function(doses, design, dm=NULL)
{
  ld  <- log(doses)
  
  if (design=="IBD"){
    desld <- dm
    nseqs <- nrow(dm)
    css2   <- 0
    for (l in 1:nseqs){
      desld[l,] <- ld[dm[l,]]
      css2 <- css2 + sum((desld[l,]-mean(desld[l,]))^2)
    }
    #print(css2)
    # next may be questionable
    desld <- as.vector(desld)
    css <- sum((desld-mean(desld))^2)
  } else {
    css <- sum((ld-mean(ld))^2)
    if (design=="crossover") css <- length(doses)*css
  }
  return(css)
}
# --------------------------------------------------------------------------
# css with ni, IBD as simple fixed effects regression
.css2 <- function(doses, design, dm=NULL, n)
{
  ld  <- log(doses)
  
  if (design=="IBD"){
    css2   <- 0
    nseqs <- nrow(dm)
    # css of one long vector of log doses
    # correct? or sum of css in sequences?
    # .css3() gives same result as long vector if omega2=0
    desld <- vector(mode="numeric")
    for (l in 1:nseqs){
      desld2    <- ld[dm[l,]]
      desld <- c(desld, rep(desld2, times=n[l]))
      css2   <- css2 + n[l]*sum((desld2-mean(desld2))^2)
    }
    desld <- as.vector(desld)
    css <- sum((desld-mean(desld))^2)
  } else {
    css <- sum(n*(ld-mean(ld))^2)
    if (design=="crossover") css <- length(doses)*css
  }
  return(css)
}

# css with ni, IBD as mixed(?) effects regression (random intercept(?))
# gives same results as fixed effects if Latin square for dm is used
# gives same results as fixed effects if omega2=0
# CSS smaller (var(beta) >) as fixed effects if omega2 > 0
# how to get a reasonable estimate of omega2?
.css3 <- function(doses, design, dm=NULL, n, s02, omega2=0.8)
{
  ld  <- log(doses)
  
  if (design=="IBD"){
    #Sethuraman et al. formula 12
    lamda <- diag(2)
    lamda[2,2] <- 0
    lamda <- lamda*omega2/s02
    per   <- ncol(dm)
    Nseq  <- nrow(dm)
    D <- matrix(0, nrow=2, ncol=2)
    for (l in 1:Nseq){
      # first row all elements=1 
      Xl     <- matrix(1, nrow=2, ncol=per)
      # second row log(doses)
      # replace index by log(doses)
      Xl[2,] <- ld[dm[l,]]
      # inverse of Xl*t(Xl)
      Ml1    <- solve(Xl%*%t(Xl))
      #why again invert?
      D      <- D + n[l]*solve(Ml1+lamda)
    }
    # and again invert?
    D  <- solve(D)
    # return equivalent of css
    return(1/D[2,2])
  } else {
    # parallel or crossover
    css <- sum(n*(ld-mean(ld))^2)
    if (design=="crossover") css <- length(doses)*css
  }
  return(css)
}
