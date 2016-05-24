RPvectors <-
function(a,m,ortho="none",distr="uniform", 
   par_unif=c(-1,1),
   par_norm=c(0,1),
   par_eq=c(-1,0,1),
   par_uneq=c(-sqrt(3),0,sqrt(3)),
   par_uneqprob=c(1/6,2/3,1/6))
{
# Generate random projection vectors

# ARGUMENTS:
#
# a           number of generated vectors (>=1)
# m           dimension of generated vectors (>=2)
#
# ortho       orthogonalization of vectors
#                "none"  no orthogonalization (default)
#                "onfly"  orthogonalization on the fly
#                "end"  orthogonalization at the end
#
# distr       distribution of generated random vector components
#                "uniform"  uniformly distributed in range par_unif (see below),
#                   default U[-1, +1]
#                "normal"  normally distributed with parameters par_norm (see below),
#                   typical N(0, 1)
#                "randeq"  random selection of values par_eq (see below) with
#                   equal probabilities,
#                   typical {-1, 0, +1}
#                "randuneq"  random selection of values par_uneq (see below) with
#                   probabilties par_uneqprob (see below)
#                   typical -(3)^0.5 with probability 1/6
#                             0      with probability 2/3
#                           +(3)^0.5 with probability 1/6
#
# par_unif       parameters for range for distr=="uniform",
#                   default to c(-1,1)
# par_norm       parameters for mean and sdev for distr=="normal",
#                   default to c(0,1)
# par_eq         values for distr=="randeq" which are replicated, 
#                   default to c(-1,0,1)
# par_uneq       values for distr=="randuneq" which are replicated with
#                   probabilties par_uneqprob
#                   default to c(-sqrt(3),0,sqrt(3))
# par_uneqprob   probabilities for distr=="randuneq" to replicate values par_uneq,
#                   default to c(1/6,2/3,1/6))
#
# VALUE:
#
# returned B       matrix with generated vectors (each in a column)
#                  m rows, a columns, 
#                  vectors are scaled to unit length
####################################################################################

genvec <- function(a,m,distr,par_unif,par_norm,par_eq,par_uneq,par_uneqprob)
# Function to generate random vectors:
{
  if (distr=="uniform"){
    B <- matrix(runif(m*a,par_unif[1],par_unif[2]),ncol=a)
  }
  else if (distr=="normal"){
    B <- matrix(rnorm(m*a,par_norm[1],par_norm[2]),ncol=a)
  }
  else if (distr=="randeq"){
    B <- matrix(sample(par_eq,m*a,replace=TRUE),ncol=a)
  }
  else if (distr=="randuneq"){
    B <- matrix(sample(par_uneq,m*a,replace=TRUE,prob=par_uneqprob),ncol=a)
  }
  else {
    stop("Method to generate random vectors not defined!")
  }
return(B)
}

# start of function

if ((a<1) | (m<2)){stop("Increase a and/or m!")}

if (a>m) { # more directions than dimensions only possible for ortho="none"
  if (ortho=="onfly" | ortho=="end"){
    a <- m
    warning("a has been set equal to m")
  }
}

if (ortho=="none"){ # no orthogonalization
  B <- genvec(a,m,distr,par_unif,par_norm,par_eq,par_uneq,par_uneqprob)
  B <- t(t(B)/sqrt(apply(B^2,2,sum))) # length 1
}
else if (ortho=="onfly"){ # orthogonalize on the fly
  B <- matrix(NA,nrow=m,ncol=a) # matrix with loadings in original space
  for (j in 1:a){  # loop over number of RP directions
    b <- genvec(a=1,m,distr,par_unif,par_norm,par_eq,par_uneq,par_uneqprob) 
    b <- b/sqrt(sum(b^2)) # norm 1
    if (j==1){B[,j] <- b}
    else {
      B[,j] <- b-B[,1:(j-1)]%*%t(B[,1:(j-1)])%*%b
      B[,j] <- B[,j]/sqrt(sum(B[,j]^2))
    }
  }
}
else if (ortho=="end"){ # orthogonalize at the end
  B <- matrix(NA,nrow=m,ncol=a) # matrix with loadings in original space
  b <- genvec(a,m,distr,par_unif,par_norm,par_eq,par_uneq,par_uneqprob) 
  b <- t(t(b)/sqrt(apply(b^2,2,sum))) # norm 1
  for (j in 1:a){  # loop over number of RP directions
    if (j==1){B[,j] <- b[,j]}
    else {
      B[,j] <- b[,j]-B[,1:(j-1)]%*%t(B[,1:(j-1)])%*%b[,j]
      B[,j] <- B[,j]/sqrt(sum(B[,j]^2))
    }
  }
}
else {
  stop("Method not defined by parameter ortho!")
}

return(B)
}

