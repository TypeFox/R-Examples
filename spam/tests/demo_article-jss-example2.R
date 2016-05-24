# This is file ../spam/tests/demo_article-jss-example2.R
# This file is part of the spam package, 
#      http://www.math.uzh.ch/furrer/software/spam/
# by Reinhard Furrer [aut, cre], Florian Gerber [ctb]
     


# INITALIZE AND FUNCTIONS:
# JSS article:
#     "spam: A Sparse Matrix R Package with Emphasis on
#            MCMC Methods for Gaussian Markov Random Fields"


# Compared to the R code in the article, here we give:
# - improved formatting
# - more comments, e.g. how to run the code using regular matrices 
# - the code to construct the figures
# - minor modifcations due to evolvement of spam
 

cat("\nThis demo contains the R code of the second example\nin the JSS article. As pointed out by Steve Geinitz\nand Andrea Riebler, the Gibbs sampler is not correct\nand contains several bugs. \n\nI'll post an updated sampler in a future release.\n\n") 


# INITALIZE AND FUNCTIONS:
require("fields", warn.conflict=FALSE)
spam.options(structurebased=TRUE)


# READ DATA:
attach(Oral)



# CONSTRUCT ADJACENCY MATRIX:
loc <- system.file("demodata/germany.adjacency", package="spam")
A <- adjacency.landkreis(loc)
n <- dim(A)[1]
# Verification that we have a symmetric matrix:
# norm(A-t(A)); display(A)


# GIBBS SETUP:
set.seed(14)

# Construct the individual block precisions
# (based on unit precision parameters kappa, denoted with k):

Q1 <- R <- diag.spam( diff(A@rowpointers)) - A   # this is R in (2)
pad(Q1) <- c(2*n,2*n)  # previously:  dim(Q1) <- c(2*n,2*n)

Q2 <-  rbind(cbind( diag.spam(n), -diag.spam(n)),
       	     cbind(-diag.spam(n),  diag.spam(n)))

# Hence the precision Q in (2) is:
# Q <- kappau*Q1 + kappav*Q2

# pre-define
diagC <- as.spam( diag.spam(c(rep(0,n),rep(1,n))))


# Recall:
# k=( kappa_u, kappa_y)'

# hyperparameters
ahyper <- c( 1, 1)
bhyper <- c( .5, .01)


# Gibbs sampler
burnin <- 50
ngibbs <- 150
totalg <- burnin+ngibbs

# Initialize parameters:
upost <- array(0, c(totalg, n))
npost <- array(0, c(totalg, n))
kpost <- array(0, c(totalg, 2)) 

# Starting values:
kpost[1,] <- c(40,500)
upost[1,] <- u <- rnorm(n,sd=.2) *1
npost[1,] <- eta <- u + rnorm(n,sd=.05)*1

uRu <- t(u) %*% (R %*% u)/2
etauetau <- t(eta-u) %*% (eta-u)/2 

postshape <- ahyper + c(n-1,n)/2

accept <- numeric(totalg)

struct <- chol(Q1 + Q2 + diag.spam(2*n),
               memory=list(nnzcolindices=5500))

# struct <- NULL        # If no update steps are wanted

# R <- as.matrix(R)     # If no spam analysis is wanted.
# Q1 <- as.matrix(Q1)
# Q2 <- as.matrix(Q2)


for (ig in 2:totalg) {

  
  kstar <- rgamma(2,postshape, bhyper + c(uRu, etauetau))	

  
  expeta0E <- exp(eta)*E
  expeta0Eeta01 <- expeta0E *(eta-1)
  diagC@entries <- expeta0E
  Q <- kstar[1]*Q1 + kstar[2]*Q2 + diagC
  b <- c( rep(0,n), Y + expeta0Eeta01)

  xstar <- rmvnorm.canonical(1,
                             # vector b:
                             b,
                             # Precision matrix
                             Q,
                             Rstruct=struct)
  

  ustar <- xstar[1:n]
  nstar <- xstar[1:n+n]

  uRustar <- t(ustar) %*% (R %*% ustar)/2
  etauetaustar <- t(nstar-ustar) %*% (nstar-ustar)/2

  
# we work on the log scale:
# logalpha <- min(0, log(ratios))=min(0, expterm+(...)log(kappa)-
  
  exptmp <- sum(expeta0Eeta01*(eta-nstar) - E*(exp(eta)-exp(nstar))) -
    sum( nstar^2*expeta0E)/2   +    sum(eta^2*expeta0E)/2 -
      kstar[1] * uRu           +    kpost[ig-1,1] * uRustar -
        kstar[2] * etauetau    +    kpost[ig-1,2] * etauetaustar
  factmp <- (postshape-1)*(log(kstar)-log(kpost[ig-1,1]))
  
  logalpha <- min(0, exptmp + sum(factmp))
  logU <- log(runif(1))

  if (logU < logalpha) { # ACCEPT draw
    upost[ig,] <- u   <- ustar
    npost[ig,] <- eta <- nstar
    kpost[ig,] <- kstar
    uRu <- uRustar
    etauetau <- etauetaustar
    accept[ig] <- 1
  } else {
    upost[ig,] <- upost[ig-1,]
    npost[ig,] <- npost[ig-1,]
    kpost[ig,] <- kpost[ig-1,]    
  }
                   
  if( (ig%%10)==0) cat('.')

}


if (FALSE) {

# POSTPROCESSING:

accept <- accept[-c(1:burnin)]
cat("\nAcceptance rate:",mean(accept),"\n")

kpost <- kpost[-c(1:burnin),]
upost <- upost[-c(1:burnin),]
npost <- npost[-c(1:burnin),]

kpostmean <- apply(kpost,2,mean)
upostmean <- apply(upost,2,mean)
npostmean <- apply(npost,2,mean)

kpostmedian <- apply(kpost,2,median)
upostmedian <- apply(upost,2,median)
npostmedian <- apply(npost,2,median)

vpost <- npost-upost
vpostmedian <- apply(vpost,2,median)


#





######################################################################
# Figures 
par(mfcol=c(1,3),mai=rep(0,4))
map.landkreis(log(Y))
map.landkreis(Y/E,zlim=c(.1,2.4))
map.landkreis(exp(upostmedian),zlim=c(.1,2.4))


par(mfcol=c(2,4),mai=c(.5,.5,.05,.1),mgp=c(2.3,.8,0))
hist(kpost[,1],main="",xlab=expression(kappa[u]),prob=TRUE)
lines(density(kpost[,1]),col=2)
tmp <- seq(0,to=max(kpost[,1]),l=500)
lines(tmp,dgamma(tmp,ahyper[1],bhyper[1]),col=4)
abline(v=kpostmedian[1],col=3)

hist(kpost[,2],main="",xlab=expression(kappa[y]),prob=TRUE)
lines(density(kpost[,2]),col=2)
tmp <- seq(0,to=max(kpost[,2]),l=500)
lines(tmp,dgamma(tmp,ahyper[2],bhyper[2]),col=4)
abline(v=kpostmedian[2],col=3)

# Trace plots:
plot(kpost[,1],ylab=expression(kappa[u]),type="l")
abline(h=kpostmedian[1],col=3)
plot(kpost[,2],ylab=expression(kappa[y]),type="l")
abline(h=kpostmedian[2],col=3)

# ACF:
acf(kpost[,1],ylab=expression(kappa[u]))
acf(kpost[,2],ylab=expression(kappa[y]))



# scatter plots
plot(kpost[,1],kpost[,2],xlab=expression(kappa[u]),ylab=expression(kappa[y]))
abline(v=kpostmedian[1],h=kpostmedian[2],col=3)


plot(accept+rnorm(ngibbs,sd=.05),pch=".",ylim=c(-1,2),yaxt="n",ylab="")
text(ngibbs/2,1/2,paste("Acceptance rate:",round(mean(accept),3)))
axis(2,at=c(0,1),label=c("Reject","Accept"))

}

detach(Oral)
######################################################################

