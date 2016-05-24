`phmclust` <-
function(x, K, method = "separate", Sdist = "weibull", cutpoint = NULL, EMstart = NA, EMoption = "classification", EMstop = 0.01, maxiter = 100)
{

# Input:
#        x ........ n x p data-matrix each column representing dwell-time of an individual
#                   session in site-area 1 to p. pages not visited should have a dwell time of 0 or NA
#        K ........ scalar with number of mixture components
#        method.... denotes the model to be fitted: possible values are "separate", "main.g", "main.p", "int.gp","main.gp", 
#                   while in method "separate" distributions of individual groups are 
#                   estimated independently, method "main.g" assumes that there is a common
#                   base-line hazard which is common to all groups, "main.p" the same for the sites. method "main.gp" 
#                   fits a main effects model whereas in model "int.gp" allows interaction effects.
#        cutpoint.. Integer value with upper bound for observed dwell time. Above this cutpoint values are regarded as censored. If NULL, no censoring.   
#        EMstart .. Vector of length n with starting values for group membership
#        EMoption.. "classification" is based on the probabilistic cluster assignment, "maximization" on deterministic assignment,
#                   "randomization" provides a posterior-based randomized cluster assignement. 
#        Sdist .... Survival distribution for WPHM model. These include "weibull", "exponential", "rayleigh". 
#        
# Output:
#        list containig the following elements
#
#
# Packages Needed: MASS, Survival

if (is.data.frame(x)) x <- as.matrix(x)
if (is.vector(x)) stop("x must be a data frame or a matrix with more than 1 columns!")
if (is.null(cutpoint)) cutpoint <- max(x, na.rm = TRUE) 

pvisit.est <- (any(is.na(x)) | any(x==0))                          # TRUE if visiting prob estimated
n <- nrow(x)                                                       # n ... number of sessions
p <- ncol(x)                                                       # p ... number of site-areas                            

d0 <- s.check(x=x,K=K,n=n,EMstart=EMstart,EMoption=EMoption,method=method,Sdist=Sdist)    #sanity checks

x <- d0$x
if (is.vector(d0$EMstart)) {
  EMstart <- d0$EMstart[1:(dim(x)[1])]
} else {
  EMstart <- t(apply(d0$EMstart, 1, function(z) z/sum(z)))    #sum 1 normalization
}

 
#EMstart <- d0$EMstart
method <- d0$method
                                                                                                    
likelihood <- numeric(maxiter)
iter <- 0                                
ConvergEM <- FALSE  


#=============================EM-estimation================================

if (EMoption == "classification") {                                  #maximization EM
    while (ConvergEM == FALSE)
    {
       iter <- iter + 1
       #print(iter)
       d1 <- Eclass(x, EMstart, K = K, method = method, Sdist = Sdist, p, cutpoint = cutpoint)	     #E-Step maximization
       d2 <- Mclass(x, d1$shape,d1$scale,d1$prior,K=K)                    #M-Step maximization
       
       likelihood[iter+1] <- d2$lik.tot                            #likeihood in the current iteration  
       
       if ((iter >= maxiter) || (abs(likelihood[iter+1]-likelihood[iter]) < EMstop)) {   
          ConvergEM <- TRUE 
       } else {
          EMstart <- d2$newgr
          }
    }
postmat <- NULL
newgr <- d2$newgr
}

#========================================================================

if (EMoption == "maximization") {                                #classification EM
    while (ConvergEM == FALSE)
    {
       iter <- iter + 1
       d1 <- Emax(x, EMstart, K=K, method=method, Sdist=Sdist,p, cutpoint = cutpoint)                          #E-Step maximization
       d2 <- Mmax(x, d1$shape,d1$scale,d1$prior,K=K)                  #M-Step maximization
       
       likelihood[iter+1] <- d2$lik.tot                            #likeihood in the current iteration  
       
       if ((iter >= maxiter) || (abs(likelihood[iter+1]-likelihood[iter]) < EMstop)) {   
          ConvergEM <- TRUE 
       } else {
          EMstart <- d2$postmat
       }
    }
postmat <- d2$postmat
newgr <- apply(postmat, 1, function(y){ind <-(1:K)[y==max(y)]})    #final group assignement                         
}

#========================================================================

if (EMoption == "randomization") {                                  #maximization EM
    while (ConvergEM == FALSE)
    {
       iter <- iter + 1
       d1 <- Eclass(x, EMstart, K=K, method=method, Sdist=Sdist,p, cutpoint = cutpoint)	     #E-Step randomization (=classification)
       d2 <- Mrandom(x, d1$shape,d1$scale,d1$prior,K=K)                    #M-Step maximization
       
       likelihood[iter+1] <- d2$lik.tot                            #likeihood in the current iteration  
       
       if ((iter >= maxiter) || (abs(likelihood[iter+1]-likelihood[iter]) < EMstop)) {   
          ConvergEM <- TRUE 
       } else {
          EMstart <- d2$newgr
       }
    }
postmat <- d2$postmat
newgr <- d2$newgr
}
#========================================================================

if (iter >= maxiter) warning("EM did not converge! Maximum iteration limit reached!")

if (pvisit.est) { 
   anzpar <- d1$anzpar + K*p       #if NA's in x --> number of estimated visiting probabilities added
} else {
   anzpar <- d1$anzpar             #no NA's in x
}


likconv <- likelihood[2:(iter+1)]
aic <- -2*(likelihood[iter+1]-anzpar)
bic <- (-2*likelihood[iter+1])+anzpar*log(n)

#-------------------- cluster means -----------------
clmean.l <- by(x,newgr,function(y) {
                          apply(y,2,function(z) mean(z[z>0]))   #compute cluster means (eliminate x==0)
                         })
clmean <- matrix(unlist(clmean.l),nrow=K,byrow=TRUE)
nobs.cl <- table(newgr)
se.clmean <- apply(x, 2, function(z1) {                        #standard errors for the cluster means (conditional on the membership)
               Smat <- tapply(z1, newgr, function(z2) {
                                           if (!is.null(z2[z2 > 0])) {
                                               sd(z2[z2>0], na.rm = TRUE)   #sample sd
                                             } else {
                                               return(NA)
                                             }
                                         })
               nobsmat <- tapply(z1, newgr, function(z3) {
                                           if (!is.null(z3[z3 > 0])) {
                                              length(z3[z3>0])   #number of observations
                                            } else {
                                              return(NA)
                                            }
                                         })
               Smat/sqrt(nobsmat)                              #s/sqrt(n)
             })
colnames(clmean) <- colnames(se.clmean) <- colnames(x)
rownames(clmean) <- rownames(se.clmean) <- paste("Cluster",1:K)
#se.clmean <- NULL

#------------------- cluster medians ----------------
clmed.l <- by(x,newgr,function(y) {
                        apply(y,2, function(z) median(z[z>0]))
                      })
clmed <- matrix(unlist(clmed.l),nrow=K,byrow=TRUE)
colnames(clmed) <- colnames(x)
rownames(clmed) <- paste("Cluster",1:K)

if (is.null(colnames(x))) {
  colnames(d1$scale) <- paste("V",1:dim(d1$scale)[2],sep="")
  colnames(d1$shape) <- paste("V",1:dim(d1$shape)[2],sep="")
} else {
  colnames(d1$scale) <- colnames(d1$shape) <- colnames(x)
}

rownames(d1$shape) <- rownames(d1$scale) <- paste("Cluster", 1:K, sep="")

#------------------- pvisit s.e. ----------------
pvisit <- d1$prior
se.pvisit <- apply(pvisit,2, function(z) sqrt(z*(1-z)/nobs.cl))
colnames(pvisit) <- colnames(se.pvisit) <- colnames(x)
rownames(pvisit) <- rownames(se.pvisit) <- paste("Cluster",1:K)

result <- list(K=K, iter=iter, method=method, Sdist=Sdist, likelihood=likconv, pvisit = pvisit, se.pvisit = se.pvisit,
               shape = d1$shape, scale = d1$scale, group = newgr, posteriors = postmat, npar = anzpar, aic = aic, bic = bic,
               clmean = clmean, se.clmean = se.clmean, clmed = clmed)
class(result) <- "mws"
result
}

