`msBIC` <-
function(x, K, method = "all", Sdist="weibull", cutpoint = NULL, EMoption="classification",
                  EMstop=0.01, maxiter=100)
{

# Input:
#        x ........ n x p data-matrix each column representing dwell-time of an individual
#                   session in site-area 1 to p. pages not visited should have a dwell time of 0 or NA
#        K ........ scalar or vector with number of mixture components
#        method...  denotes the model to be fitted: possible values are "all", "separate", "main.g", "main.p", "int.gp","main.gp",
#                   while in method "separate" distributions of individual groups are
#                   estimated independently, method "main.g" assumes that there is a common
#                   base-line hazard which is common to all groups, "main.p" the same for the sites. method "main.gp"
#                   fits a main effects model whereas in model "int.gp" allows interaction effects. "all" fits all of the models above.
#        Sdistr ... Survival distribution to be fitted; options are "weibull", "exponential", "rayleigh".
#        EMoption.. "classification" is based on the probabilistic cluster assignment, "maximization" on deterministic assignment

if (is.data.frame(x)) x <- as.matrix(x)

if ((length(method)==1) && (method=="all")) method <- c("separate", "main.g", "main.p", "int.gp","main.gp")

BICmat <- matrix(NA,length(method),length(K))

for (i in 1:length(method)) {
  mi <- method[i]
  for (j in 1:length(K)) {
    Kj <- K[j]
    resall <- phmclust(x,K=Kj,method=mi,Sdist=Sdist,cutpoint = cutpoint, EMstart=NA,EMoption=EMoption,EMstop=EMstop, maxiter=maxiter)
    BICmat[i,j] <- resall$bic
  }
}


#BIClist <- tapply(method,1:length(method), function(mx) {
#             res <- tapply(K,1:length(K), function(Kx) {
#                       phmclust(x,K=Kx,method=mx,Sdist=Sdist,EMstart=NA, 
#                                 EMoption=EMoption,EMstop=EMstop, maxiter=maxiter)$bic})})
#BICmat <- matrix(unlist(BIClist),ncol=length(K),byrow=TRUE)

result <- list(BICmat=BICmat,K=K,method=method,Sdist=Sdist)
class(result) <- "BICmat"
result
}

