Fpower2<-function(alpha=NULL, nlev=NULL,nreps=NULL, Delta=NULL, sigma=NULL)
{
##### Power Calculation for two way ANOVA ###########
# Argument list
# alpha the significance level of the test.
# nlev vector containing the number of levels of the factors. 
# nreps the number of replicates in each combination of factor levels.
# Delta the size of a practical difference in two marginal factor level means.
# sigma the standard deviation of the experimental error.
############################################################
if (is.null(alpha)|is.null(nlev)|is.null(nreps)|is.null(Delta)|is.null(sigma))
  stop("you must supply alpha, nlev, nreps, Delta and sigma")
if(length(nlev)<2)
  stop ("nlev must be a two component vecto containing levels of the 1st and 2nd factors")
a <- nlev[1]
b <- nlev[2]
cssb <- (Delta^2)/2
ncb <- a*(nreps*cssb)/(sigma^2)
cssa<-(Delta^2)/2
nca<- b*(nreps*cssa)/(sigma^2)
dfa<- a-1
dfb<- b-1
df2<-(nreps-1)*b*a
powera <- 1-pf(Fcrit(alpha,dfa,df2),dfa,df2,nca)
powerb <- 1-pf(Fcrit(alpha,dfb,df2),dfa,df2,nca)
result <-cbind(alpha,a,b,nreps,Delta,sigma,powera,powerb)
}
