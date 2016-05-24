Fpower1<-function(alpha=NULL, nlev=NULL,nreps=NULL, Delta=NULL, sigma=NULL)
{
##### Power Calculation for one way ANOVA ###########
# Argument list
# alpha the significance level of the test
# nlev the number of levels of the factor 
# nreps the number of replicates in each level of the factor
# Delta the size of a practical difference in two cell means
# sigma the standard deviation of the experimental error
#####################################################
if (is.null(alpha)|is.null(nlev)|is.null(nreps)|is.null(Delta)|is.null(sigma))
  stop("you must supply alpha, nlev, nreps, Delta and sigma")
css<-(Delta^2)/2
nc<- (nreps*css)/(sigma^2)
df1<-nlev-1
df2<-(nreps-1)*nlev
power <- 1-pf(Fcrit(alpha,df1,df2),df1,df2,nc)
result <- cbind(alpha,nlev,nreps, Delta, sigma, power)
return(result)
}
