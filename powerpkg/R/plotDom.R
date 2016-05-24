"plotDom" <-
function(lsib=1.2,alpha=0.001,power=0.85,...)
{
# Plot the effects of dominance on sample size in 
# an affected sibpair linkage study
# as a function of lambda_sib
 tb <- NULL
 beta <- 1- power
 l1.tb <- seq(1.0,lsib,0.01)
 for (l1 in l1.tb)
  tb <- rbind(tb,c(l1,nsibs(ls=lsib,lo=l1,alpha,beta)))
 tb <- data.frame(tb)
 colnames(tb) <- c("lambda.off","number.ASP")
 plot(tb$lambda.off,tb$number.ASP,type="l",
 main =paste("Required sample size for power =",power,"at alpha =",alpha),
 xlab="lambda.off",ylab="number.ASP",...)
 return(tb)
}

