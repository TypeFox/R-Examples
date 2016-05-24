"plotNsibs" <-
function(alpha=0.001,power=0.85,...)
{
# Plot sample size in an affected sibpair linkage study
# as a function of lambda_sib
 tb <- NULL
 beta <- 1- power
 ls.tb <- seq(1.1,2.0,0.01)
 for (ls in ls.tb)
  tb <- rbind(tb,c(ls,nsibs(ls,ls,alpha,beta)))
 tb <- data.frame(tb)
 colnames(tb) <- c("lambda.sib","number.ASP")
 plot(tb$lambda.sib,tb$number.ASP,type="l",xlab="lambda.off",ylab="number.ASP",
 main =paste("Required sample size for power =",power,"at alpha =",alpha),...)
 return(tb)
}

