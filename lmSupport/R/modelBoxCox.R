modelBoxCox <-
function(Model, Lambdas = seq(-2,2,by = 0.1))
#Provides plot of lambda for power transformation of response variable.  Tests if lambda <> 1
#all values of Y must be positive
#can plot smaller range of lambda as needed
#Box, G. E. P. & Cox, D. R. (1964) An analysis of transformations (with discussion).  Journal of the Royal Statistical Society, 26, 211-252.
{
  #note, changed from boxcox to boxCox to use car function
   LR <- boxCox(Model, lambda=Lambdas)  #NOTE: will crash if response variable values <= 0

   #Chi-square test to determine if log-likelihood is sig diff from best Lambda vs. 1
   Lambda1Index <- sum(LR$x < 1)
   Chi1 <- mean(c(LR$y[Lambda1Index], LR$y[Lambda1Index+1]))
   ChiLambda <- LR$y[which.max(LR$y)]
   ChiDiff <- 2 * (ChiLambda - Chi1)
   print (paste("Best Lambda=", round(LR$x[which.max(LR$y)],2)))
   print(paste("Chi-square (df=1)=", round(ChiDiff,2)))
   print(paste("p-value=", round(pchisq(ChiDiff, df=1, lower.tail=FALSE),5)))
}

