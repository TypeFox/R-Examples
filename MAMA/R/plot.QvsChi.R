#Porovnanie kvantilov Q-statistiky a Chi-square rozdelenia
plotQvsChi<-function(Q,num.studies)
{
 chisqq <- qchisq(seq(0, 0.9999, 0.001), df = num.studies-1)
 tmp <- quantile(Q, seq(0, 0.9999, 0.001))
 qqplot(chisqq, tmp, ylab = "Quantiles of Sample", pch = "*", 
 xlab = "Quantiles of Chi > square", main = "QQ Plot")
 lines(chisqq, chisqq, lty = "dotted", col = "red")
}
