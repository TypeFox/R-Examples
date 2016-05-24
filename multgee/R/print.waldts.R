print.waldts <-
function (x, ...) 
{
 cat("Goodness of Fit based on the Wald test", "\n")
 cat("\nModel under H_0: ")
 print(x$NullModel)
 cat("Model under H_1: ")
 print(x$AlternativeModel)
 cat("\nWald Statistic=",x$waldstatistic,", df=",x$df,", p-value=",ifelse(x$pvalue<=0.0001,"<0.0001",x$pvalue),"\n",sep="")
}