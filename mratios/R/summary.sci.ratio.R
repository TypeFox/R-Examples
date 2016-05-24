"summary.sci.ratio" <-
function(object, digits=4,...)
{
cat( "\n")
cat("Numerator contrast matrix:        ","\n")
print(object$Num.Contrast, digits=digits, ...)
cat( "\n")
cat("Denominator contrast matrix:        ","\n")
print(object$Den.Contrast, digits=digits, ...)

if(object$method=="Plug")
{
 cat("\n")
 cat("Estimated correlation matrix used for calculation of quantiles:","\n")
 print(object$CorrMat.est, digits=digits, ...)
}
cat("\n")

dfcrit<-cbind("Degree of freedom"=object$df, "Critical point"=object$quantile)
print(dfcrit, digits=digits, ...)

cat("\n")

if(object$method=="Unadj")

{

if(object$alternative=="two.sided")
 {
  cat("Two-sided",signif(object$conf.level*100, digits), "%", "unadjusted confidence intervals for ratios:","\n")
 }
if(object$alternative=="less")
 {
  cat("Upper",signif(object$conf.level*100, digits), "%", "unadjusted confidence limits for ratios:","\n")
 }

if(object$alternative=="greater")
 {
  cat("Lower",signif(object$conf.level*100, digits), "%", "unadjusted confidence limits for ratios:","\n")
 }

}

else

{
if(object$alternative=="two.sided")
 {
  cat("Two-sided",signif(object$conf.level*100, digits), "%", "simultaneous confidence intervals for ratios:","\n")
 }
if(object$alternative=="less")
 {
  cat("Upper",signif(object$conf.level*100, digits), "%", "simultaneous confidence limits for ratios:","\n")
 }

if(object$alternative=="greater")
 {
  cat("Lower",signif(object$conf.level*100, digits), "%", "simultaneous confidence limits for ratios:","\n")
 }
}

if(object$NSD)
{
 print(cbind(estimate=object$estimate, object$conf.int), digits=digits, ... )

 cat("                                      ","\n")
 cat("   NSD = The mean in the denominator is not significantly different from zero. ","\n")
 cat("                                      ","\n") 
}
else
{
 print( round(cbind(estimate=object$estimate, object$conf.int) ,digits=digits) ,... )
}

}

