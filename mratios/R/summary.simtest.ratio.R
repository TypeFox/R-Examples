"summary.simtest.ratio" <-
function(object, digits=4, ...)
{
cat("","\n")
cat("Numerator contrast matrix:", "\n")
print(object$Num.Contrast, digits=digits,...)

cat("","\n")
cat("Denominator contrast matrix:", "\n")
print(object$Den.Contrast, digits=digits,...)

cat("","\n")
cat("Correlation matrix under the null hypotheses:", "\n")
print(object$CorrMat, digits=digits,...)

cat("","\n")
dfcrit<-cbind("Degree of freedom"=object$df, "Critical point"=object$critical.pt)
print(dfcrit, digits=digits, ...)

cat("","\n")
if(object$alternative=="two.sided")
 {cat("Alternative hypotheses: Ratios different from margins", "\n")}
if(object$alternative!="two.sided")
 {cat("Alternative hypotheses: Ratios ",object$alternative," than margins", "\n")}

out<-cbind( object$Margin.vec, object$estimate, object$teststat, object$p.value.raw, object$p.value.adj)

rownames(out) <- object$compnames
colnames(out) <- c("margin", "estimate", "statistic", "p.value.raw", "p.value.adj")

cat("","\n")
print(round(out, digits=digits),...)
cat("","\n")
}

