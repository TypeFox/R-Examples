"print.simtest.ratio" <-
function(x, digits=4, ...)
{

cat(x$methodname)

if(x$alternative=="two.sided")
 {cat("Alternative hypotheses: Ratios different from margins", "\n")}
if(x$alternative!="two.sided")
 {cat("Alternative hypotheses: Ratios ",x$alternative," than margins", "\n", sep="")}
 cat("Degree of freedom: ",paste(signif(x$df, digits), collapse=", "),"\n", sep="")
out<-cbind( x$Margin.vec, x$estimate, x$teststat, x$p.value.raw, x$p.value.adj)

rownames(out) <- x$compnames
colnames(out) <- c("margin", "estimate", "statistic", "p.value.raw", "p.value.adj")

cat("","\n")
print(out, digits=digits, ...)
cat("","\n")

invisible(x)

}

