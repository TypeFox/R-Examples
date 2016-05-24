"print.bgtTest" <-
function(x,...)

{
args<-list(...)
if(is.null(args$digits))
 {digits<-4}
else{digits<-args$digits}

if (x$alternative=="two.sided"){alt.hyp="true proportion is not equal to"}
if (x$alternative=="less"){alt.hyp="true proportion is less than"}
if (x$alternative=="greater"){alt.hyp="true proportion is greater than"}

cat("\n")
cat(x$method, "test for one proportion in group testing\n")
cat("Alternative hypothesis:",alt.hyp,x$p.hyp,"\n", sep=" ")
cat("p-value","=",signif( x$p.value, digits),"\n", sep=" ")
cat("point estimate","=", signif( x$estimate, digits),"\n", sep=" ")
invisible(x)
}

"print.binTest" <-
function(x,...)

{
args<-list(...)
if(is.null(args$digits))
 {digits<-4}
else{digits<-args$digits}

if (x$alternative=="two.sided"){alt.hyp="true proportion is not equal to"}
if (x$alternative=="less"){alt.hyp="true proportion is less than"}
if (x$alternative=="greater"){alt.hyp="true proportion is greater than"}

cat("\n")
cat(x$method, "test for one binomial proportion\n")
cat("Alternative hypothesis:",alt.hyp,x$p.hyp,"\n", sep=" ")
cat("p-value","=",signif( x$p.value, digits),"\n", sep=" ")
cat("point estimate","=", signif( x$estimate, digits),"\n", sep=" ")
invisible(x)
}

