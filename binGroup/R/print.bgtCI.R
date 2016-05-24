"print.bgtCI" <-
function(x, ...)
{
args<-list(...)
if(is.null(args$digits))
 {digits<-4}
else{digits<-args$digits}
cat("\n")
cat(x$conf.level*100,"percent", x$method, "confidence interval:\n", sep=" ")
cat(" [",paste(signif(x$conf.int, digits), collapse=", "),"]\n", sep=" " )
cat("Point estimate:", signif(x$estimate, digits), "\n", sep=" ")
invisible(x)
}

"print.binCI" <-
function(x, ...)
{
args<-list(...)
if(is.null(args$digits))
 {digits<-4}
else{digits<-args$digits}
cat("\n")
cat(x$conf.level*100,"percent", x$method, "confidence interval\n", sep=" ")
cat(" [",paste(signif(x$conf.int, digits), collapse=", "),"]\n", sep=" " )
cat("Point estimate", signif(x$estimate, digits), "\n", sep=" ")
invisible(x)
}

