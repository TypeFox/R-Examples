print.prp.index<-function(x,digits= max(3, getOption("digits") - 3), ...){
cat("Perpendicularity index","\n\n")
cs<-x$crossing.summary
cat("Binomial estimates for crossing probability","\n")
if(!is.na(cs[1])){
rq<-structure(x$crossing.summary)
print(rq)
cat("\n")
}
ps<-x$P.summary
rq<-structure(x$P.summary)
cat("Perpendicularity summary","\n")
print(rq)
cat("\n")
invisible(x)
}
