summary.lmp<-function(object,...){
a<-object
nvar<-a$rank
N<-nrow(a$model)
rdf<-ifelse(a$knp==FALSE,N-nvar-1,N-nvar)
r<-a$residuals
Y<-a$fitted.values
p<-a$p
iter<-a$iter
mpy<-paramp(Y,p=p)$mp
rss<-sum(abs(r)^p)
resvar<-rss/rdf
sigma<-(resvar)^(1/p)
ans <- object[c("call", "terms")]
    ans$p<-p
    ans$residuals <- r
    ans$coefficients <- a$coef
    ans$sigma <- sigma
    ans$rdf <- rdf
    ans$iter <- iter
    class(ans) <- "summary.lmp"
    ans
}
  
print.summary.lmp<-function(x,...){
a<-x
cat("\nCall:\n")
cat(paste(deparse(a$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
cat("Residuals:\n")
rq <- quantile(a$residuals)
names(rq) <- c("Min", "1Q", "Median", "3Q", "Max")
print(rq,4)
cat("\nCoefficients:\n")
print(a$coef,4)
cat("\nEstimate of p\n")
print(as.symbol(a$p))
cat("\nPower deviation of order p:",
    format(signif(a$sigma,4)))
cat("\n")
if(a$iter==1) cat("\nWarning: problems in convergence in lmp\n")
invisible(a)
}

