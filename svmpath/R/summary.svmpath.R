summary.svmpath<-function(object,nsteps=5,digits=6,...){
  cat("\nCall:\n")
  cat(paste(deparse(object$call), sep = "\n", collapse = "\n"), 
        "\n\n", sep = "")
  N<-length(object$Error)
  nsupp<-apply(object$alpha,2,function(x)sum(x>1e-10))
  m<-sort(unique(c(pretty(seq(N),nsteps),1,N)))
  m<-seq(N)[match(m,seq(N),0)]
  cat("Number of steps:",N,"\n")
  cat("Selected steps:\n")
  data.frame(Lambda=round(object$lambda[m],digits),Error=object$Error[m],Size.Elbow=object$Size.Elbow[m],Support=nsupp[m],SumEps=object$SumEps[m],row.names=m)
  

}
