
summary.FuzzyStatObj<-function(object,...){
  cat(object$summary);
}

print.FuzzyStatObj<-function(x,...){
  cat(x$summary);
}