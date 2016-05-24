sort.cpg <-
function(x,decreasing,...) {
  sorted.p<-order(x$results$P.value)
  x$results<-x$results[sorted.p,]
  if(!is.factor(x$indep)) {
    x$coefficients<-x$coefficients[sorted.p,]
    }
  x$FDR.sig<-x$FDR.sig[order(x$FDR.sig$P.value),] 
  x$Holm.sig<-x$Holm.sig[order(x$Holm.sig$P.value),]
  x
  }
