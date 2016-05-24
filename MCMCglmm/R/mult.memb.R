mult.memb<-function(formula = NULL){

  if(class(formula)!="formula"){stop("formula not passed to mult.memb")}

  fac<-attr(terms(formula), "term.labels")
  Z<-model.matrix(as.formula(paste("~", fac[1], -1), env=attr(formula, ".Environment")))

  if(length(fac)>1){
    for(i in 2:length(fac)){
      Z<-Z+model.matrix(as.formula(paste("~", fac[i], -1), env=attr(formula, ".Environment")))
    }
  }
  if(any(apply(Z, 2, function(x){all(x==0)}))){
    Z<-Z[,-which(apply(Z, 2, function(x){all(x==0)})),drop=FALSE]
  }
  return(Z)
}

