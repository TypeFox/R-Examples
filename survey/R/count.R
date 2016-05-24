unwtd.count<-function(x, design,...){

  if (inherits(x, "formula")) {
    mf <- model.frame(x, model.frame(design), na.action = na.pass)
    xx <- lapply(attr(terms(x), "variables")[-1], 
                 function(tt) model.matrix(eval(bquote(~0 + .(tt))), mf)
                 )
    cols <- sapply(xx, NCOL)
    x <- matrix(nrow = NROW(xx[[1]]), ncol = sum(cols))
    scols <- c(0, cumsum(cols))
    for (i in 1:length(xx)) {
      x[, scols[i] + 1:cols[i]] <- xx[[i]]
    }
    colnames(x) <- do.call("c", lapply(xx, colnames))
  }
  else if (typeof(x) %in% c("expression", "symbol")) 
    x <- eval(x, model.frame(design))
  x <- as.matrix(x)
  out<- weights(design,"sampling")==0
  nas <- rowSums(is.na(x))

  x <- x[(nas+out) == 0, , drop = FALSE]
 
  rval<-NROW(x)
  names(rval)<-"counts"
  attr(rval,"var")<-matrix(0,1,1)
  attr(rval,"statistic")<-"counts"
  if (inherits(design,"svyrep.design"))
    class(rval)<-"svrepstat"
  else
    class(rval)<-"svystat"
  rval
  
}
