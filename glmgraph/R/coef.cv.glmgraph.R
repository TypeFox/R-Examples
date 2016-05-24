coef.cv.glmgraph <- function(object,s=c("lambda1.min","lambda1.1se"),...){
      s <- match.arg(s)
      if(s=="lambda1.min") return(object$beta.min)
      else if(s=="lambda1.1se") return(object$beta.1se)
      else stop("Invalid type")
}
