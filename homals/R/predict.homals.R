`predict.homals` <-
function(object, ...)
{ 
  #computes classification table and misclassification rate
 
  cl.table <- NULL
  nvec <- colnames(object$dframe)
  for (i in 1:length(object$loadings))                                #runs over variables
  {
    ax <- object$loadings[[i]]
    ay <- object$catscores[[i]]
    ag <- object$dframe[,i]                                           #observed response categories
    bx <- lsfit(t(ax),t(object$objscores),intercept=FALSE)$coef       #regress loadings on each of the object scores --> beta 
    ux <- crossprod(bx,ax)                                            #multiply beta with loading (for each dimension)
    d <- outer(rowSums(ux^2),rowSums(ay^2),"+")-2*tcrossprod(ux,ay)   #residuals for each category
    h <- levels(ag)[apply(d,1,which.min)]                             #pick minimum distance --> predicted response categories
    cl.table[[i]] <- table(ag,h,dnn=list("obs","pre"))                #cross-classify observed vs. predicted
    names(cl.table)[[i]] <- nvec[i]
  }
  cr.vec <- sapply(cl.table, function(x) (sum(diag(x)))/sum(x))       #
  result <- list(cl.table = cl.table, cr.vec = cr.vec)
  class(result) <- "predict.homals"
  result
}

