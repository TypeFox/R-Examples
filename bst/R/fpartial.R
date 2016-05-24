fpartial.bst <- function (object, mstop=NULL, newdata=NULL)
{   
  if(is.null(mstop))
    mstop <- object$ctrl$mstop
  else if(mstop > object$ctrl$mstop)
    stop("mstop must be equal or smaller than the one used for estimation ", object$ctrl$mstop)
  if(object$learner=="tree" && object$maxdepth > 1)
  stop("Not implemented for higher order tree\n")
  if(object$learner=="tree")
  one <- rep(1,nrow(object$x))
  x <- object$x
  if(is.null(newdata))
    newdata <- x
  if(!missing(newdata)){
    if(object$ctrl$center){
      meanx <- drop(one %*% as.matrix(x))/nrow(x)
      newdata <- scale(newdata, meanx, FALSE) # centers x
    }
  }
  ens <- object$ens
  k <- object$k
  nu <- object$ctrl$nu
  if(missing(newdata)) p <- dim(x)[1]
  else{
    newdata <- as.matrix(newdata)
    p <- dim(newdata)[1]
  }
    lp <- matrix(object$offset, nrow=p, dim(x)[2])
  if (is.matrix(newdata)) newdata <- as.data.frame(newdata)
  for(m in 1:mstop){
      if(object$learner=="tree")
      xselect <- object$ensemble[[m]]
      else xselect <- object$ensemble[m]
      if(object$learner=="tree")
        lp[,xselect] <- lp[,xselect] + nu*predict(ens[[m]], newdata = newdata)
      else if(object$learner=="sm")
        lp[,xselect] <- lp[,xselect] + nu * predict(object$ens[[m]], newdata[, object$ensemble[m]])$y
      else if(object$learner=="ls")
        lp[,xselect] <- lp[,xselect] + nu * object$coef[m] * newdata[, object$ensemble[m]]
  }
  lp 
}

