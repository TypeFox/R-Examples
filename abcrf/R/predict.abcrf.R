predict.abcrf <- function(object, obs, ntree=1000, sampsize=min(1e5, dim(sumsta)[1]), paral=FALSE, ...)
{
	modindex <- object$model.rf$y
  nmod <- nlevels(modindex)
  sumsta <- object$sumsta
	if (object$lda) nstat <- ncol(sumsta)-(nmod-1) else nstat <- ncol(sumsta)
  if (is.vector(obs)) obs <- matrix(obs,1,length(obs))
	if (is.null(colnames(obs))) {
	  colnames(obs) <- colnames(sumsta)[1:nstat]
    warning("Columns of obs have no names")
	}
	old.options <- options(); options(warn=-1)
	if (object$lda) {
	  model.lda <- lda(sumsta[,1:nstat],modindex)
		obs <- cbind(obs, predict(model.lda,obs)$x)
	}
	vote <- predict(object$model.rf, obs, type="vote", norm.votes=FALSE)
	allocation <- predict(object$model.rf, obs)
	local.error <- as.numeric(object$model.rf$predicted==modindex)
	  if (paral==TRUE) {
    ncores <- max(detectCores()-1,1) 
    cl <- makeCluster(ncores)
    registerDoParallel(cl)
    if (trunc(ntree/ncores)==ntree/ncores) ntrees <- rep(ntree/ncores, ncores) else
      ntrees <- c(rep(trunc(ntree/ncores), ncores),ntree-trunc(ntree/ncores)*ncores)
    error.rf <- foreach(ntree=ntrees, .combine= combine, .multicombine=TRUE, .packages='randomForest') %dorng% {
        randomForest(sumsta, local.error, ntree=ntree, sampsize=sampsize, ...)
      }
    stopCluster(cl)
  } else error.rf <- randomForest(sumsta, local.error, ntree=ntree, sampsize=sampsize, ...)
	options(old.options)
	tmp <- list(allocation=allocation, vote=vote, post.prob=predict(error.rf, obs))
  class(tmp) <- "abcrfpredict"
  tmp
}

summary.abcrfpredict <- function(object, ...) {
  cat("Number of affectations per model:\n")
  summary(object$allocation, ...)
}

print.abcrfpredict <- function(x, ...) {
  ret <- cbind(x$allocation, x$vote, x$post.prob)
  colnames(ret) <- c("selected model", paste("votes model",1:dim(x$vote)[2],sep=""), "post.proba")
  print(ret, ...)
}

as.matrix.abcrfpredict <- function(x, ...) {
  ret <- cbind(x$allocation, x$vote, x$post.prob)
  colnames(ret) <- c("model", "post.proba")
  ret
}

as.data.frame.abcrfpredict <- function(x, ...) {
  ret <- cbind(x$allocation, x$vote, x$post.prob)
  colnames(ret) <- c("model", "post.proba")
  as.data.frame(ret,  row.names=NULL, optional=FALSE, ...)
}

as.list.abcrfpredict <- function(x, ...) {
  list(allocation = x$allocation, vote = x$vote, post.prob = x$post.prob, ...)
}