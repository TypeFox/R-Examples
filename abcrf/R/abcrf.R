abcrf.factor <- function(modindex, sumsta, lda=TRUE, ntree=500, sampsize=min(1e5, length(modindex)), paral=FALSE, ...) 
{
	if (!is.factor(modindex)) 
    stop("modindex should be a factor containing the model indexes")
  sumsta <- as.matrix(sumsta)
	if (length(modindex)!=nrow(sumsta)) 
    stop("length of modindex differs from the number of lines in sumsta")
  if (length(modindex) == 0L)
    stop("no simulation in the reference table (modindex, sumstat)")
  if ( (!is.logical(lda)) && (length(lda) != 1L))
    stop("lda should be TRUE or FALSE")
  if (is.null(colnames(sumsta))) colnames(sumsta) <- paste("V",1:dim(sumsta)[2],sep="")
  if (lda) {
  	model.lda <- lda(sumsta, modindex)
    sumsta <- cbind(sumsta, predict(model.lda)$x)
  } else { 
    model.lda <- NULL
  }
  m <- names(match.call(expand.dots=TRUE))
  if ((!"sampsize" %in% m) && (nrow(sumsta) <= 15)) 
    sampsize <- as.integer(sampsize / 10)
  if (paral==TRUE) {
    ncores <- max(detectCores()-1,1) 
    cl <- makeCluster(ncores)
    registerDoParallel(cl)
    if (trunc(ntree/ncores)==ntree/ncores) ntrees <- rep(ntree/ncores, ncores) else
      ntrees <- c(rep(trunc(ntree/ncores), ncores),ntree-trunc(ntree/ncores)*ncores)
    model.rf <- foreach(ntree=ntrees, .combine= combine, .multicombine=TRUE, .packages='randomForest') %dorng% {
        randomForest(sumsta, modindex, ntree=ntree, sampsize=sampsize, norm.votes=FALSE, keep.inbag=TRUE, ...)
      }
    stopCluster(cl)
    tbl <- table(model.rf$predicted, model.rf$y)
    class.error = vector()

    for (i in 1:nrow(tbl)) {
      rowSum <- sum(tbl[i,])
      accurate <- diag(tbl)[i]
      error <- rowSum - accurate
      class.error[i] <- error / rowSum
    }   
    model.rf$confusion <- cbind(tbl, class.error)
    colnames(model.rf$confusion) <- c(paste(1:nlevels(modindex)),"class.error")
  } else model.rf <- randomForest(sumsta, modindex, ntree=ntree, sampsize=sampsize, norm.votes=FALSE, ...)
  cl <- match.call()
  cl[[1]] <- as.name("abcrf")
  x <- list(call=cl, lda=lda, model.rf=model.rf, model.lda=model.lda, prior.err=mean(model.rf$predicted!=model.rf$y), sumsta=sumsta)
  class(x) <- "abcrf"
  x
}

abcrf.formula <- function(formula, data=NULL, ...) {
  if (!inherits(formula, "formula"))
    stop("abcrf.formula is only for formula objects")
  m <- match.call(expand.dots=FALSE)
  names(m)[2] <- "formula"
  m[[1]] <- as.name("model.frame")
  m <- eval(m, parent.frame())
  modindex <- model.response(m)
  Terms <- attr(m, "terms")
  attr(Terms, "intercept") <- 0
  sumsta <- model.frame(terms(reformulate(attributes(Terms)$term.labels)), 
                         data.frame(m))
  ret <- abcrf(modindex, sumsta, ...)
  cl <- match.call()
  cl[[1]] <- as.name("abcrf")
  ret$call <- cl
  class(ret) <- c("abcrf.formula", "abcrf")
  return(ret)
}

abcrf.default <- function(...) {
  cl <- match.call()
  cl[[1]] <- as.name("abcrf")
  cat("call:\n")
  print(cl)
  stop("the first argument should be a formula or a factor")
}
  
abcrf <- function(...) UseMethod("abcrf")

print.abcrf <- function(x, ...) {
  cat("\nCall:\n", deparse(x$call), "\n")
  if (x$lda) cat("includes the axes of a preliminary LDA\n\n")  
  cat("Number of simulations: ", length(x$model.rf$y), "\n", sep="")
  cat("     Prior error rate: ", round(x$prior.err * 100, digits = 4), "%\n\n", sep = "")
  cat("Confusion matrix:\n")
  print(x$model.rf$confusion, ...)  
}
