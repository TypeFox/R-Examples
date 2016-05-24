LDA.format <- function (model) {
  if (!is.null(Terms<-model$terms)) {
    data <- model.frame(model)
    X <- model.matrix(delete.response(Terms),data)
    g <- model.response(data)
    xint <- match("(Intercept)",colnames(X),nomatch=0L)
    if (xint>0) {X<-X[,-xint,drop=FALSE]}
  } else {
    xname <- model$call$x
    gname <- model$call[[3L]]
    if (!is.null(sub <- model$call$subset)) {
	X <- eval.parent(parse(text=paste(deparse(xname,backtick=TRUE),
	  "[",deparse(sub,backtick=TRUE),",]")))
	g <- eval.parent(parse(text=paste(deparse(gname,backtick=TRUE),
	  "[",deparse(sub,backtick=TRUE),"]")))
    } else {
	X <- eval.parent(xname)
	g <- eval.parent(gname)
    }
    if (!is.null(nas<-model$call$na.action)) {
	df <- data.frame(g=g,X=X)
	df <- eval(call(nas,df))
	g <- df$g
	X <- df$X
    }
  }
  g.lda <- g
  means <- colMeans(model$means)
  X2 <- scale(X,center=means,scale=FALSE) %*% model$scaling
  corr <- cor(X,X2,use="pairwise")
  result <- list(x=X,grouping=g,li=X2,co=corr)
  return(result)
}
