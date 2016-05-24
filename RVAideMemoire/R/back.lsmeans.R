back.lsmeans <- function(lsm,transform=c("log","logit","sqrt","inverse"),base=exp(1),
  add=0,ord=FALSE,decreasing=TRUE) {
  transform <- match.arg(transform)
  if ("list" %in% class(summary(lsm))) {
    lsm <- summary(lsm)$lsmeans
  } else {
    lsm <- summary(lsm)
  }
  col <- which(colnames(lsm)=="lsmean")
  f <- if (col>2) {
    apply(lsm[,1:(col-1)],1,function(x) paste(x,collapse=":"))
  } else {
    lsm[,1]
  }
  moy <- lsm$lsmean
  es <- lsm$SE
  res <- data.frame(a=f,SE.inf=moy-es,Mean=moy,SE.sup=moy+es)
  colnames(res)[1] <- paste(colnames(lsm[1:(col-1)]),collapse=":")
  res[,-1] <- if (transform=="log") {
    base^res[,-1]-1
  } else if (transform=="logit") {
    exp(res[,-1])/(1+exp(res[,-1]))-add
  } else if (transform=="sqrt") {
    res[,-1]^2-add
  } else if (transform=="inverse") {
    1/res[,-1]-add
  } else {
    stop("unknown transformation")
  }
  if (transform=="inverse") {
    SE1 <- res$SE.inf
    SE2 <- res$SE.sup
    res$SE.inf <- SE2
    res$SE.sup <- SE1
  }
  if (ord) {res <- res[order(res$Mean,decreasing=decreasing),]}
  return(res)
}
