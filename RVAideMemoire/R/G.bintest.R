G.bintest <-
function(formula,data,alpha=0.05,p.method="fdr") {
  if (missing(formula)||(length(formula)!=3)) {stop("missing or incorrect formula")}
  m <- match.call()
  if (is.matrix(eval(m$data,parent.frame()))) {m$data <- as.data.frame(m$data)}
  m[[1]] <- as.name("model.frame")
  m$alpha <- m$p.method <- NULL
  mf <- eval(m,parent.frame())
  mf <- droplevels(mf[complete.cases(mf),])
  dname <- paste(names(mf)[1],paste(names(mf)[2:ncol(mf)],collapse=":"),sep=" by ")
  resp.mf <- mf[,1]
  resp <- factor(as.numeric(factor(resp.mf))-1)
  if (nlevels(resp)!=2) {stop(paste(names(mf)[1],"is not a binary variable"))}
  resp.num <- as.numeric(as.character(resp))
  fact <- interaction(mf[,2:ncol(mf)],sep=":")
  proba <- tapply(resp.num,fact,mean)
  names(proba) <- paste("proba in group ",levels(fact),sep="")
  tab.cont <- table(fact,relevel(resp,ref="1"))
  nval <- 0
  names(nval) <- "difference in probabilities"
  result <- list(data.name=dname,alternative="two.sided",null.value=nval,estimate=proba,alpha=alpha)
  test <- G.test(tab.cont)
  result$statistic <- test$statistic
  result$parameter <- test$parameter
  result$p.value <- test$p.value
  result$method.test <- "G-test"
  if (test$p.value<alpha) {
    result$p.adjust.method <- p.method
    result$p.value.multcomp <- pairwise.G.test(tab.cont,p.method=p.method)$p.value
    result$method.multcomp <- "G-tests"
  }
  class(result) <- "RVtest"
  return(result)
}
