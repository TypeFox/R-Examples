prop.bin.multcomp <-
function(formula,data,p,p.method="fdr") {
  if (missing(formula)||(length(formula)!=3)) {stop("missing or incorrect formula")}
  m <- match.call()
  if (is.matrix(eval(m$data,parent.frame()))) {m$data <- as.data.frame(m$data)}
  m[[1]] <- as.name("model.frame")
  m$p <- m$p.method <- NULL
  mf <- eval(m,parent.frame())
  mf <- droplevels(mf[complete.cases(mf),])
  dname <- paste(names(mf)[1],paste(names(mf)[2:ncol(mf)],collapse=":"),sep=" by ")
  resp.mf <- mf[,1]
  resp <- factor(as.numeric(factor(resp.mf))-1)
  if (nlevels(resp)!=2) {stop(paste(names(mf)[1],"is not a binary variable"))}
  resp.num <- as.numeric(as.character(resp))
  fact <- interaction(mf[,2:ncol(mf)],sep=":")
  tab.cont <- table(fact,relevel(resp,ref="1"))
  obs <- integer(nrow(tab.cont))
  pval <- integer(nrow(tab.cont))
  for (i in 1:nrow(tab.cont)) {
    obs[i] <- tab.cont[i,1]/sum(tab.cont[i,])
    test <- binom.test(tab.cont[i,1],sum(tab.cont[i,]),p[i])
    pval[i] <- test$p.value
  }
  p.adj <- p.adjust(pval,method=p.method)
  comp <- data.frame("observed"=obs,"expected"=p,"p-value"=p.adj," "=.psignif(p.adj),stringsAsFactors=FALSE,check.names=FALSE)
  if (!is.null(rownames(tab.cont))) {rownames(comp) <- rownames(tab.cont)}
  result <- list(method="exact binomial tests",data.name=dname,observed=obs,expected=p,p.adjust.method=p.method,p.value2=p.adj,p.value=comp)
  class(result) <- "RV.multcomp"
  return(result)
}
