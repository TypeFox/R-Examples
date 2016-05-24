chisq.bin.exp <- function(formula,data,p,graph=FALSE){
  if (missing(formula)||(length(formula)!=3)) {stop("missing or incorrect formula")}
  m <- match.call()
  if (is.matrix(eval(m$data,parent.frame()))) {m$data <- as.data.frame(m$data)}
  m[[1]] <- as.name("model.frame")
  m$p <- m$graph <- NULL
  mf <- eval(m,parent.frame())
  mf <- droplevels(mf[complete.cases(mf),])
  dname <- paste(names(mf)[1],paste(names(mf)[2:ncol(mf)],collapse=":"),sep=" by ")
  resp.mf <- mf[,1]
  resp <- factor(as.numeric(factor(resp.mf))-1)
  if (nlevels(resp)!=2) {stop(paste(names(mf)[1],"is not a binary variable"))}
  resp.num <- as.numeric(as.character(resp))
  fact <- interaction(mf[,2:ncol(mf)],sep=":")
  tab.cont <- table(fact,relevel(resp,ref="1"))
  if (length(p)!=nrow(tab.cont)){stop("number of expected probabilities and populations differ")}
  n <- integer(nrow(tab.cont))
  for (i in 1:nrow(tab.cont)) {n[i] <- sum(tab.cont[i,])}
  n.theo1 <- n*p
  n.theo2 <- n*(1-p)
  n.theo.mat <- matrix(c(n.theo1,n.theo2),nrow=nrow(tab.cont),dimnames=list(rownames(tab.cont),
    colnames(tab.cont)))
  cochran.max <- ceiling(0.8*length(tab.cont))
  cochran.min <- length(tab.cont)-cochran.max
  result <- list(p.theo=p,mat=n.theo.mat,cochran=cochran.min)
  class(result) <- "chisq.exp"
  if (graph) {mosaicplot(t(n.theo.mat),main="Expected distribution",col=TRUE)}
  return(result)
}

