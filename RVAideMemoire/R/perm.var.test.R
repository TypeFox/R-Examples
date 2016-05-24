perm.var.test <- function(x,...) {
  UseMethod("perm.var.test")
}

perm.var.test.formula <- function(formula,data,alternative=c("two.sided","less","greater"),nperm=999,...) {
  if (missing(formula)||(length(formula)!=3)) {stop("missing or incorrect formula")}
  m <- match.call()
  if (is.matrix(eval(m$data,parent.frame()))) {m$data <- as.data.frame(m$data)}
  m[[1]] <- as.name("model.frame")
  m$alternative <- m$nperm <- NULL
  mf <- eval(m,parent.frame())
  dname <- paste(paste(names(mf)[1],paste(names(mf)[2:ncol(mf)],collapse=":"),sep=" by "),"\n",nperm," permutations",sep="")
  resp <- mf[,1]
  fact <- interaction(mf[,2:ncol(mf)],sep=":")
  if (nlevels(fact)!=2) {stop(paste(paste(names(mf)[2:ncol(mf)],collapse=":")," is not a 2-levels factor",sep=""))}
  if (length(alternative)>1) {alternative <- "two.sided"}
  variance <- var(resp[fact==levels(fact)[1]],na.rm=TRUE)/var(resp[fact==levels(fact)[2]],na.rm=TRUE)
  ratio <- 1
  names(variance) <- names(ratio) <- "ratio of variances"
  F.ref <- var.test(resp~fact,ratio=ratio,alternative=alternative)$statistic
  F.perm <- numeric(nperm+1)
  F.perm[1] <- F.ref
  pb <- txtProgressBar(min=0,max=100,initial=0,style=3)
  for(i in 1:nperm) {
    F.perm[i+1] <- var.test(sample(resp)~fact,ratio=ratio,alternative=alternative)$statistic
    setTxtProgressBar(pb,round(i*100/nperm,0))
  }
  cat("\n")
  pvalue <- NULL
  if (alternative=="two.sided") {
    pvalue <- 2*min(length(which((F.perm-.Machine$double.eps/2) <= (F.ref)))/(nperm+1),length(which((F.perm+.Machine$double.eps/2) >= F.ref))/(nperm+1))
  }
  if (alternative=="less") {
    pvalue <- length(which((F.perm-.Machine$double.eps/2) <= (F.ref)))/(nperm+1)
    }
  if (alternative=="greater") {
    pvalue <- length(which((F.perm+.Machine$double.eps/2) >= F.ref))/(nperm+1)
  }
  result <- list(method="Permutational F test to compare two variances",statistic=F.ref,permutations=nperm,
    p.value=pvalue,estimate=variance,null.value=ratio,alternative=alternative,data.name=dname)
  class(result) <- "htest"
  return(result)
}

perm.var.test.default <- function(x,y,...) {
  if (!is.numeric(y)) {stop(paste(deparse(substitute(y)),"must be numeric"))}
  response <- c(x,y)
  fact <- factor(rep(LETTERS[1:2],c(length(x),length(y))))
  test <- perm.var.test(response~fact,...)
  test$data.name <- paste(deparse(substitute(x)),"and",deparse(substitute(y)))
  return(test)
}
