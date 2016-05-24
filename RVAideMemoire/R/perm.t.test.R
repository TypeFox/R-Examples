perm.t.test <- function(x,...) {
  UseMethod("perm.t.test")
}

perm.t.test.formula <- function(formula,data,alternative=c("two.sided","less","greater"),paired=FALSE,nperm=999,...) {
  if (missing(formula)||(length(formula)!=3)) {stop("missing or incorrect formula")}
  m <- match.call()
  if (is.matrix(eval(m$data,parent.frame()))) {m$data <- as.data.frame(m$data)}
  m[[1]] <- as.name("model.frame")
  m$alternative <- m$paired <- m$nperm <- NULL
  mf <- eval(m,parent.frame())
  dname <- paste(paste(names(mf)[1],paste(names(mf)[2:ncol(mf)],collapse=":"),sep=" by "),"\n",nperm," permutations",sep="")
  resp <- mf[,1]
  fact <- interaction(mf[,2:ncol(mf)],sep=":")
  if (nlevels(fact)!=2) {stop(paste(paste(names(mf)[2:ncol(mf)],collapse=":")," is not a 2-levels factor",sep=""))}
  if (paired & diff(tapply(resp,fact,length))!=0) {stop(paste("'",levels(fact)[1],"' and '",levels(fact)[2],"' lengths differ",sep=""))}
  if (length(alternative)>1) {alternative <- "two.sided"}
  method <- NULL
  moy <- NULL
  null.value <- 0
  names(null.value) <- "difference in means"
  t.ref <- t.test(resp~fact,var.equal=TRUE,alternative=alternative,paired=paired)$statistic
  t.perm <- numeric(nperm+1)
  t.perm[1] <- t.ref
  pb <- txtProgressBar(min=0,max=100,initial=0,style=3)
  if (!paired) {
    method <- "Permutational Two Sample t-test"
    moy <- tapply(resp,fact,mean)
    names(moy) <- paste("mean in group ",levels(fact),sep="")
    for(i in 1:nperm) {
	t.perm[i+1] <- t.test(sample(resp)~fact,var.equal=TRUE,alternative=alternative,paired=FALSE)$statistic
	setTxtProgressBar(pb,round(i*100/nperm,0))
    }
  } else {
    method <- "Permutational Paired t-test"
    moy <- mean(resp[fact==levels(fact)[1]]-resp[fact==levels(fact)[2]])
    names(moy) <- "mean of the differences"
    resp2 <- cbind(resp[fact==levels(fact)[1]],resp[fact==levels(fact)[2]])
    for (i in 1:nperm) {
	resp.perm <- t(apply(resp2,1,sample))
	t.perm[i+1] <- t.test(resp.perm[,1],resp.perm[,2],alternative=alternative,paired=TRUE)$statistic
	setTxtProgressBar(pb,round(i*100/nperm,0))
    }
  }
  cat("\n")
  pvalue <- NULL
  if (alternative=="two.sided") {
    pvalue <- 2*min(length(which((t.perm-.Machine$double.eps/2) <= t.ref))/(nperm+1),length(which((t.perm+.Machine$double.eps/2) >= t.ref))/(nperm+1))
  }
  if (alternative=="less") {
    pvalue <- length(which((t.perm-.Machine$double.eps/2) <= t.ref))/(nperm+1)
    }
  if (alternative=="greater") {
    pvalue <- length(which((t.perm+.Machine$double.eps/2) >= t.ref))/(nperm+1)
  }
  result <- list(statistic=t.ref,permutations=nperm,p.value=pvalue,estimate=moy,alternative=alternative,
    method=method,data.name=dname,null.value=null.value)
  class(result) <- "htest"
  return(result)
}

perm.t.test.default <- function(x,y,paired=FALSE,...) {
  if (!is.numeric(y)) {stop(paste(deparse(substitute(y)),"must be numeric"))}
  response <- c(x,y)
  fact <- factor(rep(LETTERS[1:2],c(length(x),length(y))))
  test <- perm.t.test(response~fact,paired=paired,...)
  test$data.name <- paste(deparse(substitute(x)),"and",deparse(substitute(y)))
  if (!paired) {
    names(test$estimate) <- c(paste("mean of",deparse(substitute(x))),paste("mean of",deparse(substitute(y))))
  }
  return(test)
}
