# nnet : multinom

prop.multinom.test <- function(x,p.method="fdr") {
  dname <- deparse(substitute(x))
  if (is.data.frame(x)) {x <- as.matrix(x)}
  if (is.matrix(x)) {
    if (!is.numeric(x)) {stop("incorrect 'x' format")}
    if (is.null(colnames(x))) {
	colnames(x) <- LETTERS[1:ncol(x)]
    }
    lab <- colnames(x)
  } else if (is.factor(x) | is.character(x)) {
    x <- as.factor(x)
    lab <- levels(x)
  } else {stop("incorrect 'x' format")}
  pval <- matrix(0,nrow=length(lab),ncol=length(lab),dimnames=list(lab,lab))
  z.tab <- matrix(0,nrow=length(lab),ncol=length(lab),dimnames=list(lab,lab))
  if (is.matrix(x)) {
    for (i in 1:(length(lab)-1)) {
	x.temp <- rbind(x[,c(i,0:(i-1),(i+1):ncol(x))])
	mod.temp <- nnet::multinom(x.temp~1,trace=FALSE)
	z <- as.vector(summary(mod.temp)$coefficients/summary(mod.temp)$standard.errors)
	names(z) <- lab[-i]
	z.tab[names(z),lab[i]] <- z
	p <- sapply(z,function(y) {2*min(pnorm(y),pnorm(y,lower.tail=FALSE))})
	pval[names(z),lab[i]] <- p
    }
  } else {
    for (i in 1:(length(lab)-1)) {
	x.temp <- relevel(x,ref=lab[i])
	mod.temp <- nnet::multinom(x.temp~1,trace=FALSE)
	z <- as.vector(summary(mod.temp)$coefficients/summary(mod.temp)$standard.errors)
	names(z) <- lab[-i]
	z.tab[names(z),lab[i]] <- z
	p <- sapply(z,function(y) {2*min(pnorm(y),pnorm(y,lower.tail=FALSE))})
	pval[names(z),lab[i]] <- p
    }
  }
  pval[upper.tri(pval,diag=TRUE)] <- NA
  pval[lower.tri(pval)] <- p.adjust(pval[lower.tri(pval)],method=p.method)
  pval <- pval[-1,-ncol(pval)]
  if (length(pval)==1) {names(pval) <- paste(lab,collapse="-")}
  z.tab[upper.tri(z.tab,diag=TRUE)] <- NA
  z.tab <- z.tab[-1,-ncol(z.tab)]
  if (length(z.tab)==1) {names(z.tab) <- paste(lab,collapse="-")}
  res <- list(method="Wald tests",data.name=dname,p.adjust.method=p.method,p.value=pval,
    z.tab=z.tab)
  class(res) <- "RV.multcomp"
  return(res)
}
