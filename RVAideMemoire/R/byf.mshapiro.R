byf.mshapiro <- function(formula,data) {
  if (missing(formula)||(length(formula)!=3)) {stop("missing or incorrect formula")}
  m <- match.call()
  m[[1]] <- as.name("model.frame")
  mf <- eval(m,parent.frame())
  for (i in 1:ncol(mf)) {
    if (all(is.na(suppressWarnings(as.numeric(as.character(mf[,i])))))) {
	fact1 <- i
	break
    }
  }
  resp <- mf[,1:(fact1-1)]
  fact <- interaction(mf[,fact1:ncol(mf)],sep=":")
  dname <- paste(names(mf)[1]," by ",paste(names(mf)[fact1:ncol(mf)],collapse=":"),sep="")
  nlev <- nlevels(fact)
  tab <- data.frame(W=integer(nlev),"p-value"=integer(nlev),check.names=FALSE)
  rownames(tab) <- levels(fact)
  for (i in 1:nlev) {
    test <- mshapiro.test(resp[as.numeric(fact)==i,])
    tab[i,1] <- test$statistic
    tab[i,2] <- test$p.value
  }
  result <- list(method="Multivariate Shapiro-Wilk normality tests",data.name=dname,tab=tab)
  class(result) <- "byf.test"
  return(result) 
}
