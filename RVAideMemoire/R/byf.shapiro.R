byf.shapiro <- function(formula,data) {
  if (missing(formula)||(length(formula)!=3)) {stop("missing or incorrect formula")}
  m <- match.call()
  if (is.matrix(eval(m$data,parent.frame()))) {m$data <- as.data.frame(m$data)}
  m[[1]] <- as.name("model.frame")
  mf <- eval(m,parent.frame())
  dname <- paste(names(mf)[1],paste(names(mf)[2:ncol(mf)],collapse=":"),sep=" by ")
  resp <- mf[,1]
  fact <- interaction(mf[,2:ncol(mf)],sep=":")
  nlev <- nlevels(fact)
  tab <- data.frame(W=integer(nlev),"p-value"=integer(nlev),check.names=FALSE)
  rownames(tab) <- levels(fact)
  for (i in 1:nlev) {
    test <- shapiro.test(resp[as.numeric(fact)==i])
    tab[i,1] <- test$statistic
    tab[i,2] <- test$p.value
  }
  result <- list(method="Shapiro-Wilk normality tests",data.name=dname,tab=tab)
  class(result) <- "byf.test"
  return(result) 
}
