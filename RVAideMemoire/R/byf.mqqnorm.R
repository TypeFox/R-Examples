byf.mqqnorm <- function(formula,data) {
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
  nlev <- nlevels(fact)
  opar <- par(no.readonly=TRUE)
  on.exit(par(opar))
  par(mfrow=n2mfrow(nlev))
  for (i in 1:nlev) {
    mqqnorm(resp[as.numeric(fact)==i,],main=levels(fact)[i])
  }
}
