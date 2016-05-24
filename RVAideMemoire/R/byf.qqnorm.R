byf.qqnorm <- function(formula,data,...) {
  if (missing(formula)||(length(formula)!=3)) {stop("missing or incorrect formula")}
  m <- match.call(expand.dots=FALSE)
  if (is.matrix(eval(m$data,parent.frame()))) {m$data <- as.data.frame(m$data)}
  m[[1]] <- as.name("model.frame")
  m$... <- NULL
  mf <- eval(m,parent.frame())
  resp <- mf[,1]
  fact <- interaction(mf[,2:ncol(mf)],sep=":")
  opar <- par(no.readonly=TRUE)
  on.exit(par(opar))
  par(mfrow=n2mfrow(nlevels(fact)))
  for (i in 1:nlevels(fact)) {
    qqPlot(resp[as.numeric(fact)==i],main=levels(fact)[i],lwd=1,grid=FALSE,xlab="Theoretical quantiles",
	ylab="Sample quantiles",...)
  }
}
