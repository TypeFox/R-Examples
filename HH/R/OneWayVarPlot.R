OneWayVarPlot <- function(x, data, ...,
                          main="Variability of Groups, Centers of Groups, and all Data",
                          centerFunctionName="median",
                          center=TRUE) {
  centerFunction <- get(centerFunctionName)
  lPF <- latticeParseFormula(x, data=data)
  y <- lPF$left
  group <- lPF$right
  y.name <- lPF$left.name
  group.name <- lPF$right.name

  y.center <- tapply(y, group, centerFunction)
  y.minus.center <- y - y.center[group]

  environment(x) <- environment()
  Ones <- rep("All\nData", length(y))

  panel.hovnew <- function(x, y, subscripts, groups, ...) {
    panel.abline(h=0, lty=3, col="gray20")
    panel.bwplot.superpose(x, y, groups=x,
                           subscripts=subscripts, ...)
  }

  if (!center) {
    A <- bwplot(x, data, ..., xlab=group.name, main=main,
                panel=panel.hovnew, groups=group)
    B <- bwplot(y.center ~ rep("Centers of\nGroups", length(y.center)),
                panel=panel.hovnew, groups=rep(1, length(y.center)), col="black")

    C <- bwplot(y ~ Ones, panel=panel.hovnew, groups=Ones, col="black")
  }
  else {
    A <- bwplot(y.minus.center ~ group, data, ..., xlab=group.name, main=main,
                panel=panel.hovnew, groups=group, ylab=paste("Centered", y.name))
    B <- bwplot((y.center - centerFunction(y.center)) ~ rep("Centers of\nGroups", length(y.center)),
                panel=panel.hovnew, groups=rep(1, length(y.center)), col="black")

    C <- bwplot(y-centerFunction(y) ~ Ones, panel=panel.hovnew, groups=Ones, col="black")
  }

  update(resizePanels(c(A, B, C, layout=c(3,1), y.same=TRUE), w=c(length(y.center), 1, 1)),
         between=list(x=1))
}
