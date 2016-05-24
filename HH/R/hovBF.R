hovplotBF <- function(x, data, ..., na.rm=TRUE,
                      main="Brown-Forsyth Homogeneity of Variance",
                      plotmath=TRUE) {
  lPF <- latticeParseFormula(x, data=data)
  y <- lPF$left
  group <- lPF$right
  y.name <- lPF$left.name
  group.name <- lPF$right.name

  y.median <- tapply(y, group, median, na.rm=na.rm)
  y.minus.median <- y - y.median[group]

  ## panel.hovnew ignores user-specified groups
  panel.hovnew <- function(x, y, subscripts, groups, ...) {
    panel.abline(h=0, lty=3, col="gray20")
    panel.bwplot.superpose(x, y, groups=x,
                           subscripts=subscripts, ...)
  }

  A <- bwplot(    y               ~ group, ..., panel=panel.hovnew)
  B <- bwplot(    y.minus.median  ~ group, ..., panel=panel.hovnew)
  C <- bwplot(abs(y.minus.median) ~ group, ..., panel=panel.hovnew)

  result <- update(c(y=A,
                     "y - med(y)"=B,
                     "abs(y - med(y))"=C,
                     x.same=TRUE, layout=c(3,1)),
                   xlab=group.name,
                   ylab=y.name,
                   main=main)

  if (plotmath)
    result <-
      update(result,
             strip=strip.custom(factor.levels=c(
                                  expression(y),
                                  expression(y-y^scriptstyle(symbol("\136"))),
                                  expression("| "~y-y^scriptstyle(symbol("\136"))~" |"))),
             par.strip.text=list(lines=1.5))
  result
}



hovBF <- function(x, data=sys.parent(), ..., na.rm=TRUE) {
  lPF <- latticeParseFormula(x, data=data)
  y <- lPF$left
  group <- lPF$right
  y.name <- lPF$left.name
  group.name <- lPF$right.name

  y.median <- tapply(y, group, median, na.rm=na.rm)
  y.minus.median <- y - y.median[group]

  hov.test <- summary(aov(abs(y.minus.median) ~ group))[[1]]
  result <- list(statistic=c(F=hov.test["group","F value"]),
                 parameters=hov.test[,"Df"],
                 p.value=hov.test["group","Pr(>F)"],
                 alternative="variances are not identical",
                 method="hov: Brown-Forsyth",
                 data.name=y.name)
  names(result$parameters) <- paste("df", sep=":", c(group.name, dimnames(hov.test)[[1]][2]))
  class(result) <- "htest"
  result
}
