diagplot5new <- function(linearmodel, ..., pch=19) {
  yhat <- predict(linearmodel)
  e <- resid(linearmodel)

  tmp <- na.exclude(data.frame(yhat, e))

  n <- nrow(tmp)
  f <- (1:n)/n

  xyplot(sort(yhat - mean(yhat)) + sort(e) ~ f, data=tmp,
         outer=TRUE,
         scales=list(alternating=FALSE),
         between=list(x=1),
         strip=FALSE,
         pch=pch,
         layout=c(2,1),
         ...,
         xlab="(1:n) / n",
         ## ylab=deparse(linearmodel$call[[2]][[2]]), ## y variable name
         ylab=expression("Centered Fit = " * hat(y) - bar(y)),
         ylab.right=expression("Residuals = " * y - hat(y)),
         xlab.top=
         c(expression(sort(hat(y) - bar(y))),
           expression(sort(y - hat(y)))))
}


residVSfitted <- function(linearmodel, groups=(e >= 0), ...) {
  yhat <- predict(linearmodel)
  e <- resid(linearmodel)
  xyplot(e ~ yhat, groups=groups,
         scales=list(y=list(at=-100, alternating=3)),
         panel=function(...) {
           panel.abline(h=0, lty=3, col="gray50")
           panel.xyplot(...)
           panel.axis(side='left', outside=TRUE)
           panel.axis(side='right',
                      at=-100,
                      outside=TRUE)
         },
         xlab="Fitted Values",
         xlab.top="Residuals vs Fitted",
         ylab="Residuals",
         ylab.right=" ",
         par.settings=list(clip=list(panel=FALSE)),
         ...)
}


scaleLocation <- function(linearmodel, groups=(e >= 0), ...) {
  yhat <- predict(linearmodel)
  e <- resid(linearmodel)
  ae <- abs(e)
  sae <- sqrt(ae)
  ylim <- range(0, sae)
  ylim <- ylim + diff(ylim)*c(-1, 1)*.04
  pretty.sae <- c(0, pretty(sae))
  pretty.ae <- c(0, pretty(ae))
  xyplot(sae ~ yhat, groups=groups,
         ylim=ylim,
         scales=list(y=list(at=-100)),
         panel=function(...) {
           panel.abline(h=0, lty=3, col="gray50")
           panel.xyplot(...)
           panel.axis(side='left', outside=TRUE)
           panel.axis(side='right',
                      at=sqrt(pretty.ae), labels=pretty.ae,
                      outside=TRUE)
         },
         par.settings=list(clip=list(panel=FALSE)),
         xlab="Fitted Values",
         xlab.top="Scale-Location",
         ylab=expression(sqrt(""~abs(""~Residuals~"")~"")),
         ylab.right=expression(abs(""~Residuals~"")),
         ...)
}


diagQQ <- function(lm.object, ...) {
  qqmath(resid(lm.object),
         xlab="Theoretical Quantiles",
         xlab.top="Normal Q-Q",
         ylab="Residuals",
         ylab.right="",
         pch=19,
         panel=function(..., col=col) {
           panel.qqmath( ..., col=col)
           panel.qqmathline(..., lty=3, col="gray50")
         },
         ...
         )
}
