as.rts <- function(x, ...) as.ts(x, ...)  ## compatibility with S-Plus

title.grob <- function(main=NULL, y=.99, gp=gpar(cex=1.5)) {
  grid.text(main, y=y, gp=gp, just="top")
}

title.trellis <- function(main = NULL, sub = NULL, xlab = NULL, ylab = NULL,
    line = NA, outer = FALSE, axes=NULL, ...) {
  if.R(s=title(main),
       r=title.grob(main))
}

## from ~/hh/splus.library/print.arima-bug.fix.s
as.character.arima.model <- function(x, ...) {
  if (!is.null(names(x))) x <- list(x)
  for (i in seq(along=x)) {
    mi <- x[[i]]
    mic <- paste("(",paste(mi$order, collapse=","),")",sep="")
    if (!is.null(mi$period)) mic <- paste(mic, mi$period, sep="")
    if (i == 1)
      m <- mic
    else
      m <- paste(m, "x", mic, sep="")
  }
  m
}

arima.model <- function(x) {
  result <-
    if.R(s=x$model,
         r={
           arma <- x$arma
           if (is.null(arma)) list(list(order=c(0,0,0)))
           else {
             result <-  list(list(order=arma[c(1,6,2)]),
                             list(order=arma[c(3,7,4)], period=arma[5]))
             if (all(result[[2]]$order==0)) result[[2]] <- NULL
             result
           }}
         )
  class(result) <- "arima.model"
  result
}

coef.arima.HH <- function(...)
  .Defunct("coefArimaHH", package="HH")


coefArimaHH <-
if.R(r=stats:::coef.Arima,
     s=
  function(object, ...) {
  if (!is.null(names(object$model))) object$model <- list(object$model)
  a.coef <- numeric()
  for (i in seq(along=object$model)) {
    mi <- object$model[[i]]
    a.coef <- c(a.coef, mi$ar, mi$ma)
  }
  names(a.coef) <- .arima.info.names.not.ordered(.arima.S.to.C(object$model))
  a.coef
}
     )

.arima.info.names.not.ordered <-
function(model)
{
	names <- NULL
	for(i in seq(along=model$period)) {
		n.ar <- model$order[1, i]
		if(n.ar > 0)
			names <- c(names, paste("ar(", (1:n.ar) * model$period[[
				i]], ")", sep = ""))
		n.ma <- model$order[3, i]
		if(n.ma > 0)
			names <- c(names, paste("ma(", (1:n.ma) * model$period[[
				i]], ")", sep = ""))
	}
	names
}
