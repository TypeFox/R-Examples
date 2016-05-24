dstrplotContinuous <- function(dstr, n=101, ok0=FALSE, ok1=FALSE, args=list(), ..., X=NULL)  {
  qdstr <- paste("q", dstr, sep="")
  ddstr <- paste("d", dstr, sep="")
  pp <- ppoints(n)
  pp01 <- pp[c(1, length(pp))]
  if (ok0) {
    pp01[1] <- 0
    pp <- c(0, pp)
  }
  if (ok1) {
    pp01[2] <- 1
    pp <- c(pp, 1)
  }
  qq01 <- do.call(qdstr, c(list(pp01), args))
  ## qq <- do.call(qdstr, c(list(pp), args)) ## evenly spaced in p scale
  qq <- seq(qq01[1], qq01[2], length=n) ## evenly spaced in q scale

  xlim <- list(...)$xlim
  if (!is.null(xlim)) qq <- seq(xlim[1], xlim[2], length=n) ## evenly spaced in q scale in xlim range

  if (!is.null(X)) {
    qq.left.index <- (qq <= X)
    qq <- c(qq[qq.left.index], X, qq[!qq.left.index])
  }

  plotfun <- function(x) do.call(ddstr, c(list(x), args))

  dd <- deparse(args)
  ssdd <- substring(strsplit(dd, "), .Names")[[1]][1], 16)

  yy <- plotfun(qq)

  xyplot(yy ~ qq, type="l", ...,
         main=paste(ddstr, "(x", if (length(args) > 0) ", " , ssdd, ")", sep=""),
         xlab="x", ylab=NULL,
         scales=list(tck=c(1, .3)),
         panel=function(...) {
           panel.abline(a=0, b=0, col="gray70")
           panel.xyplot(...)
         })
}

panel.dstrplotContinuousFill <-
  function(x, y, ..., X=5, col.left="skyblue3", col.right="skyblue1") {
    x.left.index <- (x <= X)
    x.left.lastindex <- sum(x.left.index)
    x.right.index <- !x.left.index
    panel.polygon(x=c(x[c(x.left.lastindex, x.left.lastindex)], x[!x.left.index],  x[c(length(x), x.left.lastindex)]),
                  y=c(0, y[x.left.lastindex],                   y[!x.left.index],  0, 0),
                  ...,
                  col=col.right, border="transparent")
    panel.polygon(x=c(x[1], x[x.left.index], x[c(x.left.lastindex, 1)]),
                  y=c(0,    y[x.left.index], 0, 0),
                  ...,
                  col=col.left, border="transparent")
  }

key.dstrplotContinuous <- function(col=c("skyblue3", "deepskyblue3"), lwd=c(3, 4), lty=c(1, 3)) {
  list(corner=c(1,1), x=.98, y=.98,
       border=TRUE,
       text=list(c("non-central", "central")),
       lines=list(lty=lty, col=col, lwd=lwd))
}


dstrplotDiscrete <- function(dstr, ok0=FALSE, ok1=FALSE, args=list(), size=10, ..., x.tick.number=5)
{
  ddstr <- paste("d", dstr, sep="")
  size <- if (missing(size)) args$size else size
  qq <- (!ok0):size
  qq01 <- qq[c(1, length(qq))]

  plotfun <- function(x) do.call(ddstr, c(list(x), args))

  yy <-  plotfun(qq)

  dd <- deparse(args)
  ssdd <- substring(strsplit(dd, "), .Names")[[1]][1], 16)

  result <-
    xyplot(yy ~ qq, horizontal=FALSE,
           panel=function(...) {
             panel.abline(a=0, b=0, col="gray10")
             panel.barchart(...)
           },
           origin=0, ...,
           xlab="x", ylab=NULL, stack=TRUE,,
           scales=list(tck=c(1, .3), x=list(tick.number=x.tick.number)),
           ylim=c(0, max(yy)) + c(-.04, .04)*(max(yy)-0),
           main=paste(ddstr, "(x", if (length(args) > 0) ", " , ssdd, ")", sep=""))

  if (size <= 10) update(result, scales=list(x=list(at=qq)))
  else result
}

panel.dstrplotDiscreteFill <-
  function(x, y, subscripts=seq(length=length(x)), ..., X=5, col.left="#3388ff", col.middle="blue", col.right="skyblue1") {
    x.left.index <- (x < X)
    x.middle.index <- (x == X)
    x.right.index <- (x > X)
    panel.barchart(x, y, subscripts=subscripts, ...,
                   groups=rep(1:3, c(sum(x.left.index), sum(x.middle.index), sum(x.right.index))),
                   col=c(col.left, col.middle, col.right))
  }
