odds.ratio <- function(...)
  .Defunct("OddsRatio", package="HH")

"OddsRatio" <-
function(x, alpha=.05) {
  p1 <- x[1,1] / (x[1,1] + x[1,2])
  p2 <- x[2,1] / (x[2,1] + x[2,2])
  omega1 <- p1 / (1-p1)
  omega2 <- p2 / (1-p2)

  if (length(alpha) != 1  || (alpha < 0 || 1 < alpha))
    stop("alpha must be a number in the range [0,1].")
  psihat <- as.numeric(omega2 / omega1)
  s.ln.psihat <- sqrt(sum(1/x))
  ci.ln.psihat <- log(psihat) + c(-1,1) * qnorm(1-alpha/2) * s.ln.psihat
  ci.psihat <- exp(ci.ln.psihat)

  prob1 <- seq(0,1,.05)
  odds1 <- prob1/(1-prob1)
  odds2 <- odds1*psihat
  ci.odds2 <- odds1 %o% ci.psihat
  prob2 <- odds2 / (odds2+1)
  prob2[is.na(prob2)] <- 1
  ci.prob2 <- ci.odds2 / (ci.odds2+1)
  ci.prob2[is.na(ci.prob2)] <- 1

  list(p1 = p1,
       p2 = p2,
       omega1 = omega1,
       omega2 = omega2,
       psihat = psihat,
       s.ln.psihat = s.ln.psihat,
       ci.ln.psihat = ci.ln.psihat,
       ci.psihat = ci.psihat,
       prob1 = prob1,
       odds1 = odds1,
       odds2 = odds2,
       ci.odds2 = ci.odds2,
       prob2 = prob2,
       ci.prob2 = ci.prob2)
}

plot.odds.ratio <- function(...)
  .Defunct("plotOddsRatio", package="HH")

"plotOddsRatio.base" <-
  function(x,
           ylab="prob(col1 | row1)",
           xlab="prob(col1 | row2)",
           alpha=c(.05, .50),
           legend.x=1.05,
           oma=c(0,0,0,5),
           ...) {
  if (missing(xlab) && missing(ylab) && length(dimnames(x))==2) {
   col1 <- dimnames(x)[[2]][1]
   row1 <- dimnames(x)[[1]][1]
   row2 <- dimnames(x)[[1]][2]
   ylab <- paste("prob(",col1,"|",row1,")")
   xlab <- paste("prob(",col1,"|",row2,")")
  }
  old.pty <- par(pty="s")
  old.oma <- par(oma=oma)
  for (i in seq(along=alpha)) {
    tmp <- OddsRatio(x, alpha[i])
    if (i == 1) {
      matplot(y=tmp$prob1, x=cbind(tmp$prob2, tmp$ci.prob2),
              type="l", lty=c(1,3,3), col=1,
              xlab=xlab,
              ylab=ylab)
      abline(a=0, b=1, lty=2)
      points(y=tmp$p1, x=tmp$p2, pch=13)
      if.R(r=old.xpd <- par(xpd=NA), s={})
      legend(x=legend.x, y=.6, lty=c(1,4,3,2), col=1,
             legend=c("  0%CI", "50%CI", "95%CI", "x=y"))
      ## if.R(s=
      ##      legend(x=legend.x, y=.4, marks=13, legend=c("MLE   "))
      ##      ,r=
           legend(x=legend.x, y=.4, pch=13,   legend=c("MLE   "))
           ## )
      if.R(r=par(old.xpd), s={})
      result <- tmp
    }
    else
      matlines(y=tmp$prob1, x=cbind(tmp$prob2, tmp$ci.prob2),
               type="l", lty=c(1,4,4), col=1)
  }
  par(old.pty)
  par(old.oma)
  invisible(result)
}

