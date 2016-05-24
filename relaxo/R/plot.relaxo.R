`plot.relaxo` <-
function(x, type= "l", lty = 1,
                        main = NULL,
                        xlab = "|beta|/max|beta| (phi=1)",
                        ylab = expression("coefficients  " * beta[j]),
                        plotphi = unique(x$phi),
                        ...)
{
  stopifnot(inherits(x, "relaxo"))
  if(any(!plotphi%in%x$phi)) stop("values of plotphi must be among the available values in object$phi")
  if(nrow(x$beta)<2) stop("plot function not useful for result of cv.relaxo selection")

  nphi <- length(plotphi)
  on.exit(par(op))
  op <-  par(mfrow = if(1 < nphi && nphi < 7) c(2, ceiling(nphi/2))  else  c(3, ceiling(nphi/3)))
  
  
  fracbeta <- apply( abs(x$beta[ x$phi==1, ]),1,sum)
  for (phi in sort(plotphi)) {
      matplot(fracbeta/max(fracbeta),
              x$beta  [x$phi == phi, ],
              type = type, lty = lty,
              main = if(is.null(main)) substitute(phi == P, list(P = format(phi))),
              xlab = xlab, ylab = ylab)
  }

}

