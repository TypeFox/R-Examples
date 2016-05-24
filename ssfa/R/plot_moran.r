plot_moran <- function (x, main, xlab, ylab, labels, listw, ...) 
{
  name <-deparse(substitute(x))
  
  if(missing(xlab))
    xlab <- paste("Residuals (",name,")")
  if(missing(ylab))
    ylab <- paste("Spatially lagged residuals (",name,")")
  if(missing(labels))
    labels <- FALSE
  if(missing(main))
    main <- "Moran scatterplot"

  moran.plot(c(residuals(x)), listw=listw, zero.policy=NULL, spChk=NULL, labels=labels,
             xlab=xlab, ylab=ylab, main=main, quiet=NULL, ...)
  }

  