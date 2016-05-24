.reg5Plot <-
function(lm.out, res.rows=NULL, pred.rows=NULL,
         scatter.coef=FALSE, X1.new=NULL,
         numeric.all, in.data.frame, c.int, p.int,
         pdf=FALSE, pdf.width=5, pdf.height=5, manage.gr=FALSE,
         scatter.3D, ...) {

  nm <- all.vars(lm.out$terms)  # names of vars in the model
  n.vars <- length(nm)
  n.pred <- n.vars - 1L
  n.obs <- nrow(lm.out$model)
  n.keep <- nrow(lm.out$model)
  if (is.null(pred.rows)) if (n.keep < 25) pred.rows <- n.keep else pred.rows <- 4 
  if (pred.rows == "all") pred.rows <- n.keep  # turn off preds with pred.rows=0


  # pdf graphics option
  if (pdf) { 
    pdf.file <- "RegScatterplot.pdf"
    if (n.pred > 1) pdf.file <- "RegScatterMatrix.pdf"
    pdf(file=pdf.file, width=pdf.width, height=pdf.height)
  }

  # keep track of the plot in this routine
  plt.i <- 0L
  plt.title  <- character(length=0)

  if (n.pred <= 1) {  # scatterplot, if one predictor variable

    #if (n.pred == 0)
      #do.predint <- FALSE
    #else {
      if ( (pred.rows==0) || !is.null(X1.new) || is.null(p.int)) 
        do.predint <- FALSE
      else 
        do.predint <- TRUE
     #}
     if (n.pred > 0) if (is.factor(lm.out$model[,nm[2]])) do.predint <- FALSE

    if (!do.predint) {
      ctitle <- "Scatterplot"
      if (n.pred > 0) if (!is.factor(lm.out$model[,nm[2]]))
        ctitle <- paste(ctitle, "and Regression Line")
      y.min <- min(lm.out$model[,nm[1]])
      y.max <- max(lm.out$model[,nm[1]])
    }
    else {
      ctitle <- "Regression Line,\nConfidence and Prediction Intervals"
      y.min <- min(p.int$lwr)
      y.max <- max(max(p.int$upr),  max(lm.out$model[,nm[1]]) )
    }

    fl <- "ls"
    if (n.pred > 0) {
      if (is.factor(lm.out$model[,nm[2]])) fl <- "none"
    }
    else
      fl <- "none"

    plt.i <- plt.i + 1L
    plt.title[plt.i] <- gsub(pattern="\n", replacement=" ", x=ctitle)

    if (n.pred > 0)
      x.values <- lm.out$model[,nm[2]]
    else
      x.values <- 1:n.obs
    y.values <- lm.out$model[,nm[1]] 
    .plt.main(x.values, y.values,
       shape=21, size=.8, xlab=nm[2], ylab=nm[1], main=ctitle,
       fit.line=fl, quiet=TRUE, ylim=c(y.min,y.max))

    if (n.pred == 0) {
      m <- lm.out$coefficients[1]  # mean of Y
      mv <- rep(m, n.obs)
      names(mv) <- NULL
      lines(x.values, mv, lwd=0.75)
    }

    if (do.predint) {
      col.ci <- getOption("col.stroke.pt")
      col.pi <- "gray30"

      lines(x.values, c.int$lwr, col=col.ci, lwd=0.75)
      lines(x.values, c.int$upr, col=col.ci, lwd=0.75)
      lines(x.values, p.int$lwr, col=col.pi, lwd=1.5)
      lines(x.values, p.int$upr, col=col.pi, lwd=1.5)
    }
  }  # end n.pred==1


  else {  # scatterplot matrix for multiple regression
    if (numeric.all && in.data.frame) {
      col.pts <- getOption("col.stroke.pt")
      col.line <- getOption("col.stroke.bar")
      col.bg=getOption("col.bg")

      panel2.smooth <- function (x, y, pch=par("pch"), cex=.9,
        col.pt=col.pts, col.smooth=col.line,
        span=2/3, iter=3, ...) 
      {
          points(x, y, pch=pch, col=col.pt, cex=cex)
          ok <- is.finite(x) & is.finite(y)
          if (any(ok)) 
            lines(lowess(x[ok], y[ok], f=span, iter=iter), col=col.smooth, ...)
      }

      plt.i <- plt.i + 1L
      plt.title[plt.i] <- "ScatterPlot Matrix"

      if (scatter.coef) {
        panel.cor <- function(x, y, digits=2, prefix="", cex.cor, ...) {
          usr <- par("usr"); on.exit(par(usr))
          par(usr=c(0, 1, 0, 1))
          r <- cor(x, y)
          txt <- format(c(r, 0.123456789), digits=digits)[1]
          txt <- paste(prefix, txt, sep="")
          if (missing(cex.cor)) cex.cor <- .9/strwidth(txt)
          cex.adj <- 2.5 - (0.18*n.pred)  # adjust size of displayed r
          text(0.5, 0.5, txt, cex=cex.adj, col=col.pts)  # or cex=cex.cor * r
        }
        pairs(lm.out$model[c(nm)],
          lower.panel=panel2.smooth, upper.panel=panel.cor)
      }
      else pairs(lm.out$model[c(nm)], panel=panel2.smooth)
    }
    else {
      cat("\n>>> No scatterplot matrix reported because not all variables are ")
      if (!in.data.frame) cat("in the data frame.\n")
      if (!numeric.all) cat("numeric.\n")
      dev.off()
    }
  }

  if (pdf) {
    dev.off()
    if (n.pred == 1)
      .showfile(pdf.file, "scatterplot")
    else
      .showfile(pdf.file, "scatterplot matrix")
    cat("\n\n")
  }

  if (scatter.3D) {  # 3d scatterplot option for 2-predictor models
    cat("\n"); stop(call.=FALSE, "\n","------\n",
      "scatter.3D option disabled\n",
      "car package no longer included because of dependencies issues\n\n",
      "If interested, directly call the scatter3d function from the car package\n",
      "First install the needed packages, then invoke the library function:\n",
      "  install.packages(\"rgl\", \"car\")\n",
      "  library(car)\n\n",
      "Example\n ",
      "  scatter3d(Ozone ~ Wind + Temp, id.method=\"identify\", data=airquality)\n\n",
      "Directions\n",
      "  Can re-size the plot window, click and drag to rotate plot\n",
      "  Press the right mouse button and drag a rectangle around any points to be\n",
      "    identified, and then release\n",
      "  To exit, right-click in a blank area of the 3d-scatterplot\n\n", sep="")

      #suppressMessages(scatter3d(lm.out$terms, id.method="identify", data=lm.out$model))  # car
  }

  # just generated plot
  invisible(list(i=plt.i, ttl=plt.title))

}
