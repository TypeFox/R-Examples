residual.plots.lattice <- function(lm.object, X=dft$x, layout=c(dim(X)[2],1),
                           par.strip.text=list(cex=.8),
                           scales.cex=.6,
                           na.action=na.pass,
                           y.relation="same",
                           ...) {
  lm.formula <- as.formula(lm.object)
  lm.data <- try(eval(lm.object$call$data), silent=TRUE)
  if (class(lm.data) == "Error" || class(lm.data)=="try-error") ## S-Plus || R
    {
      lm.data <- lm.object$x
      lm.data.y <- as.numeric(lm.object$y)
      if (is.null(lm.data) || is.null(lm.data.y))
        stop("Please recompute the 'lm.object' with 'x=TRUE, y=TRUE'.")
      lm.data <- cbind(lm.data.y, lm.data)
      names(lm.data)[1] <- as.character(lm.formula[[2]])
    }

  dft <- do.formula.trellis.xysplom(lm.formula, lm.data, na.action)

  resids <- resid(lm.object)
  yhat <- predict(lm.object, type="terms")

  X <- data.frame(X, check.names=FALSE)

  if (dim(yhat)[[2]] != dim(X)[[2]])
   stop("The model has factors or interactions.  Please use the `X=' argument.")
  partial.resids <- yhat + resids
  dimnames(partial.resids)[[2]] <- rep("part.res|X", ncol(partial.resids))

  ## data plots
  ## y against X
  y.X <- latticeresids(y ~ x | xname, # * yname,
                       data=data.frame(y=rep(dft$y[[1]], ncol(X)),
                         x=unlist(X),
                         yname=factor(names(dft$y)),
                         xname=factor(rep(names(X), each=nrow(X)), levels=names(X))),
                       main=paste(names(dft$y), " ~ x variables", sep=""),
                       par.strip.text=par.strip.text,
                       scales.cex=scales.cex,
                       y.relation=y.relation,
                       ...
                       )

  ## residual plots
  ## residuals against X  res.X <- xyresidplot(cbind(residuals=resids, X,
  res.X <- latticeresids(y ~ x | xname, # *  * yname,
                       data=data.frame(y=rep(resids, ncol(X)),
                         x=unlist(X),
                         yname=factor("residuals"),
                         ## yname=factor(dimnames(partial.resids)[[2]][1]),
                         xname=factor(rep(names(X), each=nrow(X)), levels=names(X))),
                       main="residuals ~ x variables",
                       par.strip.text=par.strip.text,
                       scales.cex=scales.cex,
                       y.relation=y.relation,
                       ...
                       )


  ## partial residuals plots
  ## partial residuals against X
  pres.X <- latticeresids(y ~ x | xname, # *  * yname,
                          data=data.frame(y=as.vector(partial.resids),
                            x=unlist(X),
                            yname=factor("partial Residuals"),
                            xname=factor(rep(names(X), each=nrow(X)), levels=names(X))),
                          main="(partial residuals of y against the other X columns) ~ x variables",
                          par.strip.text=par.strip.text,
                          scales.cex=scales.cex,
                          y.relation=y.relation,
                          ...
                          )


  ## added variable plots
  ## partial residuals against X.j
  X.res <- X.residuals(lm.object)
  names(X.res) <- paste(names(X.res), "X", sep=" | ")

  firstColumn <- function(x) {
    llx <- length(levels(x))
    if (llx==0) TRUE else c(TRUE, rep(FALSE, max(llx-2, 0)))
  }
  X.resSubscript <- X.res[ , unlist(sapply(X, firstColumn))]

  main4 <- if (length(X.resSubscript) == length(X.res))
    "partial residuals of y against the other X columns ~ residuals of x against the other X columns"
  else
     "partial residuals of y against the other X columns ~ residuals of x against the other X columns\nOnly the first dummy variable is shown for factors"

  pres.Xj <- latticeresids(y ~ x | xname, # *  * yname,
                          data=data.frame(y=as.vector(partial.resids),
                            x=unlist(X.resSubscript),
                            yname=factor("partial Residuals"),
                            xname=factor(rep(names(X.resSubscript), each=nrow(X.resSubscript)), levels=names(X.resSubscript))),
                          main="(partial residuals of y against the other X columns) ~ (residuals of x against the other X columns)",
                           par.strip.text=par.strip.text,
                           scales.cex=scales.cex,
                           y.relation=y.relation,
                           ...
                          )



  result <- list(y.X=y.X, res.X=res.X, pres.X=pres.X, pres.Xj=pres.Xj)
  class(result) <- c("latticeresids", class(result))
  result
}


latticeresids <- function(x, data,
                          main="please use an appropriate main title",
                          par.strip.text,
                          scales.cex,
                          y.relation,
                          ...) {
  LP <- xyplot(x, data, main=main,
               layout=c(length(levels(data$xname)), 1),
               panel=function(x, y, ...) {
                 panel.xyplot(x, y, ...)
                 panel.abline(lm(y ~ x))
               },
               xlab=NULL, ylab=NULL,
               par.strip.text=par.strip.text,
               between=list(x=1, y=1),
               scales=list(
                 cex=scales.cex,
                 x=list(relation="free"),
                 y=list(relation=y.relation, rot=0),
                 alternating=FALSE),
               ...
               )
  ## combineLimits(useOuterStrips(LP))
}


print.latticeresids <-
  function(x, ...,
           A321.left=0, A321.bottom=0.27,
           A4.left=0, A4.top=0.30,
           position=list(
             A321=c(A321.left,     A321.bottom, 1, 1     ),
             A4  =c(A4.left,       0,           1, A4.top)),
           panel.width=NULL,
           which=1:4) {
  yname <- strsplit(x[[1]]$main, " ~ ")[[1]][1]
  names(x) <- c(yname, "Residuals", "Partial Residuals | X", "Partial Residuals | X")
  A321 <- do.call(rbind, x[(3:1)[3:1 %in% which]])
  A321.present <- !is.null(dim(A321))
  if (A321.present)
    A321 <- combineLimits(update(A321, scales=list(relation="free")))
  A4   <- do.call(rbind, x[4[4 %in% which]]) ## this call puts the left strip label in place
  A4.present <- !is.null(dim(A4))
  if (A321.present)
    print(position=position$A321, more=A4.present, update(A321, main=NULL), panel.width=panel.width, ...)
  if (A4.present)
    print(position=position$A4, more=FALSE,
          update(A4, main=NULL, scales=list(tck=c(1,0))), panel.width=panel.width, ...)

  invisible(x)
}
