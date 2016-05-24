plot.lm2 <- function(
                     x,
                     which = 1:5,
                     caption = c("Residuals vs Fitted", "Normal Q-Q plot",
                                 "Scale-Location plot", "Cook's distance plot"),
                     panel = panel.smooth,
                     sub.caption = deparse(x$call),
                     main = "",
                     ask,
                     ...,
                     id.n = 3,
                     labels.id = names(residuals(x)),
                     cex.id = 0.75,
                     band=TRUE,
                     rug=TRUE,
                     width=1/10,
                     max.n=5000
                     )
{
  .Defunct("lmplot2", "gplots")
}


lmplot2 <- function(
                    x,
                    which = 1:5,
                    caption = c("Residuals vs Fitted", "Normal Q-Q plot",
                                "Scale-Location plot", "Cook's distance plot"),
                    panel = panel.smooth,
                    sub.caption = deparse(x$call),
                    main = "",
                    ask = interactive() && nb.fig < length(which)
                    && .Device != "postscript",
                    ...,
                    id.n = 3,
                    labels.id = names(residuals(x)),
                    cex.id = 0.75,
                    band=TRUE,
                    rug=TRUE,
                    width=1/10,
                    max.n=5000
                    )
{
  if (!inherits(x, "lm"))
    stop("Use only with 'lm' objects")
  show <- rep(FALSE, 5)
  if(!is.numeric(which) || any(which < 1) || any(which > 5))
    stop("`which' must be in 1:5")
  show[which] <- TRUE
  r <- residuals(x)
  n <- length(r)
  if(inherits(x,"glm"))
    yh <- predict(x) # != fitted() for glm
  else
    yh <- fitted(x)
  if (any(show[2:4]))
    s <- if(inherits(x, "rlm")) x$s else sqrt(deviance(x)/df.residual(x))
  if (any(show[2:3]))
    {
      ylab23 <- if(inherits(x, "glm"))
        "Std. deviance resid." else "Standardized residuals"
      hii <- lm.influence(x)$hat
      w <- weights(x)
                                        # r.w := weighted.residuals(x):
      r.w <- if(is.null(w)) r else (sqrt(w)*r)[w!=0]
      rs <- r.w/(s * sqrt(1 - hii))
    }
  if (any(show[c(1,3)]))
    l.fit <- if(inherits(x,"glm"))
      "Predicted values" else "Fitted values"
  if (is.null(id.n))
    id.n <- 0
  else {
    id.n <- as.integer(id.n)
    if(id.n < 0 || id.n > n)
      stop(paste("`id.n' must be in { 1,..,",n,"}"))
  }
  if(id.n > 0) {
    if(is.null(labels.id))
      labels.id <- paste(1:n)
    iid <- 1:id.n
    show.r <- order(-abs(r))[iid]
    if(any(show[2:3]))
      show.rs <- order(-abs(rs))[iid]
    text.id <- function(x,y, ind, adj.x = FALSE)
      text(x - if(adj.x) strwidth(" ")*cex.id else 0, y, labels.id[ind],
           cex = cex.id, xpd = TRUE, adj = if(adj.x) 1)
  }
  nb.fig <- prod(par("mfcol"))
  one.fig <- prod(par("mfcol")) == 1
  if (ask) {
    op <- par(ask = TRUE)
    on.exit(par(op))
  }

  ##---------- Do the individual plots : ----------
  if (show[1]) {
    ylim <- range(r)
    if(id.n > 0)
      ylim <- ylim + c(-1,1)* 0.08 * diff(ylim)
    plot(yh, r, xlab = l.fit, ylab = "Residuals", main = main,
         ylim = ylim, type = "n", ...)
    panel(yh, r, ...)
    if(rug)  rug(yh)                                 ## GRW 2001-06-08
    if(band) bandplot(yh,r,add=TRUE,width=width)    ## GRW 2001-06-08
    if (one.fig)
      title(sub = sub.caption, ...)
    mtext(caption[1], 3, 0.25)
    if(id.n > 0) {
      y.id <- r[show.r]
      y.id[y.id < 0] <- y.id[y.id < 0] - strheight(" ")/3
      text.id(yh[show.r], y.id, show.r, adj.x = TRUE)
    }
    abline(h = 0, lty = 3, col = "gray")
  }
  if (show[2]) {
    ylim <- range(rs)
    ylim[2] <- ylim[2] + diff(ylim) * 0.075
    qq <- qqnorm(rs, main = main, ylab = ylab23, ylim = ylim, ...)
    qqline(rs)
    if (one.fig)
      title(sub = sub.caption, ...)
    mtext(caption[2], 3, 0.25)
    if(id.n > 0)
      text.id(qq$x[show.rs], qq$y[show.rs], show.rs, adj.x = TRUE)
  }
  if (show[3]) {
    sqrtabsr <- sqrt(abs(rs))
    ylim <- c(0, max(sqrtabsr))
    yl <- as.expression(substitute(sqrt(abs(YL)), list(YL=as.name(ylab23))))
    yhn0 <- if(is.null(w)) yh else yh[w!=0]
    plot(yhn0, sqrtabsr, xlab = l.fit, ylab = yl, main = main,
         ylim = ylim, type = "n", ...)
    panel(yhn0, sqrtabsr, ...)

    abline(h=mean(sqrtabsr),lty = 3, col = "gray")
    if(rug)  rug(yh)                             ## GRW 2001-06-08
    if(band) bandplot(yhn0,sqrtabsr,add=TRUE) ## GRW 2001-06-08

    if (one.fig)
      title(sub = sub.caption, ...)
    mtext(caption[3], 3, 0.25)
    if(id.n > 0)
      text.id(yhn0[show.rs], sqrtabsr[show.rs], show.rs, adj.x = TRUE)
  }
  if (show[4]) {
    cook <- cooks.distance(x, sd=s)
    if(id.n > 0) {
      show.r <- order(-cook)[iid]# index of largest `id.n' ones
      ymx <- cook[show.r[1]] * 1.075
    } else ymx <- max(cook)
    plot(cook, type = "h", ylim = c(0, ymx), main = main,
         xlab = "Obs. number", ylab = "Cook's distance", ...)
    if (one.fig)
      title(sub = sub.caption, ...)
    mtext(caption[4], 3, 0.25)
    if(id.n > 0)
      text.id(show.r, cook[show.r] + 0.4*cex.id * strheight(" "), show.r)
  }

  if (show[5])
  {
    ## plot residuals against each predictor ##
    data  <- model.frame(x)
    for( i in 1:ncol(data) )
    {
      test <- try(
      {
        plot.default( x=data[,i], y=r,
                     xlab=names(data)[i], ylab="Residuals", type="n")
        panel( data[,i], r, ... )
        if(rug)  rug(data[,i])
        if(band) bandplot(data[,i],r,add=TRUE)
        abline(h=0,lty = 3, col = "gray")
      }
      )
    }
  }

  if (!one.fig && par("oma")[3] >= 1)
    mtext(sub.caption, outer = TRUE, cex = 1.25)
  invisible()
}
