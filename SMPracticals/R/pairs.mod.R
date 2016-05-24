"pairs.mod" <-
function(x, format="MC", labelnames=names(x), highlight=NULL, level=.9, ...)
{
  ## Arguments:
  ## x = data frame with covariates organized column-wise
  ## labelnames = names of the covariates
  ## hightlight = indexes of observations to be highlighted
  ## format = character of presentation type with
  ##   1st character for lower diagonal and 2nd character for upper diagonal.
  ##   o M = Marginal
  ##   o C = Conditional
  ##
  ## Local functions
    doaxis <- function(which, axp, visible, labelnames, srt)
    {
      if(visible && (!is.null(labelnames)))
        visible <- labelnames
      at <- if(!is.null(labelnames))
        seq(axp[1] + 1, axp[2] - 1, length = length(labelnames))
      else seq(axp[1], axp[2], length = axp[3] + 1)
      axis(which, at = round(at, digits = 2), outer =TRUE, line = -0.5, 
           labelnames = visible, srt = srt)
	}
  add.hist <- function(vector, hist.intensities, hist.breaks)
    {
    hist.out <- hist(vector, probability=TRUE, plot=FALSE, breaks=hist.breaks)
    n <- length(hist.out$density)
    xl <- hist.breaks[1]
    for(i in 1:n)
      {
      xr <- hist.breaks[i+1]
      delta <- (xr-xl)*.3
      height <- hist.out$density[i]
      if(height>0) polygon(c(xl+delta,xl+delta,xr-delta,xr-delta),
                           c(0,height,height,0), col=1)
      xl <- xr
      }
    }
  marginal.diag.plot <- function(x.mean.centered, x,
                                 highlight, highlight.char,
                                 i, j, level, format, ...)
    {
    if(format=="MM") x.mean.centered <- x
    if(!is.null(highlight))
      {
      points(as.vector(x.mean.centered[-highlight, j]),
             as.vector(x.mean.centered[-highlight, i]), ...)
      for(ihl in 1:length(highlight))
        {
        points(as.vector(x.mean.centered[highlight[ihl], j]),
               as.vector(x.mean.centered[highlight[ihl], i]),
               pch=highlight.char[ihl], ...)
        }
      }
      else
        {
        points(as.vector(x.mean.centered[, j]),
               as.vector(x.mean.centered[, i]), ...)
        }
    if(format=="MM")
      lines(ellipse(var(as.matrix(x.mean.centered[,c(j,i)])), 
                    centre=c(mean(x.mean.centered[, j]),
                      mean(x.mean.centered[, i])),
                    level=level), lty=2)
    else
      lines(ellipse(var(as.matrix(x.mean.centered[,c(j,i)])), 
                    level=level), lty=2)
    par(lty=1)
    cor.ij <- round(cor(x.mean.centered[, j],
                        x.mean.centered[, i]),2)
    if(cor.ij >= 0)
      text(min(x.mean.centered[, j]), max(x.mean.centered[, i]),
           paste(as.character(cor.ij)), col=1)
      else
        text(min(x.mean.centered[, j]), min(x.mean.centered[, i]),
             paste(as.character(cor.ij)), col=1)
    }
  conditional.diag.plot <- function(x, i, j, highlight, highlight.char,
                                    x.mean.centered, level, ...)
    {
    fit.added.ij <- lm.fit(cbind(1,as.matrix(x[,-c(i,j)])),
                              as.matrix(x[,c(i,j)]), method = "qr")
    if(!is.null(highlight))
      {
      points(as.vector(fit.added.ij$residuals[-highlight,2]),
             as.vector(fit.added.ij$residuals[-highlight,1]), ...)
      for(ihl in 1:length(highlight))
        {
        points(as.vector(fit.added.ij$residuals[highlight[ihl],2]),
               as.vector(fit.added.ij$residuals[highlight[ihl],1]),
               pch=highlight.char[ihl], ...)
        }
      }
      else
        {
        points(as.vector(fit.added.ij$residuals[,2]),
               as.vector(fit.added.ij$residuals[,1]), ...)
        }
    lines(ellipse(var(cbind(as.vector(fit.added.ij$residuals[,2]),
               as.vector(fit.added.ij$residuals[,1]))), 
                   level=level), lty=2)
    par(lty=1)
    par.cor.ij <- round(cor(as.vector(fit.added.ij$residuals[,2]),
                            as.vector(fit.added.ij$residuals[,1])),2)
    if(par.cor.ij >= 0)
      text(min(x.mean.centered[,j]), max(x.mean.centered[,i]),
           paste(as.character(par.cor.ij)), col=1)
      else
        text(min(x.mean.centered[,j]), min(x.mean.centered[,i]),
             paste(as.character(par.cor.ij)), col=1)
    }
  ## Start
  n <- nrow(x)
  p <- ncol(x)
  if(!is.null(highlight))
    {
    if(length(highlight)<=19) highlight.char <- c(1:length(highlight))-1
      else highlight.char <- rep(1, length(highlight))
    }
  oldpar <- par("oma", "mar", "cex", "tck", "mfg", "mgp", "mex", "mfrow")
	oldcex <- par("cex")
	CEX <- oldcex * max(7.7/(2 * p + 3), 0.8)
  par(mfrow=c(1,1)); plot(1,1);
	par(mfrow = c(p, p), mgp = c(2, 0.8, 0), oma = rep(3, 4), mar = rep(0.5, 4), tck = -0.03/p, pty="s")
	on.exit({par(oldpar)})
  ##par(cex = CEX, err=-1)
  x.mean.centered <- t(t(x)-apply(x,2,mean))
  xrange <- list()
  options(warn=-1)
  ## Start loop
  ##-----------
  ## 1- Get the range for the plots.
  ## Requires calculating all the residuals once
  for(i in 1:p) xrange[[i]] <- range(x.mean.centered[,i], na.rm =TRUE)
  for(i in 2:p)
    {
		for(j in 1:(i-1))
      {
      fit.added.ij <- lm.fit(cbind(1,as.matrix(x[,-c(i,j)])),
                                as.matrix(x[,c(i,j)]), method = "qr")
      xrange[[i]] <- range(c(xrange[[i]], fit.added.ij$residuals[,1]))
      xrange[[j]] <- range(c(xrange[[j]], fit.added.ij$residuals[,2]))
      }
  }
  ## 2- Now that the ranges are known, we can refit and plot them...
  for(i in 1:p)
    {
		for(j in 1:p)
      {
      if(i != j)
        {
        ## Off diagonal plots
        ## ------------------
        par(mar = rep(0.5, 4), mgp = c(2, 0.8, 0), adj=0)
        if(format!="MM")
          plot(xrange[[j]], xrange[[i]], type = "n", axes =FALSE, xlab="", ylab="",
                  ...)
        else
          plot(xrange[[j]]+mean(x[,j]), xrange[[i]]+mean(x[,i]), xlab="", ylab="",
                  type = "n", axes =FALSE, ...)
        box()
        if(i == 1)
          doaxis(3, par("xaxp"), j %% 1 == 0, NULL, 0)
        if(i == p)
          doaxis(1, par("xaxp"), j %% 1 == 0, NULL, 0)
        if(j == 1)
          doaxis(2, par("yaxp"), i %% 1 == 0, NULL, 90)
        if(j == p)
          doaxis(4, par("yaxp"), i %% 1 == 0, NULL, 90)
        if(j < i)
          {
          if((format=="MC") || (format=="MM"))
            marginal.diag.plot(x.mean.centered, x, 
                               highlight, highlight.char,
                               i, j, level, format, ...)
            else 
              conditional.diag.plot(x, i, j, highlight, highlight.char,
                                    x.mean.centered, level, ...)
          }
          else
            {
            if((format=="MC") || (format=="CC"))
              conditional.diag.plot(x, i, j, highlight, highlight.char,
                                    x.mean.centered, level, ...)
            else
              marginal.diag.plot(x.mean.centered, x,
                                 highlight, highlight.char,
                                 i, j, level, format, ...)
            }
        }
        else 
          {
          ## Diagonal plots
          ## --------------
          ##par(oma = rep(3, 4), mar = c(5,2,2,.4)+.1, tck = -0.03/p, adj=.5)
          ## 
          par(mar = c(3,0,1,0)+.1, adj=.5)
          hist.i <- hist(x[,i], probability=TRUE, plot=FALSE)
          ##
          fit.i <- lm.fit(cbind(1,as.matrix(x[,-c(i)])),
                             as.vector(x[,i]), method = "qr")
          ##
          xmc <- t(t(x)-apply(x,2,mean))
          Omega.hat <- t(xmc) %*% xmc/(nrow(x)-1)
          par.mean.var.i <- sqrt(c(var(x[,i]),
                                   var(fit.i$residuals)))
          fit.i.res <- fit.i$residuals+mean(x[,i])
          breaks <- hist.i$breaks
          breaks[1] <- min(c(breaks[1], min(fit.i.res)))
          breaks[length(breaks)] <- max(c(breaks[length(breaks)],
                                          max(fit.i.res)))
          hist.res.i <- hist(fit.i.res,
                             probability=TRUE, plot=FALSE, breaks=breaks)
          ylim <- c(0, max(c(hist.i$density, hist.res.i$density)))
          xlim <- range(breaks)
          xlab <- labelnames[i]
          if((format!="CC") && (format!="MM"))
            hist(x[,i], xlab=xlab, probability=TRUE, xlim=xlim, ylim=ylim, col=0, main="", ylab="")
            else 
              {
              plot(xlim, ylim, axes=FALSE, type="n", xlab="", ylab="")
              par(usr = c(0, 1, 0, 1))
              diag.text <- labelnames[i]
              if(format=="MM") 
                {
                diag.text <- paste(diag.text,"\n\nvar=", 
                                 as.character(signif(par.mean.var.i[1],3)))
                }
                else 
                  {
                  diag.text <- paste(diag.text,"\n\nvar=",
                                   as.character(signif(par.mean.var.i[2],3)))
                  }
              text(0.5, 0.5, diag.text, cex=1.5*CEX)
              }
          if((format != "MM") && (format != "CC"))
            add.hist(fit.i$residuals+mean(x[,i]),
                     hist.intensities=hist.i$density, hist.breaks=breaks)
          if((format=="MC"))
            {
            par(adj=0)
            text(xlim[1], ylim[2],
                 as.character(signif(par.mean.var.i[1],3)),col=1)
            par(adj=1)
            text(xlim[2], ylim[2],
                 as.character(signif(par.mean.var.i[2],3)),col=1)
            }
          if(format=="CM")
             {
             par(adj=0)
             text(xlim[1], ylim[2],
                 as.character(signif(par.mean.var.i[2],3)),col=1)
             par(adj=1)
             text(xlim[2], ylim[2],
                 as.character(signif(par.mean.var.i[1],3)),col=1)
             }
          print(par.mean.var.i)
          }
      }
    }
  }

