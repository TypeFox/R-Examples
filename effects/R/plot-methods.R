# plot.eff method for effects package, moved here from plot-summary-print-methods.R
# The plot.effpoly method remains there for now.
# 2013-10-17: Added use.splines keyword to plot.eff. Sandy
# 2013-10-17: Made ci.style="bands" default for variates; allow "bands" if multiline=TRUE
# 2013-10-29: fixed plot.eff() to handle factors with "valid" NA level. J. Fox
# 2014-03-03: modified plot.eff() to handle partial residuals. J. Fox
# 2014-09-20: fixed plot.eff() to work with partial residuals when rescale.axis=FALSE;
#             added smooth.residuals argument. J. Fox
# 2014-10-10: namespace fixes. J. Fox
# 2014-12-05: made key.args more flexible. J. Fox
# 2015-03-22: use wide columns by default only when x for legend not set. J. Fox
# 2015-03-25: use non-robust loess smooth for partial residuals for non-Gaussian families. J. Fox
# 2015-03-25: rationalized type and rescale.axis args to plot.eff(); deprecated rescale.axis arg. J. Fox
# 2015-05-28: added residuals.smooth.color argument. J. Fox
# 2015-08-28: added residuals.cex argument. J. Fox
# 2016-03-01: move computation of partial residuals to the plot.eff() method. J. Fox

# the following functions aren't exported

find.legend.columns <- function(n, target=min(4, n)){
  rem <- n %% target
  if (rem != 0 && rem < target/2) target <- target - 1
  target
}


make.ticks <- function(range, link, inverse, at, n) {
  warn <- options(warn=-1)
  on.exit(options(warn))
  link <- if (is.null(link)) 
    function(x) nlm(function(y) (inverse(y) - x)^2, 
                    mean(range))$estimate
  else link
  if (is.null(n)) n <- 5
  labels <- if (is.null(at)){
    labels <- pretty(sapply(range, inverse), n=n+1)
  }
  else at
  ticks <- sapply(labels, link)
  list(at=ticks, labels=format(labels))
}

range.adj <- function(x){
  range <- range(x, na.rm=TRUE)
  c(range[1] - .025*(range[2] - range[1]),                                              
    range[2] + .025*(range[2] - range[1]))
}

# added, modified from http://www.r-bloggers.com/confidence-bands-with-lattice-and-r/

panel.bands <- function(x, y, upper, lower, fill, col,
                        subscripts, ..., font, fontface, use.splines=FALSE)
{
  if(!missing(subscripts)) {
    upper <- upper[subscripts]
    lower <- lower[subscripts]
  }
  if (use.splines){
    up <- spline(x, upper)
    down <- spline(x, lower)
    x <- up$x
    upper <- up$y
    lower <- down$y
  }
  panel.polygon(c(x, rev(x)), c(upper, rev(lower)),
                col = fill, fill=fill, border = FALSE,
                ...)
}


# modified by Michael Friendly: added key.args:
# modified by Michael Friendly: added ci.style="bands"
# modified by Michael Friendly: added lwd= argument for llines (not used elsewhere)
# modified by Michael Friendly: added alpha.band= argument for ci.style="bands"

spline.llines <- function(x, y, ...) llines(spline(x, y), ...)

plot.eff <- function(x, x.var,
                     z.var=which.min(levels), multiline=is.null(x$se), rug=TRUE, 
                     xlab, ylab, main=paste(effect, "effect plot"),
                     colors=palette(), symbols=1:length(colors), lines=1:length(colors),
                     cex=1.5, lwd=2, ylim, xlim=NULL,
                     factor.names=TRUE, ci.style, band.transparency=0.15, band.colors=colors,
                     type=c("rescale", "response", "link"), ticks=list(at=NULL, n=5),  
                     alternating=TRUE, rotx=0, roty=0, grid=FALSE, layout, rescale.axis, 
                     transform.x=NULL, ticks.x=NULL,
                     key.args=NULL, 
                     row=1, col=1, nrow=1, ncol=1, more=FALSE, 
                     use.splines=TRUE, partial.residuals=TRUE,
                     show.fitted=FALSE,
                     residuals.color="blue", residuals.pch=1, residuals.cex=1,
                     smooth.residuals=TRUE, residuals.smooth.color=residuals.color, span=2/3, ...)
{  
  closest <- function(x, x0) apply(outer(x, x0, FUN=function(x, x0) abs(x - x0)), 1, which.min)
  .mod <- function(a, b) ifelse( (d <- a %% b) == 0, b, d)
  .modc <- function(a) .mod(a, length(colors))
  .mods <- function(a) .mod(a, length(symbols))
  .modl <- function(a) .mod(a, length(lines))
  .modb <- function(a) .mod(a, length(band.colors))
  if (is.character(partial.residuals)){
    partial.residuals <- TRUE
    warning("partial.residuals='adjusted' or 'raw' is deprecated\n",
            "partial.residuals will be set to TRUE")
  }
  ci.style <- if(missing(ci.style)) NULL else
    match.arg(ci.style, c("bars", "lines", "bands", "none")) 
  if (smooth.residuals && !is.null(x$family)){
    loess.family <- if (x$family == "gaussian") "symmetric" else "gaussian"
  }
  type <- if (!missing(rescale.axis)){
    message("NOTE: the rescale.axis argument is deprecated; use type instead")
    if (!is.logical(rescale.axis)) stop ("rescale.axis must be logical (TRUE or FALSE)")
    if (rescale.axis) "rescale" else "response"
  }
  else match.arg(type)
  switch(type,
         rescale = {
           type <- "response"
           rescale.axis <- TRUE
         },
         response = {
           type <- "response"
           rescale.axis <- FALSE
         },
         link = {
           type <- "link"
           rescale.axis <- TRUE
         }
  )
  levels <- sapply(x$variables, function(z) length(as.vector(z[["levels"]])))
  thresholds <- x$thresholds
  has.thresholds <- !is.null(thresholds)
  effect.llines <- llines
  if (missing(ylab)){
    ylab <- if (has.thresholds) paste(x$response, ": ", paste(x$y.levels, collapse=", "), sep="")
    else x$response
  }     
  if (has.thresholds){ 
    threshold.labels <- abbreviate(x$y.levels, minlength=1)
    threshold.labels <- paste(" ", 
                              paste(threshold.labels[-length(threshold.labels)], threshold.labels[-1], sep=" - "),
                              " ", sep="")
  }
  trans.link <- x$transformation$link
  trans.inverse <- x$transformation$inverse
  residuals <- if (partial.residuals) x$residuals else NULL
  partial.residuals.range <- x$partial.residuals.range
  if (!is.null(residuals) && !rescale.axis) {
    residuals <- trans.inverse(residuals)
  }
  if (!rescale.axis){
    x$lower[!is.na(x$lower)] <- trans.inverse(x$lower[!is.na(x$lower)])
    x$upper[!is.na(x$upper)] <- trans.inverse(x$upper[!is.na(x$upper)])
    x$fit[!is.na(x$fit)] <- trans.inverse(x$fit)[!is.na(x$fit)]
    trans.link <- trans.inverse <- I
  }
  x.all <- x$x.all
  split <- c(col, row, ncol, nrow)
  ylab # force evaluation
  if (missing(x.var)) x.var <- x$x.var
  if (!is.null(x.var) && is.numeric(x.var)) x.var <- names(x.var)
  x.data <- x$data
  effect <- paste(sapply(x$variables, "[[", "name"), collapse="*")
  vars <- x$variables
  x <- as.data.frame(x, transform=I)
  for (i in 1:length(vars)){
    if (!(vars[[i]]$is.factor)) next
    x[,i] <- factor(x[,i], levels=vars[[i]]$levels, exclude=NULL)
  }
  has.se <- !is.null(x$se)
  n.predictors <- ncol(x) - 1 - 3*has.se
  if (n.predictors == 1){
    predictor <- names(x)[1]
    ### factor no other predictors
    if (is.factor(x[,1])){
      ci.style <- if(is.null(ci.style)) "bars" else ci.style
      range <- if(has.se & ci.style!="none")
        range(c(x$lower, x$upper), na.rm=TRUE) else range(x$fit, na.rm=TRUE)
      ylim <- if (!missing(ylim)) ylim else c(range[1] - .025*(range[2] - range[1]),                                              
                                              range[2] + .025*(range[2] - range[1]))
      tickmarks <- if (type == "response") make.ticks(ylim, 
                                                      link=trans.link, inverse=trans.inverse, at=ticks$at, n=ticks$n)
      else make.ticks(ylim, link=I, inverse=I, at=ticks$at, n=ticks$n)
      levs <- levels(x[,1])  
      plot <- xyplot(eval(parse(
        text=paste("fit ~ as.numeric(", names(x)[1], ")"))), 
        strip=function(...) strip.default(..., strip.names=c(factor.names, TRUE)),
        panel=function(x, y, lower, upper, has.se, ...){
          if (grid) panel.grid()
          good <- !is.na(y)
          if (has.se){ 
            if (ci.style == "bars"){
              larrows(x0=x[good], y0=lower[good], x1=x[good], y1=upper[good], angle=90, 
                      code=3, col=colors[.modc(2)], length=0.125*cex/1.5)
            }
            else if(ci.style == "lines") {
              effect.llines(x[good], lower[good], lty=2, col=colors[.modc(2)])
              effect.llines(x[good], upper[good], lty=2, col=colors[.modc(2)])
            }
            else{ if(ci.style == "bands") { 
              panel.bands(x[good], y[good], upper[good], lower[good], fill=band.colors[1],
                          alpha=band.transparency, use.splines=FALSE)
            }}
          }
          effect.llines(x[good], y[good], lwd=lwd, col=colors[1], type='b', pch=19, cex=cex, ...)
          if (has.thresholds){
            panel.abline(h=thresholds, lty=3)
            panel.text(rep(current.panel.limits()$xlim[1], length(thresholds)), 
                       thresholds, threshold.labels, adj=c(0,0), cex=0.75)
            panel.text(rep(current.panel.limits()$xlim[2], length(thresholds)), 
                       thresholds, threshold.labels, adj=c(1,0), cex=0.75)
          }
        },
        ylim=ylim,
        ylab=ylab,
        xlab=if (missing(xlab)) names(x)[1] else xlab,
        scales=list(x=list(at=1:length(levs), labels=levs, rot=rotx), 
                    y=list(at=tickmarks$at, labels=tickmarks$labels, rot=roty),
                    alternating=alternating, y=roty),
        main=main,
        lower=x$lower, upper=x$upper, has.se=has.se, data=x, ...)
      result <- update(plot, layout = if (missing(layout)) c(0, prod(dim(plot))) 
                       else layout)
      result$split <- split
      result$more <- more
      class(result) <- c("plot.eff", class(result))
    }  
    ### variate, no other predictors  ***
    else {
      effect.llines <- if(use.splines) spline.llines else effect.llines
      ci.style <- if(is.null(ci.style)) "bands" else ci.style
      range <- if(has.se & ci.style!="none")
        range(c(x$lower, x$upper), na.rm=TRUE) else range(x$fit, na.rm=TRUE)
      
      ylim <- if (!missing(ylim)) ylim
      else if (is.null(residuals)) c(range[1] - .025*(range[2] - range[1]), range[2] + .025*(range[2] - range[1]))
      else c(min(partial.residuals.range[1], range[1] - .025*(range[2] - range[1])), 
             max(partial.residuals.range[2], range[2] + .025*(range[2] - range[1])))
      
      tickmarks <- if (type == "response") make.ticks(ylim, 
                                                      link=trans.link, inverse=trans.inverse, at=ticks$at, n=ticks$n)
      else make.ticks(ylim, link=I, inverse=I, at=ticks$at, n=ticks$n)
      nm <- names(x)[1]
      x.vals <- x.data[, nm]   
      if (nm %in% names(ticks.x)){
        at <- ticks.x[[nm]]$at
        n <- ticks.x[[nm]]$n
      }
      else{
        at <- NULL
        n <- 5
      }
      xlm <- if (nm %in% names(xlim)){
        xlim[[nm]]
      }
      else range.adj(x[nm]) # range(x.vals)
      tickmarks.x <- if ((nm %in% names(transform.x)) && !(is.null(transform.x))){
        trans <- transform.x[[nm]]$trans
        make.ticks(trans(xlm), link=transform.x[[nm]]$trans, inverse=transform.x[[nm]]$inverse, at=at, n=n)
      }
      else {
        trans <- I
        make.ticks(xlm, link=I, inverse=I, at=at, n=n)
      }
      
      if (is.null(x.var)){
        if (!is.null(residuals)){
          x.var <- names(x)[1]
        }
        else x.var <-  which.max(levels)
      }
      if (!is.null(residuals)) x.fit <- x.data[, predictor]
      if (is.numeric(x.var)) x.var <- predictor
      plot <- xyplot(eval(parse(
        text=paste("fit ~ trans(", x.var, ")"))),
        strip=function(...) strip.default(..., strip.names=c(factor.names, TRUE)),
        panel=function(x, y, x.vals, rug, lower, upper, has.se, ...){
          if (grid) panel.grid()
          good <- !is.na(y)
          axis.length <- diff(range(x))
          effect.llines(x[good], y[good], lwd=lwd, col=colors[1], ...)
          if (rug && is.null(residuals)) lrug(trans(x.vals))
          if (has.se){  
            if (ci.style == "bars"){
              larrows(x0=x[good], y0=lower[good], 
                      x1=x[good], y1=upper[good], 
                      angle=90, code=3, col=eval(colors[.modc(2)]),
                      length=.125*cex/1.5)
            }
            else if(ci.style == "lines") {
              effect.llines(x[good], lower[good], lty=2, col=colors[.modc(2)])
              effect.llines(x[good], upper[good], lty=2, col=colors[.modc(2)])
            }
            else{ if(ci.style == "bands") {
              panel.bands(x[good], y[good], upper[good], lower[good], fill=band.colors[1],
                          alpha=band.transparency, use.splines=use.splines)
            }}
          }
          if (has.thresholds){
            panel.abline(h=thresholds, lty=3)
            panel.text(rep(current.panel.limits()$xlim[1], length(thresholds)), 
                       thresholds, threshold.labels, adj=c(0,0), cex=0.75)
            panel.text(rep(current.panel.limits()$xlim[2], length(thresholds)), 
                       thresholds, threshold.labels, adj=c(1,0), cex=0.75)
          }
          if (!is.null(residuals)){ 
            fitted <- y[good][closest(x.fit, x[good])]
            partial.res <- fitted + residuals
            lpoints(trans(x.fit), partial.res, col=residuals.color, pch=residuals.pch, cex=residuals.cex)
            if (show.fitted) lpoints(trans(x.fit), fitted, pch=16, col=residuals.color)  # REMOVE ME
            if (smooth.residuals){
              llines(loess.smooth(trans(x.fit), partial.res, span=span, family=loess.family), lwd=lwd, lty=2, col=residuals.smooth.color)
            }
          }
          
        },
        ylim=ylim,
        xlim=suppressWarnings(trans(xlm)),
        ylab=ylab,
        xlab=if (missing(xlab)) names(x)[1] else xlab,
        x.vals=x.vals, rug=rug,
        main=main,
        lower=x$lower, upper=x$upper, has.se=has.se, data=x, 
        scales=list(y=list(at=tickmarks$at, labels=tickmarks$labels, rot=roty),
                    x=list(at=tickmarks.x$at, labels=tickmarks.x$labels, rot=rotx), alternating=alternating), ...)
      result <- update(plot, layout = if (missing(layout)) c(0, prod(dim(plot))) 
                       else layout)
      result$split <- split
      result$more <- more
      class(result) <- c("plot.eff", class(result))
    }
    return(result)
  }
  ###  more than one variate
  predictors <- names(x)[1:n.predictors]
  levels <- sapply(apply(x[,predictors], 2, unique), length)
  if (is.null(x.var)){
    if (!is.null(residuals)){
      x.var <- names(x)[1]
    }
    else x.var <-  which.max(levels)
  }
  if (!is.null(residuals)) x.fit <- x.data[, x.var]
  if (is.character(x.var)) {
    which.x <- which(x.var == predictors)
    if (length(which.x) == 0) stop(paste("x.var = '", x.var, "' is not in the effect.", sep=""))
    x.var <- which.x
  }
  if (is.character(z.var)) {
    which.z <- which(z.var == predictors)
    if (length(which.z) == 0) stop(paste("z.var = '", z.var, "' is not in the effect.", sep=""))
    z.var <- which.z
  }    
  if (x.var == z.var) z.var <- z.var + 1
  ### multiline
  if (multiline){ 
    ci.style <- if(is.null(ci.style)) "none" else ci.style
    if(ci.style == "lines") { 
      cat("Confidence interval style 'lines' changed to 'bars'\n")
      ci.style <- "bars"}
    range <- if (has.se && ci.style !="none")
      range(c(x$lower, x$upper), na.rm=TRUE) else range(x$fit, na.rm=TRUE)
    ylim <- if (!missing(ylim)) ylim else c(range[1] - .025*(range[2] - range[1]),                                              
                                            range[2] + .025*(range[2] - range[1]))
    tickmarks <- if (type == "response") make.ticks(ylim, link=trans.link, 
                                                    inverse=trans.inverse, at=ticks$at, n=ticks$n)
    else make.ticks(ylim, link=I, inverse=I, at=ticks$at, n=ticks$n)
    zvals <- unique(x[, z.var])
    ### multiline factor
    if (is.factor(x[,x.var])){
      levs <- levels(x[,x.var])
      key <- list(title=predictors[z.var], cex.title=1, border=TRUE,
                  text=list(as.character(zvals)), 
                  lines=list(col=colors[.modc(1:length(zvals))], lty=lines[.modl(1:length(zvals))], lwd=lwd),
                  points=list(col=colors[.modc(1:length(zvals))], pch=symbols[.mods(1:length(zvals))]),
                  columns = if ("x" %in% names(key.args)) 1 else find.legend.columns(length(zvals)))
      for (k in names(key.args)) key[k] <- key.args[k]
      plot <- xyplot(eval(parse( 
        text=paste("fit ~ as.numeric(", predictors[x.var], ")",
                   if (n.predictors > 2) paste(" |", 
                                               paste(predictors[-c(x.var, z.var)])), collapse="*"))),
        strip=function(...) strip.default(..., strip.names=c(factor.names, TRUE)),
        panel=function(x, y, subscripts, z, lower, upper, show.se, ...){
          if (grid) panel.grid()
          for (i in 1:length(zvals)){
            sub <- z[subscripts] == zvals[i]
            good <- !is.na(y[sub])
            os <- if(show.se)
              (i - (length(zvals) + 1)/2) * (2/(length(zvals)-1)) * 
              .01 * (length(zvals) - 1) else 0
            effect.llines(x[sub][good]+os, y[sub][good], lwd=lwd, type='b', col=colors[.modc(i)],
                          pch=symbols[.mods(i)], lty=lines[.modl(i)], cex=cex, ...)
            if (show.se){
              larrows(x0=x[sub][good]+os, y0=lower[subscripts][sub][good], 
                      x1=x[sub][good]+os, y1=upper[subscripts][sub][good], 
                      angle=90, code=3, col=eval(colors[.modc(i)]),
                      length=.125*cex/1.5)
            }
          } 
          if (has.thresholds){
            panel.abline(h=thresholds, lty=3)
            panel.text(rep(current.panel.limits()$xlim[1], length(thresholds)), 
                       thresholds, threshold.labels, adj=c(0,0), cex=0.75)
            panel.text(rep(current.panel.limits()$xlim[2], length(thresholds)), 
                       thresholds, threshold.labels, adj=c(1,0), cex=0.75)
          }
        },        
        ylim=ylim,
        ylab=ylab,
        xlab=if (missing(xlab)) predictors[x.var] else xlab,
        z=x[,z.var],
        scales=list(x=list(at=1:length(levs), labels=levs, rot=rotx), 
                    y=list(at=tickmarks$at, labels=tickmarks$labels, rot=roty),
                    alternating=alternating),
        zvals=zvals,
        main=main,
        key=key,
        lower=x$lower, upper=x$upper, 
        show.se=has.se && ci.style=="bars", 
        data=x, ...)
      result <- update(plot, layout = if (missing(layout)) 
        c(0, prod(dim(plot))) else layout)
      result$split <- split
      result$more <- more
      class(result) <- c("plot.eff", class(result))
    } 
    ### multiline variate   
    else{
      effect.llines <- if(use.splines) spline.llines else effect.llines
      nm <- names(x)[x.var]
      x.vals <- x.data[, nm]   
      if (nm %in% names(ticks.x)){
        at <- ticks.x[[nm]]$at
        n <- ticks.x[[nm]]$n
      }
      else{
        at <- NULL
        n <- 5
      }
      xlm <- if (nm %in% names(xlim)){
        xlim[[nm]]
      }
      else range.adj(x[nm]) 
      tickmarks.x <- if ((nm %in% names(transform.x)) && !(is.null(transform.x))){
        trans <- transform.x[[nm]]$trans
        make.ticks(trans(xlm), link=transform.x[[nm]]$trans, inverse=transform.x[[nm]]$inverse, at=at, n=n)
      }
      else {
        trans <- I
        make.ticks(xlm, link=I, inverse=I, at=at, n=n)
      }
      key <- list(title=predictors[z.var], cex.title=1, border=TRUE,
                  text=list(as.character(zvals)), 
                  lines=list(col=colors[.modc(1:length(zvals))], lty=lines[.modl(1:length(zvals))], lwd=lwd),
                  columns = if ("x" %in% names(key.args)) 1 else find.legend.columns(length(zvals)))
      for (k in names(key.args)) key[k] <- key.args[k]
      plot <- xyplot(eval(parse( 
        text=paste("fit ~trans(", predictors[x.var], ")", 
                   if (n.predictors > 2) paste(" |", 
                                               paste(predictors[-c(x.var, z.var)])), collapse="*"))),
        strip=function(...) strip.default(..., strip.names=c(factor.names, TRUE)),
        panel=function(x, y, subscripts, x.vals, rug, z, lower, upper, show.se, ...){
          if (grid) panel.grid()
          if (rug && is.null(residuals)) lrug(trans(x.vals))
          axis.length <- diff(range(x))
          for (i in 1:length(zvals)){
            sub <- z[subscripts] == zvals[i]
            good <- !is.na(y[sub])
            effect.llines(x[sub][good], y[sub][good], lwd=lwd, type='l',
                          col=colors[.modc(i)], lty=lines[.modl(i)], cex=cex, ...)
            if(show.se){ 
              if(ci.style == "bars"){
                os <- (i - (length(zvals) + 1)/2) * (2/(length(zvals)-1)) * 
                  .01 * axis.length
                larrows(x0=x[sub][good]+os, y0=lower[subscripts][sub][good], 
                        x1=x[sub][good]+os, y1=upper[subscripts][sub][good], 
                        angle=90, code=3, col=eval(colors[.modc(i)]),
                        length=.125*cex/1.5)
              }
              if(ci.style == "bands"){ 
                panel.bands(x[sub][good], y[sub][good], 
                            upper[subscripts][sub][good], lower[subscripts][sub][good], 
                            fill=eval(band.colors[.modb(i)]),
                            alpha=band.transparency, use.splines=use.splines)
              }
            }
          }
          if (has.thresholds){
            panel.abline(h=thresholds, lty=3)
            panel.text(rep(current.panel.limits()$xlim[1], length(thresholds)), 
                       thresholds, threshold.labels, adj=c(0,0), cex=0.75)
            panel.text(rep(current.panel.limits()$xlim[2], length(thresholds)), 
                       thresholds, threshold.labels, adj=c(1,0), cex=0.75)
          }
        },
        ylim=ylim,
        xlim=suppressWarnings(trans(xlm)), 
        ylab=ylab,
        xlab=if (missing(xlab)) predictors[x.var] else xlab,
        x.vals=x.vals, rug=rug,
        z=x[,z.var],
        zvals=zvals,
        main=main,
        key=key, 
        #
        lower=x$lower, upper=x$upper, 
        show.se=has.se && ci.style %in% c("bars", "bands"),
        #
        data=x, scales=list(y=list(at=tickmarks$at, labels=tickmarks$labels),
                            rot=roty, x=list(at=tickmarks.x$at, labels=tickmarks.x$labels, rot=rotx), 
                            alternating=alternating),  ...)
      result <- update(plot, layout = if (missing(layout)) c(0, prod(dim(plot))) 
                       else layout)
      result$split <- split
      result$more <- more
      class(result) <- c("plot.eff", class(result))
    }
    return(result)
  } 
  # multiplot
  ci.style <- if(is.null(ci.style)){
    if(is.factor(x[, x.var])) "bars" else "bands"} else ci.style
  range <- if (has.se && ci.style !="none")
    range(c(x$lower, x$upper), na.rm=TRUE) else range(x$fit, na.rm=TRUE)
  # multiplot factor
  if (is.factor(x[,x.var])){
    
    ylim <- if (!missing(ylim)) ylim else c(range[1] - .025*(range[2] - range[1]),                                              
                                            range[2] + .025*(range[2] - range[1]))
    tickmarks <- if (type == "response") make.ticks(ylim, link=trans.link, 
                                                    inverse=trans.inverse, at=ticks$at, n=ticks$n)
    else make.ticks(ylim, link=I, inverse=I, at=ticks$at, n=ticks$n)  
    
    levs <- levels(x[,x.var])
    plot <- xyplot(eval(parse( 
      text=paste("fit ~ as.numeric(", predictors[x.var], ") |", 
                 paste(predictors[-x.var], collapse="*")))),
      strip=function(...) strip.default(..., strip.names=c(factor.names, TRUE)),
      panel=function(x, y, subscripts, lower, upper, has.se, ...){  
        if (grid) panel.grid()
        good <- !is.na(y)
        if (has.se){
          if (ci.style == "bars"){
            larrows(x0=x[good], y0=lower[subscripts][good], x1=x[good], y1=upper[subscripts][good], 
                    angle=90, code=3, col=colors[.modc(2)], length=0.125*cex/1.5)
          }
          else if(ci.style == "lines") {
            effect.llines(x[good], lower[subscripts][good], lty=2, col=colors[.modc(2)])
            effect.llines(x[good], upper[subscripts][good], lty=2, col=colors[.modc(2)])
          }
          else{ if(ci.style == "bands") { 
            panel.bands(x[good], y[good], upper[subscripts][good], lower[subscripts][good], 
                        fill=band.colors[1], alpha=band.transparency, use.splines=FALSE)
          }}
        }
        effect.llines(x[good], y[good], lwd=lwd, type='b', col=colors[1], pch=19, cex=cex, ...)
        if (has.thresholds){
          panel.abline(h=thresholds, lty=3)
          panel.text(rep(current.panel.limits()$xlim[1], length(thresholds)), 
                     thresholds, threshold.labels, adj=c(0,0), cex=0.75)
          panel.text(rep(current.panel.limits()$xlim[2], length(thresholds)), 
                     thresholds, threshold.labels, adj=c(1,0), cex=0.75)
        }
      },
      ylim=ylim,
      ylab=ylab,
      xlab=if (missing(xlab)) predictors[x.var] else xlab,
      scales=list(x=list(at=1:length(levs), labels=levs, rot=rotx), 
                  y=list(at=tickmarks$at, labels=tickmarks$labels, rot=roty),
                  alternating=alternating),
      main=main,
      lower=x$lower, upper=x$upper, has.se=has.se, data=x, ...)
    result <- update(plot, layout = if (missing(layout)) c(0, prod(dim(plot))) else layout)
    result$split <- split
    result$more <- more
    class(result) <- c("plot.eff", class(result))
  } 
  ### multiplot variate   ***
  else{
    effect.llines <- if(use.splines) spline.llines else effect.llines
    nm <- names(x)[x.var]
    x.vals <- x.data[, nm]   
    if (nm %in% names(ticks.x)){
      at <- ticks.x[[nm]]$at
      n <- ticks.x[[nm]]$n
    }
    else{
      at <- NULL
      n <- 5
    }
    xlm <- if (nm %in% names(xlim)){
      xlim[[nm]]
    }
    else range.adj(x[nm]) 
    tickmarks.x <- if ((nm %in% names(transform.x)) && !(is.null(transform.x))){
      trans <- transform.x[[nm]]$trans
      make.ticks(trans(xlm), link=transform.x[[nm]]$trans, inverse=transform.x[[nm]]$inverse, at=at, n=n)
    }
    else {
      trans <- I
      make.ticks(xlm, link=I, inverse=I, at=at, n=n)
    }
    ylim <- if (!missing(ylim)) ylim
    else if (is.null(residuals)) c(range[1] - .025*(range[2] - range[1]), range[2] + .025*(range[2] - range[1]))
    else c(min(partial.residuals.range[1], range[1] - .025*(range[2] - range[1])), 
           max(partial.residuals.range[2], range[2] + .025*(range[2] - range[1])))
    tickmarks <- if (type == "response") make.ticks(ylim, link=trans.link, 
                                                    inverse=trans.inverse, at=ticks$at, n=ticks$n)
    else make.ticks(ylim, link=I, inverse=I, at=ticks$at, n=ticks$n)  
    x.fit <- x.data[, predictors[x.var]]
    use <- rep(TRUE, length(residuals))
    xx <- x[, predictors[-x.var], drop=FALSE]
    plot <- xyplot(eval(parse( 
      text=paste("fit ~ trans(", predictors[x.var], ") |", 
                 paste(predictors[-x.var], collapse="*")))),
      strip=function(...) strip.default(..., strip.names=c(factor.names, TRUE)),
      panel=function(x, y, subscripts, x.vals, rug, lower, upper, has.se, ...){
        if (grid) panel.grid()
        good <- !is.na(y)
        effect.llines(x[good], y[good], lwd=lwd, col=colors[1], ...)
        if (rug && is.null(residuals)) lrug(trans(x.vals))
        if (has.se){  
          if (ci.style == "bars"){ 
            larrows(x0=x[good], y0=lower[subscripts][good], 
                    x1=x[good], y1=upper[subscripts][good], 
                    angle=90, code=3, col=eval(colors[.modc(2)]),
                    length=.125*cex/1.5)
          }
          else if(ci.style == "lines") {
            effect.llines(x[good], lower[subscripts][good], lty=2, col=colors[.modc(2)])
            effect.llines(x[good], upper[subscripts][good], lty=2, col=colors[.modc(2)])
          }
          else if(ci.style == "bands") { 
            panel.bands(x[good], y[good], upper[subscripts][good], lower[subscripts][good], 
                        fill=band.colors[1], alpha=band.transparency, use.splines=use.splines)
          }
          if (!is.null(residuals)){
            predictors <- predictors[-x.var]
            factors <- sapply(xx, is.factor)
            for (predictor in predictors){
              use <- use & if(factors[predictor]) x.all[, predictor] == xx[subscripts[1], predictor]
              else x.all[, predictor] == xx[subscripts[1], predictor]
            }
            n.in.panel <- sum(use)
            if (n.in.panel > 0){
              fitted <- y[good][closest(x.fit[use], x[good])]
              partial.res <- fitted + residuals[use]
              lpoints(trans(x.fit[use]), partial.res, col=residuals.color, pch=residuals.pch, cex=residuals.cex)
              if (show.fitted) lpoints(trans(x.fit[use]), fitted, pch=16, col=residuals.color)  # REMOVE ME
              if (smooth.residuals && n.in.panel >= 10) {
                llines(loess.smooth(x.fit[use], partial.res, span=span, family=loess.family), 
                       lwd=lwd, lty=2, col=residuals.smooth.color)
              }
            }
          }
        }
        if (has.thresholds){
          panel.abline(h=thresholds, lty=3)
          panel.text(rep(current.panel.limits()$xlim[1], length(thresholds)), 
                     thresholds, threshold.labels, adj=c(0,0), cex=0.75)
          panel.text(rep(current.panel.limits()$xlim[2], length(thresholds)), 
                     thresholds, threshold.labels, adj=c(1,0), cex=0.75)
        }
      },
      ylim=ylim,
      xlim=suppressWarnings(trans(xlm)),
      ylab=ylab,
      xlab=if (missing(xlab)) predictors[x.var] else xlab,
      x.vals=x.vals, rug=rug,
      main=main,
      lower=x$lower, upper=x$upper, has.se=has.se, data=x, 
      scales=list(y=list(at=tickmarks$at, labels=tickmarks$labels, rot=roty),
                  x=list(at=tickmarks.x$at, labels=tickmarks.x$labels, rot=rotx), 
                  alternating=alternating), ...)
    result <- update(plot, layout = if (missing(layout)) c(0, prod(dim(plot))) else layout)
    result$split <- split
    result$more <- more
    class(result) <- c("plot.eff", class(result))
  }
  return(result)
}

print.plot.eff <- function(x, ...){
  NextMethod(split=x$split, more=x$more, ...)
  invisible(x)
}

plot.efflist <- function(x, selection, rows, cols, ask=FALSE, graphics=TRUE, ...){
  if (!missing(selection)){
    if (is.character(selection)) selection <- gsub(" ", "", selection)
    return(plot(x[[selection]], ...))
  }
  effects <- gsub(":", "*", names(x))
  if (ask){
    repeat {
      selection <- menu(effects, graphics=graphics, title="Select Term to Plot")
      if (selection == 0) break
      else print(plot(x[[selection]], ...))
    }
  }
  else {
    neffects <- length(x)
    mfrow <- mfrow(neffects)
    if (missing(rows) || missing(cols)){
      rows <- mfrow[1]
      cols <- mfrow[2]
    }
    for (i in 1:rows) {
      for (j in 1:cols){
        if ((i-1)*cols + j > neffects) break
        more <- !((i-1)*cols + j == neffects)
        print(plot(x[[(i-1)*cols + j]], row=i, col=j, nrow=rows, ncol=cols, more=more, ...))
      }
    }
  }
}

