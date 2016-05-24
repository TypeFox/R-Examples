panel.xy <- function(x, y, type = "p,r,r.pred,r=,N=,grid",
                     r.min=0.5, level=0.8, slope=0.0, intercept=0.0,
                     unicolor=FALSE, ...) {
  #---------- nested function panel.lminterval: ----------
  panel.lminterval = function(x, y, 
                              lmplot=TRUE, lmconfplot=FALSE, 
                              lmpredplot=TRUE, lqsplot=FALSE,
                              unicolor=FALSE, verbose=FALSE, level, ...) {
    # Cannot be used directly because of missing missing value handling.
    rl = trellis.par.get("regression.line")
    if (!(lmplot|lmconfplot|lmpredplot|lqsplot)) return() # nothing to do
    #--- lm line:
    lm.fit = lm(y~x, data=data.frame(list(x=x, y=y)))
    coefs = coef(lm.fit)
    if (lmplot) {
      if (unicolor)
        panel.abline(a=coefs[1], b=coefs[2], ...)
      else
        panel.abline(a=coefs[1], b=coefs[2], 
                     col=rl$col.lm, lty=rl$lty.lm, lwd=rl$lwd.lm)
    }
    if (verbose) {
      cat(paste(c("   lm:  y =", "+"), signif(coefs[1:2], 3)), "x\n")
    }
    #--- conf and pred lines:
    N = length(x)
    s <- sum(lm.fit$residuals^2)/(N-2)  # variance of the fit
    xmean <- mean(x)                 # mean of the x values
    xsumsq <- sum((x-xmean)^2)   # sum of squares of x values
    if (lmconfplot || lmpredplot) {
      xrange <- grid::current.viewport()$xscale  # undocumented!!!
      # or use: convertX(unit(0:1, "npc"), "native", valueOnly=TRUE)
      xpoints <- seq(xrange[1], xrange[2], length.out=50)
      ypoints <- coefs[2]*xpoints + coefs[1]
      tconf <- qt(1-(1-level)/2, df=N-2)
    }
    if (lmconfplot) {        #-- lm confidence interval
      yint <- tconf * sqrt(s) * sqrt(1/N + (xpoints-xmean)^2 / xsumsq)
      if (unicolor) {
        llines(xpoints, ypoints + yint, ...)
        llines(xpoints, ypoints - yint, ...)
      } else {
        llines(xpoints, ypoints + yint, 
               col=rl$col.conf, lwd=rl$lwd.conf, lty=rl$lty.conf)
        llines(xpoints, ypoints - yint,
               col=rl$col.conf, lwd=rl$lwd.conf, lty=rl$lty.conf)
      }
    }
    if (lmpredplot) {        #-- lm prediction interval
      yint <- tconf * sqrt(s) * sqrt(1 + 1/N + (xpoints-xmean)^2 / xsumsq)
      if (unicolor) {
        llines(xpoints, ypoints + yint, ...)
        llines(xpoints, ypoints - yint, ...)
      } else {
        llines(xpoints, ypoints + yint, 
               col=rl$col.pred, lwd=rl$lwd.pred, lty=rl$lty.pred)
        llines(xpoints, ypoints - yint,
               col=rl$col.pred, lwd=rl$lwd.pred, lty=rl$lty.pred)
      }
    }
    #--- lqs line:
    if (lqsplot) {
      lqs.fit = MASS::lqs(y~x, data=data.frame(list(x=x, y=y)))
      #lqs.fit = MASS::lqs(x, y, method="lqs", 
      #                    control=list(nsamp="best",adjust=TRUE), ...)
      lqs.intercept = lqs.fit$coef[1]
      lqs.slope = lqs.fit$coef[2]
      if (unicolor)
        panel.abline(b=lqs.slope, a=lqs.intercept, ...)
      else {
        panel.abline(b=lqs.slope, a=lqs.intercept,
                     col.line=rl$col.lqs, lty=rl$lty.lqs, lwd=rl$lwd.lqs)
      }
      if (verbose) {
        cat(paste(c("  lqs:  y =", "+"), signif(c(intercept,slope), 3)), "x\n")
      }
    }
  }
  #---------- end of nested function panel.lminterval ----------
  ok = !is.na(x) & !is.na(y) # remove all pairs having missing values
  x = x[ok]  
  y = y[ok]
  N = length(x)
  if (length(type) == 1) type = strsplit(type, ",")[[1]]
  if (length(intersect(c("r","r.pred","r.conf","r=","v"),type)) > 0)
    corr = cor(x, y)
  rl = trellis.par.get("regression.line")
  #--- type == g or grid:
  if ("grid" %in% type || "g" %in% type) {
    panel.grid(h=-1, v=-1, col="grey50", lty=3)
    type = setdiff(type, c("g","grid"))
  }
  #--- type == abline:
  if ("abline" %in% type) {
    # caller must provide slope= and intercept=
    a.l = trellis.par.get("ab.line") # nonstandard extension for panel.xy
    panel.abline(a=intercept, b=slope, col=a.l$col, lty=a.l$lty, lwd=a.l$lwd)  
    type = setdiff(type, "abline")
  }
  #--- type == rug:
  if ("rug" %in% type) {
    panel.rug(x, y, ...)
    type = setdiff(type, "rug")
  }
  #--- type == r, r.conf, r.pred, or lqs:
  regressions = c("r", "r.conf", "r.pred", "lqs")
  if (length(intersect(regressions, type))) {
    do.lm = N >= 8 && r.min <= 1 && abs(corr) >= r.min
    panel.lminterval(x, y,
                     lmplot=(do.lm &
                             length(intersect(c("r","r.conf","r.pred"),type))),
                     lmconfplot=(do.lm & ("r.conf" %in% type)),
                     lmpredplot=(do.lm & ("r.pred" %in% type)),
                     lqsplot=("lqs" %in% type),
                     level=level,
                     unicolor=unicolor, verbose=("v" %in% type), ...)
    type = setdiff(type, regressions)
  }
  #--- type == smooth, loess, or lowess:
  loesses = c("smooth", "loess", "lowess")
  if (length(intersect(loesses, type))) {
    if (N >= 8) {
      if (unicolor)
        panel.loess(x, y, ...)
      else
        panel.loess(x, y,
                    col=rl$col.loess, lty=rl$lty.loess, lwd=rl$lwd.loess, ...)
    }
    type = setdiff(type, loesses)
  }
  #--- type == v:
  if ("v" %in% type) {
    cat(paste(c("median x=", " y=", " mean x=", " y=",
                " sdev x=", " y=", " cor(x,y)="),
              signif(c(median(x), median(y), mean(x), mean(y), 
                       sqrt(var(x)), sqrt(var(y)), corr), 3), sep=""), 
        "\n", sep="")
    type = setdiff(type, "v")
  }
  #--- type == 'r=' or 'N=':
  texts = c("r=", "N=")
  if (length(intersect(texts, type))) {
    result = character(0)
    if ("r=" %in% type)
      result = c(result, paste("r=", signif(corr, 2), sep=""))
    if ("N=" %in% type)
      result = c(result, paste("N=", length(x), sep=""))
    grid.text(label=paste(paste(result, collapse=", "), " \n\n"),
              x = unit(1, "npc"), y = unit(0, "npc"), just="right")
    type = setdiff(type, texts)
  }
  #--- type == p, l, b, o, s, S, h:
  panel.xyplot(x, y, type, ...)
  type = setdiff(type, c("p","l", "b", "o", "s", "S", "h"))
  #--- type == anything else:
  if (length(type) > 0)
    warning("panel.xy: non-standard values in type/type: '",
            paste(type,collapse=","), "'")
}
