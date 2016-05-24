plot.HOF <- function (
		x, 
		marginal =  c('bar', 'rug', 'hist', 'points', 'n'), 
		boxp = TRUE, 
		las.h = 1, 
		yl, 
		main, 
		model, 
		test = "AICc", 
		modeltypes, 
		onlybest = TRUE, 
		penal, 
		para = FALSE, 
		gam.se = FALSE, 
		color, 
		newdata = NULL, 
		lwd=1, 
		leg = TRUE, 
		add=FALSE, 
		xlabel, 
		...)    {
  resp <- x
  ow <- options("warn")
  options(warn = -1)
  yaxt = TRUE
  independ <- resp$x 
  depend <- resp$y 
  if(missing(marginal)) {
    dit <- ifelse(para, 'n', if(length(independ) > 200) 'hist' else  'bar') 
    if(resp$family != 'binomial') dit <- 'points'
    } else dit <-  match.arg(marginal)
  if(missing(modeltypes)) modeltypes <- eHOF.modelnames else onlybest = FALSE
  if(missing(penal)) penal <- 'df'
  if(missing(model)) {
	    model <- pick.model(resp, modeltypes = modeltypes, penal = penal, test = test, gam=FALSE, ...)
  }
  cols <- if(missing(color)) c("black", "red", "green", "blue", "sienna", "violet", "pink") else color
#  cols <- c("black", "red", "green", "blue", "sienna", "violet", "pink")
  #if(length(cols) != length(modeltypes)) cols <- rep(cols, length.out=length(modeltypes))
  if(missing(xlabel)) xlabel <- resp$x.name
  if (missing(main)) main <- resp$y.name
  if (missing(yl)) yl <- c(0, 1)
  if(missing(yaxt)) yaxt <- 's'
  if (is.null(newdata)) newdata <- seq(min(resp$x), max(resp$x), length.out=10000)
    
  logi.box <- function(independ, depend, col.box = "grey", x.lab = xlabel, yrange = yl, ylab= NULL, las = las.h, multi=FALSE, ...) {
        if(is.null(ylab)) ylab <- if(resp$M == 1) 'Predicted probability' else 'Response'
	plot(independ, depend, ylim = c(yrange[1] - diff(yrange)/10, yrange[2] + diff(yrange)/10), yaxp = c(0, yrange[2], 5), type = "n", ylab = ylab, xlab = x.lab, yaxt = 'n', ...)
	ax <- axTicks(2)
 	if(yaxt!='n') axis(2, at=ax, labels= eval(ax*resp$M), ...)
      indep.1 <- independ[depend > 0]
      indep.0 <- independ[depend == 0]
      boxplot(indep.1, horizontal = TRUE, add = TRUE, at = (yrange[2] + diff(yrange)/20), boxwex = diff(yrange)/10, col = col.box, notch = TRUE, axes = FALSE, ...)
      axis(1, mgp=c(2,.5,0))
        boxplot(indep.0, horizontal = TRUE, add = TRUE, at = 0 - diff(yrange)/20, boxwex = diff(yrange)/10, col = col.box, notch = TRUE, axes = FALSE, ...)
       if(para & leg) title(main, adj = 0, ...) else title(main, ...)
    }

  logi.scatter <- function(independ, depend, x.lab = xlabel, yrange = yl, las = las.h, ylab=NULL, ...) {
        if(is.null(ylab)) ylab <- if(resp$M == 1) 'Predicted probability' else 'Response'
        plot(independ, depend / resp$M, type = 'n', ylim = yrange, yaxp = c(0, yrange[2], 5), ylab = ylab, xlab = x.lab, las = las, yaxt='n', ...)
	ax <- axTicks(2)
 	if(yaxt!='n') axis(2, at=ax, labels= eval(ax*resp$M))
        if(para) title(main, adj = 0, ...) else title(main, ...)
   }

  logi.rug <- function(independ, depend, cex.rug = 1, yrange = yl, ...) {
    # p <- ifelse(depend > 0, 2, 1)
   	points(independ, (depend > 0 )* yrange[2], pch = '|', cex = cex.rug, ...)
   }

  logi.points <- function(independ, depend, pch.p = c(1, 19), cex.rug = .7, yrange = yl, ...) {
    p <- ifelse(depend > 0, 2, 1) #pch = pch.p[p],
    points(independ, depend / resp$M, cex = cex.rug, pch = pch.p[p], ...)
  }
  
  logi.bar <- function(independ, depend, incre = 0.02, cex.p = .5, pch.rug = c(1, 19), yrange = yl, ...) {
     p <- depend > 0 ; p[p] <- 2; p[!p] <- 1
     points(independ, (depend > 0 )* yrange[2], pch = pch.rug[p], cex = cex.p, ...)
     indep.0 <- independ[depend == 0]
     indep.1 <- independ[depend > 0]
     uni.plot.0 <- function(x) length(which(indep.0 == x))
     uni.plot.1 <- function(x) length(which(indep.1 == x))
     cosa.0 <- apply(as.matrix(unique(indep.0)), 1, uni.plot.0)
     cosa.1 <- apply(as.matrix(unique(indep.1)), 1, uni.plot.1)
     for (i in 1:max(cosa.0)) {
         for (j in 1:i) points(unique(indep.0)[which(cosa.0 == i + 1)], 
               rep(0 + incre * j, length(which(cosa.0 == i + 1))), pch = 1, cex = cex.p)
     }
     for (i in 1:max(cosa.1)) {
         for (j in 1:i) points(unique(indep.1)[which(cosa.1 == i + 1)], 
               rep(1 - incre * j, length(which(cosa.1 == i + 1))) - (1 - yrange[2]), pch = 19, cex = cex.p)
     }
}  

  logi.hist <- function(independ, depend, scale.hist = 5, yrange = yl, col.hist = gray(0.9), bord.hist = gray(0.8), count.hist = FALSE,  interval = 0, las.h1 = las.h, ...) {
      h.br <- hist(independ, plot = FALSE)$br
      if (interval > 0) 
          h.br <- seq(from = range(h.br)[1], to = range(h.br)[2], y = interval)
      h.x <- hist(independ[depend == min(depend)], breaks = h.br, plot = FALSE)$mid
      h.y0 <- hist(independ[depend == min(depend)], breaks = h.br, plot = FALSE)$counts
      h.y1 <- hist(independ[depend > min(depend)], breaks = h.br, plot = FALSE)$counts
      h.y0n <- h.y0/(max(c(h.y0, h.y1)) * scale.hist)
      h.y1n <- 1 - h.y1/(max(c(h.y0, h.y1)) * scale.hist)
      for (i in 1:length(h.y0n)) {
          if (h.y0n[i] > 0) 
              polygon(c(rep(h.br[i], 2), rep(h.br[i + 1], 2)), 
                c(0, rep(h.y0n[i], 2), 0), col = col.hist, border = bord.hist)
      }
      for (i in 1:length(h.y1n)) {
          if (h.y1n[i] < 1) 
            polygon(c(rep(h.br[i], 2), rep(h.br[i + 1], 2)), 
              c(yrange[2] - (1 - h.y1n[i]), yrange[2], #resp$M
                yrange[2], yrange[2] - (1 - h.y1n[i])), #resp$M
              col = col.hist, border = bord.hist)
      }
      if (count.hist == TRUE) 
        for (i in 1:length(h.x)) {
         text(h.x[i], yrange[2] - (1 - h.y1n[i]), h.y1[i], cex = 1, pos = 1, col='grey')
         text(h.x[i], h.y0n[i], h.y0[i], cex = 1, pos = 3, col='grey') #
        }
      axis.hist <- function(h.y0, h.y1, scale.hist, las = las.h1) {
          tope.0 <- max(h.y0)
          tope.1 <- max(h.y1)
          label.down<- c(0, (ceiling(tope.0/10)) * 5, (ceiling(tope.0/10)) * 10)
          label.up  <- c((ceiling(tope.1/10)) * 10,  0)
          at.down <- label.down/(tope.0 * scale.hist)
          at.up <- max(yl) - label.up/(tope.0 * scale.hist)
     if(onlybest) { axis(side = 4, at = c(at.down, at.up), col.lab = bord.hist, col.axis = bord.hist, labels = c(label.down, label.up), las = las)
          mtext("Frequency", col = bord.hist, side = 4, line = 2, ...)
          } else {
          axis(side = 4, at = at.down, col.lab = bord.hist, col.axis = bord.hist, labels = label.down, las = las)
          } 
      }
      axis.hist(h.y0, h.y1, scale.hist)
  }

 logi.curve <- function(resp, cex.l = .8, ...) {
   x <- sort(newdata) 
   if (onlybest) {
     cols <- cols[match(model, eHOF.modelnames)]
     fv <- predict(resp, model, newdata=x)
     linewd <- lwd
     lines(x, fv, col =  cols[1], lwd = linewd, ...)
   } else {
     fv <- matrix(ncol=length(modeltypes), nrow=length(newdata))
     for(i in 1:length(modeltypes)) 
         fv[,i] <- predict(resp, model=modeltypes[i], newdata=x)
     # fv <- sweep(fv, 1, resp$M, "/")

	  if(length(cols) == length(eHOF.modelnames)) cols <- cols[match(modeltypes, eHOF.modelnames)]
	  if(length(lwd)==1) linewd <- rep(lwd, length(modeltypes))
	  matlines(x, fv, lty = 1, col = cols, lwd = linewd, ...)
	}
  par(xpd = TRUE)
  if(leg) {
       legend(par("usr")[2] + 0.05, par("usr")[4], if(onlybest) model else rev(modeltypes), ncol = 1, bty = "n", col = if(onlybest) cols else rev(cols[match(modeltypes, eHOF.modelnames)]), lty = 1, lwd = rev(linewd), title="Model", cex = cex.l)
       if (para) legend((par("usr")[2]), par("usr")[4], legend=c("optimum", "expectancy value", "raw mean", "inflection points", "outer niche", "central niche"), title='Response parameters', ncol=2, 
                        pch=c(1,NA,NA,1,15,15), 
                        lty = c(0,2,1,0,0,0), col = c(4, 2, 1, 'orange', 'lightgrey', 'grey'), bg = "transparent", xjust = 1, yjust = 0, lwd = c(1,2,2,1,NA,NA), cex = cex.l)
      }
  invisible()
  }

para.fun <- function(resp, cex.pl = .8, ...) {
    if(dit=='hist') warning("It is not recommended to use marginal='hist' together with para='TRUE'.")
    p <- Para(resp, model, ...)
    tl <- yl[2] - yl[2]/15

## Niche
  x <- c(p$outerBorder[1],rep(p$outerBorder[2], 2),p$outerBorder[1])
  y <- rep(c(yl[2] - yl[2]/35, yl[2]),each=2)
  polygon(x,y, col=grey(.9), border = NA)

  x <- c(p$centralBorder[1],rep(p$centralBorder[2], 2),p$centralBorder[1])
  y <- rep(c(yl[2] - yl[2]/20, yl[2]),each=2)
  polygon(x,y, col=grey(.85), border = NA)

  if(model %in% c('VI','VII')) {
    x <- c(p$outerBorder[3],rep(p$outerBorder[4], 2),p$outerBorder[3])
    y <- rep(c(yl[2] - yl[2]/35, yl[2]),each=2)
    polygon(x,y, col=grey(.9), border = NA)
    x <- c(p$centralBorder[3],rep(p$centralBorder[4], 2),p$centralBorder[3])
    y <- rep(c(yl[2] - yl[2]/20, yl[2]),each=2)
    polygon(x,y, col=grey(.85), border = NA)
  }
## Top
    lines(c(par('usr')[1], p$range[1]+diff(p$range)/100), rep(p$top[1], length.out=2))

## Range
    for(i in 1:length(p$expect)) lines(rep.int(p$expect[i], 2), c(yl[2], tl), col = 2, lty=2, lwd = 2)

## Raw mean
    lines(rep.int(p$raw.mean, 2), c(yl[2] + yl[2]/15, yl[2] - yl[2]/15), col = 1, lwd = 2)

## Optima
   if(model %in% c('VI','VII'))  {
    	lines(rep(p$opt['opt1'],2), c(tl, yl[2]), col = 4, lwd = 2)
    	lines(rep(p$opt['opt2'],2), c(tl, yl[2]), col = 4, lwd = 2)
    	points(p$opt, predict(resp, model, newdata=p$opt), col='blue')
    } else { 
      if(length(p$opt)==1) lines(c(p$opt, p$opt), c(tl, yl[2]), col = 4, lwd = 2) else
    	lines(p$opt, c(yl[2], yl[2]), col='blue', lwd = 2)
      points(p$opt, predict(resp, model, newdata=p$opt), col='blue')
    }
## Pessimum
#  points(p$pess, predict(resp, model, newdata=p$pess), col='blue')

## Inflection
  if(length(p$inflection) > 0)
     for(n in 1:length(p$inflection)) {
       y <- predict(resp, model, p$inflection[n], ...)
       points(p$inflection, predict(resp, model, newdata=p$inflection), col='orange')
     }
}
  
 gam.conf <- function(independ, depend, bs = 'cr', k = 4, ...) {
    if (!is.null(newdata)) warning('GAM confidence intervalls only calculated for original predictor values!')
    fit <- mgcv::gam(depend ~ s(independ, bs=bs, k=k),family=get(resp$family), scale=0, ...)
    pred.fit <- mgcv::predict.gam(fit, type='response', se.fit=TRUE)
    o <- order(independ)
    ll <- pred.fit$fit - 2*pred.fit$se.fit; ll[ll<0] <- 0
    ul <- pred.fit$fit + 2*pred.fit$se.fit; ul[ul>resp$M] <- resp$M
    lines(independ[o], ll[o]/resp$M , lty=2)
    lines(independ[o], ul[o]/resp$M , lty=2)
    lines(independ[o], pred.fit$fit[o]/resp$M , lty=3)
 }

## Call plot functions
  old.mar <- par()$mar
  if(leg) {
  	par(mar=par()$mar + c(0,0,2))
  	if(para) par(mar=par()$mar + c(0,0,0,3))
    }
  if(!add) {
    if (boxp == TRUE) logi.box(independ, depend, xaxt='n', ...) else  logi.scatter(independ, depend, ...)
    if (para)  para.fun(resp, xaxt='n', ...)
    switch(dit,
       	hist =logi.hist(independ, depend, ...),
       	bar = logi.bar(independ, depend, ...),
       	rug = logi.rug(independ, depend, ...),
        points = logi.points(independ, depend, ...)
        )
  }
  if(gam.se) gam.conf(independ, depend, ...)
  logi.curve(resp, xaxt='n', ...)
  par(mar = old.mar)
  options(ow) # reset
}
