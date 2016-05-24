plot.HOF.list <- function (
	x, 
	plottype = c('layout', "lattice", "all"), 
	xlabel = NULL, 
	test = "AICc",
	modeltypes,
	border.top = 0.1,
	color,
	yl,
	leg = FALSE,
	...) {
# ow <- options('warn')  
# options(warn = 10)
ncol = 4
plottype <- match.arg(plottype) 
cols <- if(missing(color)) c("black", "red", "green", "blue", "sienna", "violet", "pink") else color

if(missing('modeltypes')) modeltypes <- eHOF.modelnames
xlabel <- if(is.null(xlabel)) x[[1]]$x.name else xlabel
mods <- x
N <- length(mods)
nobs <- length(mods[[1]]$x)

fitfun <- function(x, test, modeltypes,...) fitted(x, model = pick.model(x, test = test, modeltypes=modeltypes, gam=FALSE, ...))/x$M

   if (plottype == "layout") {
  	M <- mods[[1]]$M
  	minresp <- if(missing(yl)) 0 else yl[1]
    maxresp <- 
  	if(missing(yl)) max(sapply(mods, function(x) max(x$models[[pick.model(x, test = test, modeltypes, gam=FALSE)]]$fitted )/M))  else yl[2]
    layoutfun <- function(mods, N, mar=NULL, ...) {
  	  if(is.null(mar)) mar <- if(N < 30) c(2,2,2,0) else c(1,0,0,0)
  	  autolayout(N)
  	  par(mar=mar)
  	  for(i in 1:N) plot(mods[[i]], test=test, leg=leg, yl=c(minresp,maxresp), color=cols,  ...)
  	}
    layoutfun(mods, N=N, ...)
  }

  if (plottype == "lattice") {
      Response <- unlist(lapply(mods, function(x) x$y/x$M))
      Gradient <- unlist(lapply(mods, function(x) x$x))
      Species <- rep(names(mods), each = nobs)
      Fit <- unlist(lapply(mods, fitfun, test, modeltypes))
      mod <- sapply(mods, pick.model, test = test, modeltypes, gam=FALSE,  ...)
      mod <- rep(mod, each = nobs)
      fit.panel <- function(x, y, subscripts, Fit, ...) {
	  panel.xyplot(x, y, ...)
	  i <- order(x)
	  fv <- Fit[subscripts]
	  sp <- unique(cbind(x[i], fv[i]))
	  panel.xyplot(sp[, 1], sp[, 2], type = "l", lwd = 4, col = cols[match(mod[min(subscripts)], modeltypes)], ...)
      }
      mykey <- list(text = list(text = modeltypes), lines = list(lty = 1, 
	  col = cols[match(modeltypes, eHOF.modelnames)]), columns = length(modeltypes))
      out <- xyplot(Response ~ Gradient | Species, xlab = mods[[1]]$x.name, 
	  Fit = Fit, key = mykey, panel = fit.panel)
      return(out)
  }
 if (plottype == "all") {
      lplot <- function( ..., xlim = c(min(grad), max(grad)), ylim=c(0,m), ylab = "Predicted probability", xlab = xlabel, type = "n", para) {
        if(!missing(para)) message('Option "para" is available only for plottype "layout".')
        plot(...)
      } 
      m <- max(sapply(mods, function(x) fitted(x, model=pick.model(x, gam=FALSE, test = test, ...)))/mods[[1]]$M)
      grad <- mods[[1]]$x
#     plot(x=0, xlim = c(min(grad), max(grad)), ylim=c(0,m), ylab = "Predicted probability", xlab = xlabel, type = "n")
    lplot(mods[[1]]$x, mods[[1]]$y, ...)
      at <- order(grad)
      models <- pick.model(object = mods, test = test, modeltypes = modeltypes, gam=FALSE, ...)
      for (i in 1:N) {
#        par(xpd = TRUE)
         lines(grad[at], fitfun(mods[[i]], test, modeltypes)[at], lty=(1:N)[i], col = cols[match(models[i], eHOF.modelnames)])
        }
        legtext <- paste(names(mods), "(", sapply(mods, function(x, ...) pick.model(x,modeltypes, gam=FALSE, test = test, ...)), ")")
       if(leg) legend(par("usr")[1], par("usr")[4] + border.top, legtext, ncol = ncol, bty = "n", fill = rainbow(N))
    }
# options(ow)
}
