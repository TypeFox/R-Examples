"run.hist.demo" <-
function(x) {

  if(!requireNamespace('tcltk', quietly=TRUE)){stop('The tcltk package is needed')}

  pr <- pretty(x)
  xr <- range(pr)
  xr[1] <- 4*xr[1] - 3*min(x)
  xr[2] <- 4*xr[2] - 3*max(x)

  hist.refresh <- function(...) {

    hist(x,seq( slider(no=2), slider(no=3), length=slider(no=1)+1),
         xlim=xr)

  }

  slider(hist.refresh, c('Number of bins','Minimum','Maximum'),
         c(1, xr[1], max(x)),
         c(length(x),min(x),xr[2]),
         c(1, (min(x)-xr[1])/50, (xr[2]-max(x))/50),
         c(nclass.Sturges(x),min(pr),max(pr)),
         title="Histogram Demo")

}


## # create new version using tkrplot changing min, max, and nbins, include rug
##
## hist.demo <- function(x,xmin,xmax,n,xlab=deparse(substitute(x))) {
##     br <- seq(xmin,xmax, length.out=n+1)
##  print(range(x))
##     hist(x,br, xlab=xlab, main='')
## }
##
##
## run.hist.demo <- function(xx,...) {
##   if(!require(tkrplot)) { stop('The tkrplot package is needed')}
##
##   xlab <- deparse(substitute(xx))
##
##   pr <- pretty(xx)
##   h1 <- hist(xx, plot=FALSE)
##   plist <- list( xmin =  list('spinbox', from=min(pr), to=min(x),
##                  increment=(min(x)-min(pr))/10 ),
##                 xmax = list('spinbox', from=max(x), to=max(pr),
##                  increment=(max(pr)-max(x))/10, init=max(pr) ),
##                 n = list('slider', from=1, to=length(x), resolution=1,
##                  init=length(h1$breaks)-1))
##   tkexamp( hist.demo(xx,xlab=xlab),  plist )
## }
##
##
