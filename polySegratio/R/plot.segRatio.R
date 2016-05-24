`plot.segRatio` <-
function(x, main=deparse(substitute(x)), xlab="",
                           xlab.segRatio="Segregation ratio",
                           xlab.nobs="Number of dominant markers",
                           xlab.miss="Number of missing markers per individual",
                           NCLASS=100,
                           type=c("seg.ratio","all","no","missing"), ...)
{

  ## Description: Plots object of class segRatio

  ## Arguments:
  ## seg.ratio : object of class segRatio containing counts and
  ##             segregation ratios
  ## main: title for plots
  ## xlab.segRatio: x axis label when plotting segregation ratios
  ## xlab.nobs: x axis label when plotting no. of 1's
  ## xlab.miss:  x axis label when plotting number of missing
  ##             individuals per marker
  ## NCLASS: number of classes for histograms
  ##         Default: 100
  ## type: type of plot "all","seg.ratio","no","missing"
  ##       Deafult: "seg.ratio"

  ## NB: could rewrite this to get a plot of missing markers for
  ## individual as well as "missing" which is individuals per marker
  ## but this would mean adding axtra stuff to clas 'segRatio'
  
  ## Values:
  ## Used for its side-effects
  
  type <- match.arg(type)

  if (type == "no") hist(x$r, breaks=NCLASS, main = main,
        xlab=xlab.nobs, ...)
  
  if (type == "seg.ratio") {
    hist(x$seg.ratio, breaks=NCLASS, main=main, xlab=xlab.segRatio, ...)
  }
  
  if (type == "missing") {
    hist(x$n.individuals - x$n, main=main, xlab=xlab.miss,...)
  }

  if (type == "all") {
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))
    layout(matrix(c(1,2,3,4), nrow=2, ncol=2))
    hist(x$r, breaks=NCLASS, main = main, xlab=xlab.nobs, ...)
    hist(x$seg.ratio, breaks=NCLASS, main=main, xlab=xlab.segRatio, ...)
    hist(x$n.individuals - x$n, main=main, xlab=xlab.miss,...)
    
  }
}

