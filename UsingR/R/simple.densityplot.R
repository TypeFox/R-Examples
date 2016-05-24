## courtesy of boxplot from base
## plot empirical densities for datasets
## useful for comparing distributions. See also violinplot

## "x" can be a dataframe, a formula x ~ f, or a sequence of variable
## names as in boxplot.

##' make a simple density plot
##'
##' @param x data
##' @param ... passed on
##' @return NULL
##'
##' @rdname simple.densityplot
##' @export
simple.densityplot <- function(x, ...) UseMethod("simple.densityplot")

##' default method for simple densityplot
##'
##' @param x data
##' @param ... passed on
##' @param bw bandwidth
##' @param do.legend do legend?
##' @param names names
##' @param pars pars
##' @return NULL
##'
##' @rdname simple.densityplot
##' @export
simple.densityplot.default <-
  function(x,...,
           bw = "nrd0",
           do.legend = "auto",
           names=FALSE,
           pars = NULL
           ) 
  {
    args <- list(x, ...)
    namedargs <-
      if(!is.null(attributes(args)$names))
        attributes(args)$names != ""
      else
        rep(FALSE, length = length(args))
    pars <- c(args[namedargs], pars)
    groups <- if(is.list(x)) x else args[!namedargs]
    if(0 == (n <- length(groups)))
      stop("invalid first argument")
    if(length(class(groups)))
	groups <- unclass(groups)
    if(!missing(names))
      attr(groups, "names") <- names
    else {
      if(is.null(attr(groups, "names")))
        attr(groups, "names") <- 1:n
      names <- attr(groups, "names")
    }

    ## work on the gropu by group level
    xvals <- matrix(0,nrow=512,ncol=n)
    yvals <- matrix(0,nrow=512,ncol=n)
    for(i in 1:n) {
      tmp.dens <- density(groups[[i]],bw=bw)
      xvals[,i] <- tmp.dens$x
      yvals[,i] <- tmp.dens$y
    }

    ## Now plot
    xrange <- range(xvals)
    yrange <- range(yvals)              #all of them

    plot.new()
    plot.window(xlim = xrange, ylim = yrange)

    for (i in 1:n)
      dnstyplt(xvals[,i],yvals[,i],
               lty=i,
               ...)
    
    axis(1)
    axis(2)

    ## add a legend
    if (do.legend != FALSE) {
      if (do.legend == "auto") {
        lambda = c(1,3)/4
        x.where = lambda %*% xrange;y.where = lambda %*% yrange
        legend(x.where,y.where,names,lty=1:n)
      } else {
        ## need to click
        print("click in graph to place legend")
        names
        legend(locator(1),names,lty=1:n)
      }
    }
}                                     # end of default

##' simple density plot with a formula interface
##' 
##' @param formula formula
##' @param data data frame
##' @param ... passed on
##' @param subset subset data
##'
##' @rdname simple.densityplot
##' @export
simple.densityplot.formula <- function(formula, data = NULL, ..., subset)
{
    if(missing(formula) || (length(formula) != 3))
        stop("formula missing or incorrect")
    m <- match.call(expand.dots = FALSE)
    if(is.matrix(eval(m$data, parent.frame())))
        m$data <- as.data.frame(data)
    m$... <- NULL
    m[[1]] <- as.name("model.frame")
    mf <- eval(m, parent.frame())
    response <- attr(attr(mf, "terms"), "response")
    simple.densityplot(split(mf[[response]], mf[-response]), ...)
}

## Make a simple density plot call from densityplot. values are x,y to plot
"dnstyplt" <-
  function(x,y,center,
           add=TRUE,
           ...) {

    lines(x,y,...)

  }

