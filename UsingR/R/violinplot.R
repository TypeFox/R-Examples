## courtesy of boxplot from base

##' violinplot
##'
##' See the \pkg{vioplot} package
##' @param x data
##' @param ... passed on
##'
##' @rdname violinplot
##' @export
violinplot <- function(x, ...) UseMethod("violinplot")


##' violinplot
##'
##' @param x data
##' @param ... passed on
##' @param orientation orientation
##' @param bw bandwidth
##' @param names names
##' @param pars pars
##' @return NULL
##'
##' @rdname violinplot
##' @export
violinplot.default <-
  function(x,...,
           orientation = "vertical", bw = "nrd0",
           names = NULL,
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

    ## grab the global parameters
    ## n- number of plots, N[i] number in each sample,
    ## 


    ## work on the gropu by group level
    xvals <- matrix(0,nrow=512,ncol=n)
    yvals <- matrix(0,nrow=512,ncol=n)
    center<- 1:n                        # where are they centered
    for(i in 1:n) {
      tmp.dens <- density(groups[[i]],bw=bw)
      xvals[,i] <- tmp.dens$x
      yvals.needtoscale <- tmp.dens$y
      ## we scale so largest size is less than 1/2
      yvals.scaled <- 7/16 * yvals.needtoscale/max(yvals.needtoscale)
      yvals[,i] <- yvals.scaled
    }

    ## now plot
    ## need to first make the plot range
    ## this depends on horizontal or vertical
    if(orientation == "vertical") {
      xrange <- c(1/2,n+1/2)            # each gets 1 unit centered on integers
      yrange <- range(xvals)
    } else {                            #horizontal
      xrange <- range(xvals)
      yrange <- c(min(yvals),max(yvals))
    }
    plot.new()
    plot.window(xlim = xrange, ylim = yrange)

    for (i in 1:n)
      vlnplt(xvals[,i],yvals[,i],center[i],
             bordercolor = rainbow(i),
             bgcolor = rainbow(n-i),
             orientation = orientation, ...)
    
    axis(1,at=1:n,labels = names)
    axis(2)


}                                     # end of default

##' formula interface for violinplot
##'
##' @param formula formula
##' @param data data frame
##' @param ... passed on
##' @param subset for subsetting data
##' @return NULL
##'
##' @rdname violinplot
##' @export
violinplot.formula <- function(formula, data = NULL, ..., subset)
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
    violinplot(split(mf[[response]], mf[-response]), ...)
}

## Make a simple violin plot call from violinplot. values are x,y to plot
"vlnplt" <-
  function(x,y,center,
           add=TRUE,
           orientation="horizontal",
           bgcolor=NA , bordercolor='red',
           ...) {

    ## double up first
    x <- c(x,rev(x))
    y <- c(y,-rev(y))
    y <- y + center

    
    if (orientation == "vertical") {
      # swtich x and y
      tmp=x;x=y;y=tmp;
    }

    if(add) {
      polygon(x,y,...)
    }

  }

