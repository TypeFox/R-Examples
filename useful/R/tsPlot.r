#' vplayout
#' 
#' Viewport
#' 
#' Creates viewport for pushing ggplot objects to parts of a console.
#' 
#' @author Jared P. Lander
#' @aliases vplayout
#' @export vplayout
#' @return An R object of class viewport.
#' @param x The x cell of the viewport to push into.
#' @param y The y cell of the viewport to push into.
#' @examples
#' 
#' library(ggplot2)
#' library(grid)
#' 
vplayout <- function(x, y)
{
    grid::viewport(layout.pos.row=x, layout.pos.col=y)
}



#' fortify.ts
#' 
#' Fortify a ts object.
#' 
#' Prepares a ts object for plotting with ggplot.
#' 
#' @author Jared P. Lander
#' @aliases fortify.ts
#' @export
#' @return \code{\link{data.frame}} for plotting with ggplot.
#' @param model A \code{\link{ts}} object.
#' @param data A vector of the same length of \code{x} that specifies the time component of each element of \code{x}.
#' @param name Character specifying the name of x if it is to be different that the variable being inputed.
#' @param \dots Further arguments.
#' @examples
#' 
#' fortify(sunspot.year)
#' 
fortify.ts <- function(model, data=NULL, name=as.character(m[[2]]), ...)
{
    m <- match.call()
    
    # if time is provided use that as the x values
    if(!is.null(data))
    {
        theX <- data
    }else
        # otherwise use the built in attributes
    {
        theTime <- attr(model, which="tsp")
        theX <- seq(from=theTime[1], to=theTime[2], by=1/theTime[3])
        rm(theTime)
    }
    
    data <- data.frame(theX, model)
    names(data) <- c("Time", name)
    return(data)
}



#' ts.plotter
#' 
#' Plot a ts object
#' 
#' Fortifies, then plots a \code{\link{ts}} object.
#' 
#' @export ts.plotter
#' @aliases ts.plotter
#' @author Jared P. Lander
#' @import ggplot2
#' @return A ggplot object
#' @param data A \code{\link{ts}} object to be plotted.
#' @param time A vector of the same length of \code{data} that specifies the time component of each element of \code{data}.
#' @param title Title of plot.
#' @param xlab X-axis label.
#' @param ylab Y-axis label.
#' @examples
#' 
#' ts.plotter(sunspot.year)
#' 
ts.plotter <- function(data, time=NULL, 
                       title="Series Plot", xlab="Time", ylab="Rate")
{
    # grab the name of the ts that was provided    
    # fortify the ts so it is usable in ggplot
    data <- fortify.ts(data, data=time, name=as.character(match.call()[[2]]))
    
    # figure out the names returned by fortifying
    x <- names(data)[1]
    y <- names(data)[2]
    
    # build the plot
    ggplot(data, aes_string(x=x, y=y)) + geom_line(aes(group=1)) +
        labs(title=title, x=xlab, y=ylab)
}


#' fortify.acf
#'
#' Fortify an acf/pacf object
#' 
#' Prepares acf (and pacf) objects for plotting with ggplot.
#' 
#' @author Jared P. Lander
#' @aliases fortify.acf
#' @export
#' @return \code{\link{data.frame}} for plotting with ggplot. 
#' @param model An \code{\link{acf}} object.
#' @param data Not used.  Just for consistency with the fortify method.
#' @param \dots Other arguments
#' @examples
#' 
#' fortify(acf(sunspot.year, plot=FALSE))
#' fortify(pacf(sunspot.year, plot=FALSE))
#' 
fortify.acf <- function(model, data=NULL, ...)
{
    # the different tpe of acf objects
    theNames <- c(correlation="ACF", covariance="ACF", partial="Partial.ACF")
    
    # build a data.frame consisting the lag number and the acf value
    data <- data.frame(model$lag, model$acf)
    
    # name the data "Lag" and the appropriate type of acf
    names(data) <- c("Lag", theNames[model$type])
    
    return(data)
}


#' plot.acf
#' 
#' Plot acf objects
#' 
#' Plot acf (and pacf) objects.
#' 
#' @author Jared P. Lander
#' @aliases plot.acf
#' @export plot.acf
#' @export
#' @return A ggplot object.
#' @param x An \code{\link{acf}} object.
#' @param xlab X-axis label.
#' @param ylab y-axis label.
#' @param title Graph title.
#' @param \dots Further arguments.
#' @examples
#' 
#' plot(acf(sunspot.year, plot=FALSE))
#' plot(pacf(sunspot.year, plot=FALSE))
#'
plot.acf <- function(x, 
                     xlab=x, ylab=y,
                     #xlab=x, ylab=sub("\\.", " ", y), 
                     title=sprintf("%s Plot", y), ...
                     #title=sprintf("%s Plot", sub("\\.", " ", y))
                     )
{
    # fortify the acf object
    data <- fortify.acf(x)
    
    # get the names we are using
    x <-names(data)[1]
    y <- names(data)[2]
    
    # build plot
    ggplot(data, aes_string(x=x)) + 
        geom_linerange(aes_string(ymin=pmin(y, 0), ymax=pmax(y, 0))) +
        labs(title=title, x=xlab, y=ylab)
}

#' @title plot
#' @description Overwitten plot generic so that plot.acf can be defined in this package
#' @details Overwitten plot generic so that plot.acf can be defined in this package
#' @author Jared P. Lander
#' @param x Object to be plotted
#' @param \dots Further arguments
#' @return A plot
#' 
plot <- function(x, ...)
{
    UseMethod('plot')
}

#' plotTimesSeries
#' 
#' Plot ts object
#' 
#' Plot a ts object and, if desired, it's acf and pacf.
#' 
#' @aliases plot.times.series
#' @author Jared P. Lander
#' @export plotTimesSeries
#' @importFrom stats na.fail
# @S3method plot ts
# @method plot ts
#' @return A ggplot object if \code{acf} is \code{FALSE}, otherwise \code{TRUE} indicating success.
#' @param x a \code{\link{ts}} object.
#' @param time A vector of the same length of \code{x} that specifies the time component of each element of \code{x}.
#' @param acf Logical indicating if the acf and pacf should be plotted.
#' @param lag.max maximum lag at which to calculate the acf. Default is 10*log10(N/m) where N is the number of observations and m the number of series. Will be automatically limited to one less than the number of observations in the series.
#' @param na.action function to be called to handle missing values. na.pass can be used.
#' @param demean logical. Should the covariances be about the sample means?
#' @param xlab X-axis label.
#' @param ylab Y-axis label.
#' @param title Graph title.
#' @param \dots Further arguments.
#' @seealso ts.plotter plot.acf fortify.ts
#' @examples
#' 
#' plot(sunspot.year)
#' plot(sunspot.year, acf=TRUE)
#' 
plotTimesSeries <- function(x, time=NULL, acf=FALSE,
                    lag.max=NULL, na.action=na.fail, demean=TRUE, 
                    title=sprintf("%s Plot", name), xlab="Time", ylab=name, ...)
{
    # get real name of x
    name <- as.character(match.call()[[2]])
    
    # build ts plot
    tsPlot <- ts.plotter(data=x, time=time, title=title, xlab=xlab, ylab=ylab)
    
    # if we're not doing acf just return the ts plot
    if(!acf)
    {
        return(tsPlot)
    }
    
    # calculate the acf/pacf
    theAcf <- acf(x=x, lag.max=lag.max, na.action=na.action, demean=demean, type="correlation", plot=FALSE)
    thePacf <- acf(x=x, lag.max=lag.max, na.action=na.action, demean=demean, type="partial", plot=FALSE)
    
    # build plots for acf/pacf
    acfPlot <- plot.acf(theAcf, title=NULL)
    pacfPlot <- plot.acf(thePacf, title=NULL)
    
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(layout=grid::grid.layout(nrow=2, ncol=2)))
    
    print(tsPlot, vp=vplayout(1, 1:2))
    print(acfPlot, vp=vplayout(2, 1))
    print(pacfPlot, vp=vplayout(2, 2))
    
    #list(tsPlot, acfPlot, pacfPlot)
    invisible(TRUE)
}
