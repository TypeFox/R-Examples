###########################################################################
## PROGRAM: plot.pred.int.R
## PURPOSE: Defines a plot function for an object of type "pred.int".
## INPUT: 
##   x: Object of type pred.int containing the PIPs data for plotting
##   conf.int: Print confidence interval for observed data?  TRUE (default)/FALSE
##   vline: Vector of x-values for vertical lines to print on the graph. These may
##    represent superiority/inferiority bounds or other x-values of interest
##   which: Only create graphs for some of the comparisons
##   axes:  TRUE/FALSE. Print axes on graph?  Default is TRUE.  Probably
##    you should only suppress axes now if you will add them later (for example 
##    if you don't like the default axes)
##   pi.col.fun: An optional one parameter function that takes a number between
##    0 and 100 and returns a color.  This can be used to color the predicted 
##    intervals different colors.  The input argument is the position on the 
##    vertical axis of the graph.  Default coloring is three shades of gray. 
##    (0-10) and (90-100) are light gray, (10-25) and (75-90) are darker, and 
##    (25-75) is darkest
##   ci.col: Color for effect estimate and confidence interval for the observed 
##    data.  Default is 2 (second color in palette)
##   main: Main title of graph.  If blank, a default will be used.  Can be either
##    a single title or a vector of titles.  If a vector, the first will be used
##    for the first graph, the second for the second graph, etc...  If any title
##    contains the string "#BY#", this will be replaced with the name of the
##    comparison (i.e. "B vs A")
##   xlab: Label on xaxis.  If blank, a default will be used.
##   ylab: Label on xaxis.  If blank, a default will be used.
##   xlim: Limits of xaxis (as vector of length 2).  Default is large enough to 
##    contain the predicted intervals and confidence interval.  Limits narrower than
##    the defaults will be ignored.
##   ...: Other options will be passed through to the plot.default function. For
##    example, you can set the color of the main title to blue by passing
##    col.main="blue"
## OUTPUT: Draws a figure
## MACROS USED:   None 
## CALLED BY:     None
## AUTHOR:        Daniel Muenz / Ray Griner
## CREATION DATE: November, 2011 (started in 2010)
## NOTES:         
## MODIFICATIONS:  
## [RG20120107] Ray Griner standardized stop message
## [RG20120201] Ray Griner added support for ci.col and pi.col.fun
## [RG20120830] Ray Griner standardized header
###########################################################################
plot.pred.int <- function(x, conf.int=TRUE, vline=NA, which=NULL, axes=TRUE,
    pi.col.fun=NULL, ci.col=2, main, xlab, ylab, xlim, ...)
{
    ## If the study group doesn't have a name, give it a temporary
    ## made-up one.
    if (fake.name <- is.null(names(x$pi))) {
        names(x$pi) <- "FAKE"
        names(x$ci) <- "FAKE"
    }

    ## Location of default vertical line.
    if (x$data.type %in% c("t.test","binary"))
        if (length(x$obs.n) == 1)
            def.line <- NA
        else
            def.line <- 0
    else
        def.line <- 1

    ## Which predicted intervals should we plot?  Default is all.
    if (is.null(which))
        plots <- names(x$pi)
    else {
        if (is.character(which)) plots <- which
        else plots <- names(x$pi)[which]
        if (any(bad <- which==0) || any(bad <- !(plots %in% names(x$pi))))
            stop("illegal value(s) in 'which': ",
                 paste(which[bad], collapse=", "))
    }

    ## Use default title if none passed in input.
    if (missing(main) && ( (x$conf.level==x$obs.conf.level) || (conf.int==FALSE)) )
       { mymain<-paste(x$conf.level*100, "% Predicted Interval Plot", if (!fake.name) " for #BY#", sep="") }
    else if (missing(main))
       { mymain<-paste("Predicted Interval Plot", if (!fake.name) " for #BY#", sep="") }
    else mymain<-main

    mymain <- rep(mymain, length.out=length(plots))
    names(mymain) <- plots

    ## Loop through all the comparison groups.  Each one gets a plot.
    for(comp in plots) {

        point <- x$pi[[comp]][,"point"]
        lower <- x$pi[[comp]][,"lower"]
        upper <- x$pi[[comp]][,"upper"]
        ## Determine the ranges for both axes.
        myxlim <- c(min(lower), max(upper))
        ylim <- c(0, 100)

        if (conf.int)
            myxlim <- c(min(myxlim[1], x$ci[[comp]]["lower"]),
                        max(myxlim[2], x$ci[[comp]]["upper"]))

        myxlim <- c(min(myxlim[1], min(vline), na.rm=TRUE),
                    max(myxlim[2], max(vline), na.rm=TRUE))
        if (!missing(xlim)) {
          if (length(xlim)!=2) stop ("'xlim' should have length=2") ## [RG20120107]
          myxlim <- c(min(myxlim[1],xlim[1]),max(myxlim[2],xlim[2]))
        }

        if (length(x$obs.n) > 1) {
            if (x$data.type %in% c("t.test","binary"))
                comp2 <- gsub(" vs ", " - ", comp)
            else
                comp2 <- gsub(" vs ", " / ", comp)
            comp2 <- paste("(",comp2,")",sep="")
        }
        else
            comp2 <- ""
        
        xlab <- if (missing(xlab)) paste("Effect Size", comp2)
                else xlab
        ylab <- if (missing(ylab)) "Percentile of Point Estimate Distribution"
                else ylab
        thismain<-gsub("#BY#",comp,mymain[[comp]],ignore.case=TRUE)

        ## Start a new plot and add axes and a frame
        ## plot.new()
        plot.default(NULL, xlab=xlab, ylab=ylab, xlim=myxlim, ylim=ylim, 
                     main=thismain, axes=FALSE, ...)
        if (axes) {
          axis(1)
          axis(2, at=c(0,10,25,50,75,90,100), las=1)
        }
        box(bty="l")
        ##axis(1, at=seq(from=-1,to=1,by=.05), las=1)
        ##axis(2, at=seq(from=0,to=100,by=5), las=1)


        ## Get cumulative distribution function of the point estimates
        ## But the ecdf will return the same value when the point
        ## estimate is the same, so do we really want to use this as the
        ## y-axis in the plot?  If we do, it means we might have fewer lines
        ## than simulations, because some are plotted right on top of each other
        ################################################################
        cidist    <- function(y) { 100*ecdf(point)(y) }
        negcidist <- function(y) { 100-100*ecdf(-point)(-y) }
        dists  <- seq(from=0, to=100, length.out=length(point))

        ## ecdf(point)(y) is the proportion of predicted interval estimates <= y.
        ## but what if the ecdf of points with value 1 is from 10 pct - 90 pct. 
        ## Then really we want ecdf(1) to be 50 pct, not 90 pct.  This is the 
        ## reason we take the average of cidist and negcidist later
 
        ## Only use three shades of gray
        ## mygrays<-(cidist(point)/100)
        ## mycols<-mygrays
        ## mygrays[dists<10 | dists>90]<-gray(.8)
        ## mygrays[dists>=10 & dists<=90]<-gray(.55)
        ## mygrays[dists>=25 & dists<=75]<-gray(.3)
        ## [RG20120201]S Add support for pi.col.fun
        if (is.null(pi.col.fun)) {
          pi.col.fun <- function (x) {
            if (x<10 | x>90)   { retval<-gray(.8) }
            if (x>=10 & x<=90) { retval<-gray(.55) }
            if (x>=25 & x<=75) { retval<-gray(.3) }
            retval
          }  
        }
        ## [RG20120201]E 

        ## Predicted intervals
        segments(lower, dists, upper, dists, col=sapply(dists, pi.col.fun)) # [RG20120201] 
        points(x=point, y=dists)

        ## Confidence interval from observed data
        if (conf.int) {
            ycidist<-( cidist(x$ci[[comp]]["point"]) + negcidist(x$ci[[comp]]["point"]))/2
            segments(x$ci[[comp]]["lower"], ycidist,
                     x$ci[[comp]]["upper"], ycidist,
                     col=ci.col, lwd=2)  # [RG20120201] 
            points(x=x$ci[[comp]]["point"], ycidist,
                   col=ci.col, lwd=2, cex=3) # [RG20120201]
        }

        ## Vertical reference line(s)
        if (!is.null(vline)) {
          abline(v=vline, lty=2, lwd=2)
        }
        ## abline(h=c(10,25,75,90))  ## useful for debugging

        ## Always have a line through 0 or 1
        if (!is.null(def.line)) abline(v=def.line)
    }

}
