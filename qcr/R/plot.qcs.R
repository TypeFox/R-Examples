#-------------------------------------------------------------------------
# plot.qcs
#-------------------------------------------------------------------------
##' function to create a plotting 'qcs' object
##' 
##' Generic function for plotting Shewhart charts of object of class 'qcs' to perform statistical 
##' quality control.
##' 
##' @method plot qcs
##' @param x  Object qcs (Quality Control Statical)
##' @param title an overall title for the plot
##' @param subtitle a sub title for the plot
##' @param xlab a title for the x axis
##' @param ylab a title for the y axis
##' @param ylim the y limits of the plot
##' @param center.nominal a value specifying the 
##' center of group statistics or the "target" value 
##' of the process
##' @param limits.specification a two-values vector 
##' specifying control limits.
##' @param limits.alert a two-values vector 
##' specifying control alert limits.
##' @param label.index logical. If TRUE label index is plotted
##' @param type.data  a string specifying el type de data.
##' @param conf.nsigma  a numeric value used to compute control limits, specifying the
##' number of standard deviations (if \code{conf.nsigma} > 1) or the confidence level (if 0
##' < \code{conf.nsigma} < 1).
##' @param conf.nsigma.alert  a numeric value used to compute control limits, specifying the
##' number of standard deviations (if \code{conf.nsigma} > 1) or the confidence level (if 0
##' < \code{conf.nsigma} < 1).
##' @param ...  arguments to be passed to or from methods.
##' @export
##' 
plot.qcs <- function(x, title, subtitle, xlab, ylab, ylim, center.nominal = NULL, 
                     limits.specification = NULL, limits.alert = NULL, label.index = NULL, 
                     type.data  =  c("continuous", "atributte", "dependence"), ...)
#.........................................................................                     
{
  
  type.data <- match.arg(type.data)
  
  oldpar <- par(mar = c(5, 4, 4, 3) + 0.1)
  
  limits <- x$limits

  if (is.null(ylim)){
    if (!is.null(limits.specification)){
      dist <- max(abs(x$center - max(x$statistics[[1]])),
                  abs(x$center - min(x$statistics[[1]])),
                  abs(x$center - limits.specification[1]),
                  abs(x$center - limits.specification[2]),
                  abs(x$center - limits[1]),
                  abs(x$center - limits[2]))        
    } else{      
      dist <- max(abs(x$center - max(x$statistics[[1]])),
                  abs(x$center - min(x$statistics[[1]])),
                  abs(x$center - limits[1]),
                  abs(x$center - limits[2]))        
    }
    ylim <- x$center + c(-dist,dist)
  }

  
  sample <- rownames(x$statistics)
  if (inherits(x, "qcs.cusum")){
    plot(x$pos ~ sample, type =  "n",pch  =  16, axes  =  FALSE, 
       main  =  title, sub  =  subtitle, xlab  = xlab, 
       ylab  =  ylab, ylim  =  ylim)
  } else {
    plot(x$statistics[[1]] ~ sample, type =  "n",pch  =  16, axes  =  FALSE, 
         main  =  title, sub  =  subtitle, xlab  = xlab, 
         ylab  =  ylab, ylim  =  ylim)
    
  }
  axis(1, at  =  sample, cex.axis  =  0.7)
  axis(2, cex.axis  =  0.7)

  if (inherits(x, "qcs.cusum")){
  axis(4, at = c(0, max(x$limits), min(x$limits)), 
       labels = c("CL","LCL","UCL"), adj = 0, las = 1)
  } else{
    axis(4, at = c(x$center, max(x$limits), min(x$limits)), 
         labels = c("CL","LCL","UCL"), adj = 0, las = 1)
    
  }
    
  rect(par("usr")[1],
       par("usr")[3],
       par("usr")[2],
       par("usr")[4],
       col  =  "#CCCCCC")
  box(col  =  "#CCCCCC")
  grid(col  =  "#EEEEEE")
 
  
if (type.data == "dependence"){  
  if (inherits(x, "qcs.ewma")){  
    points(sample, x$statistics[[1]], pch = 3, cex = 0.8)
    lines(x$y ~ sample, type = "o", pch=20)
    abline(h  =  x$center, lwd  =  2, col  =  "red")
  } else {
    lines(x$pos ~ sample, type = "o", pch=20)
    lines(x$neg ~ sample, type = "o", pch=20)  
    abline(h = 0, lwd = 2)
  }
  
} else  {  
    lines (x$statistics[[1]] ~ sample, type  =  "b",pch  =  16)
    abline(h  =  x$center, lwd  =  2, col  =  "red")
  }   
  
  if (length(x$limits) == 2){
    lcl <- x$limits[1]
    ucl <- x$limits[2]
    
    abline(h  =  lcl, lwd  =  2, col  =  "red", lty = 2)
    abline(h  =  ucl, lwd  =  2, col  =  "red",lty = 2) 
    
  } else{
    lcl <- x$limits[,1]
    ucl <- x$limits[,2]
    
    lines(lcl, lwd  =  2, col  =  "red", lty = 2, type="s")
    lines(ucl, lwd  =  2, col  =  "red",lty = 2, type="s") 
    
  }
  
  
  
  
  m<-length(x$statistics)
  if (!is.null(label.index) & (m>1)){
    etiqueta <- x$statistics[[label.index+1]]
    text(label  =  etiqueta, x  =  sample, y  =  x$statistics[[1]],
         pos  =  4, cex = 0.7)
  }
  
  
  if (type.data != "dependence") {
    beyond.limits <- x$violations[[1]]
    runs.limits <- x$violations[[2]]
    
    points(x$statistics[[1]][beyond.limits]~sample[beyond.limits],
           col = "red",
           pch  =  19)
    
    points(x$statistics[[1]][runs.limits] ~ sample[runs.limits],
           col = "orange",
           pch  =  19)

    if (!is.null(limits.alert)){
      lcl.alert <- limits.alert[1]
      abline(h  =  lcl.alert, lwd  =  2, col  =  "yellow", lty = 2)
      
      ucl.alert <- limits.alert[2]
      abline(h  =  ucl.alert, lwd  =  2, col  =  "yellow",lty = 2) 
    }
    
    if (!is.null(center.nominal)){  
      abline(h  = center.nominal, lwd  =  2, col  =  "green")
    }
    
    if (!is.null(limits.specification)){  
      lel <- limits.specification[1]
      abline(h  =  lel, lwd  =  2, col  =  "green", lty = 2)
      
      uel <- limits.specification[2]
      abline(h  =  uel, lwd  =  2, col  =  "green",lty = 2) 
    }
    
    
  }  
  
  
  par(oldpar)
  #.........................................................................
} # plot.qcs
#.........................................................................

#-------------------------------------------------------------------------
# plot.qcs.xbar
#-------------------------------------------------------------------------
##' @rdname  plot.qcs
##' @method plot qcs.xbar
##' @export
##' 
plot.qcs.xbar <- function(x, title  =  NULL, subtitle  =  NULL, xlab  =  NULL,
                          ylab  =  NULL, ylim  =  NULL, conf.nsigma.alert  =  NULL,
                          center.nominal  =  NULL, limits.specification  =  NULL,
                          label.index  =  NULL, ...)
#.........................................................................                     
  {
  
  m<-length(x$statistics)
  
  if (!is.null(label.index))
    if (m-1 < label.index) stop("the covariable index is out of range")
  
  if(!is.null(conf.nsigma.alert)){
    limits.alert <- limits.xbar(x$center, x$std.dev, x$sizes,
                                conf.nsigma.alert)
  } else {
    limits.alert <- NULL
  }
  
  if (is.null(title)) title <- x$data.name
  
  if (is.null(subtitle)) subtitle <- expression(paste("Chart of control ", bar(x)," "))
  if (is.null(xlab)) xlab <- "Sample"
  if (is.null(ylab)) ylab <- expression(bar(x))
  
  plot.qcs(x, title , subtitle, xlab, ylab, ylim, center.nominal,
           limits.specification, limits.alert, label.index)
} #plot.qcs.xbar
#.........................................................................                     

#-------------------------------------------------------------------------
# plot.qcs.S
#-------------------------------------------------------------------------
##' @rdname  plot.qcs
##' @method plot qcs.S
##' @export
##' 
plot.qcs.S <- function(x, title  =  NULL, subtitle  =  NULL, xlab  =  NULL,
                       ylab  =  NULL, ylim  =  NULL, conf.nsigma.alert  =  NULL,
                       center.nominal  =  NULL, limits.specification  =  NULL,
                       label.index  =  NULL, ...)
#.........................................................................                       
  {
  
  m<-length(x$statistics)
  
  if (!is.null(label.index))
    if (m-1 < label.index) stop("the covariable index is out of range")
  
  if(!is.null(conf.nsigma.alert)){
    limits.alert <- limits.S(x$center, x$std.dev, x$sizes,
                             conf.nsigma.alert)
  } else {
    limits.alert <- NULL
  }
  
  if (is.null(title)) title <- x$data.name
  
  if (is.null(subtitle)) subtitle <- "Chart of control S"
  if (is.null(xlab)) xlab <- "Sample"
  if (is.null(ylab)) ylab <- "S"
  
  plot.qcs(x, title , subtitle, xlab, ylab, ylim, center.nominal,
           limits.specification, limits.alert, label.index)
} #plot.qcs.S
#.........................................................................                     


#-------------------------------------------------------------------------
# plot.qcs.R
#-------------------------------------------------------------------------
##' @rdname  plot.qcs
##' @method plot qcs.R
## @param conf.nsigma  a numeric value used to compute control limits, specifying the
##' @export
##' 
plot.qcs.R <- function(x, title  =  NULL, subtitle  =  NULL, xlab  =  NULL,
                       ylab  =  NULL, ylim  =  NULL, conf.nsigma.alert  =  NULL,
                       center.nominal  =  NULL, limits.specification  =  NULL,
                       label.index  =  NULL, ...)
#.........................................................................                     
  {
  
  m<-length(x$statistics)
  
  if (!is.null(label.index))
    if (m-1 < label.index) stop("the covariable index is out of range")
  
  if(!is.null(conf.nsigma.alert)){
    limits.alert <- limits.R(x$center, x$std.dev, x$sizes,
                             conf.nsigma.alert)
  } else {
    limits.alert <- NULL
  }
  
  if (is.null(title)) title <- x$data.name
  
  if (is.null(subtitle)) subtitle <- "Chart of control R"
  if (is.null(xlab)) xlab <- "Sample"
  if (is.null(ylab)) ylab <- "R"
  
  plot.qcs(x, title , subtitle, xlab, ylab, ylim, center.nominal,
           limits.specification, limits.alert, label.index)
} #plot.qcs.R
#.........................................................................                     

#-------------------------------------------------------------------------
# plot.qcs.one
#-------------------------------------------------------------------------
##' @rdname  plot.qcs
##' @method plot qcs.one
##' @export
##' 
plot.qcs.one <- function(x, title  =  NULL, subtitle  =  NULL, xlab  =  NULL, 
                         ylab  =  NULL, ylim  =  NULL, conf.nsigma.alert  =  NULL, 
                         center.nominal  =  NULL, limits.specification  =  NULL, 
                         label.index  =  NULL, ...)
#.........................................................................                     
  {
  
  m<-length(x$statistics)
  
  if (!is.null(label.index))
    if (m-1 < label.index) stop("the covariable index is out of range")
  
  if(!is.null(conf.nsigma.alert)){
    limits.alert <- limits.xbar.one(x$center, x$std.dev, x$sizes, 
                                    conf.nsigma.alert)
  } else {                                
    limits.alert <- NULL
  }
  
  if (is.null(title)) title <- x$data.name
  
  if (is.null(subtitle)) subtitle <- "Chart of control Xbar.one"
  if (is.null(xlab)) xlab <- "Sample"
  if (is.null(ylab)) ylab <- "Xbar.one"
  
  plot.qcs(x, title , subtitle, xlab, ylab, ylim, center.nominal, 
           limits.specification, limits.alert, label.index)
} #plot.qcs.one
#.........................................................................                     

#-------------------------------------------------------------------------
# plot.qcs.p
#-------------------------------------------------------------------------
##' @rdname  plot.qcs
##' @method plot qcs.p
##' @export
##' 
plot.qcs.p <- function(x, title  =  NULL, subtitle  =  NULL, xlab  =  NULL,
                       ylab  =  NULL, ylim  =  NULL, conf.nsigma.alert  =  NULL,
                       center.nominal  =  NULL, limits.specification  =  NULL,
                       label.index  =  NULL, ...)
#.........................................................................                     
  {
  
  m<-length(x$statistics)
  
  if (!is.null(label.index))
    if (m-1 < label.index) stop("the covariable index is out of range")
  
  if(!is.null(conf.nsigma.alert)){
    limits.alert <- limits.xbar(x$center, x$std.dev, x$sizes,
                                conf.nsigma.alert)
  } else {
    limits.alert <- NULL
  }
  
  if (is.null(title)) title <- x$data.name
  
  if (is.null(subtitle)) subtitle <- expression(paste("Chart of control p"," "))
  if (is.null(xlab)) xlab <- "Sample"
  if (is.null(ylab)) ylab <- "p"
  
  plot.qcs(x, title , subtitle, xlab, ylab, ylim, center.nominal,
           limits.specification, limits.alert, label.index)
} #plot.qcs.p
#.........................................................................                     


#-------------------------------------------------------------------------
# plot.qcs.np
#-------------------------------------------------------------------------
##' @rdname  plot.qcs
##' @method plot qcs.np
##' @export
##' 
plot.qcs.np <- function(x, title  =  NULL, subtitle  =  NULL, xlab  =  NULL,
                        ylab  =  NULL, ylim  =  NULL, conf.nsigma.alert  =  NULL,
                        center.nominal  =  NULL, limits.specification  =  NULL,
                        label.index  =  NULL, ...)
#.........................................................................                     
  {
  
  m<-length(x$statistics)
  
  if (!is.null(label.index))
    if (m-1 < label.index) stop("the covariable index is out of range")
  
  if(!is.null(conf.nsigma.alert)){
    limits.alert <- limits.xbar(x$center, x$std.dev, x$sizes,
                                conf.nsigma.alert)
  } else {
    limits.alert <- NULL
  }
  
  if (is.null(title)) title <- x$data.name
  
  if (is.null(subtitle)) subtitle <- expression(paste("Chart of control np"," "))
  if (is.null(xlab)) xlab <- "Sample"
  if (is.null(ylab)) ylab <- "np"
  
  plot.qcs(x, title , subtitle, xlab, ylab, ylim, center.nominal,
           limits.specification, limits.alert, label.index)
} #plot.qcs.np
#.........................................................................                     

#-------------------------------------------------------------------------
# plot.qcs.c
#-------------------------------------------------------------------------
##' @rdname  plot.qcs
##' @method plot qcs.c
##' @export
##' 
plot.qcs.c <- function(x, title  =  NULL, subtitle  =  NULL, xlab  =  NULL,
                       ylab  =  NULL, ylim  =  NULL, conf.nsigma.alert  =  NULL,
                       center.nominal  =  NULL, limits.specification  =  NULL,
                       label.index  =  NULL, ...)
#.........................................................................                     
  {
  
  m<-length(x$statistics)
  
  if (!is.null(label.index))
    if (m-1 < label.index) stop("the covariable index is out of range")
  
  if(!is.null(conf.nsigma.alert)){
    limits.alert <- limits.xbar(x$center, x$std.dev, x$sizes,
                                conf.nsigma.alert)
  } else {
    limits.alert <- NULL
  }
  
  if (is.null(title)) title <- x$data.name
  
  if (is.null(subtitle)) subtitle <- expression(paste("Chart of control c"," "))
  if (is.null(xlab)) xlab <- "Sample"
  if (is.null(ylab)) ylab <- "c"
  
  plot.qcs(x, title , subtitle, xlab, ylab, ylim, center.nominal,
           limits.specification, limits.alert, label.index)
} #plot.qcs.c
#.........................................................................                     

#-------------------------------------------------------------------------
# plot.qcs.u
#-------------------------------------------------------------------------
##' @rdname  plot.qcs
##' @method plot qcs.u
##' @export
##' 
plot.qcs.u <- function(x, title  =  NULL, subtitle  =  NULL, xlab  =  NULL,
                       ylab  =  NULL, ylim  =  NULL, conf.nsigma.alert  =  NULL,
                       center.nominal  =  NULL, limits.specification  =  NULL,
                       label.index  =  NULL, ...)
#.........................................................................                     
  {
  
  m<-length(x$statistics)
  
  if (!is.null(label.index))
    if (m-1 < label.index) stop("the covariable index is out of range")
  
  if(!is.null(conf.nsigma.alert)){
    limits.alert <- limits.xbar(x$center, x$std.dev, x$sizes,
                                conf.nsigma.alert)
  } else {
    limits.alert <- NULL
  }
  
  if (is.null(title)) title <- x$data.name
  
  if (is.null(subtitle)) subtitle <- expression(paste("Chart of control u"," "))
  if (is.null(xlab)) xlab <- "Sample"
  if (is.null(ylab)) ylab <- "u"
  
  plot.qcs(x, title , subtitle, xlab, ylab, ylim, center.nominal,
           limits.specification, limits.alert, label.index)
} #plot.qcs.u
#.........................................................................                     

#-------------------------------------------------------------------------
# plot.qcs.ewma
#-------------------------------------------------------------------------
##' @rdname  plot.qcs
##' @method plot qcs.ewma
##' @export
##' 
plot.qcs.ewma <- function(x, title  =  NULL, subtitle  =  NULL, xlab  =  NULL,
                          ylab  =  NULL, ylim  =  NULL, 
                          label.index  =  NULL, ...)
#.........................................................................                     
  {
  
  
  
  
  m<-length(x$statistics)
  
  if(is.null(ylim)) 
    ylim <-  range(x$statistics, x$limits)
  
  
  if (!is.null(label.index))
    if (m-1 < label.index) stop("the covariable index is out of range")
  
  if (is.null(title)) title <- x$data.name
  
  if (is.null(subtitle)) subtitle <- expression(paste("Chart of control ewma"," "))
  if (is.null(xlab)) xlab <- "Sample"
  if (is.null(ylab)) ylab <- "ewma"
  
  plot.qcs(x, title , subtitle, xlab, ylab, ylim, label.index, type.data = "dependence")
} #plot.qcs.ewma
#.........................................................................                     


#-------------------------------------------------------------------------
# plot.qcs.cusum
#-------------------------------------------------------------------------
##' @rdname  plot.qcs
##' @method plot qcs.cusum
##' @export
##' 
plot.qcs.cusum <- function(x, title  =  NULL, subtitle  =  NULL, xlab  =  NULL,
                          ylab  =  NULL, ylim  =  NULL, 
                          label.index  =  NULL, ...)
  #.........................................................................                     
{
  
  m<-length(x$statistics)
  
  if(is.null(ylim)) 
    ylim <-  range(x$pos, x$neg, x$limits)
  
  
  if (!is.null(label.index))
    if (m-1 < label.index) stop("the covariable index is out of range")
  
  if (is.null(title)) title <- x$data.name
  
  if (is.null(subtitle)) subtitle <- expression(paste("Chart of control cusum"," "))
  if (is.null(xlab)) xlab <- "Sample"
  if (is.null(ylab)) ylab <- "Cusum"
  
  plot.qcs(x, title , subtitle, xlab, ylab, ylim, label.index, type.data = "dependence")
} #plot.qcs.cusum
#.........................................................................                     


