#' nos.legend Generates a text-based legend for plots.
#' @description Generates a text-based legend that relates plot characters to functional values or outcomes.
#' 
#' @param legend A character vector of the plot symbols.
#' @param nams A vector that specifies what each plot symbol represents.
#' @param width A numeric value that specifies the number of characters to be printed before a line-break.
#' @author Hien D. Nguyen
#' @seealso \code{\link[txtplot]{txtplot}}.
#' @examples
#' ## Generates a legend that relates the plot symbols c('a','b','c') to 
#' ## the functional values c(1,2,3) with line-breaks every 14 characters.
#' nos.legend(c(1,2,3),c('a','b','c'),14)
# nos.legend -----------------------------------------------------------
nos.legend <- function(nams,legend=1:length(nams),width){
  insEOL <- function(x, width){
    lenx <- length(x)
    lenins <- floor(lenx/width)
    if(lenins > 0)
      out <- character(lenx + lenins+1)
    else
      return(c(x, "\n"))
    ind <- setdiff(1:(length(out)-1), (1:lenins)*(width+1))
    out[ind] <- x
    ind <- c((1:lenins)*(width+1), lenx + lenins+1)
    out[ind] <- "\n"
    out
  }
  leg <- paste(legend, nams, sep=" ~ ")
  leg <- paste(leg, collapse = "; ")
  leg <- insEOL(strsplit(leg, NULL)[[1]], width)
  cat("Legend: ")
  if(length(leg)>width)
    cat("\n")
  cat(leg, sep="")
}

#' nos.ecdf text-based plot of an empirical CDF.
#' @description Plots a text-based empirical cumulative distribution function.
#' 
#' @param data A numeric vector containing the values to be plotted. \code{NA} and \code{NaN} are allowed, but are removed for the plot. Infinities will cause an error.
#' @param pch A single-character plot symbol.
#' @param width Width of the plot in points.
#' @param height Height of the plot in points.
#' @param xlab Label of the x-axis of the plot.
#' @param ratio Coefficient that controls the aspect ratio of the plot.
#' @note
#' Due to rounding to a relatively crude grid results can only be approximate. The equally spaced axis ticks, for example, may be non-equally spaced in the plot. Further, due to the crude grid also there might be several points per character. The function uses the same plotting symbol no matter how many points coincide on one character position.
#' @author Hien D. Nguyen
#' @seealso \code{\link[stats]{ecdf}}, \code{\link[graphics]{plot}} and \code{\link[txtplot]{txtplot}}.
#' @examples
#' ## Plot the empirical CDF of 10 random standard normal points with 'o' shaped points.
#' data <- rnorm(10)
#' nos.ecdf(data)
#' 
#' ## Plot the empirical CDF of 100 random stanard normal points with '*' shaped points
#' data <- rnorm(100)
#' nos.ecdf(data,pch='*')
# nos.ecdf --------------------------------------------------------------
nos.ecdf <- function(data,
                     xlab = NULL,
                     ratio = 0.25,
                     width = round(options()$width*0.8),
                     height = round(ratio*width),
                     pch = 'o') {
  ndat <- length(data)
  sortdat <- sort(data)
  res <- max(width*2,height*2,min(ndat,width*10,height*10))
  points <- c(rep('~',res),rep(pch,ndat))
  app <- approx(sortdat,1:ndat/ndat,method='constant',n=res,ties=max)
  txtplot(c(app$x,sortdat),c(app$y,1:ndat/ndat),xlab=xlab,ylab='ECDF',pch=points,width=width,height=height)
}

#' nos.xyplot text-based scatter plot.
#' @description Plots a text-based scatter plot, with the option of having a regression line, or a linear interpolation of the points.
#' 
#' @param x A numeric vector containing the x-values to be plotted. \code{NA} and \code{NaN} are allowed, but are removed for the plot. Infinities will cause an error.
#' @param y A numeric vector containing the y-values to be plotted. \code{NA} and \code{NaN} are allowed, but are removed for the plot. Infinities will cause an error.
#' @param pch A two dimensional vector of single-character symbols. The first symbol is for the x and y coordinate points, and the second symbol is for the interpolation or regression line.
#' @param width Width of the plot in points.
#' @param height Height of the plot in points.
#' @param ratio Coefficient that controls the aspect ratio of the plot.
#' @param type One of the values in the set \code{c('l','p','r','lp','pr')}. Type \code{'l'} plots a linear interpolation, type \code{'p'} plots the x and y coordinates as points, type \code{'r'} plots a linear regression line, type \code{'lp'} plots x and y coordinates along with a linear interpolation, and type \code{'pr'} plots x and y coordinates along with a linear regression line.
#' @param xlab Label of the x-axis of the plot.
#' @param ylab Label of the y-axis of the plot.
#' @note
#' Due to rounding to a relatively crude grid results can only be approximate. The equally spaced axis ticks, for example, may be non-equally spaced in the plot. Further, due to the crude grid also there might be several points per character. The function uses the same plotting symbol no matter how many points coincide on one character position.
#' @author Hien D. Nguyen
#' @seealso \code{\link[graphics]{plot}} and \code{\link[txtplot]{txtplot}}.
#' @examples
#' ## Plot 10 correlated points
#' x <- 10*runif(10)
#' y <- x + rnorm(10)
#' nos.xyplot(x,y,type='p',xlab='x',ylab='y')
#' 
#' ## Plot 10 correlated points with a regression line
#' x <- 10*runif(10)
#' y <- x + rnorm(10)
#' nos.xyplot(x,y,type='pr',xlab='x',ylab='y')
#' 
#' ## Plot 10 correlated points with a linear interpolation
#' x <- 10*runif(10)
#' y <- x + rnorm(10)
#' nos.xyplot(x,y,type='lp',xlab='x',ylab='y')
# nos.xyplot ------------------------------------------------------------
nos.xyplot <- function(x=1:length(y),
                       y,
                       xlab = NULL,
                       ylab = NULL,
                       ratio = 0.25,
                       width = round(options()$width*0.8),
                       height = round(ratio*width),
                       pch = c('o','~'),
                       type='p') {
  ndat <- length(x)
  res <- max(width*2,height*2,min(ndat,width*10,height*10))
  if (type=='p') {
    txtplot(x,y,xlab=xlab,ylab=ylab,pch=pch[1],width=width,height=height)
  } 
  if (type=='l') {
    app <- approx(x,y,n=res)
    txtplot(app$x,app$y,xlab=xlab,ylab=ylab,pch=pch[2],width=width,height=height)
  } 
  if (type=='lp') {
    points=c(rep(pch[2],res),rep(pch[1],ndat))
    app <- approx(x,y,n=res)
    txtplot(c(app$x,x),c(app$y,y),xlab=xlab,ylab=ylab,pch=points,width=width,height=height)
  }
  if (type=='r') {
    linmod <- lm(y~x)
    pred <- predict(linmod)
    app <- approx(x,pred,n=res)
    txtplot(app$x,app$y,xlab=xlab,ylab=ylab,pch=pch[2],width=width,height=height)
  }
  if (type=='pr') {
    linmod <- lm(y~x)
    points=c(rep(pch[2],res),rep(pch[1],ndat))
    pred <- predict(linmod)
    app <- approx(x,pred,n=res)
    txtplot(c(app$x,x),c(app$y,y),xlab=xlab,ylab=ylab,pch=points,width=width,height=height)
  }
  if ( is.na(match(type,c('l','p','r','lp','pr'))) ) {
    stop('Type is required to be l, p, r, lp, or pr.')
  }
}

#' nos.density text-based plot of a kernel density function.
#' @description Plots a text-based of a kernel density function, with the option of plotting the location of the data points along the x-axis.
#' 
#' @param data A numeric vector containing the values to be plotted. \code{NA} and \code{NaN} are allowed, but are removed for the plot. Infinities will cause an error.
#' @param pch A two dimensional vector of single-character symbols. The first symbol is for data point, and the second symbol is for the density curve.
#' @param width Width of the plot in points.
#' @param height Height of the plot in points.
#' @param ratio Coefficient that controls the aspect ratio of the plot.
#' @param locations If \code{TRUE}, the location of the data points are plotted along the x-axis. If \code{FALSE}, then the locations of the data points are not plotted.
#' @param xlab Label of the x-axis of the plot.
#' @param bw A numerical value or character string to specify the bandwidth used for kernel density estimation. The default setting is \code{'nrd0'}, which is the default option in \code{\link[stats]{density}}; see \code{\link[stats]{density}} for other options.
#' @param kernel A character string to specify the type of smoothing kernel used. The default setting is \code{'gaussian'}, which is the default option in \code{\link[stats]{density}}; see \code{\link[stats]{density}} for other options.
#' @note
#' Due to rounding to a relatively crude grid results can only be approximate. The equally spaced axis ticks, for example, may be non-equally spaced in the plot. Further, due to the crude grid also there might be several points per character. The function uses the same plotting symbol no matter how many points coincide on one character position.
#' @author Hien D. Nguyen
#' @seealso \code{\link[stats]{density}}, \code{\link[graphics]{plot}}, \code{\link[txtplot]{txtplot}}, and \code{\link[txtplot]{txtdensity}}.
#' @examples
#' ## Plot a kernel density function of 10 random standard normal points with 
#' ## a Gaussian kernel and with the locations of the data plotted along the x-axis.
#' data <- rnorm(10)
#' nos.density(data)
#' 
#' ## Plot a kernel density function of 100 random stanard normal points with 
#' ## a triangular kernel and without the locations of the data plotted along the x-axis.
#' data <- rnorm(100)
#' nos.density(data,kernel='triangular',location=FALSE)
# nos.density -----------------------------------------------------------
nos.density <- function(data,
                        xlab = NULL,
                        ratio = 0.25,
                        bw = "nrd0",
                        kernel = 'gaussian',
                        locations = T,
                        width = round(options()$width*0.8),
                        height = round(ratio*width),
                        pch = c('o','~') ) {
  dense <- density(data,bw=bw,kernel=kernel)
  ndat <- length(data)
  if (locations == T) {
    xax <- c(data,dense$x)
    yax <- c(rep(0,ndat),dense$y)
    points <- c(rep(pch[1],ndat),rep(pch[2],length(dense$y)))
    txtplot(xax,yax,xlab=xlab,ylab='Density',pch=points,width=width,height=height)
  } else {
    txtplot(dense$x,dense$y,xlab=xlab,ylab='Density',pch=pch[2],width=width,height=height)
  }
}

#' nos.qqnorm text-based normal quantile-quantile plots
#' @description Produces a text-based normal quantile-quantile plot where the theoretical quantiles are plotted along the x-axis and the sample quantiles are plotted along the y-axis.
#' 
#' @param data A numeric vector containing the values to be plotted. \code{NA} and \code{NaN} are allowed, but are removed for the plot. Infinities will cause an error.
#' @param pch A two dimensional vector of single-character symbols. The first symbol is for data point, and the second symbol is for the line of theoretical quantile equality.
#' @param width Width of the plot in points.
#' @param height Height of the plot in points.
#' @param ratio Coefficient that controls the aspect ratio of the plot.
#' @param line If \code{TRUE}, the line of theoretical quantile equality is plotted. If \code{FALSE}, then the line of theoretical quantile equality is not plotted.
#' @note
#' Due to rounding to a relatively crude grid results can only be approximate. The equally spaced axis ticks, for example, may be non-equally spaced in the plot. Further, due to the crude grid also there might be several points per character. The function uses the same plotting symbol no matter how many points coincide on one character position.
#' @author Hien D. Nguyen
#' @seealso \code{\link[stats]{qqnorm}}, \code{\link[stats]{qqline}} and \code{\link[txtplot]{txtplot}}.
#' @examples
#' ## Produce a normal quantile-quantile plot of 10 random standard normal points, 
#' ## without the line of theoretical quantile equality.
#' data <- rnorm(10)
#' nos.qqnorm(data,line=FALSE)
#' 
#' ## Produce a normal quantile-quantile plot of 100 random chi-squared(3) points, 
#' ## with the line of theoretical quantile equality.
#' data <- rchisq(100,3)
#' nos.qqnorm(data)
# nos.qqnorm ------------------------------------------------------------
nos.qqnorm <- function(data,
                       line = T,
                       ratio = 0.25,
                       width = round(options()$width*0.8),
                       height = round(ratio*width),
                       pch = c('o','~') ) {
  ndat <- length(data)
  res <- max(width*2,height*2,min(ndat,width*10,height*10))
  if (line == T) {
    lin <- seq(min(data),max(data),length.out=res)
    xax <- (c(lin,qnorm(ppoints(data),mean(data),sd(data)))-mean(data))/sd(data)
    yax <- c(lin,sort(data))
    points <- c(rep(pch[2],length(lin)),rep(pch[1],ndat))
    txtplot(xax,yax,xlab='Theoretical Qs',ylab='Sample Qs',pch=points,width=width,height=height)    
  } else {
    txtplot(qnorm(ppoints(data)),sort(data),xlab='Theoretical Qs',ylab='Sample Qs',pch=pch[1],width=width,height=height)    
  }
}

#' nos.qqplot text-based generic quantile-quantile plots
#' @description Produces a text-based quantile-quantile plot between two equal size samples, with the option of plotting a line of theoretical quantile equality.
#' 
#' @param x A numeric vector containing the values to be plotted along the x-axis. \code{NA} and \code{NaN} are allowed, but are removed for the plot. Infinities will cause an error.
#' @param y A numeric vector containing the values to be plotted along the y-axis. \code{NA} and \code{NaN} are allowed, but are removed for the plot. Infinities will cause an error.
#' @param pch A two dimensional vector of single-character symbols. The first symbol is for data point, and the second symbol is for the line of theoretical quantile equality.
#' @param width Width of the plot in points.
#' @param height Height of the plot in points.
#' @param xlab Label of the x-axis of the plot.
#' @param ylab Label of the y-axis of the plot.
#' @param ratio Coefficient that controls the aspect ratio of the plot.
#' @param line If \code{TRUE}, the line of theoretical quantile equality is plotted. If \code{FALSE}, then the line of theoretical quantile equality is not plotted.
#' @note
#' Due to rounding to a relatively crude grid results can only be approximate. The equally spaced axis ticks, for example, may be non-equally spaced in the plot. Further, due to the crude grid also there might be several points per character. The function uses the same plotting symbol no matter how many points coincide on one character position.
#' @author Hien D. Nguyen
#' @seealso \code{\link[stats]{qqplot}}, \code{\link[stats]{qqline}}, \code{\link[stats]{ppoints}}, and \code{\link[txtplot]{txtplot}}.
#' @examples
#' ## Produce a quantile-quantile plot between two samples of 10 random standard normal points, 
#' ## without the line of theoretical quantile equality.
#' x <- rnorm(10)
#' y <- rnorm(10)
#' nos.qqplot(x,y,line=FALSE)
#' 
#' ## Produce a quantile-quantile plot of 100 random chi-squared(3) points against the 
#' ## true theoretical distribution, with the line of theoretical quantile equality.
#' y <- rchisq(100,3)
#' x <- qchisq(ppoints(100),3)
#' nos.qqplot(x,y,xlab='Theoretical Qs',ylab='Sample Qs')
# nos.qqplot ------------------------------------------------------------
nos.qqplot <- function(x,
                       y,
                       xlab=NULL,
                       ylab=NULL,
                       line = T,
                       ratio = 0.25,
                       width = round(options()$width*0.8),
                       height = round(ratio*width),
                       pch = c('o','~') ) {
  ndat <- length(x)
  res <- max(width*2,height*2,min(ndat,width*10,height*10))
  if (line == T) {
    lin <- seq(min(x),max(x),length.out=res)
    xax <- c(lin,sort(x))
    yax <- c(lin,sort(y))
    points <- c(rep(pch[2],length(lin)),rep(pch[1],ndat))
    txtplot(xax,yax,xlab=xlab,ylab=ylab,pch=points,width=width,height=height)    
  } else {
    txtplot(sort(x),sort(y),xlab=xlab,ylab=ylab,pch=pch[1],width=width,height=height)    
  }
}

#' nos.hist text-based plot of a histogram.
#' @description Plots a text-based of a histogram, with the option to plot densities or frequencies.
#' 
#' @param data A numeric vector containing the values to be plotted. \code{NA} and \code{NaN} are allowed, but are removed for the plot. Infinities will cause an error.
#' @param pch A single-character plot symbol.
#' @param width Width of the plot in points.
#' @param height Height of the plot in points.
#' @param ratio Coefficient that controls the aspect ratio of the plot.
#' @param xlab Label of the x-axis of the plot.
#' @param breaks Either a numerical value or a character string to specify the number of breaks to used. The default setting is \code{'Sturges'}, which is the default option in \code{\link[graphics]{hist}}; see \code{\link[graphics]{hist}} for other options.
#' @param freq If \code{TRUE}, the y-coordinate will display frequencies. If \code{FALSE}, then the y-coordinate will display densities.
#' @note
#' Due to rounding to a relatively crude grid results can only be approximate. The equally spaced axis ticks, for example, may be non-equally spaced in the plot. Further, due to the crude grid also there might be several points per character. The function uses the same plotting symbol no matter how many points coincide on one character position. Histogram columns are plotted at the midpoint of each bin.
#' @author Hien D. Nguyen
#' @seealso \code{\link[graphics]{hist}}, and \code{\link[txtplot]{txtplot}}.
#' @examples
#' ## Plot a histogram for the frequencies of 100 random standard normal points 
#' ## using 'Sturges' breaks and plot symbol 'o'.
#' data <- rnorm(100)
#' nos.hist(data)
#' 
#' ## Plot a histogram for the densities of 1000 random chi-squared(3) points 
#' ## using 'FD' breaks and plot symbol '#'.
#' data <- rchisq(1000,3)
#' nos.hist(data,breaks='FD',freq=FALSE,pch='#')
# nos.hist --------------------------------------------------------------
nos.hist <- function(data,
                     breaks = "Sturges", 
                     freq = T,
                     xlab = NULL,
                     ratio = 0.25,
                     width = round(options()$width*0.8),
                     height = round(ratio*width),
                     pch = 'o') {
  ndat <- length(data)
  histo <- hist(data,breaks = breaks, plot=F)
  xax <- histo$mids
  res <- max(width*2,height*2,min(ndat,width*10,height*10))
  if (freq == T) {
    yax <- histo$counts
    ydum <- xdum <- c()
    for (ii in 1:length(yax)) {
      ydum <- c(ydum,seq(0,yax[ii],length.out=res))
      xdum <- c(xdum,rep(xax[ii],res))
    }
    txtplot(xdum,ydum,xlab=xlab,ylab='Frequency',pch=pch,width=width,height=height)
  } else {
    yax <- histo$density
    ydum <- xdum <- c()
    for (ii in 1:length(yax)) {
      ydum <- c(ydum,seq(0,yax[ii],length.out=res))
      xdum <- c(xdum,rep(xax[ii],res))
    }
    txtplot(xdum,ydum,xlab=xlab,ylab='Density',pch=pch,width=width,height=height)
  }
}

#' nos.image text-based image plot.
#' @description Produces a text-based image plot for the visualization of three-dimensional data that are stored in a matrix.
#' 
#' @param data A numeric matrix of data to be plotted. \code{NA} and \code{NaN} are allowed, but are removed for the plot. Infinities will cause an error.
#' @param pch A vector of single-character symbols. The length of the vector determines the number of bins that the data is partitioned into. From left to right, the symbols represent an increasing order of binned.
#' @param xmin A numeric value indicating the smallest x-axis value, which corresponds to the first row of data.
#' @param xmax A numeric value indicating the greatest x-axis value, which corresponds to the last row of data.
#' @param ymin A numeric value indicating the smallest y-axis value, which corresponds to the first column of data.
#' @param ymax A numeric value indicating the greatest y-axis value, which corresponds to the last column of data.
#' @param width Width of the plot in points.
#' @param height Height of the plot in points.
#' @param ratio Coefficient that controls the aspect ratio of the plot.
#' @param xlab Label of the x-axis of the plot.
#' @param ylab Label of the y-axis of the plot.
#' @note
#' Due to rounding to a relatively crude grid results can only be approximate. The equally spaced axis ticks, for example, may be non-equally spaced in the plot. Further, due to the crude grid also there might be several points per character. The function uses the same plotting symbol no matter how many points coincide on one character position. The bins are equally spaced between the minimum and the maximum value of the data matrix, and are produced using the \code{\link[graphics]{hist}} function. The legend reports the smallest and largest value of each bin associated with each plot symbol.
#' @author Hien D. Nguyen
#' @seealso \code{\link[graphics]{hist}}, \code{\link[graphics]{image}}, \code{\link[datasets]{volcano}}, and \code{\link[txtplot]{txtplot}}.
#' @examples
#' ## Produce the image plot of the volcano dataset, using the default plotting symbols.
#' library(datasets)
#' nos.image(volcano)
#' 
#' ## Produce the image plot of the volcano dataset, using the plotting symbols 1:9.
#' library(datasets)
#' nos.image(volcano,pch=1:9)
# nos.image -------------------------------------------------------------
nos.image <- function(data,
                      xmin = 1,
                      xmax = dim(data)[1],
                      ymin = 1,
                      ymax = dim(data)[2],
                      xlab = NULL,
                      ylab = NULL,
                      ratio = 0.35,
                      width = round(options()$width*0.8),
                      height = round(ratio*width),
                      pch = c('.','o','x','X','#') ) {
  na_status <- !is.na(data)
  longform <- cbind(c(data)[na_status],which(na_status,arr.ind = T))
  longform[,2] <- (longform[,2]-1)/(dim(data)[1]-1)*(xmax-xmin)+xmin
  longform[,3] <- (longform[,3]-1)/(dim(data)[2]-1)*(xmax-xmin)+xmin
  histo <- hist(longform[,1],breaks = seq(min(longform[,1]),max(longform[,1]),length.out=(length(pch)+1)), plot=FALSE)
  alloc <- sapply(longform[,1],function(x) which.min(abs(x - histo$mids)))
  points <- pch[alloc]
  txtplot(longform[,2],longform[,3],xlab=xlab,ylab=ylab,pch=points,width=width,height=height)
  up_bound <- histo$breaks[-1]
  low_bound <- histo$breaks[-length(histo$breaks)]
  nam <- paste('(',signif(low_bound,3),',',signif(up_bound,3),')',sep='')
  nos.legend(nam,pch,width)
}

#' nos.image text-based contour plot.
#' @description Produces a text-based contour plot for the visualization of three-dimensional data that are stored in a matrix.
#' 
#' @param data A numeric matrix of data to be plotted. \code{NA} and \code{NaN} are allowed, but are removed for the plot. Infinities will cause an error.
#' @param pch A vector of single-character symbols. The length of the vector determines the number of contours to be plotted. From left to right, the symbols represent an increasing order of contour values.
#' @param xmin A numeric value indicating the smallest x-axis value, which corresponds to the first row of data.
#' @param xmax A numeric value indicating the greatest x-axis value, which corresponds to the last row of data.
#' @param ymin A numeric value indicating the smallest y-axis value, which corresponds to the first column of data.
#' @param ymax A numeric value indicating the greatest y-axis value, which corresponds to the last column of data.
#' @param width Width of the plot in points.
#' @param height Height of the plot in points.
#' @param ratio Coefficient that controls the aspect ratio of the plot.
#' @param xlab Label of the x-axis of the plot.
#' @param ylab Label of the y-axis of the plot.
#' @note
#' Due to rounding to a relatively crude grid results can only be approximate. The equally spaced axis ticks, for example, may be non-equally spaced in the plot. Further, due to the crude grid also there might be several points per character. The function uses the same plotting symbol no matter how many points coincide on one character position. The contour values are equally spaced quantiles between the 5\% and 95\% level, and are produced using the \code{\link[stats]{quantile}} function. The legend reports the value associated with each contour.
#' @author Hien D. Nguyen
#' @seealso \code{\link[stats]{quantile}}, \code{\link[graphics]{contour}}, \code{\link[grDevices]{contourLines}}, \code{\link[datasets]{volcano}}, and \code{\link[txtplot]{txtplot}}.
#' @examples
#' ## Produce the contour plot of the volcano dataset, using the default plotting symbols.
#' library(datasets)
#' nos.contour(volcano)
#' 
#' ## Produce the contour plot of the volcano dataset, using the plotting symbols letters[1:3].
#' library(datasets)
#' nos.contour(volcano,pch=letters[1:3])
# nos.contour -------------------------------------------------------------
nos.contour <- function(data,
                        xmin = 1,
                        xmax = dim(data)[1],
                        ymin = 1,
                        ymax = dim(data)[2],
                        xlab = NULL,
                        ylab = NULL,
                        ratio = 0.35,
                        width = round(options()$width*0.8),
                        height = round(ratio*width),
                        pch = 1:5 ) {
  nlev <- length(pch)
  con_line <- contourLines(x=seq(xmin,xmax,length.out=nrow(data)),y=seq(ymin,ymax,length.out=ncol(data)),data,nlevels=nlev,levels=quantile(c(data),seq(0.05,0.95,length.out=nlev),na.rm=T))
  xax <- unlist(sapply(con_line,'[[',2))
  yax <- unlist(sapply(con_line,'[[',3))
  con_length <- sapply(sapply(con_line,'[[',2),length)
  levs <- rep(unlist(sapply(con_line,'[[',1)),con_length)
  unique_levs <- sort(unique(levs))
  match_unique <- match(levs,unique_levs)
  points <- pch[match_unique]
  txtplot(xax,yax,xlab=xlab,ylab=ylab,pch=points,width=width,height=height)
  nam <- signif(unique_levs,3)
  nos.legend(nam,pch,width)
}