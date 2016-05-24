#' Calculate statistics for Bland-Altman-Plot
#' 
#' Does the computation for Bland Altman plots. This will usually be called from
#' graphic functions like \code{bland.altman.plot} but will be usefull for 
#' customized plot (see examples for color coded BA plot). Offers symmetric 
#' confidence intervalls for bias and upper and lower limits. 
#'  
#' @param group1 vector of numerics to be compared to group2
#' @param group2 vector of numerics to be compared to group1
#' @param two numeric defines how many standard deviations from mean are to be 
#'        computed, defaults to 1.96 as this gives proper 95 percent CI. However, in the
#'        original publication a factor of 2 is used.
#' @param mode if 1 then difference group1 minus group2 is used, if 2 then
#'        group2 minus group1 is used. Defaults to 1.
#' @param conf.int usefull
#' @return \code{means} vector of means, i. e. data for the x axis
#' @return \code{diffs} vector of differences, i. e. data for the y axis 
#' @return \code{groups} data.frame containing pairwise complete cases of group1 and 
#' group2. NAs are removed.
#' @return \code{based.on} count of pairwise complete cases in groups
#' @return \code{lower.limit} lower limit for BA plot
#' @return \code{mean.diffs} mean of differences, also called 'bias'
#' @return \code{upper.limit} upper limit for BA plot
#' @return \code{lines} vector containing y values where to draw horizontal 
#' lines, i. e. mean of differences minus "two" standard deviations, mean of 
#' differences and mean of differences plus "two" standard deviations (i. e.
#' \code{c(lower.limit, mean.diffs, upper.limit}). This is convenient for 
#' printing.
#' @return \code{CI.lines} vector of confidence intervalls for the values of 
#' lines (based on the assumption of normal distribution of differences 
#' \code{diffs}).
#' @return \code{two} the argument 'two'
#' @return \code{critical.diff} critical difference, i. e. 'two' times standard 
#' deviation of differences, equals half the difference of lower.limit and 
#' upper.limit
#' @author Bernhard Lehnert <bernhard.lehnert@@uni-greifswald.de> 
#' @seealso \code{\link{bland.altman.plot}}
#' @examples 
#' # simple calculation of stats:
#' a <- rnorm(20)
#' b <- jitter(a)
#' print(bland.altman.stats(a, b))
#' print(bland.altman.stats(a, b)$critical.diff)
#' 
#' # drawing Bland-Altman-Plot with color coding sex:
#' example.data <- data.frame(sex = gl(2,6,labels=c("f","m")),
#'                  m1 = c(16,10,14,18,16,15,18,19,14,11,11,17),
#'                  m2 = c(18, 9,15,19,19,13,19,20,14,11,13,17))
#' ba <- bland.altman.stats(example.data$m1, example.data$m2)
#' plot(ba$means, ba$diffs, col=example.data$sex, ylim=c(-4,4))
#' abline(h=ba$lines, lty=2)
#'               
#' # compute 95%-CIs for the bias and upper and lower limits of PEFR data as 
#' # in Bland&Altman 1986
#' bland.altman.stats(bland.altman.PEFR[,1],bland.altman.PEFR[,3])$CI.lines
#' # apparently wrong results? CAVE: Bland&Altman are using two=2, thus
#' bland.altman.stats(bland.altman.PEFR[,1],bland.altman.PEFR[,3], two=2)$CI.lines
#' @export
#' @importFrom graphics abline plot sunflowerplot
#' @importFrom stats na.omit qt sd
bland.altman.stats <- function(group1, group2, two=1.96, mode=1, conf.int=.95){
    if(length(group1) != length(group2)) 
        stop("Error in bland.altman.stats: groups differ in length.")
    if(!is.numeric(group1))
        stop("Error in bland.altman.stats: group1 is not numeric.")
    if(!is.numeric(group2))
        stop("Error in bland.altman.stats: group2 is not numeric.")
    if(two<=0)
        stop("Error in bland.altman.stats: inproper value of two.")
    if(mode!=1 & mode !=2)
        stop("Error in bland.altman.stats: mode must be either 1 oder 2.")
    
    dfr <- data.frame(group1 = group1, group2 = group2, check.names=FALSE)
    dfr <- na.omit(dfr)
    called.with <- length(group1)
    based.on <- length(dfr[[1]])
    if(based.on < 2)
        warning("Warning in bland.altman.stats:less than 2 data pairs after deleting NAs.", 
                call.=FALSE)
    if(mode==1) 
        diffs <- dfr[[1]]-dfr[[2]]
    if(mode==2) 
        diffs <- dfr[[2]]-dfr[[1]]
    means <- (dfr[[1]]+dfr[[2]])/2
    critical.diff <- two*sd(diffs)
    mean.diffs <- mean(diffs)  # aka 'bias'
    lower.limit <- mean.diffs-critical.diff
    upper.limit <- mean.diffs+critical.diff
    lines <- c(lower.limit = lower.limit, 
               mean.diffs = mean.diffs, 
               upper.limit = upper.limit)
    #confidence intervals
    t1 <- qt((1-conf.int)/2, df = based.on-1)
    t2 <- qt((conf.int+1)/2, df = based.on-1)
    CI.lines <- c(lower.limit.ci.lower=lower.limit + t1*sqrt(sd(diffs)^2*3/based.on),
                  lower.limit.ci.upper=lower.limit + t2*sqrt(sd(diffs)^2*3/based.on),
                  mean.diff.ci.lower=mean.diffs+t1*sd(diffs)/sqrt(based.on),
                  mean.diff.ci.upper=mean.diffs+t2*sd(diffs)/sqrt(based.on),
                  upper.limit.ci.lower=upper.limit + t1*sqrt(sd(diffs)^2*3/based.on),
                  upper.limit.ci.upper=upper.limit + t2*sqrt(sd(diffs)^2*3/based.on)
                  )
    
    return(list(
        means = means,
        diffs = diffs,
        groups = dfr,
        based.on = based.on,
        lower.limit = lower.limit,
        mean.diffs = mean.diffs,
        upper.limit = upper.limit,
        lines = lines,
        CI.lines = CI.lines,
        two = two,
        critical.diff = critical.diff))
}
    
#' Produce Bland-Altman Plot
#'  
#' Bland-AltmanPlots for assessing agreement between two measuring methods or
#' repeatability (test-retest agreement) of measurements. Using either base graphics
#' or ggplot2.
#' @param group1 Measurements with first method or first measurement
#' @param group2 Measurements with second method or second measurement
#' @param two Lines are drawn "two" standard deviations from mean differences.
#' This defaults to 1.96 for proper 95 percent confidence interval estimation
#' but can be set to 2.0 for better agreement with e. g. the Bland Altman publication.
#' @param mode if 1 then difference group1 minus group2 is used, if 2 then
#' group2 minus group1 is used. Defaults to 1.
#' @param graph.sys Graphing system within R. This defaults to "base" but can be
#' one out of \code{c("base", "ggplot2")}, providing ggplot2 is installed.
#' @param conf.int Defaults to 0 which draws the usual Bland Altman plot which 
#' contains no confidence intervalls. Change to .95 for 95 percent confidence
#' intervalls to be drawn.
#' @param silent logical. If graph.sys=="base" and silent==TRUE then no return value. 
#' If graph.sys=="base" and silent==FALSE then returns statistics.
#' @param sunflower logical. If TRUE, the plot will be based on a sunflower plot
#' and ties will be marked accordingly. Try with data with ties. Works only with
#' \code{graph.sys=="base"}.
#' @param geom_count logical. If TRUE, the dots will get larger the more frequent
#' given pair is. Use in presence of ties. Works only with
#' \code{graph.sys=="ggplot2"} version >=2.0.0.
#' @param ... passed on to graphics functions if \code{graph.sys=="base"}
#' @author Bernhard Lehnert <bernhard.lehnert@@uni-greifswald.de>
#' @return Depends on graphic system chosen. In case of "base" depending on whether
#' silent==TRUE. If silent==TRUE then no returns. If silent==FALSE than returns
#' list of statistics as returned by \code{bland.altman.stats()}. In case the 
#' graphics system is "ggplot2" than the graphic object is returned so that it 
#' can be printed or altered.
#' @seealso \code{\link{bland.altman.stats}}
#' @export
#' @examples
#' bland.altman.plot(rnorm(20), rnorm(20), xlab="mean measurement", 
#'                   ylab="differences", main="Example plot")
#'                   
#' bland.altman.plot(rnorm(20), 2+.8*rnorm(20), xlab="mean measurement", 
#'                   ylab="differences", conf.int=.95)
#'                   
#' bland.altman.plot(rnorm(200), 2+.8*rnorm(200), xlab="mean measurement", 
#'                   ylab="differences", conf.int=.95)
#'                   
#' # this is what fig.2 in Bland&Altman1986 would have looked like
#' PEFR1 <- bland.altman.PEFR[,1]
#' PEFR2 <- bland.altman.PEFR[,3]
#' bland.altman.plot(PEFR1, PEFR2, silent=TRUE, xlim=c(0,800),
#'                   xlab="Average PEFR by two meters",
#'                   ylab="Difference in PEFR (large-mini)")
#' 
#' # and this is the same but with additional 95 percent CIs
#' data(bland.altman.PEFR)
#' bland.altman.plot(PEFR1, PEFR2, silent=TRUE, conf.int=.95, xlim=c(0,800))
#'                   
#' # an example with many ties and the 'sunflower'-option
#' a <- rep(c(1,1,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,4,4,4,5,6,6),2)
#' b <- rep(c(1,1,1,2,2,2,3,1,4,2,5,3,3,3,3,3),3)
#' bland.altman.plot(a,b,sunflower=TRUE, xlab="Mean",ylab="Difference",
#'                   main="discrete values lead to ties")
#'                   
#' library(ggplot2)
#' a <- bland.altman.plot(rnorm(20), rnorm(20), graph.sys="ggplot2", conf.int=.9)
#' print(a + xlab("you can change this later") + ggtitle("Title goes here"))                  
bland.altman.plot <- function(group1, group2, two=1.96, mode=1,
                              graph.sys="base", conf.int=0, silent=TRUE,
                              sunflower = FALSE, geom_count=FALSE, ...){
    if(graph.sys=="base")
        return(bland.altman.base(group1, group2, two, mode, conf.int, 
                                 silent, sunflower, geom_count, ...))
    if(graph.sys=="ggplot2")
        return(bland.altman.ggplot2(group1, group2, two, mode, conf.int, 
                                 silent, sunflower, geom_count, ...))
    stop("ERROR: graph.sys must be element of c('base','ggplot2').")
} 
 
  
bland.altman.base <- function(group1, group2, two, mode, conf.int, 
                               silent, sunflower, geom_count, ...){
    if(geom_count)
        warning("bland.altman.base: No geom_count option for base graphics implemented.")
    
    ba <- bland.altman.stats(group1 = group1, group2 = group2, two=two, 
                             mode = mode, conf.int=conf.int)
    if(sunflower==TRUE){
        plotfun=sunflowerplot
    }
    else{
        plotfun=plot
    }
    if(conf.int==0){ 
        ymin <- min(c(ba$diffs, ba$lines))
        ymax <- max(c(ba$diffs, ba$lines))
        plotfun(ba$means, ba$diffs, ylim=c(ymin,ymax),...)
        abline(h=ba$lines, lty=2)
    }
    else{
        ymin <- min(c(ba$diffs, ba$CI.lines))
        ymax <- max(c(ba$diffs, ba$CI.lines))
        plotfun(ba$means, ba$diffs, ylim=c(ymin,ymax),...)
        abline(h=ba$lines, lty=2, col=c("blue","red","blue"),lwd=3)
        abline(h=ba$CI.lines[c(1,2,5,6)], lty=2, col="darkblue", lwd=2)
        abline(h=ba$CI.lines[c(3,4)], lty=2, col="darkred", lwd=2)
    }
   
    if(silent==FALSE)
        return(ba)
    else
        return()
}

bland.altman.ggplot2 <- function(group1, group2, two, mode, conf.int, 
                              silent, sunflower, geom_count, ...){
    if(sunflower)
        warning("bland.altman.plot: No sunflower option for ggplot2 implemented.")
    
    ba <- bland.altman.stats(group1 = group1, group2 = group2, two=two, 
                             mode = mode, conf.int=conf.int)
    if(!requireNamespace("ggplot2"))
        stop("Could not load ggplot2. Sorry.")
    values <- data.frame(m = ba$means, d = ba$diffs)
    geom_my <- if(geom_count) ggplot2::geom_count
               else ggplot2::geom_point
    m <- NULL; d <- NULL #this is useless but helps with CRAN tests ;-(
    p <- ggplot2::ggplot(values, ggplot2::aes(x=m, y=d))+
         geom_my()+
         ggplot2::geom_hline(yintercept=ba$lines, linetype=2, size=1.0)+
         ggplot2::xlab("mean of measurements")+
         ggplot2::ylab("difference")
    if(conf.int>0){
        p <- p+ggplot2::geom_hline(yintercept=ba$CI.lines, linetype=2, size=0.7)
    }
    #print(p)   # deleted after version 0.2.1
    return(p)
}