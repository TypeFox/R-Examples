#<<BEGIN>>
plot.mc <- function(x, prec=0.001, stat = c("median","mean"), lim = c(0.025, 0.25, 0.75, 0.975), na.rm=TRUE, griddim = NULL, xlab = NULL, ylab = "Fn(x)", main = "", draw = TRUE, paint=TRUE, xlim=NULL,ylim=NULL,...)
#TITLE Plots Results of a Monte Carlo Simulation
#DESCRIPTION
# Plots the empirical cumulative distribution function of a \samp{mcnode} or a \samp{mc} object ("0" and "V" nodes) or the
#empirical cumulative distribution function of the estimate of a \samp{mcnode} or \samp{mc} object ("U" and "VU" nodes).
#KEYWORDS hplot
#INPUTS
# {x}<<a \samp{mcnode} or a \samp{mc} objects>>
#[INPUTS]
#{prec}<<the precision of the plot. 0.001 will
#provide an ecdf from the 0.000, 0.001, .002,  ..., 1.000 quantiles.>>
#{stat}<<the function used for estimates (2D \samp{mc} or \samp{mcnode}). By default the median.>>
#{lim}<<a vector of numbers (between 0 and 1) indicating the enveloppe (2D \samp{mc} or \samp{mcnode}) . Maybe \samp{NULL} or empty.>>
#{na.rm}<<Should NA values be discarded>>
#{griddim}<<a vector of two integers, indicating the size of the grid of the graph. If \samp{NULL}, the grid is calculated to produce a "nice" graph.>>
#{xlab}<<vector of labels for the x-axis. If \samp{NULL}, use the name of the node.>>
#{ylab}<<vector of labels for the y-axis.>>
#{main}<<vector of main titles of the graph.>>
#{draw}<<Should the plot be drawn?>>
#{paint}<<Should the enveloppes be filled?>>
#{xlim}<<x coordinate range. \samp{xlim} is either a vector of length 2, used for each graph, or a list of vectors of length 2, 
#whose ith element is used for the ith graph. By default, the data range is used as \samp{xlim}.>>
#{ylim}<<y coordinate range. \samp{ylim} is either a vector of length 2, used for each graph, or a list of vectors of length 2, 
#whose ith element is used for the ith graph. By default, the data range is 0-1.>>
#{\dots}<<further arguments to be passed to \samp{plot.stepfun}.>>
#DETAILS
#\samp{plot.mcnode} is a user-friendly function that send the \samp{mcnode} to \samp{plot.mc}.</>
#For \samp{"VU"} and \samp{"U"} \samp{mcnode}s, quantiles are calculated using \code{\link{quantile.mc}}
#within each of the \samp{nsu} simulations (i.e. by columns of each \samp{mcnode}). The medians (but may be
#the means using \samp{stat="mean"}) calculated from the \samp{nsu} values are plotted. The 0.025 and 0.975 quantiles, 
#and the 0.25 and 0.75 quantiles (default values of \samp{lim}) of these quantiles are used as the enveloppe.
#REFERENCE
#Cullen AC and Frey HC (1999) Probabilistic techniques in exposure assessment. Plenum Press, USA, pp. 81-155. 
#VALUE
#A \samp{plot.mc} object, list of the quantiles used to plot the draw.
#SEE ALSO
#\code{\link{ecdf}}, \code{\link{plot}}, \code{\link{quantile.mc}}
#EXAMPLE
#data(total)

#plot(xVUM3)
### only one enveloppe corresponding to quantiles 0.025 and 0.975
#plot(xVUM3,lim=c(0.025,0.975)) 
### only one enveloppe not painted
#plot(xVUM3,lim=c(0.025,0.975),paint=FALSE) 

#def.par <- par(no.readonly = TRUE)
#par(mar=c(4,4,1,1))
#plot(total)
#par(def.par)

#CREATED 07-08-01
#REVISED 10-02-10
#--------------------------------------------
#
{
  if(inherits(x,"mc")==TRUE) {
    x <- quantile.mc(x, probs=seq(0,1,prec),lim = lim, na.rm=na.rm, lnames=xlab)
    }

  if(draw) {
	y <- x                           # for a correct return
    stat <- match.arg(stat)

	 beau <- function(n){
		nc <- round(sqrt(n))
		nr <- ceiling(n/nc)
		c(nc,nr)}

   noms <- names(rapply(x,function(x) 1))    #moche mais efficace
   if(is.null(xlab)) xlab <- noms
   n <- length(noms)

   if(!is.null(ylim) & ((is.list(ylim) & length(ylim)!= n)|(is.vector(ylim) & length(ylim)!= 2))) stop("ylim should be NULL, a vector of 2 elements or a list of length the number of nodes") 
   if(!is.null(xlim) & ((is.list(xlim) & length(xlim)!= n)|(is.vector(xlim) & length(xlim)!= 2))) stop("xlim should be NULL, a vector of 2 elements or a list of length the number of nodes") 

	 main <- rep(main,n)
	 xlab <- rep(xlab,n)
	 ylab <- rep(ylab,n)

  if(is.null(griddim)) griddim <- beau(n)
  if(prod(griddim) < n) op <- par(mfrow=griddim,ask=TRUE,mar=c(5,4,.2,.2))
     else op <- par(mfrow=griddim, mar=c(5,4,.2,.2))

  try({   #to restore par in case of error

  i <- 1
  env <- environment()
  
  LEPLOT <- function(y,...){
      if(nrow(y) != 1) {
        if(stat=="median") y <- y[-2,,drop=FALSE]
        else y <- y[-1,,drop=FALSE]}                                              #Retrait median or mean
  		nr <- nrow(y)
      i <- get("i",envir=env)
      xlima <- if(is.null(xlim)) range(y,na.rm=TRUE) else 
		xlima <- if(is.list(xlim)) xlim[[i]] else xlim
  	  if(xlima[1]==xlima[2]) xlima <- xlima + c(-0.5,0.5)
      ylima <- if(is.null(ylim)) c(0,1) else 
		ylima <- if(is.list(ylim)) ylim[[i]] else ylim
      x50 <- plot.stepfun(y[1,], main=main[i], xlim=xlima, ylim=ylima, ylab = ylab[i], verticals = TRUE, do.points = FALSE, xlab=xlab[i], lwd=3, ...)
      abline(h = c(0, 1), col =  "gray70", lty = 3)

      # Points for the polygon used to fill the enveloppe
      if(paint){
        ti.l <- x50$t[-length(x50$t)]
        ti.r <- x50$t[-1L]
        y50 <- x50$y
        thex50 <- rev(as.vector(rbind(ti.l,ti.r)))
        they50 <- rev(as.vector(rbind(y50, y50)))
      }
      
      if(nr > 1){
        rankplot <- 1 + order(-abs(lim-0.5))      # in order to draw over in the good order
        for(j in rankplot) {
          xp <- plot.stepfun(y[j,], verticals=TRUE, do.points=FALSE, add= TRUE, lty=3 ,col="gray30",...)
          if(paint){
            ti.lp <- xp$t[-length(xp$t)]
            ti.rp <- xp$t[-1L]
            yp <- xp$y
            thexp <- as.vector(rbind(ti.lp,ti.rp))
            theyp <- as.vector(rbind(yp, yp))
            color <- grey(abs(lim[j-1]-.5)+.25)
            polygon(c(thexp,thex50), c(theyp,they50), col= color)
          }
        }
      }
    assign("i",i+1,envir=env) }
    
  rapply(y,LEPLOT)
    
  })
	par(op)
  }
  class(x) <- "plotmc"
  return(invisible(x))}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<BEGIN>>
plot.mcnode <- function(x, ...)
#ISALIAS plot.mc
#--------------------------------------------
{ nom <- deparse(substitute(x))
  x <- list(x)
  names(x) <- nom
  class(x) <- "mc"
  x <- plot.mc(x, ... )
  return(invisible(x))}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<BEGIN>>
plot.plotmc <- function(x, ...)
#ISALIAS plot.mc
#--------------------------------------------
{ x <- plot.mc(x, ... )
  return(invisible(x))}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

