#' Standard 2x2 plot for sample data
#' 
#' Standard 2x2 plot for sample data
#'
#' @param object   An object of ecd class.
#' @param ts       The xts object for the timeseries.
#' @param EPS      Logical, indicating whether to save the plot to EPS, default = FALSE
#' @param eps_file File name for eps output
#'
#' @keywords plot
#'
#' @export plot_2x2
#'
#' @importFrom grDevices dev.off postscript
#'
#' @examples
#' \dontrun{
#' plot_2x2(d, ts)
#' }
### <======================================================================>
"plot_2x2.ecd" <- function (object, ts, EPS = FALSE, eps_file = NA)
{
	if ( EPS ) {
		postscript(file= eps_file,
			paper="special", width=8, height=8, horizontal=FALSE)
	}
	attr <- xtsAttributes(ts)
    if (! "tail" %in% names(attr)) {
        stop("ts does not contain tail stats")
    }

    main1 <- "ECD PDF"
    if ("symbol" %in% names(attr)) {
        main1 <- paste("PDF of symbol:", attr$symbol)
    }
    object <- quantilize(object)
    
	par(mfcol=c(2,2)); # 4 plots
	par(oma = c(1, 1, 4, 1))

	ecd.plot_pdf    (object, ts, main = main1)
	ecd.plot_logpdf (object, ts)
	ecd.plot_cdf    (object, ts)
	ecd.plot_qq     (object, ts)

	if ( EPS ) dev.off()

}
### <---------------------------------------------------------------------->
#' @rdname plot_2x2.ecd
setGeneric("plot_2x2", function(object, ts, EPS = FALSE, eps_file = NA) standardGeneric("plot_2x2"))
#' @rdname plot_2x2.ecd
setMethod("plot_2x2", signature("ecd"), plot_2x2.ecd)
### <---------------------------------------------------------------------->
### Standard plot facility
# TODO split the following to separate files
### <---------------------------------------------------------------------->

ecd.plot_pdf <- function (object, ts, xlab="$\\log(r)$", main="ECD PDF" )
{
    if (length(object@stats)==0) {
        stop("stats is not computed in ecd object")
    }

    hist <- xtsAttributes(ts)$hist
    m1 <- object@stats$m1
    s <- object@stats$stdev
    xe <- ellipticity(object)

    Q <- 0.005
    xmin <- qec(Q, object)
    xmax <- qec(1-Q, object)

    N <- 400
    logr  <- seq(xmin, xmax, length.out = N)
    pdf <- dec(logr, object)
    ymax <- max(max(pdf),max(hist$density))
    
	par(mar=c(4,4,2,1)) 

	plot( hist$mids, hist$density, 
		pch = "o", col="red",
		xlim=c(xmin, xmax),
		ylim=c(0, ymax),
		ylab="PDF", xlab=xlab, main=main)

	lines(logr, pdf, pch = ".", col="blue", lwd=2)
    abline(v=xe$xe1)
    abline(v=xe$xe2)
    
	sc = 0.8
	yseg <- 0.20*max(hist$density)
  
    # std dev indicator
	segments( m1-s, yseg, 
              m1+s, yseg, lwd=3 ) 

    text( m1, yseg*0.5, "stdev", cex=sc )
    
    yr1 <- format(min(zoo::index(ts)), "%Y")
    yr2 <- format(max(zoo::index(ts)), "%Y")
    text( xmax*0.7, ymax*0.8, yr1, cex=sc )
    text( xmax*0.7, ymax*0.7, "to", cex=sc )
    text( xmax*0.7, ymax*0.6, yr2, cex=sc )
    text( xe[1], ymax*0.4, pos=2, paste("$x_e$=", sprintf("%.4f", xe$xe1)), cex=sc)
    text( xe[2], ymax*0.4, pos=4, paste("$x_e$=", sprintf("%.4f", xe$xe2)), cex=sc)
    
	# legends
	legend(xmin, ymax, c("data","fit"), 
           cex=sc, 
		   col=c("red","blue"), 
           pch=c("o",NA),
           lwd=c(1,2),
           lty=c(NA,1));
}
# ------------------------------------------------------------------
ecd.plot_logpdf <- function (object, ts, xlab="$\\log(r)$", main="ECD Log PDF" )
{
    if (length(object@stats)==0) {
        stop("stats is not computed in ecd object")
    }

    hist <- xtsAttributes(ts)$hist
    m1 <- object@stats$m1
    s <- object@stats$stdev
    xe <- ellipticity(object)
    
    par(mar=c(4,4,2,1)) 

    fnt <- log(hist$density) [ which(is.finite(log(hist$density))) ]

    xmin <- min(hist$mids)
    xmax <- max(hist$mids)
    
    N <- 400
    logr  <- seq(xmin, xmax, length.out = N)
    pdf <- dec(logr, object)
    
    plot( hist$mids, log(hist$density), 
          pch = "o", col="red", 
        ylim=c( min(max(fnt)-8,min(fnt)), max(fnt)+0.5 ),
        xlim=c( xmin, xmax ),
        ylab="$\\log$(PDF)", xlab=xlab, main=main)
  
    lines(logr, log(pdf), pch = ".", col="blue", lwd=2)
    abline(v=xe$xe1)
    abline(v=xe$xe2)
    
    # legends

    # parameters
    xpos <- 0.7 * if ( abs(xmin) > abs(xmax) ) xmin else xmax 
    ypos <- max(log(hist$density))+0.5
    sc <- 0.8
    sp <- 0.8

    text(xpos,ypos-1*sp,labels=c("ecd fit"),cex=sc)
    text(xpos,ypos-2*sp,labels=paste("$\\alpha$=", sprintf("%.4f", object@alpha)), cex=sc)
    text(xpos,ypos-3*sp,labels=paste("$\\gamma$=", sprintf("%.4f", object@gamma)), cex=sc)
    text(xpos,ypos-4*sp,labels=paste("$\\sigma$=", sprintf("%.6f", object@sigma)), cex=sc)
    text(xpos,ypos-5*sp,labels=paste("$\\beta$=",  sprintf("%.4f", object@beta)), cex=sc)
    text(xpos,ypos-6*sp,labels=paste("$\\mu$=",    sprintf("%.6f", object@mu)), cex=sc)
}

# ------------------------------------------------------------------
ecd.plot_cdf <- function (object, ts, xlab="$\\log(r)$", main="ECD CDF" )
{
    if (length(object@stats)==0) {
        stop("stats is not computed in ecd object")
    }

    attr <- xtsAttributes(ts)
    hist <- attr$hist
    
    stats <- object@stats
    m1 <- stats$m1
    s <- stats$stdev

    Q <- 0.005
    xmin <- qec(Q, object)
    xmax <- qec(1-Q, object)
    
    N <- 400
    logr  <- seq(xmin, xmax, length.out = N)
    cdf <- pec(logr, object)
        
    par(mar=c(4,4,2,1)) 
    
    hist_cdf <- cumsum(hist$counts) / sum(hist$counts)
    plot(hist$mids, hist_cdf, 
         lwd=2, type="s", col="red", 
         xlim=c(xmin, xmax),
         ylab="CDF", xlab=xlab, main=main)
    
    dr <- logr[2]-logr[1]
    # lines(hist$mids, hist_cdf, col="red") 
    lines(logr, cdf, pch = ".", col="blue", lwd=2)
    
    # statistics
    sc <- 0.8    
    ell <- ellipticity(object)$avg
    
    q <- attr$tail_quantile
    asd <- attr$tail$asymp_stdev_0 # observed
    ask <- attr$tail$asymp_skewness_0 # observed
    aku <- attr$tail$asymp_kurt_0 # observed

    d_astats <- ecd.asymp_stats(object,q)[[1]] # theoretical
    d_asd <- d_astats$stdev # theoretical
    d_ask <- d_astats$skewness # theoretical
    d_aku <- d_astats$kurtosis # theoretical

    text(xmin*0.5,0.95,labels=c("data stats"),cex=sc)
    text(xmin*0.5,0.85,labels=paste("mean", sprintf("%.6f",attr$mean)),cex=sc)
    text(xmin*0.5,0.77,labels=paste("stdev",  sprintf("%.4f",attr$stdev)),cex=sc)
    text(xmin*0.5,0.69,labels=paste("skew", sprintf("%.4f",attr$skewness)),cex=sc)
    text(xmin*0.5,0.61,labels=paste("kurt", sprintf("%.4f",attr$kurtosis)),cex=sc)
    text(xmin*0.5,0.53,labels=paste("asymp stdev", sprintf("%.4f",asd)),cex=sc)
    text(xmin*0.5,0.45,labels=paste("asymp skew", sprintf("%.2f",ask)),cex=sc)
    text(xmin*0.5,0.37,labels=paste("asymp kurt", sprintf("%.1f",aku)),cex=sc)
    text(xmin*0.5,0.29,labels=paste("tail quant $e^\\wedge$", sprintf("%.1f",log(q))),cex=sc)
    
    text(xmax*0.5,0.75,labels=c("ecd fit"),cex=sc)
    text(xmax*0.5,0.60,labels=paste("mean", sprintf("%.6f",stats$mean)),cex=sc)
    text(xmax*0.5,0.53,labels=paste("stdev",  sprintf("%.4f",stats$stdev)),cex=sc)
    text(xmax*0.5,0.46,labels=paste("skew", sprintf("%.4f",stats$skewness)),cex=sc)
    text(xmax*0.5,0.39,labels=paste("kurt", sprintf("%.4f",stats$kurtosis)),cex=sc)
    text(xmax*0.5,0.32,labels=paste("asymp stdev", sprintf("%.4f",d_asd)),cex=sc)
    text(xmax*0.5,0.25,labels=paste("asymp skew", sprintf("%.2f",d_ask)),cex=sc)
    text(xmax*0.5,0.18,labels=paste("asymp kurt", sprintf("%.1f",d_aku)),cex=sc)
    text(xmax*0.5,0.11,labels=paste("ellipticity", sprintf("%.4f",ell)),cex=sc)
    
}
# ------------------------------------------------------------------
ecd.plot_qq <- function (object, ts, main="ECD QQ-Plot" ) {

    htu <- xtsAttributes(ts)$histuple
    merge_tails <- xtsAttributes(ts)$merge_tails

    x <- htu$hx
    hq <- cumsum(htu$hy) / (sum(htu$hy)+1)
    y <- qec(hq, object)
    
    par(mar=c(4,4,2,1))
 
    plot(x, y, 
         ylim=c(min(y)*1.2, max(y)*1.2),
         pch="o", col="red", 
         main=main, 
         xlab="Observed Quantile", 
         ylab="Theoretical Quantile")
    
    abline(0,1,col="black", lwd=2)
    abline(h=0,col="yellow")
    abline(v=0,col="yellow")
  
    lines(x, abs(y-x), pch = ".", col="green")

    # legends
    sc <- 0.8
    xmin <- min(x)
    ymax <- max(y)
    legend(xmin, ymax, 
           c("qq data","45 degree","error"), 
           cex=sc, 
           col=c("red","black","green"), 
           pch=c("o",NA,NA),
           lwd=c(1,2,1),
           lty=c(NA,1,1));

    xmax <- max(x)
    ymin <- min(y)
    if (sum(merge_tails)>0) {
        text(xmax*0.5, ymin*0.4, "tail dropped:",cex=sc);
        text(xmax*0.5, ymin*0.6, paste("left",  merge_tails[1]),cex=sc);
        text(xmax*0.5, ymin*0.8, paste("right", merge_tails[2]),cex=sc);
    }
}

