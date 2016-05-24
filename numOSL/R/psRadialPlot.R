#####
psRadialPlot <- 
function(EDdata, addsigma=0, dose=NULL, zmin=NULL, zmax=NULL, ntick=6, 
         digits=2, pcolor="blue", psize=1, rg=2, zlabel="De (Gy)") {
    UseMethod("psRadialPlot")
} #
### 2015.05.03.
psRadialPlot.default <- 
function(EDdata, addsigma=0, dose=NULL, zmin=NULL, zmax=NULL, ntick=6, 
         digits=2, pcolor="blue", psize=1, rg=2, zlabel="De (Gy)") {
    ### Stop if not.
    stopifnot(ncol(EDdata)==2L, nrow(EDdata)>=5L,
              all(EDdata[,1L,drop=TRUE]>0), all(EDdata[,2L,drop=TRUE]>0),
              length(addsigma)==1L, is.numeric(addsigma), addsigma>=0,
              is.null(dose) || (is.numeric(dose) && all(dose>0)),
              is.null(zmin) || (length(zmin)==1L && is.numeric(zmin) && zmin>0), 
              is.null(zmax) || (length(zmax)==1L && is.numeric(zmax) && zmax>0), 
              is.numeric(ntick), length(ntick)==1L, ntick %in% (2L:10L),
              is.numeric(digits), length(digits)==1L, digits %in% c(0L:5L),
              is.character(pcolor), length(pcolor)==1L, 
              is.numeric(psize), length(psize)==1L, 
              is.numeric(rg), length(rg)==1L, rg %in% (0L:2L),
              is.character(zlabel), length(zlabel)==1L) 
    ###
    if(!is.null(zmin) && !is.null(zmax)) {
        if(zmin>=zmax) {
            stop("Error: zmin should be smaller than zmax!")
        } # end if.
    } # end if.
    if(!is.null(dose) && !is.null(zmin)) {
        if(any(dose<zmin)) {
            stop("Error: all dose values should be larger than zmin!")
        } # end if.
    } # end if.
    if (!is.null(dose) && !is.null(zmax)) {
        if(any(dose>zmax)) {
            stop("Error: all dose values should be smaller than zmax!")
        } # end if. 
    } # end if.
    ###
    ###
    ndat <- nrow(EDdata)
    ed1 <- as.numeric(EDdata[,1L,drop=TRUE])
    sed1 <- as.numeric(EDdata[,2L,drop=TRUE])
    ###
    ### Transform to log-scale.
    x <- sqrt( (sed1/ed1)^2L + addsigma^2L)
    y <- log(ed1)
    centralDose <- sum(y/x^2L)/sum(1.0/x^2L)
    ###
    xx <- 1.0/x
    yy <- (y-centralDose)/x
    ###
    ### Range of precision.
    minPrecision <- 0
    maxPrecision <- max(xx)*1.3
    ###
    ### Range of z-axis.
    if (is.null(zmin)) {
        zmin <- min(ed1)*0.9
    } # end if.
    if (is.null(zmax)) {
        zmax <- max(ed1)*1.1
    } # end if.
    ###
    ### transform z-scale values to y-scale values.
    miny <- (log(zmin)-centralDose)*maxPrecision
    maxy <- (log(zmax)-centralDose)*maxPrecision
    ###
    ### Set plot margin.
    par("mar"=c(4.1, 4.1, 2.1, 5.1))
    ### 
    ###
    plot(xx, yy, xaxt="n", yaxt="n", bty="n", xlim=c(minPrecision, maxPrecision), 
         xaxs="i", yaxs="i", ylim=c(miny, maxy), xlab="", ylab="", cex.lab=1.0, type="n")
    locx <- axTicks(side=1L)
    locx <- seq(from=min(locx), to=max(locx), by=(max(locx)-min(locx))/5L)
    locy <- seq(from=miny, to=maxy, by=(maxy-miny)/(ntick-1L))
    ###
    ### Set relative standard error.
    axis(side=1L, at=locx[-1L], lwd=2.0, labels=round(100.0/locx[-1L],1L), line=0, tck=0.02, padj=-4)
    mtext(text="Relative standard error", side=1L, cex=1.0, line=-3)
    ### Set precision (x-axis).
    axis(side=1L, at=locx, lwd=2.0, labels=locx, line=0, padj=0)
    mtext(text="Precision", side=1L, cex=1.0, line=2)
    ###
    ### Set y-axis.
    axis(side=2L, at=c(-2L,0L,2L), lwd=2.0, labels=FALSE)
    mtext(text="-2", side=2L, at=-2.0, cex=0.8, line=1)
    mtext(text="0",side=2L, at=0.0, cex=0.8, line=1)
    mtext(text="2", side=2L, at=2.0, cex=0.8, line=1)
    mtext(text="Standardised Estimate", side=2L, at=0.0, cex=1.0, line=2)
    ###
    ### Set z-axis.
    axis(side=4L, at=locy, labels=round(exp(locy/maxPrecision+centralDose),digits=digits),lwd=2.0)
    abline(v=maxPrecision, lwd=2.0)
    mtext(text=zlabel, side=4L, cex=1.0, line=3) 
    ###
    ###
    if (!is.null(dose))  {
        ### Plot polygons.
        for (i in seq(dose)) {      
            y1 <- (log(dose[i])-centralDose)*maxPrecision
            if (rg %in% c(1L, 2L))  {
                polygon(x=rep(c(minPrecision, maxPrecision), each=2L),
                        y=c(-rg, rg, y1+rg, y1-rg), 
                        col="grey90", lty="blank")
            } # end if.
        } # end for.
        ###
        ### Plot lines.
        for (i in seq(dose)) {      
            y1 <- (log(dose[i])-centralDose)*maxPrecision
            lines(c(minPrecision, maxPrecision), c(0.0, y1), lwd=1.0)
            if (rg %in% c(1L, 2L))  {
                lines(c(minPrecision, maxPrecision), c(-rg, y1-rg), lwd=0.3, lty="dashed")
                lines(c(minPrecision, maxPrecision), c(rg, y1+rg), lwd=0.3, lty="dashed")
            } # end if.
        } # end for.
    } # end if.
    ###
    ###
    points(xx, yy, pch=21, col="black", bg=pcolor, cex=psize)
    ###
    par("mar"=c(5.1, 4.1, 4.1, 2.1))
} # end function psRadialPlot.
#####
