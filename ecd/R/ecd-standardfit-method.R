# ------------------------------------------
# Fitting Functions
# ------------------------------------------
# optimx: http://cran.r-project.org/web/packages/optimx/optimx.pdf
# spg: http://cran.r-project.org/web/packages/BB/BB.pdf

ecd.standardfit <- function(object, ts, 
                            trace = 0, iter = 5000, 
                            plotqq = 1, weights=list()) 
{

    print(paste("standardfit:",
                "iter=",iter,"plotqq=",plotqq,
                "weights= (",paste(weights, collapse = ','),")"))

    # vectorize ecd object with mu first, and sigma second
    # both normalized with stdev of ts to make them signoficant enough
    attr <- xtsAttributes(ts)
    dv <- attr$stdev
    ecdc <- c(object@mu/dv*1000, object@sigma/dv*10, object@alpha, object@gamma, object@beta)

    opt.out <- optimx::optimx(ecdc, ecd.standardfit_fn, 
                      method = c("spg"), 
                      itnmax = iter, control = list(trace=trace,eps=1e-5), 
                      ecd_ts = ts, 
                      plotqq = plotqq, 
                      weights = weights)

    delta <- opt.out$value
    msg <- opt.out$message
    p <- opt.out
    print(paste("Finished fit: detla=", delta, "par alpha= ", p$p3))
    
    dist <- ecd(mu=p$p1*dv/1000, sigma=p$p2*dv/10, alpha=p$p3, gamma=p$p4, beta=p$p5)
    list(dist = dist, opt.out = opt.out)
}

# ------------------------------------------
# TODO to be documented and formalized
# ------------------------------------------
ecd.standardfit_qa <- function(object, ts, plotqq=1, weights=list(), debug=1)
{
    attr <- xtsAttributes(ts)
    dv <- attr$stdev
  	ecdc <- c(object@mu/dv*1000, object@sigma/dv*10, object@alpha, object@gamma, object@beta)
  	print("Entering ecd.standardfit_fn")
    ecd.standardfit_fn(ecdc, ts, plotqq=plotqq, weights=weights, debug=debug)
}


# This function calculates delta by comparing the data to the theoretical ecd
# Delta is the sum of the following metrics with user specified weights
#   m1, m2, skewness, kurtosis
#   diff on pdf, diff on QQ-line

ecd.standardfit_fn <- function(ecdc, ecd_ts, 
	                           plotqq=1, weights=list(), debug=0)
{

    tm0 <- unclass(Sys.time())

    # construct ecd object from vector
    ts <- ecd_ts # rename it to the simple ts
    attr <- xtsAttributes(ts)
    hist <- attr$hist
    
    dv <- attr$stdev
    d <- ecd(mu= ecdc[1]*dv/1000, sigma= ecdc[2]*dv/10, alpha= ecdc[3], gamma= ecdc[4], beta= ecdc[5] )
    d <- quantilize(d)
    
    # forbidden area
    if ( abs(d@beta) >= 10.0 ) {
        delta <- if ( debug == 1 ) list( delta= abs(d@beta)*1000 ) else abs(d@beta)*1000
        return(delta)
    } 

    # default weights
    if (is.null(weights$NS_mmx)) weights$NS_mmx <- 1
    if (is.null(weights$NS_pdf)) weights$NS_pdf <- 2

    if (is.null(weights$m1)) weights$m1 <- 100
    if (is.null(weights$m2)) weights$m2 <- 5
    if (is.null(weights$m3)) weights$m3 <- 5
    if (is.null(weights$m4)) weights$m4 <- 5
    if (is.null(weights$mmx)) weights$mmx <- 1 # diff on peak pdf
    if (is.null(weights$pdf_df)) weights$pdf_df <- 0.5 # diff on log(pdf)
    if (is.null(weights$qq_df))  weights$qq_df  <- 0.5 # diff on qq
    if (is.null(weights$discr_bump)) weights$discr_bump <- 1
  
    tm1 <- unclass(Sys.time()) - tm0
    # --------------------------------

    if (length(d@stats)==0) {
        stop("stats is not computed in ecd object (d)")
    }


    st <- d@stats
    NS <- weights$NS_mmx # number of stdev for xs range for mmx
    xs <- seq( st$mean-NS*st$stdev, st$mean+NS*st$stdev, by=2*NS*st$stdev/200 )
    ds <-ecd.pdf(d,xs) # pdf of ecd

    max_hist_density <- max(hist$density)
    mmx <- abs( max(ds)-max_hist_density ) / max_hist_density

    tm2 <- unclass(Sys.time()) - tm0
    # --------------------------------
    
    # calculate pdf diff within NS*sigma
    NS <- weights$NS_pdf # number of stdev for log pdf_df
    pdf_df <- 0.0
    xpdf <- c() 
    tpdf <- c()
    j <- 1
    for ( i in 1:length(hist$mids) ) {
        x <- hist$mids[i]
        if ( abs(x-st$mean) <= NS*st$stdev ) {
            xpdf[j] <- x
            tpdf[j] <- ecd.pdf(d, x)
            df <- (log(hist$density[i]) - log(tpdf[j]))^2
            pdf_df <- pdf_df + if ( abs(df) <= 1000.0 ) df else 0
            j <- j + 1
        }
    }
    pdf_df <- sqrt(abs(pdf_df))
    pdf_df <- if ( pdf_df == Inf ) 10000.0 else pdf_df

    tm3 <- unclass(Sys.time()) - tm0
    # --------------------------------
    
    qqp_pct <- 0.10 # 10%
    htu <- attr$histuple

    # ------------------------------------------
    # -- qqp object --
    # x: data's x
    # xq: data's cdf
    # y: x of the theoretical distribution
    # yq: cdf(y) of the theoretical distribution
    # The QQ-plot is y vs x
    
    qqp <- list (
        x = htu$hx,
        xq = cumsum(htu$hy) / (sum(htu$hy)+1)
    )
    qqp$y <- qec(qqp$xq, d)
    qqp$yq <- ecd.cdf(d, qqp$y)

    Nq <- length(qqp$x)
    qqp_rng <- seq(floor(Nq*qqp_pct)+1, floor(Nq*(1-qqp_pct)))
    qq_df <- sqrt(sum(abs(qqp$y[qqp_rng]-qqp$x[qqp_rng])^2))/max(hist$mids)

    tm4 <- unclass(Sys.time()) - tm0
    # --------------------------------
    
    discr_bump <- ifelse(discr(d,no.validate=TRUE) <= 0 & d@alpha>=0 & d@gamma<=0,
                         (abs(d@alpha) + abs(d@gamma))*1000, 0)
                         
    m1 <- abs( st$mean - attr$mean ) / attr$stdev # stdev pct
    m2 <- abs( st$stdev - attr$stdev ) / attr$stdev # stdev pct
    m3 <- abs( st$skewness - attr$skewness )
    m4 <- abs( st$kurtosis - attr$kurtosis ) / attr$kurtosis # kurt pct

    delta <- 0
    delta <- delta + (m1 * weights$m1)^2
    delta <- delta + (m2 * weights$m2)^2
    delta <- delta + (m3 * weights$m3)^2
    delta <- delta + (m4 * weights$m4)^2
    delta <- delta + (discr_bump * weights$discr_bump)^2
    delta <- delta + (mmx * weights$mmx)^2
    delta <- delta + (qq_df * weights$qq_df)^2
    delta <- delta + (pdf_df * weights$pdf_df)^2
    delta <- sqrt(delta)


    if ( plotqq == 1 ) {
        par(mfcol=c(2,2)); # 4 plots
        par(oma = c(1, 1, 4, 1))

        # --------------------
        par(mar=c(4,4,2,2)) 
        
        xmin <- min(c(min(xs), min(xpdf)*1.2))
        xmax <- max(c(max(xs), max(xpdf)*1.2))
        
        plot(hist$mids, hist$density, 
             type="o", pch=22, lty=2, col="red", 
             main="PDF", xlab="r", ylab="PDF",
             xlim=c(xmin,xmax),
             ylim=c(0, max(hist$density)*1.3))
        
        lines(xs, ds, pch = ".", col="black", lwd=4)
        lines(xpdf, tpdf, pch = ".", col="green")

        # --------------------
        par(mar=c(4,4,2,2)) 
        plot(hist$mids, log(hist$density), 
             type="o", pch=22, lty=2, col="red", 
             main="Log PDF", xlab="log(r)", ylab="Log PDF")
        
        lines(xpdf, log(tpdf), pch = ".", col="black", lwd=4)
        lines(hist$mids, log(ecd.pdf(d,hist$mids)), pch = ".", col="green")
        
        xmin <- min(hist$mids)
        xmax <- max(hist$mids)
        xpos <- xmin*0.5
        ypos <- max(log(hist$density))
        sc <- 0.8
        text(xmin*0.5,ypos,labels=c("ecd fit"),cex=sc)
        text(xmin*0.5,ypos-1,labels=paste("mean", sprintf("%.6f",st$mean)),cex=sc)
        text(xmin*0.5,ypos-2,labels=paste("std",  sprintf("%.6f",st$stdev)),cex=sc)
        text(xmin*0.5,ypos-3,labels=paste("skew", sprintf("%.6f",st$skewness)),cex=sc)
        text(xmin*0.5,ypos-4,labels=paste("kurt", sprintf("%.6f",st$kurtosis)),cex=sc)
        text(xmin*0.5,ypos-5,labels=paste("discr", sprintf("%.1f",discr(d,no.validate=T))),cex=sc)
        text(xmin*0.5,ypos-6,labels=paste("j-inv", sprintf("%.1f",jinv(d,no.validate=T))),cex=sc)
        
        text(xmax*0.5,ypos,labels=c("data"),cex=sc)
        text(xmax*0.5,ypos-1,labels=paste("mean", sprintf("%.6f",attr$mean)),cex=sc)
        text(xmax*0.5,ypos-2,labels=paste("stdev",  sprintf("%.4f",attr$stdev)),cex=sc)
        text(xmax*0.5,ypos-3,labels=paste("skew", sprintf("%.4f",attr$skewness)),cex=sc)
        text(xmax*0.5,ypos-4,labels=paste("kurt", sprintf("%.4f",attr$kurtosis)),cex=sc)
        
        # ----------------------------------
        ymin <- min(qqp$y) - 0.1*abs(min(qqp$y))
        ymax <- max(qqp$y)*1.5
        
        par(mar=c(4,4,2,2))
        plot(qqp$x, qqp$y, pch = "o", col="red", 
            main=paste("QQ Plot",Sys.time()), 
            ylim=c(ymin, ymax),
            xlab="Observed Quantile", ylab="Theoretical Quantile")
    
        abline(h=0,col="black")
        abline(v=0,col="black")
    
        lines(qqp$x[qqp_rng], abs(qqp$y[qqp_rng]-qqp$x[qqp_rng])*4, 
              pch = ".", col="green",lwd=3)
        abline(0,1,col="blue")

        tm9 <- unclass(Sys.time()) - tm0

        xmin <- min(qqp$x)
        xpos <- xmin*0.5
        ypos <- max(qqp$y)*1.5
        ydf <- (max(qqp$y)-min(qqp$y))/10
        sc <- 0.8

        text(xpos,ypos-1*ydf,labels=c("timing"),cex=sc)
        text(xpos,ypos-2*ydf,labels=paste("tm1",   sprintf("%.2f",tm1)),cex=sc)
        text(xpos,ypos-3*ydf,labels=paste("tm2",   sprintf("%.2f",tm2)),cex=sc)
        text(xpos,ypos-4*ydf,labels=paste("tm3",   sprintf("%.2f",tm3)),cex=sc)
        text(xpos,ypos-5*ydf,labels=paste("tm4",   sprintf("%.2f",tm4)),cex=sc)
        text(xpos,ypos-6*ydf,labels=paste("tm9",   sprintf("%.2f",tm9)),cex=sc)

        # -----------------------------
        par(mar=c(4,4,2,2)) 
        plot(qqp$x, qqp$y, pch = "o", col="green", 
            main=paste("QQ Data",Sys.time()), 
            ylim=c( min(qqp$y), max(qqp$y)*1.5 ),
            xlab="Observed Quantile", ylab="Theoretical Quantile")

        text(xpos,ypos-1*ydf,labels=c("ecd fit"),cex=sc)
        text(xpos,ypos-2*ydf,labels=paste("delta",  sprintf("%.8f",delta)),cex=sc)
        text(xpos,ypos-4*ydf,labels=paste("alpha",  sprintf("%.8f",d@alpha)),cex=sc)
        text(xpos,ypos-5*ydf,labels=paste("gamma",  sprintf("%.8f",d@gamma)),cex=sc)
        text(xpos,ypos-6*ydf,labels=paste("sigma",  sprintf("%.8f",d@sigma)),cex=sc)
        text(xpos,ypos-7*ydf,labels=paste("beta",   sprintf("%.8f",d@beta)),cex=sc)
        text(xpos,ypos-8*ydf,labels=paste("mu",     sprintf("%.8f",d@mu)),cex=sc)
        text(xpos,ypos-9*ydf,labels=paste("ecdc1",  sprintf("%.6f",ecdc[1])),cex=sc)
        text(xpos,ypos-10*ydf,labels=paste("ecdc2", sprintf("%.6f",ecdc[2])),cex=sc)
        
        xmax <- max(qqp$x)
        xpos <- xmax*0.4
        text(xpos,ypos-1*ydf,labels=c("delta parts"),cex=sc)
        text(xpos,ypos-2*ydf,labels=paste("max_pdf", sprintf("%.6f",max_hist_density)),cex=sc)
        text(xpos,ypos-3*ydf,labels=paste("mmx", sprintf("%.6f %.2f",mmx, weights$mmx)),cex=sc)
        text(xpos,ypos-4*ydf,labels=paste("m1",  sprintf("%.6f %.2f",m1, weights$m1)),cex=sc)
        text(xpos,ypos-5*ydf,labels=paste("m2",  sprintf("%.6f %.2f",m2, weights$m2)),cex=sc)
        text(xpos,ypos-6*ydf,labels=paste("m3",  sprintf("%.6f %.2f",m3, weights$m3)),cex=sc)
        text(xpos,ypos-7*ydf,labels=paste("m4",  sprintf("%.6f %.2f",m4, weights$m4)),cex=sc)
        text(xpos,ypos-8*ydf,labels=paste("discr",  sprintf("%.6f %.2f",discr_bump, weights$discr_bump)),cex=sc)
        text(xpos,ypos-9*ydf,labels=paste("pdf_df", sprintf("%.6f %.2f",pdf_df, weights$pdf_df)),cex=sc)
        text(xpos,ypos-10*ydf,labels=paste("qq_df",  sprintf("%.8f %.2f",qq_df, weights$qq_df)),cex=sc)

    } # end of qq plot
  
    # ------------------------------------------
    # return the delta as a measure of deviation
    if ( debug == 1 ) list( 
        delta= delta,
        mmx= mmx, m2=m2, m3=m3, m4=m4, 
        discr_bump=discr_bump,
        max_density=c(max(ds), max_hist_density), 
        qq_df= qq_df ) 
    else delta
}

