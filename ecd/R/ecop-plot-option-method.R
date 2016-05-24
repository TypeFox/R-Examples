#' Plot option chain charts using conf from option sample data
#'
#' This utility produces standardized plots of 3.
#' The first plot is the option state price and fits.
#' The second plot is the log-slope of option state prices and fits.
#' The thrid plot is the implied volatility and fits.
#'
#' @param object an ecop object with conf
#' @param otype option type
#' @param simulate logic, if \code{TRUE}, simulate according to lambda transformation and lambda distribution.
#'
#' @return The \code{ecop.opt} object
#'
#' @keywords plot
#'
#' @author Stephen H-T. Lihn
#'
#' @export
#'
#' @examples
#' \dontrun{
#'     op <- ecop.from_symbol_conf("spx2_1d")
#'     par(mfcol=c(3,2))
#'     ecop.plot_option(op, otype="c")
#'     ecop.plot_option(op, otype="p")
#' }
### <======================================================================>
"ecop.plot_option" <- function(object, otype, simulate=TRUE) {

    if (!(otype %in% c("c","p"))) {
        stop(paste("Unknown option type:", otype))
    }

    if (class(object) != "ecop") {
        stop(paste("Expect object to be ecop class, but got:", class(object)))
    }

    key <- object@key
    DAYS <- object@days
    DT <- object@datadate
    r <- object@int_rate
    T <- object@ttm
    y <- object@div_yield
  
    opd <- if (otype=="c") object@call_data else object@put_data 
    
    S <- opd@S
    K <- opd@strike
    V <- opd@V
    V_bid <- opd@V_bid
    V_ask <- opd@V_ask
    
    k <- opd@k

    imp_vol <- ecop.bs_implied_volatility(V, K, S, ttm=T, div_yield=y, otype=otype)
    imp_vol_bid <- ecop.bs_implied_volatility(V_bid, K, S, ttm=T, div_yield=y, otype=otype)
    imp_vol_ask <- ecop.bs_implied_volatility(V_ask, K, S, ttm=T, div_yield=y, otype=otype)

    # 
    momen <- opd@momentum
    epsilon <- opd@epsilon
    k_cusp <- opd@k_cusp
    ld1 <- opd@ecldOrEcd
    mu_D <- if (class(ld1)=="ecld") ecld.mu_D(ld1) else ecd.mu_D(ld1)

    # use a higher density for theoretical calculation
    k2 <- seq(min(k), max(k), length.out=100)
    V_fit <- ecop.polyfit_option(k, V, k_cusp, k.new=k2)

    V_ecop <- NULL
    V0_ecop <- NULL
    V2_ecop <- NULL
    ivol_ecop <- NULL
    ivol0_ecop <- NULL
    
    L  <- NULL # L is with epsilon's effect
    L0 <- NULL # L0 is without epsilon
    if (simulate) {
        if (class(ld1)=="ecld") {
            L0 <- ecld.ogf(ld1, k2, otype=otype, RN=FALSE)
        } else {
            L0 <- ecd.ogf(ld1, k2, otype=otype)
        }
        L <- L0 + epsilon
        
        V_ecop <- L*S
        V0_ecop <- L0*S
        ivol_ecop <- 1/sqrt(2*T) * ecld.op_V(L, k2, otype=otype)
        ivol0_ecop <- 1/sqrt(2*T) * ecld.op_V(L0, k2, otype=otype)
        
        # apply momentum operator 
        V2_ecop <- ecop.bs_option_price(ivol_ecop, exp(k2)*S*(1+momen), S, T,
                                        int_rate=r, div_yield=y, otype=otype)
    }

    # ready to plot
    oname <- if (otype=="c") "Call" else "Put"
    sc = 0.8
    
    vertical_1 <- function(x, segment_y=NaN, lty) {
        if (is.na(segment_y)) {
            abline(v=x, lty=lty)
        } else {
            segments(x, 0, x1=x, y1=segment_y, lty=lty)
        }
    }
    verticals <- function(segment_y=NaN) {
        vertical_1(x=ecd.mp2f(mu_D), segment_y, lty=3) # dotted
        vertical_1(x=ecd.mp2f(ld1@mu), segment_y, lty=2) # dashed
        vertical_1(x=ecd.mp2f(ld1@mu+log(1+momen)), segment_y, lty=1) # solid
    }

    # -----------------------------------------------------------
    # log-price chart
    plot(k, log(V), xlab="$k$", ylab=paste("$\\log(", toupper(otype), "(k))$", sep=""),
        main=paste(oname, "Option on", DT, sprintf("/ %dd", DAYS)))
    
    px <- c(V, V_bid, V_ask)
    px <- px[px > 0]
    sml_px <- sort(unique(px))[1:2] # smallest price quote
    abline(h=log(sml_px[1]))
    abline(h=log(sml_px[2]))
    verticals()
    
    lines(k, log(V_bid), type="o", pch=4, col="green")
    lines(k, log(V_ask), type="o", pch=2, col="green")

    if (simulate) {
        lines(k2, log(V_ecop), col="red", lwd=1)
        lines(k2, log(V0_ecop), col="red", lwd=1, lty=3) # dotted
        lines(k2+log(1+momen), log(V2_ecop), col="blue", lwd=3)
    }
    
    xmax <- max(k)
    xmin <- min(k)
    ymax <- max(log(V))
    dy <- (max(log(V)) - min(log(V)))/10
    
    ypos <- ymax
    if (otype=="p") ypos <- ymax - dy*2
    xpos <- xmax * 0.25
    text( xpos, ypos-dy,   paste(sprintf("S= %.2f", S)), cex=sc, pos=4 )
    text( xpos, ypos-dy*2, paste(sprintf("days= %d", DAYS)), cex=sc, pos=4 )
    text( xpos, ypos-dy*3, paste(sprintf("T= %d/365", T*365)), cex=sc, pos=4 )
    text( xpos, ypos-dy*4, paste(sprintf("$r_f$= %.2f$\\%%$", r*100)), cex=sc, pos=4 )
    text( xpos, ypos-dy*5, paste(sprintf("div yld= %.2f$\\%%$", y*100)), cex=sc, pos=4 )
    
    # legends
    xpos <- xmin*0.9
    if (otype=="c") {
        legend(xpos, ymax-dy*2, bty="n",
            c("bid","ask", "mid", 
              "local w/o $\\epsilon$", 
              "local w/ $\\epsilon$", 
              "$\\lambda$ trfm"),
            cex=sc,
            col=c("green", "green", "black", "red", "red", "blue"),
            pch=c(4, 2, 1, NA, NA, NA),
            lwd=c(NA, NA, NA, 1, 1, 3),
            lty=c(NA, NA, NA, 3, 1, 1)
        )
    }
    if (otype=="p") {
        text( xpos, ymax-dy, "vertical", cex=sc, pos=4 )
        legend(xpos, ymax-dy, bty="n",
            c("$\\mu_D$","$\\mu$", "$\\mu+r_M$"),
            cex=sc,
            col=c("black", "black", "black"),
            pch=c(NA, NA, NA),
            lwd=c(1, 1, 1),
            lty=c(3, 2, 1)
        )
    }

    # -----------------------------------------------------------
    sd1 <- if (class(ld1)=="ecld") ecld.sd(ld1) else ecd.sd(ld1)

    dlogV_fit   <- ecld.op_U_lag(V_fit, k2, sd1, 1)
    dlogV       <- ecld.op_U_lag(V, k, sd1, 2)
    dlogV_ecop  <- ecld.op_U_lag(V_ecop,  k2, sd1, 1)
    dlogV0_ecop <- ecld.op_U_lag(V0_ecop,  k2, sd1, 1)
    dlogV2_ecop <- ecld.op_U_lag(V2_ecop, k2, sd1, 1)
    all_dlogV   <- c(dlogV_fit, dlogV, dlogV_ecop, dlogV2_ecop)
    
    lVmin <- if (otype=="c") min(all_dlogV, na.rm=TRUE) else 0
    lVmax <- if (otype=="c") 0 else max(all_dlogV, na.rm=TRUE)

    # log-slope of price chart
    plot(k2, dlogV_fit, type="l", col="black", lwd=2, lty=3,
         xlab="$k$", ylab=paste("$U_", otype, "(k)$", sep=""),
         main=paste("Log-slope of Option Prices"), 
         ylim=c(lVmin, lVmax))
    points(k, dlogV, col="black")
    
    verticals()

    if (simulate) {
        lines(k2, dlogV_ecop, col="red", lwd=1)
        lines(k2, dlogV0_ecop, col="red", lwd=1, lty=3)
        lines(k2+log(1+momen), dlogV2_ecop, col="blue", lwd=3)
    }
    
    # legends
    xpos <- if (otype=="c") min(k)*0.9 else max(k)*0.4
    ypos <- if (otype=="c") lVmin*2/3 else lVmax*2/3

    legend(xpos, ypos, bty="n",
        c("poly fit"),
        cex=sc,
        col=c("black"),
        pch=c(NA),
        lwd=c(2),
        lty=c(3)
    )
    ypos <- if (otype=="c") lVmin/2 else lVmax*0.8
    text( xpos, ypos, paste(sprintf("$k_{cusp}$= %.5f", k_cusp)), cex=sc, pos=4 )


    # -----------------------------------------------------------
    # implied volatility
    min_ivol <- min(ivol_ecop, na.rm=TRUE)
    ymax <- max(c(imp_vol, imp_vol_bid, imp_vol_ask), na.rm=TRUE)
    ymin <- min(c(imp_vol, imp_vol_bid, imp_vol_ask), na.rm=TRUE)
    if (ymax > 0.4) {
        if (ymax < 0.8) ymax <- 0.8
    }
    if (ymin > 0.0) ymin <- 0.0
    
    plot(k, imp_vol, ylim=c(ymin, ymax), xlab="$k$", ylab="iVol",
         main=paste("Implied Volatility"))
    lines(k, imp_vol_bid, type="o", pch=4, col="green")
    lines(k, imp_vol_ask, type="o", pch=2, col="green")
    abline(h=0)
    verticals(segment_y=0.2)

    if (simulate) {
        lines(k2, ivol_ecop, col="red", lwd=1)
        lines(k2, ivol0_ecop, col="red", lwd=1, lty=3)
        lines(k2+log(1+momen), ivol_ecop, col="blue", lwd=3)
    }
    
    xpos <- -0.005
    if (otype=="p") xpos <- 0.3*min(k)
    dy <- (ymax-ymin)/12
    sged_str <- if (class(ld1)=="ecld") { if (ld1@is.sged) "(SGED)" else "" } else ""
    text( xpos, ymax-dy,   paste(sprintf("$r_M$= %.5f", momen)), cex=sc, pos=4 )
    text( xpos, ymax-dy*2, paste(sprintf("$\\mu$= %.5f", ld1@mu)), cex=sc, pos=4 )
    text( xpos, ymax-dy*3, paste(sprintf("$\\mu_D$= %.5f", mu_D)), cex=sc, pos=4 )
    
    if (class(ld1)=="ecld") {
        text( xpos, ymax-dy*4, paste(sprintf("$\\lambda$= %.3f", ld1@lambda)), cex=sc, pos=4 )
    } else {
        text( xpos, ymax-dy*4, paste(sprintf("$deg/R$= %.1f/ %.3f", ld1@theta/pi*180, ld1@R)), cex=sc, pos=4 )
    }
    text( xpos, ymax-dy*5, paste(sprintf("$\\beta$= %.3f %s", ld1@beta, sged_str)), cex=sc, pos=4 )
    text( xpos, ymax-dy*6, paste(sprintf("$\\sigma$= %.5f", ld1@sigma)), cex=sc, pos=4 )
    text( xpos, ymax-dy*7, paste(sprintf("$\\epsilon$= %.5f", epsilon)), cex=sc, pos=4 )

    text(0.8*min(k), min_ivol, paste(sprintf("min(iVol)= %.3f", min_ivol)), cex=sc, pos=4 )
    
    # -----------------------------------------------------------
    invisible(opd)
    
}
