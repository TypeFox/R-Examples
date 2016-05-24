aiDynFit <- function(w, dum.dif = FALSE, AR1 = FALSE, 
  rho.sel = c("all", "mean"), ...)
{
    if (!inherits(w, "aiStaFit")) {
      stop("Please provide an object from 'aiStaFit'.\n")}
    y <- w$y
    # lose two observations by diff and lag
    hShare <- bsLag(h = diff(y[, w$share]), lag = 1, prefix = "diff.") 
    hExpen <- diff(y[, w$expen])
    hPrice <- diff(y[, w$price])
    colnames(hPrice) <- paste("diff.", colnames(hPrice), sep="")
    
    # residuals from static model loses one obs if it has AR correction
    if (w$AR1) {begin <- start(y) + c(0,1)} else {begin <- start(y)}   
    resid <- ts(residuals(w$est), start = begin, frequency = tsp(y)[3])
    colnames(resid) <- paste("resi.", w$share[-w$nOmit], sep="")
    hResid <- bsLag(h = resid, lag = 1, include.orig = FALSE)   
              
    hShift <- y[, w$shift]
    if (!is.null(w$shift)) {shift <- w$shift} else {shift <- NULL}
    if (dum.dif) {
       if (!is.null(w$shift)) {
           hShift <- diff(hShift)
           shift <- paste("diff.", w$shift, sep="")}
    }
                           
    if (is.null(w$shift)) {
        comb <- ts.union(hShare, hResid, hExpen, hPrice)
        colnames(comb) <- c(colnames(hShare), colnames(hResid), 
            "diff.expen", colnames(hPrice))
    } else {                                        
        comb <- ts.union(hShare, hResid, hExpen, hPrice, hShift)
        colnames(comb) <- c(colnames(hShare), colnames(hResid), 
            "diff.expen", colnames(hPrice), shift)
    }   
    vaa <- list(hShare, hResid, hExpen, hPrice, hShift)
    sta <- c(tsp(hShare)[1], tsp(hResid)[1], tsp(hExpen)[1],
             tsp(hPrice)[1], tsp(hShift)[1])      
    loc <- which(sta == max(sta))  # there may be multiple matches      
    beg <- start(vaa[[loc[1]]])
    daDyn <- window(comb, start = beg, end = end(y), frequency = tsp(y)[3])

    share.d  <- paste("diff.", w$share, ".t_0", sep="")
    share.dl <- paste("diff.", w$share, ".t_1", sep="")
    resid.d  <- colnames(hResid)  
    expen.d  <- "diff.expen"
    price.d  <- colnames(hPrice)
    nOmit    <- w$nOmit
    omit     <- share.d[nOmit]
    
    nShare <- length(share.d)
    nExoge <- 0 + 3 + length(shift)
    nParam <- nExoge + nShare
    nTotal <- (nShare - 1) * nParam
    m.h <- matrix(0, (nShare-1), nTotal)
    for (i in 1:(nShare - 1)) {
        for (j in 1:nShare) { m.h[i, (i-1)*nParam + nExoge + j] <- 1 }
    }   
    ns  <- (nShare - 1) * (nShare - 2)/2 
    k   <- 0
    m.s <- matrix(0, ns, nTotal)
    for ( i in 1:(nShare - 2) ) {
        for ( j in (i+1):(nShare-1) ) {
            k <- k + 1
            m.s[k, (i-1)*nParam + nExoge + j] <-  1
            m.s[k, (j-1)*nParam + nExoge + i] <- -1
        }
    }
    m.hs <- rbind(m.h,m.s)
    r.s  <- rep(0, nrow(m.s))    
    r.h  <- rep(0, nrow(m.h))
    r.hs <- rep(0,nrow(m.hs))
    if (!w$hom & !w$sym) {aa <- NULL; bb <- NULL}
    if ( w$hom & !w$sym) {aa <- m.h ; bb <- r.h }
    if (!w$hom &  w$sym) {aa <- m.s ; bb <- r.s }
    if ( w$hom &  w$sym) {aa <- m.hs; bb <- r.hs} 
    
    sa <- list()
    for (i in 1:(length(share.d) - 1)) {
       sa[[i]] <- bsFormu(intercept = FALSE, 
           name.y = share.d[-nOmit][i], 
           name.x = c(share.dl[-nOmit][i], resid.d[i], shift, expen.d, price.d))
    }
    est <- systemfitAR(formula = sa, data = data.frame(daDyn), method="SUR",
       restrict.matrix = aa, restrict.rhs = bb, AR1 = AR1, rho.sel = rho.sel,
       model="dynamic")
       
    ## Omitted equation: coefficients
#    co <- matrix(0, nrow = nShare - 1, ncol = nParam)
#    for (i in 1:(nShare - 1)) {
#      co[i, ] <- coef(est)[((i-1) * nParam + 1):(i * nParam)]
#    }
#    co.omit <- rep(0, nParam) - colSums(co) # diff in dynamic model so all sum 0
#   
#    # Omitted equation: variance for gamma_last and exoge variables
#    cof <- coef(est);  vco <- vcov(est)
#    gama.vc <- 0
#    for (i in 1:(nShare-1)) { 
#      for (j in 1:(nShare-1)) {
#        gama.vc <- gama.vc + vco[nParam*i, nParam*j]
#      }
#    }
#    ex.vc <- rep(0, nExoge)
#    for(k in 1:nExoge) {
#      for (i in 1:(nShare-1)) { 
#        for (j in 1:(nShare-1)) {
#          ex.vc[k] <- ex.vc[k] + vco[nParam*(i-1) + k, nParam*(j-1) + k]
#        }
#      }
#    }
#    
#    # Omitted equation: combined
#    df <- df.residual(est)
#    c.ex <- co.omit[c(1:nExoge, nParam)]
#    e.ex <- sqrt(c(ex.vc, gama.vc))
#    t.ex <- c.ex / e.ex
#    p.ex <- 2 * ( 1- pt(abs(t.ex), df) )
#    ex <- data.frame(cbind(c.ex, e.ex, t.ex, p.ex))
    
    result <- listn(w, y=w$y, dum.dif, daDyn,
        share=share.d, price=price.d, expen=expen.d, 
        shift, omit, nOmit, hom=w$hom, sym=w$sym, nShare, nExoge, 
        nParam, nTotal, formula=sa, res.matrix=aa, res.rhs=bb, est=est)
    class(result) <- c("aiDynFit", "aiFit")
    return(result)
}