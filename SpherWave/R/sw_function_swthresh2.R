"swthresh" <-
function (swd, policy, by.level, type, nthresh, value = 0.1, Q = 0.05) 
{
    if (class(swd) != "swd") 
        stop("Input must be swd class.")
        
    if (type != "soft" && type != "hard" && type != "Lorentz")
        stop("Thresholding type must be soft, hard or Lorentz.")
    if (by.level)
        if (policy != "universal" && policy != "Lorentz"&& policy != "fdr")
            stop("Only universal, Lorentz and fdr are provided for level-dependent thresholding.")
    else 
        if (policy != "universal" && policy != "probability" && policy != "fdr" && policy != "Lorentz" && policy != "sure")
            stop("Only universal, probability, Lorentz, fdr and sure are provided for global thresholding.")
                
    if ((policy == "probability" && value < 0) || (policy == "probability" && value > 1) || (policy == "probability" && is.null(value))) 
        stop("Value of quantile must be provided properly.")
    if ((policy == "fdr" && Q < 0) || (policy == "fdr" && Q > 1) || (policy == "fdr" && is.null(Q))) 
        stop("A false discovery rate must be provided properly.")
 
    if (policy == "sure" & type != "soft")     
        stop("Can only do soft thresholding with sure policy.")  
        
    site <- cbind(swd$latlon[, 1] * 0.5 * pi/90, swd$latlon[, 2] * pi/180)
    
    grid <- NULL
    grid$theta <- swd$gridlon * pi/180
    grid$phi <- swd$gridlat * pi/180  
    
    rangedetail <- NULL
    for (i in 1:(swd$nlevels - 1))
        rangedetail <- rbind(rangedetail, range(swd$detail[[i]]))
        
    if (by.level) {
        thresh.detail <- mrsfield.comp.thresh.level(grid, swd$coeff, site, swd$netlab, swd$eta, nthresh, policy, Q, type)
    } else {
        lam <- lambda.global(swd, policy, nthresh, value, Q)$lam
        thresh.detail <- mrsfield.comp.thresh.global(grid, swd$coeff, site, swd$netlab, swd$eta, lam, nthresh, type)
    }

    for (i in 1:(swd$nlevels - 1)) 
        swd$detail[[i]] <- t(matrix(thresh.detail[, i], nrow=length(swd$gridlat))) 
    
    thresh.info <- list(policy=policy, by.level=by.level, type=type, nthresh=nthresh, rangedetail=rangedetail)
    if (policy == "probability") thresh.info <- c(thresh.info, list(value=value))
    if (policy == "fdr") thresh.info <- c(thresh.info, list(Q=Q))
                        
    out <- list(obs=swd$obs, latlon=swd$latlon,
        netlab=swd$netlab, eta=swd$eta, method=swd$method, approx=swd$approx, grid.size=swd$grid.size, lambda=swd$lambda, p0=swd$p0,
        gridlat=swd$gridlat, gridlon=swd$gridlon, nlevels=swd$nlevels, coeff=swd$coef, field=swd$field, density1=swd$density,
        latlim=swd$latlim, lonlim=swd$lonlim, global=swd$global, density=swd$density, detail=swd$detail, swcoeff=swd$swcoeff,
        thresh.info=thresh.info)
     
    class(out) <- "swd"  
    out        
}

"mrsfield.comp.thresh.level" <-
function (grid, coef, site, netlab, eta, K, policy, Q, type) 
{
    dfield <- NULL
    alpha <- mrs.comp.thresh.level(coef, site, netlab, eta, K, policy, Q, type)$talpha
    for (l in 1:ncol(alpha)) {
        temp <- msbf.comp(grid, site, alpha[, l], netlab, eta, l - 1)
        dfield <- cbind(dfield, temp$field)
    }
    dfield
}

mrs.comp.thresh.level <-
function (coef, site, netlab, eta, K, policy, Q, type) 
{
    J <- length(netlab)
    L <- length(eta)
    alp <- NULL
    blp <- NULL
    a1 <- mcov.comp(site, netlab, eta)
    beta1 <- coef
    mrs <- coefmatrix(beta1, a1, netlab, 1)
    norm <- sqrt(diag(mrs$norm))
    thr <- thresh.level(mrs$gamma1, norm, policy, Q, type)
    talpha <- t(mrs$wcoef) %*% thr$tgamma
    beta2 <- mrs$beta2
    if (L == 2)
        talpha <- as.matrix(talpha)
    else {
        if (K == (L - 1))
            for (l in 2:K) {
                temp <- coefmatrix(beta2, a1, netlab, l)
                norm <- sqrt(diag(temp$norm))
                thr <- thresh.level(temp$gamma1, norm, policy, Q, type)
                talpha1 <- t(temp$wcoef) %*% thr$tgamma
                beta2 <- temp$beta2
                temp <- rep(NA, J)
                temp[netlab > (l - 1)] <- talpha1
                alp <- cbind(alp, temp)
            }
        else {
            for (l in 2:K) {
                temp <- coefmatrix(beta2, a1, netlab, l)
                norm <- sqrt(diag(temp$norm))
                thr <- thresh.level(temp$gamma1, norm, policy, Q, type)
                talpha1 <- t(temp$wcoef) %*% thr$tgamma
                beta2 <- temp$beta2
                temp <- rep(NA, J)
                temp[netlab > (l - 1)] <- talpha1
                alp <- cbind(alp, temp)
            }
            for (l in (K + 1):(L - 1)) {
                tmp <- coefmatrix(beta2, a1, netlab, l)
                talpha1 <- tmp$alpha1
                beta2 <- tmp$beta2
                tmp <- rep(NA, J)
                tmp[netlab > (l - 1)] <- talpha1
                blp <- cbind(blp, tmp)
            }
        }
    }
    list(talpha = cbind(talpha, alp, blp))
}

#######################################
#thresh.level.gamma(x, norm)  -> thresh.level(x, norm, policy="univ", Q=NULL, type) 
#fdr.level.gamma(x, norm, al) -> thresh.level(x, norm, policy="fdr", Q, type)
#Lorentz.level.gamma(x, norm) -> thresh.level(x, norm, policy="Lorentz", Q=NULL, type) 

thresh.level <- function(x, norm, policy, Q, type) {

    if (type != "soft" && type != "hard" && type != "Lorentz")
        stop("Thresholding type must be soft, hard or Lorentz.")
        
    if (policy != "universal" && policy != "fdr" && policy != "Lorentz")
        stop("Policy must be universal, fdr or Lorentz.")
       
    if (policy == "fdr" && is.null(Q))        
        stop("Provide the parameter Q for fdr.")
                
    if (policy == "universal") {
        tmp <- as.vector(na.omit(x)) * norm
        ml <- length(tmp)
        sigma <-  mad(tmp)
        lam <- sqrt(2 * log(ml)) * sigma
        if (type == "soft")
            x <- sign(x) * (abs(x) - lam/norm) * (abs(x) > lam/norm)
        if (type == "hard")
            x <- x * (abs(x) > lam/norm)
        if (type == "Lorentz")
            x <- sign(x) * sqrt(ifelse(abs(x) > lam/norm, x^2, (lam/norm)^2) - (lam/norm)^2) 
    }
    if (policy == "fdr") {
        ### Adapt threshold.wd in wavethresh4_4.0-2
        d <- as.vector(na.omit(x)) * norm
        m <- length(d)
        noise.level <- mad(d)
        thinit <- qnorm(1 - Q/2) * noise.level
        dinit <- d[abs(d) >= thinit]
        minit <- length(dinit)
        if (minit == 0) 
            lam <- max(abs(d)) * 1.0001
        else {
            if (noise.level > 0) {
                p <- (2 - 2 * pnorm(abs(dinit)/noise.level))
                index <- order(p)
                j <- seq(1, minit)
                m0 <- j[p[index] <= (Q * j)/m]
                if (length(m0) != 0)
                    m0 <- max(m0)
                if (length(m0) != 0  && m0 < minit) 
                    lam <- abs(dinit[index[m0]])
                else {
                    if (length(m0) == 0) 
                        lam <- max(abs(dinit)) * 1.0001
                    else lam <- 0
                }          
            }
            else lam <- 0  
        }
        if (type == "soft")
            x <- sign(x) * (abs(x) - lam/norm) * (abs(x) > lam/norm) 
        if (type == "hard")
            x <- x * (abs(x) > lam/norm)
        if (type == "Lorentz")
            x <- sign(x) * sqrt(ifelse(abs(x) > lam/norm, x^2, (lam/norm)^2) - (lam/norm)^2) 
    }
    if (policy == "Lorentz") {
        tmp <- as.vector(na.omit(x)) * norm
        ml <- length(tmp)
        lam <- sqrt(sum(tmp * tmp))/sqrt(ml)
        if (type == "soft")
            x <- sign(x) * (abs(x) - lam/norm) * (abs(x) > lam/norm) 
        if (type == "hard")
            x <- x * (abs(x) > lam/norm)
        if (type == "Lorentz")
            x <- sign(x) * sqrt(ifelse(abs(x) > lam/norm, x^2, (lam/norm)^2) - (lam/norm)^2) 
    }
        
    list(tgamma = x)    
   
}

#######################################
#thresh.total.mrs.comp(coef, site, netlab, eta, K) -> mrs.comp.total(coef, site, netlab, eta, K, lam=NULL, policy="univ", type="soft")
#fdr.total.mrs.comp(coef, site, netlab, eta, lam, K) -> mrs.comp.total(coef, site, netlab, eta, K, lam=NULL, policy="fdr", type="soft")
#Lorentz.total.mrs.comp(coef, site, netlab, eta, lam, K) -> mrs.comp.total(coef, site, netlab, eta, K, lam=NULL, policy="Lorentz", type="Lorentz")
#cv.total.mrs.comp(coef, site, netlab, eta, lam, K) -> mrs.comp.total(coef, site, netlab, eta, K, lam=NULL, policy="sure", type="soft")

lambda.global <- function(swd, policy, nthresh, value, Q) 
{
    coeff <- NULL
    for (i in 1:nthresh)
        coeff <- c(coeff, swd$swcoeff[[i]])

    len <- length(coeff)
    n <- length(swd$obs)
    
    if (policy == "universal") {
        sigma <- mad(swd$swcoeff[[1]])
        lam <- sqrt(2 * log(n)) * sigma    
    }        
    if (policy == "probability") {
        if (length(value) != 1) 
            stop("Length of value should be 1")
        lam <- quantile(abs(coeff), prob = value)
    }
    if (policy == "fdr") {    
        ### Adapt threshold.wd in wavethresh4_4.0-2
        sigma <- mad(swd$swcoeff[[1]])
        minit <- length(coeff)
        dinit <- coeff
        thinit <- qnorm(1 - Q/2) * sigma
        if (log(n, 2) > 12) 
            ninit <- 3
        else {
            if (log(n, 2) > 10) 
              ninit <- 2
            else ninit <- 1
        }
        for (k in seq(1, ninit)) {
            dinit1 <- dinit[abs(dinit) >= thinit]
            minit <- length(dinit1)
            if (minit == 0) 
                lam <- max(abs(coeff)) * 1.0001
            else {
                thinit <- qnorm(1 - (Q * minit)/(2 * n)) * sigma
                minit1 <- length(dinit1[abs(dinit1) >= thinit])
                if (minit1 == minit || minit1 == 0) 
                    break
                dinit <- dinit1
            }
        }
        if (sigma > 0) {
            m <- length(coeff)
            minit <- length(dinit)
            p <- (2 - 2 * pnorm(abs(dinit)/sigma))
            index <- order(p)
            j <- seq(1, minit)
            m0 <- j[p[index] <= (Q * j)/m]
            if (length(m0) != 0)
                m0 <- max(m0)
            if (length(m0) != 0 && m0 < minit) 
                lam <- abs(dinit[index[m0]])
            else {
                if (length(m0) == 0) 
                    lam <- max(abs(dinit)) * 1.0001
                else lam <- 0
            }
        }
        else lam <- 0
    }
    if (policy == "Lorentz")
        lam <- sqrt(sum(coeff * coeff))/sqrt(len)
    if (policy == "sure") {
        coeff <- coeff/mad(coeff)
        coeff <- abs(coeff)
        y <- sort(coeff)
        cy <- cumsum(y^2)
        cy <- c(0, cy[1:(len - 1)])
        ans <- len - 2 * 1:len + cy + len:1 * y^2
        lam <- y[ans==min(ans)] * mad(coeff)
    }              
    
    list(lam = lam)
}

mrsfield.comp.thresh.global <-
function (grid, coef, site, netlab, eta, lam, K, type) 
{
    dfield <- NULL
    alpha <- mrs.comp.thresh.global(coef, site, netlab, eta, K, lam, type)$talpha
    for (l in 1:ncol(alpha)) {
        temp <- msbf.comp(grid, site, alpha[, l], netlab, eta, l - 1)
        dfield <- cbind(dfield, temp$field)
    }
    dfield
}

mrs.comp.thresh.global <-
function (coef, site, netlab, eta, K, lam, type) 
{
    if (type != "soft" && type != "hard" && type != "Lorentz")
        stop("Thresholding type must be soft, hard or Lorentz.")
    
    J <- length(netlab)
    L <- length(eta)
    alp <- NULL
    blp <- NULL
    a1 <- mcov.comp(site, netlab, eta)
    beta1 <- coef
    mrs <- coefmatrix(beta1, a1, netlab, 1)
    norm <- sqrt(diag(mrs$norm))
           
    thr <- thresh.global(mrs$gamma1, lam/norm, type=type)
    talpha <- t(mrs$wcoef) %*% thr$tgamma
    beta2 <- mrs$beta2
    if (L == 2) 
        talpha <- as.matrix(talpha)
    else {
        if (K == (L - 1))
            for (l in 2:K) {
                temp <- coefmatrix(beta2, a1, netlab, l)
                norm <- sqrt(diag(temp$norm))
                thr <- thresh.global(temp$gamma1, lam/norm, type=type)
                talpha1 <- t(temp$wcoef) %*% thr$tgamma
                beta2 <- temp$beta2
                temp <- rep(NA, J)
                temp[netlab > (l - 1)] <- talpha1
                alp <- cbind(alp, temp)
            }
        else {
            for (l in 2:K) {
                temp <- coefmatrix(beta2, a1, netlab, l)
                norm <- sqrt(diag(temp$norm))
                thr <- thresh.global(temp$gamma1, lam/norm, type=type)
                talpha1 <- t(temp$wcoef) %*% thr$tgamma
                beta2 <- temp$beta2
                temp <- rep(NA, J)
                temp[netlab > (l - 1)] <- talpha1
                alp <- cbind(alp, temp)
            }
            for (l in (K + 1):(L - 1)) {
                tmp <- coefmatrix(beta2, a1, netlab, l)
                talpha1 <- tmp$alpha1                               
                beta2 <- tmp$beta2
                tmp <- rep(NA, J)
                tmp[netlab > (l - 1)] <- talpha1
                blp <- cbind(blp, tmp)
            }
        }
    }
    
    list(talpha = cbind(talpha, alp, blp))
}

#######################################
#thresh.total.gamma(x, lam) -> thresh.global(x, lam, type="soft")
#thresh.Lorentz.gamma(x, lam) -> thresh.global(x, lam, type="Lorentz")
#thresh.hard.gamma(x, lam) -> thresh.global(x, lam, type="hard")

thresh.global <-
function(x, lam, type)
{
    if (type != "soft" && type != "hard" && type != "Lorentz")
        stop("Thresholding type must be soft, hard or Lorentz.")
        
    if (type == "soft")
        x <- sign(x) * (abs(x) - lam) * (abs(x) > lam)
    if (type == "hard")
        x <- x * (abs(x) > lam)
    if (type == "Lorentz")
        x <- sign(x) * sqrt(ifelse(abs(x) > lam/norm, x^2, (lam/norm)^2) - (lam/norm)^2)             
    
    list(tgamma = x)
}

"swr" <- 
function (swd) 
{
    if (class(swd) != "swd") 
        stop("Input must be swd class.")   
    
    recon <- swd$global[[swd$nlevels - 1]]
    
    for(i in 1:(swd$nlevels - 1)) 
          recon <- recon + swd$detail[[i]] 
}
