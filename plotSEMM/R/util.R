examplePlot <- function(){
    cex <- 2.5
    plot(c(0,100), c(0,80), axes=FALSE, frame.plot=FALSE, type='n', xlab='', ylab='', asp=1)
    
    #latent variables
    draw.ellipse(20, 55, 10,8)
    text(20, 55, expression(eta[1]), cex=cex)
    draw.ellipse(68, 55, 10,8)
    text(68, 55, expression(eta[2]), cex=cex)
    
    #beta[21] and arrow
    arrows(30,55, 58,55, length=.1)
    text(45, 59, expression(beta[21]["[k]"]), cex=cex)
    
    
    #psi[11] 
    draw.arc(7, 55, 3, deg1=320, deg2=40)
    arrows(8.34, 57.94, 9.25, 57.26,  length=.08)
    arrows(7.66, 52.03, 9.25, 53.16, length=.08)
    text(5, 62, expression(psi[11]["[k]"]), cex=cex)
    #curvedarrow(from=c(9.2,57.4), to=c(9.2,52.8), lwd=1, curve=.75, 
    #arr.type="simple", arr.length=.2, arr.pos=1)
    #curvedarrow(from=c(9.2,52.8), to=c(9.2,57.4), lwd=1, curve=-.75, 
    #arr.type="simple", arr.length=.2, arr.pos=1)
    
    
    
    #psi[22]
    draw.circle(89, 55, 5)
    text(89, 55, expression(zeta), cex=cex)
    draw.arc(97, 55, 3, deg1=220, deg2=500)
    text(96, 62, expression(psi[22]["[k]"]), cex=cex)
    arrows(95.43,57.71, 94.74, 57.03, length=.08)
    arrows(95.43,52.26, 94.52, 53.17, length=.08)
    arrows(84,55, 78,55, length=.1)
    
    
    #latent means
    polygon(c(15,25,20), c(70,70,80))
    text(20,74, '1', cex=cex)
    polygon(c(15,25,20)+48, c(70,70,80))
    text(68,74, '1', cex=cex)
    arrows(20,70,20,63, length=.1)
    text(14,66, expression(alpha[1]["[k]"]), cex=cex)
    arrows(68,70,68,63, length=.1)
    text(74,66, expression(alpha[2]["[k]"]), cex=cex)
    
    #measured variables
    polygon(c(9,19,19,9)+3, c(34,34,24,24))
    polygon(c(5,15,15,5)-9+3, c(34,34,24,24))
    polygon(c(5,15,15,5)+25+3, c(34,34,24,24))
    
    polygon(c(9,19,19,9)+48+3, c(34,34,24,24))
    polygon(c(5,15,15,5)-9+48+3, c(34,34,24,24))
    polygon(c(5,15,15,5)+25+48+3, c(34,34,24,24))
    
    text(14+3,30-1,expression(italic(y)[2]), cex=cex)
    text(2+3,30-1,expression(italic(y)[1]), cex=cex)
    text(24.5+3,30-1, expression(...), cex=3)
    text(35+3,30-1,expression(italic(y)[italic(i)]), cex=cex)
    
    text(14+48+3,30-1,expression(italic(y)[italic(i)+3]), cex=cex)
    text(2+48+3,30-1,expression(italic(y)[italic(i)+1]), cex=cex)
    text(24.5+48+3,30-1, expression(...), cex=3)
    text(35+48+3,30-1,expression(italic(y)[italic(n)]), cex=cex)
    
    arrows(20,47, 1+3, 35.2-1, length=.1)
    arrows(20,47, 14+3, 35-1, length=.1)
    arrows(20,47, 35+3, 35.2-1, length=.1)
    
    arrows(20+48,47, 1+48+3, 35.2-1, length=.1)
    arrows(20+48,47, 14+48+3, 35-1, length=.1)
    arrows(20+48,47, 35+48+3, 35.2-1, length=.1)
    
    #residual of measured variables
    draw.circle(1+3, 10+1, 5)
    draw.circle(14+3, 10+1, 5)
    draw.circle(35+3, 10+1, 5)
    
    draw.circle(1+48+3, 10+1, 5)
    draw.circle(14+48+3, 10+1, 5)
    draw.circle(35+48+3, 10+1, 5)
    
    text(1+3, 10+1, expression(epsilon[1]), cex=cex)
    text(14+3, 10+1, expression(epsilon[2]), cex=cex)
    text(35+3, 10+1, expression(epsilon[italic(i)]), cex=cex)
    
    text(1+48+3, 10+1, expression(epsilon[1+italic(i)]), cex=cex)
    text(14+48+3, 10+1, expression(epsilon[2+italic(i)]), cex=cex)
    text(35+48+3, 10+1, expression(epsilon[italic(n)]), cex=cex)
    
    arrows(1+3, 24, 1+3,16, length=.1)
    arrows(14+3, 24, 14+3,16, length=.1)
    arrows(35+3, 24, 35+3,16, length=.1)
    
    arrows(1+48+3, 24, 1+48+3,16, length=.1)
    arrows(14+48+3, 24, 14+48+3,16, length=.1)
    arrows(35+48+3, 24, 35+48+3,16, length=.1)
    
    draw.arc(1+3, 2+1, 3, deg1=135, deg2=405)
    draw.arc(14+3, 2+1, 3, deg1=135, deg2=405)
    draw.arc(35+3, 2+1, 3, deg1=135, deg2=405)
    
    draw.arc(1+48+3, 2+1, 3, deg1=135, deg2=405)
    draw.arc(14+48+3, 2+1, 3, deg1=135, deg2=405)
    draw.arc(35+48+3, 2+1, 3, deg1=135, deg2=405)
    
    arrows(1.07, 3.14+1, 1.75, 4.05+1, length=.08)
    arrows(6.75, 2.91+1, 6.07, 4.05+1, length=.08)
    arrows(1.07+13, 3.14+1, 1.75+13, 4.05+1, length=.08)
    arrows(6.75+13, 2.91+1, 6.07+13, 4.05+1, length=.08)
    arrows(1.07+34, 3.14+1, 1.75+34, 4.05+1, length=.08)
    arrows(6.75+34, 2.91+1, 6.07+34, 4.05+1, length=.08)
    
    arrows(1.07+48, 3.14+1, 1.75+48, 4.05+1, length=.08)
    arrows(6.75+48, 2.91+1, 6.07+48, 4.05+1, length=.08)
    arrows(1.07+13+48, 3.14+1, 1.75+13+48, 4.05+1, length=.08)
    arrows(6.75+13+48, 2.91+1, 6.07+13+48, 4.05+1, length=.08)
    arrows(1.07+34+48, 3.14+1, 1.75+34+48, 4.05+1, length=.08)
    arrows(6.75+34+48, 2.91+1, 6.07+34+48, 4.05+1, length=.08)
}

# Bootstrap ellipse function. Be very careful not to allocate huge objects to avoid memory issues
bs.CE <- function(read, x, alpha, boot = FALSE){
    ACOV <- read$tech3$paramCov.savedata
    nclass <- length(read$tech1[[1L]]) - 1L
    omitpars <- c()
    for(i in 1L:nclass){
        tmp <- read$tech1[[1L]][[i]]
        omitpars <- c(omitpars, unique(na.omit(as.numeric(tmp$theta))), 
            unique(na.omit(as.numeric(tmp$lambda))), unique(na.omit(as.numeric(tmp$nu))))
    }
    nomitpars <- sum(unique(omitpars) != 0) + nclass # also to remove latent residual vars count
    draws <- if(boot) 1000 else pickNdraws(ncol(ACOV) - nomitpars, alpha=alpha)
    iter <- 0
    is.variance <- logical(ncol(ACOV))
    cholL <- chol(ACOV)
    points <- 250
    bs.yall <- NULL
    
    #construct spcification model with real parameters as list
    spec <- read$tech1$parameterSpecification
    samplepars <- spec
    pars <- read$parameters$unstandardized
    opars <- colnames(spec[[1L]]$theta)
    lpars <- unique(pars[, 'param'][!(pars[, 'param'] %in% opars)])
    exo <- pars[pars[, 'param'] %in% lpars & pars[, 'paramHeader'] == 'Means' & 
                     pars[, 'LatentClass'] == 1, 'param']
    endo <- colnames(spec[[1L]]$lambda)
    endo <- endo[endo != exo]
    for(i in 1L:nclass){
        cls <- pars[pars[,'LatentClass'] == i, ]
        samplepars[[i]]$nu <- cls[cls[, 'paramHeader'] == 'Intercepts' & 
                                      cls[, 'param'] %in% opars, 'est']
        tmp <- cls[grepl('\\.BY', cls[, 'paramHeader']), ]
        endoind <- grepl(endo, tmp[,'paramHeader'])
        lambdas <- samplepars[[i]]$lambda
        tmp <- cls[grepl('\\.BY', cls[, 'paramHeader']) & 
                       cls[, 'param'] %in% opars, 'est']
        lambdas[endoind, colnames(lambdas) == endo] <- tmp[endoind]
        lambdas[!endoind, colnames(lambdas) != endo] <- tmp[!endoind]
        samplepars[[i]]$lambda <- lambdas
        diag(samplepars[[i]]$theta) <- cls[cls[, 'paramHeader'] == 'Residual.Variances' & 
                                         cls[, 'param'] %in% opars, 'est']
        samplepars[[i]]$alpha <- c(cls[cls[, 'paramHeader'] == 'Means'
                                       & cls[, 'param'] == exo, 'est'],
                                   cls[cls[, 'paramHeader'] == 'Intercepts'
                                       & cls[, 'param'] == endo, 'est'])
        samplepars[[i]]$beta[2,1] <- cls[grepl('\\.ON', cls[, 'paramHeader']), 'est']
        diag(samplepars[[i]]$psi) <- c(cls[cls[, 'paramHeader'] == 'Variances'
                                         & cls[, 'param'] == exo, 'est'],
                                       cls[cls[, 'paramHeader'] == 'Residual.Variances'
                                           & cls[, 'param'] == endo, 'est'])
    }
    samplepars[[nclass + 1L]]$alpha.c <- c(pars[!(pars[,'param'] %in% 
                                                      c(opars, exo, endo)), 'est'], 0)
    
    while(TRUE){
        
        jitter <- rnorm(ncol(ACOV)) %*% cholL
        
        #load the jittered pars
        tmpmod <- loadMplusJitter(samplepars, spec, jitter)
        redraw <- FALSE
        for(i in 1L:nclass){
            tmp <- tmpmod[[i]]$psi
            tmp[is.na(tmp)] <- 0
            if(is(try(chol(tmp), silent=TRUE), 'try-error'))
                redraw <- TRUE
        }
        if(redraw) next
        
        #perform computations
        bs.pi <- exp(tmpmod[[nclass + 1]][[1]]) / sum(exp(tmpmod[[nclass + 1]][[1]]))
        bs.phi <- matrix(NA, length(x), nclass)
        for(i in 1L:nclass)
            bs.phi[,i] <- dnorm(x, tmpmod[[i]]$alpha[1L], sqrt(tmpmod[[i]]$psi[1L,1L]))
        tmp <- t(bs.pi * t(bs.phi))
        bs.pi_ <- tmp / rowSums(tmp)
        bs.y <- numeric(length(x))
        for(i in 1L:nclass)
            bs.y <- bs.y + bs.pi_[,i] * (tmpmod[[i]]$alpha[2L] + 
                                            tmpmod[[i]]$beta[2L,1L] * x)
        if(iter == 0L){
            lb.CE <- ub.CE <- bs.y
            if(boot) bs.yall <- matrix(bs.y, draws, length(bs.y), byrow=TRUE)
        } else {
            lb.CE[lb.CE > bs.y] <- bs.y[lb.CE > bs.y]
            ub.CE[ub.CE < bs.y] <- bs.y[ub.CE < bs.y]
            if(boot) bs.yall[iter + 1L, ] <- bs.y
        }
        
        #increment and break
        iter <- iter + 1
        if(iter == draws) break
    }
    
    ret <- list(lb = lb.CE, ub = ub.CE, bs.yall = bs.yall)
    return(ret)
}

pickNdraws <- function(npars, alpha){ 
    #95% only for now
    if(alpha == .025){
        draws <- switch(as.character(npars), 
                        "1" = 25,
                        "2" = 85,
                        "3" = 232,
                        "4" = 577,
                        "5" = 1353,
                        "6" = 3047,
                        "7" = 6665,
                        "8" = 14255,
                        "9" = 29945,
                        "10" = 61976,
                        "11" = 126671,
                        "12" = 256130,
                        "13" = 513079,
                        "14" = 1019382,
                        "15" = 2010562,
                        "16" = 3939636,
                        "17" = 7674248,
                        "18" = 14868770,
                        "19" = 28667700,
                        "20" = 55024760,
                        "21" = 105117600,
                        "22" = 200272600,
                        "23" = 379989400,
                        "24" = 718584600,
                        "25" = 1354673000,
                        "26" = 2546395000,
                        "27" = 4773392000,
                        "28" = 8924985000,
                        "29" = 16646750000,
                        "default" = stop('bootstrap CIs only supported in models with 8 to 29 parameters'))
    } else {
        draws <- switch(as.character(npars), 
                        "1" = 12,
                        "2" = 39,
                        "3" = 98,
                        "4" = 228,
                        "5" = 504,
                        "6" = 1077,
                        "7" = 2244,   
                        "8" = 4586,
                        "9" = 9232,
                        "10" = 18352,
                        "11" = 36095,
                        "12" = 70351,
                        "13" = 136042,
                        "14" = 261257,
                        "15" = 498651,
                        "16" = 946541,
                        "17" = 1787852,
                        "18" = 3361809,
                        "19" = 6295567,
                        "20" = 11745310,
                        "21" = 21836900,
                        "22" = 40469420,
                        "23" = 74777710,
                        "24" = 137789200,
                        "25" = 253241500,
                        "26" = 464303700,
                        "27" = 849338400,
                        "28" = 1550345000,
                        "29" = 2824215000,
                        "default" = stop('bootstrap CIs only supported in models with 8 to 29 parameters'))
    }
        
    draws
}
    
loadMplusJitter <- function(samplepars, spec, jitter){
    for(g in 1:length(samplepars)){
        for(i in 1:length(samplepars[[g]])){
            tmp <- samplepars[[g]][[i]] 
            pick <- spec[[g]][[i]] != 0 & !is.na(spec[[g]][[i]])
            tmp[pick] <- tmp[pick] + jitter[spec[[g]][[i]][pick]]
            samplepars[[g]][[i]] <- tmp
        }
    }
    samplepars
}