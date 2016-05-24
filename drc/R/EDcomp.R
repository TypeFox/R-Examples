"EDcomp" <-
function(object, percVec, percMat = NULL, compMatch = NULL, od = FALSE, vcov. = vcov, reverse = FALSE, 
interval = c("none", "delta", "fieller", "fls"), level = ifelse(!(interval == "none"), 0.95, NULL), 
reference = c("control", "upper"), type = c("relative", "absolute"),
display = TRUE, pool = TRUE, logBase = NULL, ...)
{
     ## Matching the argument 'method'
     interval <- match.arg(interval)
     reference <- match.arg(reference)
     type <- match.arg(type)     

    if ( (is.null(logBase)) && (interval == "fls") )
    {
        stop("Argument 'logBase' not specified for interval = 'fls'")
    }
#    sifct <- createsifct(object$"fct"$"edfct", logBase, identical(interval, "fls"), object$"indexMat", length(coef(object)))

    ## Checking contain of percVec vector ... should be numbers between 0 and 100
    if ( (type == "relative") && any(percVec<=0 | percVec>=100) ) 
    {
        stop("Percentages outside the interval [0, 100] not allowed")
    }

    if (missing(compMatch)) {matchNames <- FALSE} else {matchNames <- TRUE}

    lenPV <- length(percVec)

    ## Retrieving relevant quantities
    indexMat <- object$"indexMat"
#    lenEB <- ncol(indexMat)    
    parmMat <- object$"parmMat"
#    strParm <- colnames(parmMat)    
#    varMat <- vcov(object, od = od, pool = pool)
##    compNames <- as.character(strParm)  # converting a factor


    curveNames <- colnames(object$"parmMat")
    options(warn = -1)  # switching off warnings caused by coercion in the if statement
    if (any(is.na(as.numeric(curveNames))))
    {
        curveOrder <- order(curveNames)
    } else { # if names are numbers then skip re-ordering
        curveOrder <- 1:length(curveNames)
    }
    options(warn = 0)  # normalizing behaviour of warnings
    
    strParm0 <- curveNames[curveOrder]
    indexMat <- indexMat[, curveOrder, drop = FALSE]
    lenEB <- ncol(indexMat) 
    sifct <- createsifct(object$"fct"$"edfct", logBase, identical(interval, "fls"), indexMat, length(coef(object)))    
    
    parmMat <- parmMat[, curveOrder, drop = FALSE]
    
    strParm <- strParm0
#    varMat <- vcov(object, od = od, pool = pool)
    varMat <- vcov.(object)    

    ## Calculating SI values
    numComp <- (lenPV*(lenPV-1)/2)*(lenEB * (lenEB - 1) / 2)
    
#    if (!identical(interval, "none"))
#    {
#        siMat <- matrix(0, numComp, 3)
#        cNames <- c("Estimate", "Lower", "Upper")
#    
#    } else {
#        siMat <- matrix(0, numComp, 4)
#        cNames <- c("Estimate", "Std. Error", "t-value", "p-value")
#    }
    matchVec <- rep(TRUE, numComp)
    rNames <- rep("", numComp)
    oriMat <- matrix(0, numComp, 2)    
    degfree <- df.residual(object)  
    rowIndex <- 1
    
#    require(gtools, quietly = TRUE)
    
    pairsMat <- combinations(lenEB, 2)  # canonical "2" as pairs are considered
    if (is.null(percMat))
    {
        percMat <- combinations(lenPV, 2)  # canonical "2" as pairs are considered
    }
    if (reverse)
    {
        pairsMat <- pairsMat[, 2:1, drop = FALSE]
        percMat <- percMat[, 2:1, drop = FALSE]
#        strParm <- rev(strParm)
    }
#    print(strParm)
#    print(pairsMat)

    appFct1 <- function(percVal)
    {
        apply(pairsMat, 1, siInner, pVec = percVec[percVal], compMatch = compMatch, object = object, indexMat = indexMat, parmMat = parmMat, 
        varMat = varMat, level = level, reference = reference, type = type, sifct = sifct, interval = interval, degfree = degfree, logBase = logBase)
    }
    SImat0 <- matrix(apply(percMat, 1, appFct1), nrow = nrow(pairsMat) * nrow(percMat), byrow = TRUE)
    SImat <- SImat0[, 1:4, drop = FALSE]
#    print(SImat)
    dSImat <- SImat0[, 5:ncol(SImat0), drop = FALSE] 
    
#    matchVec[rowIndex] <- (is.null(compMatch) || all(c(strParm[j], strParm[k]) %in% compMatch))
    
#    strParm0 <- sort(colnames(object$"parmMat"))
    appFct2 <- function(percVal)
    {
        apply(pairsMat, 1, 
        function(indPair, percVal) 
        {
#            paste(strParm0[indPair[1]], "/", strParm0[indPair[2]], ":", percVec[percVal[1]], "/", percVec[percVal[2]], sep = "")
            paste(strParm[indPair[1]], "/", strParm[indPair[2]], ":", percVec[percVal[1]], "/", percVec[percVal[2]], sep = "")
        }, percVal = percVal)
    }    
    rownames(SImat) <- apply(percMat, 1, appFct2) 

    appFct3 <- function(percVal)
    {
        apply(pairsMat, 1, 
        function(indPair, percVal) 
        {
#            (is.null(compMatch) || all(c(strParm0[indPair[1]], strParm0[indPair[2]]) %in% compMatch))
            (is.null(compMatch) || all(c(strParm[indPair[1]], strParm[indPair[2]]) %in% compMatch))
        })
    }
    SImat <- SImat[as.vector(apply(percMat, 1, appFct3)), , drop = FALSE]
    
    if (!identical(interval, "none"))
    {
 #       siMat <- matrix(0, numComp, 3)
        SImat <- SImat[, -4, drop = FALSE]
        cNames <- c("Estimate", "Lower", "Upper")
    
    } else {
#        siMat <- matrix(0, numComp, 4)
        cNames <- c("Estimate", "Std. Error", "t-value", "p-value")
    }    
    colnames(SImat) <- cNames 
#    print(SImat)
    
    ciLabel <- switch(interval,
    "delta" = "Delta method",
    "tfls" = "To and from log scale",
    "fls" = "From log scale",
    "fieller" = "Fieller")
    
#    resPrint(SImat, "Estimated ratios of effect doses\n", interval, ciLabel, display = display) 
    resPrint(SImat, "Estimated ratios of effect doses", interval, ciLabel, display = display)
    
#    invisible(SImat) 
    ## require(multcomp, quietly = TRUE)
    invisible(list(SIdisplay = SImat, SImultcomp = list(SImat[, 1], dSImat %*% varMat %*% t(dSImat))))  
#    invisible(list(SImat, SImultcomp = list(EDest = EDmat[, 1], EDvcov = dEDmat %*% vcMat %*% t(dEDmat))))      

#    
#    if (FALSE)
#    {
#    for (i in 1:lenPV)
#    {
#        for (ii in 1:lenPV)
#        {
#            if (i>=ii) {next}
#            pVec <- percVec[c(i, ii)]
#
#            for (j in 1:lenEB)
#            {
#                for (k in 1:lenEB)
#                {
#                    if (j>=k) {next}
#                    matchVec[rowIndex] <- (is.null(compMatch) || all(c(strParm[j], strParm[k]) %in% compMatch))  
# 
#                    jInd <- j
#                    kInd <- k
#                    if (reverse) 
#                    {
#                        jInd <- k; kInd <- j; pVec <- pVec[c(2, 1)]
#                    }
#                   
#                    parmInd1 <- indexMat[, jInd]
#                    parmInd2 <- indexMat[, kInd]
#                    
#                    splInd <- splitInd(parmInd1, parmInd2)
#                    
#                    parmChosen1 <- parmMat[, jInd]
#                    parmChosen2 <- parmMat[, kInd]
#
#                    SIeval <- 
#                    sifct(parmChosen1, parmChosen2, pVec, 
#                    splInd[[1]][, 1], splInd[[2]][, 1], splInd[[3]][, 1], splInd[[3]][, 2], reference, type, jInd, kInd, ...)
#
#                    indInOrder <- c(splInd[[1]][, 2], splInd[[2]][, 2], splInd[[3]][, 3])
#                                       
#                    SIval <- SIeval$"val"  # SIeval[[1]]
#                    dSIval <- SIeval$"der"  # SIeval[[2]]
##                    print(dSIval)
# 
#                    oriMat[rowIndex, 1] <- SIval
##                    oriMat[rowIndex, 2] <- sqrt(dSIval %*% varMat[indInOrder, indInOrder] %*% dSIval)  # sqrt(dSIval%*%varCov%*%dSIval)
#                    oriMat[rowIndex, 2] <- sqrt(t(dSIval) %*% varMat %*% dSIval)
#
#                    siMat[rowIndex, 1] <- SIval
#                    rNames[rowIndex] <- paste(strParm[jInd], "/", strParm[kInd], ":", pVec[1], "/", pVec[2], sep="")
#
#                    ## Using t-distribution for continuous data
#                    ##  only under the normality assumption
#                    if (identical(object$"type", "continuous"))
#                    {
#                        qFct <- function(x) {qt(x, degfree)}
#                        pFct <- function(x) {pt(x, degfree)}
#                    } else {
#                        qFct <- qnorm
#                        pFct <- pnorm
#                    }
#
#                    if (identical(interval, "none"))
#                    {
#                        siMat[rowIndex, 2] <- oriMat[rowIndex, 2]  # sqrt(dSIval%*%varCov%*%dSIval)
#
#                        ## Testing SI equal to 1
#                        tempStat <- (siMat[rowIndex, 1] - 1)/siMat[rowIndex, 2]
#                        siMat[rowIndex, 3] <- tempStat
#                        siMat[rowIndex, 4] <- pFct(-abs(tempStat)) + (1 - pFct(abs(tempStat)))
#                    }
#                    if ( (identical(interval, "delta")) || (identical(interval, "fls")) )
#                    {
#                        stErr <- oriMat[rowIndex, 2]  # sqrt(derEval%*%varCov%*%derEval)
#                        tquan <- qFct(1 - (1 - level)/2)
#                        
#                        siMat[rowIndex, 2] <- siMat[rowIndex, 1] - tquan * stErr
#                        siMat[rowIndex, 3] <- siMat[rowIndex, 1] + tquan * stErr
#                        ciLabel <- "Delta method"
#                    }
#                    if (identical(interval, "tfls"))
#                    {
#                        lsVal <- log(oriMat[rowIndex, 1])
#                        lsdVal <- oriMat[rowIndex, 2]/oriMat[rowIndex, 1]
#                        tquan <- qFct(1 - (1 - level)/2)
#                        
#                        siMat[rowIndex, 2] <- exp(lsVal - tquan * lsdVal)
#                        siMat[rowIndex, 3] <- exp(lsVal + tquan * lsdVal)
#                        ciLabel <- "To and from log scale"
#                    }
#                    if ((!is.null(logBase)) && (identical(interval, "fls")))
#                    {
#                        siMat[rowIndex, 1] <- logBase^(siMat[rowIndex, 1])
#                        siMat[rowIndex, 2] <- logBase^(siMat[rowIndex, 2])
#                        siMat[rowIndex, 3] <- logBase^(siMat[rowIndex, 3])
#                        ciLabel <- "From log scale"
#                    }
#                    if (identical(interval, "fieller"))  # using t-distribution
#                    {
#                        vcMat <- matrix(NA, 2, 2)
#                        vcMat[1, 1] <- SIeval$"der1"%*%varMat[parmInd1, parmInd1]%*%SIeval$"der1"
#                        vcMat[2, 2] <- SIeval$"der2"%*%varMat[parmInd2, parmInd2]%*%SIeval$"der2"
#                        vcMat[1, 2] <- SIeval$"der1"%*%varMat[parmInd1, parmInd2]%*%SIeval$"der2"
#                        vcMat[2, 1] <- vcMat[1, 2]
#                        muVec <- c(SIeval$"valnum", SIeval$"valden")
#                        
#                        siMat[rowIndex, 2:3] <- fieller(muVec, degfree, vcMat, level = level)  
#                        ciLabel <- "Fieller"
#                    }
#
#
#                    rowIndex <- rowIndex+1
#                }
#            }
#        }
#    }
#    dimnames(siMat) <- list(rNames, cNames)
#    siMat <- siMat[matchVec, , drop = FALSE]
#
#    resPrint(siMat, "Estimated ratios of effect doses\n", interval, ciLabel, display = display)
    
#    if (display)
#    { 
#        cat("\n")
#        cat("Estimated ratios of effect doses\n")
#        if (!(ci == "none")) 
#        {
#            ciText <- paste("(", ciLabel, "-based confidence interval(s))\n", sep = "")
#            cat(ciText)
#        } 
#        cat("\n")
#        printCoefmat(siMat)
#    }
#    invisible(siMat)   
}

SI <- EDcomp


"fieller" <-
function(mu, df, vcMat, level = 0.95, finney = FALSE, resVar)
{
    tper <- qt(1-0.5*(1-level), df)^2 
    
    if (!finney)
    {
        ## Based on the entry on Fieller's theorem 
        ##  in Encyclopedia of Statistical Sciences Vol. 3 (1983), p. 86 
        ##  essentially same formula as in Finney (see below)

        mup <- prod(mu)
    
        fVec0 <- mup - tper*vcMat[1,2]
        y2 <- mu[2]^2    
        fVec <- (fVec0)^2 - (mu[1]^2 - tper*vcMat[1,1])*(y2 - tper*vcMat[2,2])
    
        denom <- y2 - tper*vcMat[2,2]
        lowerL <- (fVec0 - sqrt(fVec))/denom
        upperL <- (fVec0 + sqrt(fVec))/denom
        
    } else {
    
        ## Using the formula 
        ##  in Finney: Statistical Method in Biological Assay p. 81 (3rd edition, 1978)
        ##  OOPS: uses the estimated residual variance
        fac <- sqrt(tper)*sqrt(resVar)/mu[2]
        g <- tper*vcMat[2,2]/(mu[2]^2)
        if (g >= 1) {stop("Fieller's theorem not useful!")} 
        ratio <- mu[1]/mu[2]
    
        v11 <- vcMat[1,1]/(resVar)
        v12 <- vcMat[1,2]/(resVar)
        v22 <- vcMat[2,2]/(resVar)
        innerBr <- g*(v11 - (v12^2)/v22)
        inBr <- v11 - 2*ratio*v12 + (ratio^2)*v22 - innerBr

        firstTerm <- ratio - g*vcMat[1,2]/vcMat[2,2]
        secondTerm <- fac*sqrt(inBr)
        denom <- 1 - g
        lowerL <- (firstTerm - secondTerm)/denom
        upperL <- (firstTerm + secondTerm)/denom
    }
    return(c(lowerL, upperL))
}

"splitInd"  <- function(ind1, ind2)
{
    matchVec1 <- ind1%in%ind2
    matchVec2 <- ind2%in%ind1
    
#    inCommon <- list(pos1 = (1:length(ind1))[matchVec1], pos2 = (1:length(ind2))[matchVec2], val = ind1[matchVec1])
#
    lmv1 <- sum(matchVec1)
    if (lmv1 > 0.01)
    {
        inCommon <- matrix( c( (1:length(ind1))[matchVec1], (1:length(ind2))[matchVec2], ind1[matchVec1]), lmv1, 3)
    } else {
        inCommon <- NULL
    }
    
#    only1 <- list(pos = (1:length(ind1))[!matchVec1], val = ind1[!matchVec1])
    only1 <- matrix( c( (1:length(ind1))[!matchVec1], ind1[!matchVec1] ), sum(!matchVec1), 2)
    
#    only2 <- list(pos = (1:length(ind2))[!matchVec2], val = ind2[!matchVec2])
    only2 <- matrix( c( (1:length(ind2))[!matchVec2], ind2[!matchVec2] ), sum(!matchVec2), 2)

    return(list(only1, only2, inCommon))
}

createsifct <- function(edfct, logBase = NULL, fls = FALSE, indexMat, lenCoef)
{
    if (is.null(edfct)) 
    {
        stop("SI values cannot be calculated")
    } else {
        
        if (!fls)
        {
            if (is.null(logBase))  # this clause has been updated October 12 2010
            {
                "sifct" <- function(parm1, parm2, pair, jInd, kInd, reference, type, ...)
                ## ind1, ind2, cmonInd1, cmonInd2 not used
                {
#                    print(parm1)
#                    print(parm2)
#                    print(indexMat)

                    ED1 <- edfct(parm1, pair[1], reference = reference, type = type, ...)
                    ED1v <- ED1[[1]]
                    ED1d <- rep(0, lenCoef)
#                    print(indexMat[, jInd])
#                    print(ED1[[2]])
                    ED1d[indexMat[, jInd]] <- ED1[[2]]        
#                    print(ED1v)
#                    print(ED1d)
        
                    ED2 <- edfct(parm2, pair[2], reference = reference, type = type, ...)
                    ED2v <- ED2[[1]]
                    ED2d <- rep(0, lenCoef)
                    ED2d[indexMat[, kInd]] <- ED2[[2]]        
#                    print(ED2v)
#                    print(ED2d)

        
                    SIpair <- ED1v / ED2v  # calculating the SI value        
                    SIder <- (ED1d - SIpair * ED2d) / ED2v  # calculating the derivative of SI

#                    SIder1 <- ED1d/ED2v
#                    SIder2 <- (-ED2d/ED2v)*SIpair
#                    SIder12 <- commonParm(SIder1, SIder2, cmonInd1, cmonInd2)
#                    SIder12 <- ED1d/ED2v - (ED2d/ED2v)*SIpair
                     
#                    return(list(val = SIpair, der = c(SIder1[ind1], SIder2[ind2], SIder12),                    
                    return(list(val = SIpair, der = SIder,
                    der1 = ED1d, der2 = ED2d, valnum = ED1v, valden = ED2v))
                }
            } else {
        
                "sifct" <- function(parm1, parm2, pair, jInd, kInd, reference, type, ...)
                {
                    ED1 <- edfct(parm1, pair[1], reference = reference, type = type, ...)
                    ED1v <- ED1[[1]]
                    ED1d <- rep(0, lenCoef)
                    ED1d[indexMat[, jInd]] <- ED1[[2]]        
        
                    ED2 <- edfct(parm2, pair[2], reference = reference, type = type, ...)
                    ED2v <- ED2[[1]]
                    ED2d <- rep(0, lenCoef)
                    ED2d[indexMat[, kInd]] <- ED2[[2]]        

#                    ED1 <- edfct(parm1, pair[1], reference = reference, type = type, ...)
#                    ED1v <- ED1[[1]]
#                    ED1d <- ED1[[2]]       
#                    ED2 <- edfct(parm2, pair[2], reference = reference, type = type, ...)
#                    ED2v <- ED2[[1]]
#                    ED2d <- ED2[[2]]
        
                    SIpair <- logBase^(ED1v - ED2v)  # calculating the SI value
                    SIder <- SIpair * log(logBase) * (ED1d - ED2d)
        
#                    SIder1 <- SIpair*log(logBase)*ED1d
#                    SIder2 <- SIpair*log(logBase)*(-ED2d)
#                    SIder12 <- commonParm(SIder1, SIder2, cmonInd1, cmonInd2)

#                    return(list(val = SIpair, der = c(SIder1[ind1], SIder2[ind2], SIder12),
                    return(list(val = SIpair, der = SIder,
                    der1 = (log(logBase)*logBase^ED1v)*ED1d, der2 = (log(logBase)*logBase^ED2v)*ED2d, 
                    valnum = logBase^ED1v, valden = logBase^ED2v))
                }        
            }
        } else {
            
            "sifct" <- function(parm1, parm2, pair, jInd, kInd, reference, type, ...)
            {
                ED1 <- edfct(parm1, pair[1], reference = reference, type = type, ...)
                ED1v <- ED1[[1]]
                ED1d <- rep(0, lenCoef)
                ED1d[indexMat[, jInd]] <- ED1[[2]]        
        
                ED2 <- edfct(parm2, pair[2], reference = reference, type = type, ...)
                ED2v <- ED2[[1]]
                ED2d <- rep(0, lenCoef)
                ED2d[indexMat[, kInd]] <- ED2[[2]]        

#                ED1 <- edfct(parm1, pair[1], reference = reference, type = type, ...)
#                ED1v <- ED1[[1]]
#                ED1d <- ED1[[2]]       
#                ED2 <- edfct(parm2, pair[2], reference = reference, type = type, ...)
#                ED2v <- ED2[[1]]
#                ED2d <- ED2[[2]]        
            
                SIpair <- ED1v - ED2v  # calculating the log SI value
                SIder <- ED1d - ED2d
        
#                SIder1 <- ED1d
#                SIder2 <- -ED2d
#                SIder12 <- commonParm(SIder1, SIder2, cmonInd1, cmonInd2)

#                return(list(val = SIpair, der = c(SIder1[ind1], SIder2[ind2], SIder12),
                return(list(val = SIpair, der = SIder,
                der1 = ED1d, der2 = ED2d, valnum = ED1v, valden = ED2v))
            }
        }        
        return(sifct)
    }
}

#commonParm <- function(SIder1, SIder2, cmonInd1, cmonInd2)
#{
#    lind1 <- length(cmonInd1)
#    retVec <- rep(NA, lind1)
#    for (i in 1:lind1)
#    {
#        retVec[i] <- SIder1[cmonInd1[i]] + SIder2[cmonInd2[i]]
#    
#    }
#    return(retVec)
#}
#

