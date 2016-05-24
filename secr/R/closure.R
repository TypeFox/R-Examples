############################################################################################
## package 'secr'
## closure.R
## Tests for closure of population sampled by capture-recapture
## At present does not fully allow for removals (losses on capture)
## last changed 2010 03 26
############################################################################################

#######################################
## Model Mt closed population estimate
Mtnegloglik <- function (N, ni, Mt1, nocc) {
    ## Otis et al 1978 p106-7
    nocc * N * log(N) -
        sum (ni * log(ni) + (N-ni) * log(N-ni)) -
        sum( log ( (N - Mt1 + 1) : N) )
}
N.Mt <- function (ni, Mt1) {
    optimize (f = Mtnegloglik, lower = Mt1, upper = 1e6, ni = ni,
        Mt1 = Mt1, nocc = length(ni))$minimum
}
#######################################

closure.test <- function (object, SB = FALSE, min.expected = 2) {

    if (!inherits(object, 'capthist'))
        stop ("requires 'capthist' object")

    if (inherits(object,'list')) {
        ## multiple sessions
        lapply(object, closure.test, SB = SB, min.expected = min.expected)
    }
    else {

        chisq <- function (x) {
            test <- suppressWarnings(chisq.test(matrix(x, ncol = 2), correct=F))
            if (any(test$expected<min.expected)) NA else test$statistic
        }

        ## single session
        nocc <- ncol(object)
        
        ## 2015-10-14
        if (nocc < 3) {
            cstat <- NA
            warning ("too few occasions for Otis et al. closure test")
        }
     
        cts <- as.matrix(summary(object)$counts)
        di <- apply(object, 2, function(x) sum(x<0))
        if (sum(di[-nocc])>0)
            warning ("estimates not adjusted for losses on capture")

        #######################################
        ## closure test of Otis et al. 1978
        ## robust to individual heterogeneity
        ## not robust to behavioural response etc.

        object <- abs(object)   ## ignore deads for now
        if (length(dim(object))>2)
            object <- apply(object,1:2, sum)
        object <- object>0
        mint <- apply(object, 1, function(x) min((1:nocc)[x]))
        maxt <- apply(object, 1, function(x) max((1:nocc)[x]))
        f <- apply(object, 1, sum)
        qtemp <- tapply (maxt-mint, f, sum)
        q <- numeric (nocc); names(q) <- 1:nocc
        q[names(qtemp)] <- qtemp
        f <- tabulate (f, nbins = nocc)
        k <- 2:(nocc-1)
        fik <- f[k]
        qik <- q[k]
        OK <- fik>0
        fik <- fik[OK]; qik <- qik[OK]; k <- k[OK]
        sumqi <- sum(qik/fik - (k-1) * (nocc+1) / (k+1))
        denom <- sum( 2*(nocc-k)*(k-1)*(nocc+1) / (k+2)/(k+1)^2/fik)
        cstat <- sumqi/sqrt(denom)
        if (sum(fik)<10) warning ("sum(fik) < 10")

        if (!SB) {
            c(statistic = cstat, p = pnorm(cstat))
        }
        else {
            #########################################
            ## closure test of Stanley & Burnham (1978) Env & Ecol Stat 6: 197-209
            ## see also Stanley & Richards (2005) WSB 33:782-785
            ## not robust to behavioural response etc.

            ni <- cts['n',1:nocc]
            ui <- cts['u',1:nocc]
            uistar <- rev(cumsum(rev(ui)))
            mi <- (ni-ui)
            ri <- (ni - tabulate(maxt, nbins=nocc) - di)
            zi <- numeric(nocc)
            for (i in 1:(nocc-2)) zi[i+1] <- zi[i] + ri[i] - mi[i+1] ## Seber 1982 p201
            Ti <- mi + zi
            Tistar <- c(0,Ti[-1]) + uistar

            Ri <- (ni - di)                       ## released
            Qi <- c(0,cumsum (Ri - mi)[-nocc])

            #######################################
            ## No-recruitment model
            stats <- cbind(mi, Ti-mi, ui, uistar-ui)
            XNRi <- apply(stats[-c(1,nocc),], 1, chisq)   ## i = 2..k-1
            XNR  <- sum(XNRi, na.rm=T)
            dfNR <- sum(!is.na(XNRi))
            #######################################
            ## No-mortality model
            stats <- cbind(zi, Qi-mi-zi,ri,Ri-ri)
            XNMi <- apply(stats[-c(1,nocc),], 1, chisq)   ## i = 2..k-1
            XNM  <- sum(XNMi, na.rm=T)
            dfNM <- sum(!is.na(XNMi))
            #######################################

            if (any(ni==0)) {
                warning ("one or more ni = 0")
                XtNRi <- rep(NA, nocc-1)
                XtNR <- NA
                dftNR <- NA
            }
            else {
                pi <- ni / N.Mt (ni, sum(ui))   ## simple form for no losses
                qi <- 1 - pi
                li <- c(1 - rev(cumprod(rev(qi[-1]))), NA)
                taui <- pi / (1 - rev(cumprod(rev(qi))))
                XtNRi <- (ri - Ri*li)^2 * (1/(Ri*li) + 1/(Ri*(1-li))) +
                    (ni - Tistar*taui)^2 * (1/(Tistar*taui) + 1/(Tistar*(1-taui)))
                XtNRi <- XtNRi[-nocc]
                XtNR <- sum(XtNRi, na.rm=T)
                dftNR <- nocc-2 - sum (is.na(XtNRi) | ((li*taui)[-nocc] == 0))
            }

            ## Following code is not reliable
            ## Note Stanley & Burnham 1999:203 calculate this by difference (see below)
            ## XtNMi <- ((ui - uistar * taui)^2 * (1/(uistar * taui) + 1/(uistar * (1-taui))))[1:(nocc-1)] +
            ##     ((mi - Qi * pi)^2 * (1/(Qi * pi) + 1/(Qi * (1-pi))))[2:nocc]
            ## XtNM <- sum(XtNMi, na.rm=T)
            ## dftNM <- sum (!is.na(XtNMi) & (taui[-nocc] != 0))

            XtNM <- XNR + XtNR - XNM
            dftNM <- dfNR + dftNR - dfNM

            #######################################
            ## Stanley & Burnham overall statistic
            Xc <- sum(XNRi, na.rm=T) + sum(XtNRi, na.rm=T)
            dfXc <- dfNR + dftNR
            #######################################

            list(
                 Otis   = data.frame(statistic = cstat, p = pnorm(cstat), row.names=''),
                 Xc     = data.frame(statistic = Xc,  df = dfXc, p = 1-pchisq(Xc, dfXc), row.names=''),
                 NRvsJS = data.frame(statistic = XNR, df = dfNR, p = 1-pchisq(XNR, dfNR), row.names=''),
                 NMvsJS = data.frame(statistic = XNM, df = dfNM, p = 1-pchisq(XNM, dfNM), row.names=''),
                 MtvsNR = data.frame(statistic = XtNR, df = dftNR, p = 1-pchisq(XtNR, dftNR), row.names=''),
                 MtvsNM = data.frame(statistic = XtNM, df = dftNM, p = 1-pchisq(XtNM, dftNM), row.names=''),
                 compNRvsJS = data.frame(Occasion=2:(nocc-1), Chisquare=XNRi, df=ifelse(is.na(XNRi),NA,1),
                     p = pchisq(XNRi, 1, lower.tail=FALSE), row.names=1:(nocc-2)),
                 compNMvsJS = data.frame(Occasion=2:(nocc-1), Chisquare=XNMi, df=ifelse(is.na(XNMi),NA,1),
                     p = pchisq(XNMi, 1, lower.tail=FALSE), row.names=1:(nocc-2))
            )
        }
    }
}
############################################################################################


