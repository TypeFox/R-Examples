############################################################################################
## package 'secr'
## dbar.R, ARL.R, MMDM.R, RPSV.R
## all in homerange.R 2014-09-01
## 2014-09-01 modified for userdist
## 2014-09-10 does NOT work when userdist fn requires mask covariates..
############################################################################################

dbar <- function (capthist, userdist = NULL, mask = NULL) {
    if (inherits (capthist, 'list')) {
        lapply(capthist, dbar, userdist, mask)   ## recursive
    }
    else {
        ## 2014-09-01
        dbarx    <- function (x) {
            x <- abs(unlist(x))
            ## sqrt(diff(traps$x[x])^2 + diff(traps$y[x])^2)
            distmat[cbind(x[-length(x)], x[-1])]  ## vector
        }
        dbarxy    <- function (xy) {
            sqrt(diff(xy$x)^2 + diff(xy$y)^2)
        }

        traps <- traps(capthist)
        ## 2014-09-01
        ## NOT USING PARAMETERS noneuc ETC
        distmat <- valid.userdist(userdist, detector(traps), traps, traps, mask )
        if (!(detector(traps) %in% .localstuff$individualdetectors))
            stop ("require individual detector type for dbar")
        
        if (detector(traps) %in% .localstuff$polydetectors) {
            if (is.null(xy(capthist)))
                temp <- NA
            else {
                lxy <- split (xy(capthist), animalID(capthist))
                d <- try(lapply(lxy,dbarxy), silent = TRUE)
                if (inherits(d, 'try-error'))
                    d <- NA
                mean(unlist(d), na.rm=T)
            }
        }
        else {
            w <- split(trap(capthist, names=F), animalID(capthist))
            d <- try(unlist(lapply(w,dbarx)), silent = TRUE)
            if (inherits(d, 'try-error'))
                d <- NA
            mean(d, na.rm=T)
        }
    }
}
############################################################################################

moves <- function (capthist, userdist = NULL, mask = NULL) {
    if (inherits (capthist, 'list')) {
        lapply(capthist, moves, userdist, mask)   ## recursive
    }
    else {
        movex    <- function (x) {
            x <- abs(unlist(x))
            distmat[cbind(x[-length(x)], x[-1])]  ## vector
            ## sqrt(diff(traps$x[x])^2 + diff(traps$y[x])^2)
        }
        movexy    <- function (xy) {
            sqrt(diff(xy$x)^2 + diff(xy$y)^2)
        }
        traps <- traps(capthist)
        distmat <- valid.userdist(userdist, detector(traps), traps, traps, mask)
        if (!(detector(traps) %in% .localstuff$individualdetectors))
            stop ("require individual detector type for moves")
        if (detector(traps) %in% .localstuff$polydetectors) {
            if (is.null(xy(capthist)))
                NA
            else {
                lxy <- split (xy(capthist), animalID(capthist))
                lapply (lxy, movexy)
            }
        }
        else {
            w <- split(trap(capthist, names=F), animalID(capthist))
            lapply(w,movex)
        }
    }
}
############################################################################################


ARL <- function (capthist, min.recapt = 1, plt = FALSE, full = FALSE, userdist = NULL,
                 mask = NULL) {
    if (inherits (capthist, 'list')) {
        lapply(capthist, ARL, plt = plt, full = full, userdist = userdist, mask)   ## recursive
    }
    else {
        MMDMx <- function (cx) {
            cx <- abs(cx)  # no special trt for deads
            if (sum(cx>0, na.rm=T) == 1) NA
            else {
              ## x <- traps$x[cx]
              ## y <- traps$y[cx]
              ## max(dist(cbind(x,y)))
              as.numeric(max(distmat[cx, cx]))
            }
        }
        MMDMxy <- function (xy) {
            max(dist(cbind(xy$x, xy$y)))
        }
        traps <- traps(capthist)
        if (!(detector(traps) %in% .localstuff$individualdetectors))
            stop ("require individual detector type for ARL")
        distmat <- valid.userdist(userdist, detector(traps), traps, traps, mask )
        prox  <- length(dim(capthist)) > 2
        if (detector(traps) %in% .localstuff$polydetectors) {
            if (is.null(xy(capthist)))
                stop("no xy coordinates")
            else {
                lxy <- split (xy(capthist), animalID(capthist))
                maxd <- unlist(lapply (lxy, MMDMxy))
                n <- unlist(lapply (lxy, nrow))
            }
        }
        else {
            w <- split(trap(capthist, names=F), animalID(capthist))
            maxd <- unlist(lapply(w, MMDMx))
            n <- unlist(lapply(w, length))
        }
        maxd <- maxd[n>min.recapt]
        n <- n[n>min.recapt]

        temp <- try(coef(nls (maxd ~ aa * (1 - exp(bb * (n-1))),
            start= list (aa = max(maxd)*1.2, bb = -0.4))))

        if (inherits(temp, "try-error")) {
            warning ("nls failure")
            aa <- NA
            bb <- NA
        }
        else {
            aa <- temp[1]
            bb <- temp[2]
            if (plt) {
                plot (jitter(n, amount=0.1), maxd,
                    xlim=c(0,max(c(n,ncol(capthist)))),
                    xlab='Number of captures', ylab='ARL')
                xv <- seq(2,max(n),0.01)
                lines(xv, aa * (1 - exp(bb * (xv-1))))
            }
        }
        attr(aa,'names') <- NULL
        attr(bb,'names') <- NULL

        if (!full) aa
        else list (ARL = aa, b = bb, data = data.frame(maxd = maxd, n=n))
    }
}
############################################################################################

MMDM <- function (capthist, min.recapt = 1, full = FALSE, userdist = NULL, mask = NULL) {
    if (inherits (capthist, 'list')) {
        lapply(capthist, MMDM, full = full, userdist = userdist, mask = mask)   ## recursive
    }
    else {
        MMDMx <- function (cx) {
              cx <- abs(cx)  # no special trt for deads
              if (sum(cx>0, na.rm=T) == 1) NA
              else {
                ## x <- traps$x[cx]
                ## y <- traps$y[cx]
                ## max(dist(cbind(x,y)))
                as.numeric(max(distmat[cx, cx]))
              }
        }
        MMDMxy    <- function (xy) {
                max(dist(cbind(xy$x, xy$y)))
        }
        traps <- traps(capthist)
        distmat <- valid.userdist(userdist, detector(traps), traps, traps, mask )
        if (!(detector(traps) %in% .localstuff$individualdetectors))
            stop ("require individual detector type for MMDM")
        if (detector(traps) %in% .localstuff$polydetectors) {
            if (is.null(xy(capthist)))
                stop ("no xy coordinates")
            else {
                lxy <- split (xy(capthist), animalID(capthist))
                maxd <- unlist(lapply (lxy, MMDMxy))
                n <- unlist(lapply (lxy, nrow))
            }
        }
        else {
            ## streamlined 2010 03 30
            w <- split(trap(capthist, names=F), animalID(capthist))
            maxd <- unlist(lapply( w, MMDMx))
            n <- unlist(lapply(w, length))
        }
        maxd <- maxd[n>min.recapt]
        n <- n[n>min.recapt]

        temp <- mean (maxd, na.rm = TRUE)

        if (!full) temp
        else {
            SE <- function(x) sqrt(var(x, na.rm=T)/sum(!is.na(x)))
            summaryd <- data.frame (Ncapt = names(table(n)),
                           n = as.numeric(table(n)),
                           mean = tapply(maxd, n, mean, na.rm=T),
                           se = tapply(maxd, n, SE))
            summaryd$mean[is.na(summaryd$mean)] <- NA  ## tidier
            summaryd <- rbind(summaryd, data.frame(Ncapt='Total', n=sum(!is.na(maxd)),
                mean=temp, se=SE(maxd)))
            list (MMDM = temp, data = data.frame(maxd = maxd, n=n), summary = summaryd )
        }
    }
}
############################################################################################

RPSV <- function (capthist, CC = FALSE)
{
    if (inherits (capthist, 'list')) {
        lapply(capthist, RPSV, CC)   ## recursive
    }
    else {
        RPSVx <- function (cx) {
            cx <- abs(cx)  # no special trt for deads
            x <- traps$x[cx]
            y <- traps$y[cx]
            n <- length(x)
            c(n = n-1, ssx = sum(x^2) - (sum(x))^2/n, ssy = sum(y^2) - (sum(y))^2/n)
        }
        RPSVxy <- function (xy) {
            x <- xy$x
            y <- xy$y
            n <- length(x)
            c(n = n-1, ssx = sum(x^2) - (sum(x))^2/n, ssy = sum(y^2) - (sum(y))^2/n)
        }
        ## 2014-12-04
        if (nrow(capthist) < 1) return(NA)
        traps <- traps(capthist)
        if (!(detector(traps) %in% .localstuff$individualdetectors))
            stop ("require individual detector type for RPSV")
        if (detector(traps) %in% .localstuff$polydetectors) {
  
            if (is.null(xy(capthist)))
                temp <- NA
            else {
                lxy <- split ( xy(capthist), animalID(capthist))
                temp <- lapply (lxy, RPSVxy)
            }
        }
        else {
            w <- split(trap(capthist, names=F), animalID(capthist))
            temp <- lapply(w,RPSVx)
        }
        temp <- matrix(unlist(temp), nrow = 3)
        temp <- apply(temp,1,sum, na.rm=T)
        if (CC)
            temp <- sqrt((temp[2]+temp[3]) / (2 * temp[1]))
        else
            temp <- sqrt((temp[2]+temp[3]) / (temp[1]-1))
        attr(temp,'names') <- NULL
        temp
    }
}

############################################################################################
## source ('d:\\density secr 1.3\\secr\\r\\mmdm.R')
## data(Peromyscus)

## MMDM(Peromyscus.WSG, full=T)$summary
##   Ncapt  n     mean        se
## 1     1  9       NA        NA
## 2     2  9 28.32839  9.434631
## 3     3 10 24.05921  9.335062
## 4     4  8 33.87949  5.471227
## 5     5  8 52.37655 15.470420
## 6     6  7 34.24929  5.961350
## 7 Total 42 33.93669  4.495646

## MMDM(Peromyscus.WSG, full=T)$summary[,3:4]/15.2
##       mean        se
## 1       NA        NA
## 2 1.863710 0.6206994
## 3 1.582843 0.6141488
## 4 2.228914 0.3599492
## 5 3.445826 1.0177908
## 6 2.253243 0.3921941
## 7 2.232677 0.2957662   <<<< 0.575?

## cf Otis et al 1978 p 87 Fig 23a

##################################################
## source ('d:\\density secr 1.3\\secr\\r\\ARL.R')
##  data(Peromyscus)
##  ARL(Peromyscus.WSG, full=T)
