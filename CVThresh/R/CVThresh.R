cvwavelet <- function(
    y=y, ywd=ywd,
    cv.optlevel, cv.bsize=1, cv.kfold, cv.random=TRUE, cv.tol=0.1^3, cv.maxiter=100,
    impute.vscale="independent", impute.tol=0.1^3, impute.maxiter=100,
    filter.number=10, family="DaubLeAsymm", thresh.type ="soft", ll=3)
{
    if (impute.vscale != "independent" && impute.vscale != "level") stop("impute.vscale must be one of independent or level.")

    cv.index <- cvtype(n=length(y), cv.bsize=cv.bsize, cv.kfold=cv.kfold, cv.random=cv.random)$cv.index

    yimpute <- cvimpute.by.wavelet(y=y, impute.index=cv.index, impute.tol=impute.tol,
                        impute.maxiter=impute.maxiter, impute.vscale=impute.vscale,
                        filter.number=10, family="DaubLeAsymm", ll=3)$yimpute

    tmpout <- cvwavelet.after.impute(y=y, ywd=ywd, yimpute=yimpute, cv.index=cv.index,
                     cv.optlevel=cv.optlevel, cv.tol=cv.tol, cv.maxiter=cv.maxiter,
                     filter.number=10, family="DaubLeAsymm", thresh.type="soft", ll=3)

    return(list(y=y, yimpute=yimpute, yc=tmpout$yc, cvthresh=tmpout$cvthresh))
}


cvtype <- function(n, cv.bsize=1, cv.kfold, cv.random=TRUE) {
    if(n < cv.bsize * cv.kfold) stop("Block size or no. of fold is too large.")
    if(n <= 0) stop("The number of data must be greater than 0.")    
    
    cv.index <- tmp.index <- NULL

    if(cv.random) {
        cv.nblock <- trunc(n / (cv.bsize * cv.kfold))
        for (k in 1:cv.kfold)
            tmp.index <- rbind(tmp.index,
                sort(sample(seq(sample(1:cv.bsize, 1), n-cv.bsize+1, by=cv.bsize), cv.nblock)))
    } else {
        for (k in 1:cv.kfold)
            tmp.index <- rbind(tmp.index,
                seq(1, n-cv.bsize*(cv.kfold-1), by=cv.bsize*cv.kfold)+ cv.bsize * (k-1))
    }

    if(cv.bsize > 1) {
        cv.index <- tmp.index
        for (i in 1:(cv.bsize-1))
            cv.index <- cbind(cv.index, tmp.index + i)
        cv.index <- t(apply(cv.index, 1, sort))
    } else {
        cv.index <- tmp.index
    }

    list(cv.index=cv.index)
}


cvimpute.by.wavelet <- function(
    y, impute.index, impute.tol=0.1^3, impute.maxiter=100, impute.vscale="independent",
    filter.number=10, family="DaubLeAsymm", ll=3)
{
    if (impute.vscale != "independent" && impute.vscale != "level") stop("impute.vscale must be one of independent or level.")

    n <- length(y)
    sdy <- sd(y)

    impute.nrep <- nrow(impute.index)
    fracmiss <- ncol(impute.index) / n

    yimpute <- y
    slevel <- log2(n) - ll
    nlevelm1 <- log2(n) - 1
    ll.nlevelm1 <- ll:nlevelm1

    for (k in 1:impute.nrep) {

        tmpyimpute <- y
        tmpyimputewd <- wd(tmpyimpute, filter.number=filter.number, family=family)

        if(impute.vscale=="independent") {
            sdev0 <- mad(accessD(tmpyimputewd, nlevelm1))
            tmpyimpute[impute.index[k,]] <-
                wr(ebayesthresh.wavelet(tmpyimputewd, vscale=sdev0, smooth.levels=slevel))[impute.index[k,]]
        }

        if(impute.vscale=="level") {
            sdev0 <- NULL
            for (i in ll.nlevelm1) {
                tmpsdev0 <- mad(accessD(tmpyimputewd, i))
                tmpyimputewd <- putD(tmpyimputewd, i, ebayesthresh(accessD(tmpyimputewd, i), sdev=tmpsdev0))
                sdev0 <- c(sdev0, tmpsdev0)
            }
            tmpyimpute[impute.index[k,]] <- wr(tmpyimputewd)[impute.index[k,]]
        }

        j <- 0
        repeat{
            tmpyimputewd <- wd(tmpyimpute, filter.number=filter.number, family=family)

            if(impute.vscale=="independent") {
                sdev1 <- mad(accessD(tmpyimputewd, nlevelm1))
                sdev1 <- sqrt(sdev1^2 + fracmiss * sdev0^2)
                tmpyimputewr <- wr(ebayesthresh.wavelet(tmpyimputewd, vscale=sdev1, smooth.levels=slevel))
            }

            if(impute.vscale=="level") {
                sdev1 <- NULL
                for (i in ll.nlevelm1) {
                    tmpsdev1 <- mad(accessD(tmpyimputewd, i))
                    tmpsdev1 <- sqrt(tmpsdev1^2 + fracmiss * sdev0[i-ll+1]^2)
                    tmpyimputewd <- putD(tmpyimputewd, i, ebayesthresh(accessD(tmpyimputewd, i), sdev=tmpsdev1))
                    sdev1 <- c(sdev1, tmpsdev1)
                }
                tmpyimputewr <- wr(tmpyimputewd)
            }

            ydiff <- mean(abs(tmpyimpute[impute.index[k,]] - tmpyimputewr[impute.index[k,]])) / sdy
            if (ydiff < impute.tol || j > impute.maxiter) {
                tmpyimpute[impute.index[k,]] <- tmpyimputewr[impute.index[k,]]
                break
            }

            tmpyimpute[impute.index[k,]] <- tmpyimputewr[impute.index[k,]]
            sdev0 <- sdev1
            j <- j + 1
        } # End of repeat loop
        yimpute[impute.index[k,]] <- tmpyimpute[impute.index[k,]]
    } # End of k-th impute loop

    return(list(yimpute=yimpute))
}


cvwavelet.after.impute <- function(
    y, ywd, yimpute,
    cv.index, cv.optlevel, cv.tol=0.1^3, cv.maxiter=100,
    filter.number=10, family="DaubLeAsymm", thresh.type="soft", ll=3)
{ ### We adapt grid search algorithm of R function "WaveletCV" in WaveThresh3 of Nason (1998).
  ### This algorithm is a special form of grid search algorithm called golden section search.
  ### When bisectioning the interval containing solution, the golden number is used to keep
  ### the ratio of two sectioned subinterval fixed.

    if (length(y) != length(yimpute)) stop("The length of data and imputed values must be the same.")
    
    n <- length(y)
    cd.kfold <- nrow(cv.index)
    nlevel <- log2(n)
    ll.nlevelm1 <- ll:(nlevel - 1)

    cv.ndim <- length(cv.optlevel)

    R <- (sqrt(5)-1)/2 # R is the golden number 0.61803399000000003 
    C <- 1 - R

    lambda <- matrix(0, 4, cv.ndim)

    if(cv.ndim != 1) {
        ethresh <- NULL
        for (i in ll.nlevelm1)
            ethresh <- c(ethresh, ebayesthresh(accessD(ywd, i), verbose=TRUE)$threshold.origscale)
        lambda.range <- 2 * ethresh[cv.optlevel + diff(c(cv.optlevel, nlevel))-ll]
    }
    else {
        lambda.range <- threshold(ywd, policy="universal", type=thresh.type, return.threshold=TRUE, lev=ll.nlevelm1)[1]
    }

    lambda[1, ] <- rep(0, cv.ndim)
    lambda[4, ] <- lambda.range
    lambda[2, ] <- lambda[4, ]/2
    lambda[3, ] <- lambda[2, ] + C * (lambda[4, ] - lambda[2, ])
    optlambda <- lambda[3, ]

    tmpyimputewr <- y

    for (k in 1:cd.kfold) {
        yimpute2 <- y
        yimpute2[cv.index[k,]] <- yimpute[cv.index[k,]]
        tmpyimputewr[cv.index[k,]] <- wr(threshold(wd(yimpute2, filter.number=filter.number, family=family),
                                policy="manual", value=rep(optlambda, diff(c(cv.optlevel, nlevel))),
                                by.level=TRUE, type=thresh.type, lev=ll.nlevelm1))[cv.index[k,]]
    }
    perr <- mean((tmpyimputewr - y)^2)

    j <- 0
    repeat{

        for (i in 1:cv.ndim) {

            optlambda[i] <- lambda[2, i]
            for (k in 1:cd.kfold) {
                yimpute2 <- y
                yimpute2[cv.index[k,]] <- yimpute[cv.index[k,]]
                tmpyimputewr[cv.index[k,]] <- wr(threshold(wd(yimpute2, filter.number=filter.number, family=family),
                                    policy="manual", value=rep(optlambda, diff(c(cv.optlevel, nlevel))),
                                    by.level=TRUE, type=thresh.type, lev=ll.nlevelm1))[cv.index[k,]]
            }
            f2 <- mean((tmpyimputewr - y)^2)

            optlambda[i] <- lambda[3, i]
            for (k in 1:cd.kfold) {
                yimpute2 <- y
                yimpute2[cv.index[k,]] <- yimpute[cv.index[k,]]
                tmpyimputewr[cv.index[k,]] <- wr(threshold(wd(yimpute2, filter.number=filter.number, family=family),
                                    policy="manual", value=rep(optlambda, diff(c(cv.optlevel, nlevel))),
                                    by.level=TRUE, type=thresh.type, lev=ll.nlevelm1))[cv.index[k,]]
            }
            f3 <- mean((tmpyimputewr - y)^2)

            if(f3 < f2) {
                optlambda[i] <- lambda[3, i]
                optf <- f3
                lambda[1, i] <- lambda[2, i]
                lambda[2, i] <- lambda[3, i]
                lambda[3, i] <- R * lambda[2, i] + C * lambda[4, i]

            }
            else {
                optlambda[i] <- lambda[2, i]
                optf <- f2
                lambda[4, i] <- lambda[3, i]
                lambda[3, i] <- lambda[2, i]
                lambda[2, i] <- R * lambda[3, i] + C * lambda[1, i]
            }
            perr <- c(perr, optf)
        }

        stopping <- NULL
        for (i in 1:cv.ndim)
            stopping <- c(stopping, abs(lambda[4, i] - lambda[1, i]) / (abs(lambda[2, i]) + abs(lambda[3, i])))
        if (all(stopping < cv.tol) || abs(perr[cv.ndim*(j+1)+1]-perr[cv.ndim*j+1])/perr[cv.ndim*j+1] < cv.tol || j > cv.maxiter) break
        j <- j + 1
    }

    yc <- wr(threshold(ywd, policy="manual", value=rep(optlambda, diff(c(cv.optlevel, nlevel))),
                        by.level=TRUE, type=thresh.type, lev=ll.nlevelm1))
    if(cv.ndim !=1)
        list(yc=yc, cvthresh=rep(optlambda, diff(c(cv.optlevel, nlevel))))
    else
        list(yc=yc, cvthresh=optlambda)
}


cvwavelet.image <- function(
    images, imagewd,
    cv.optlevel, cv.bsize=c(1,1), cv.kfold, cv.tol=0.1^3, cv.maxiter=100,
    impute.tol=0.1^3, impute.maxiter=100,
    filter.number=2, ll=3)
{
    if(!is.matrix(images)) stop("images must be a matrix format.")

    tmp <- cvtype.image(n=c(nrow(images), ncol(images)), cv.bsize=cv.bsize, cv.kfold=cv.kfold)
    cv.index1 <- tmp$cv.index1
    cv.index2 <- tmp$cv.index2

    imageimpute <- cvimpute.image.by.wavelet(images=images, impute.index1=cv.index1, impute.index2=cv.index2,
                            impute.tol=impute.tol, impute.maxiter=impute.maxiter,
                            filter.number=filter.number, ll=ll)$imageimpute

    tmpout <- cvwavelet.image.after.impute(images=images, imagewd=imagewd, imageimpute=imageimpute,
                     cv.index1=cv.index1, cv.index2=cv.index2,
                     cv.optlevel=cv.optlevel, cv.tol=cv.tol, cv.maxiter=cv.maxiter,
                     filter.number=2, ll=ll)

    list(imagecv=tmpout$imagecv, cvthresh=tmpout$optlambda)
}


cvtype.image <- function(n, cv.bsize=c(1,1), cv.kfold) {
    if(length(n) != 2 || length(cv.bsize) != 2) stop("Two dimension")
    if(n[1]*n[2] < cv.bsize[1] * cv.bsize[2] * cv.kfold) stop("Block size or no. of fold is too large.")

    cv.perc <- 1 / cv.kfold
    cv.nblock <- trunc(n * sqrt(cv.perc) / cv.bsize)

    for (j in 1:2) {
        cv.index <- tmp.index <- NULL
        for (k in 1:cv.kfold)
            tmp.index <- rbind(tmp.index,
                sort(sample(seq(sample(1:cv.bsize[j], 1), n[j]-cv.bsize[j]+1, by=cv.bsize[j]), cv.nblock[j])))

        if(cv.bsize[j] > 1) {
            cv.index <- tmp.index
            for (i in 1:(cv.bsize[j]-1))
                cv.index <- cbind(cv.index, tmp.index + i)
            cv.index <- t(apply(cv.index, 1, sort))
        } else {
            cv.index <- tmp.index
        }
        #assign(paste("cv.index", j, sep=""), cv.index)
        if(j==1) cv.index1 <- cv.index
        if(j==2) cv.index2 <- cv.index        
    }
    list(cv.index1=cv.index1, cv.index2=cv.index2)
}


cvimpute.image.by.wavelet <- function(
   images, impute.index1, impute.index2, impute.tol=0.1^3, impute.maxiter=100,
   filter.number=2, ll=3)
{
    if(!is.matrix(images)) stop("images must be a matrix format.")

    fracmiss <- (ncol(impute.index1) * ncol(impute.index2)) / (nrow(images) * ncol(images))
    sdimage <- sd(as.numeric(images))

    tmp <- c(8,9,10)
    slevel <- NULL
    for (k in 1:(log2(nrow(images))-ll))
        slevel <- c(slevel, tmp+4*(k-1))
    slevel <- slevel[length(slevel):1]

    yimpute <- images

    for (k in 1:nrow(impute.index1)) {
        ym <- images
        ymwd <- imwd(ym, filter.number=filter.number)
        sdev0 <- mad(c(ymwd[[8]], ymwd[[9]], ymwd[[10]]))
        for (i in 1:length(slevel))
            ymwd[[slevel[i]]] <- ebayesthresh(ymwd[[slevel[i]]], sdev=sdev0)

        ym[impute.index1[k, ], impute.index2[k, ]] <- imwr(ymwd)[impute.index1[k, ], impute.index2[k, ]]

        j <- 0
        repeat{
            ymwd <- imwd(ym, filter.number=filter.number)
            sdev1 <- mad(c(ymwd[[8]], ymwd[[9]], ymwd[[10]]))
            sdev1 <- sqrt(sdev1^2 + fracmiss * sdev0^2)
            for (i in 1:length(slevel))
                ymwd[[slevel[i]]] <- ebayesthresh(ymwd[[slevel[i]]], sdev=sdev1)

            tmpymwr <- imwr(ymwd)

            ydiff <- mean(abs(ym[impute.index1[k, ], impute.index2[k, ]] -
                              tmpymwr[impute.index1[k, ], impute.index2[k, ]])) / sdimage
            if (ydiff < impute.tol || j > impute.maxiter) {
                ym[impute.index1[k, ], impute.index2[k, ]] <- tmpymwr[impute.index1[k, ], impute.index2[k, ]]
                break
            }

            sdev0 <- sdev1
            ym[impute.index1[k, ], impute.index2[k, ]] <- tmpymwr[impute.index1[k, ], impute.index2[k, ]]
            j <- j + 1
        }
        yimpute[impute.index1[k, ], impute.index2[k, ]] <- ym[impute.index1[k, ], impute.index2[k, ]]
    }
    list(imageimpute=yimpute)
}


cvwavelet.image.after.impute <- function(
   images, imagewd, imageimpute,
   cv.index1=cv.index1, cv.index2=cv.index2,
   cv.optlevel=cv.optlevel, cv.tol=cv.tol, cv.maxiter=cv.maxiter,
   filter.number=2, ll=3)
{ ### We adapt grid search algorithm of R function "WaveletCV" in WaveThresh3 of Nason (1998).
  ### This algorithm is a special form of grid search algorithm called golden section search.
  ### When bisectioning the interval containing solution, the golden number is used to keep
  ### the ratio of two sectioned subinterval fixed.

    if(!is.matrix(images)) stop("images must be a matrix format.")
  
    sdimage <- sd(as.numeric(images))
    cv.kfold <- nrow(cv.index1)
    tmp <- c(8,9,10)
    slevel <- NULL

    for (i in 1:(diff(range(cv.optlevel))+1))
        slevel <- c(slevel, tmp+4*(i-1))
    slevel <- slevel[length(slevel):1]

    sdev <- mad(c(imagewd[[8]], imagewd[[9]], imagewd[[10]]))

    lambda.range <- NULL
    for (i in slevel)
        lambda.range <- c(lambda.range, ebayesthresh(imagewd[[i]], sdev=sdev, verbose=TRUE)$threshold.origscale)

    R <- (sqrt(5)-1)/2 # R is the golden number 0.61803399000000003 
    C <- 1 - R

    ndim <- (diff(range(cv.optlevel))+1) * 3
    lambda <- matrix(0, 4, ndim)

    lambda[1, ] <- rep(0, ndim)
    lambda[4, ] <- 2.0 * lambda.range
    lambda[2, ] <- lambda[4, ]/2
    lambda[3, ] <- lambda[2, ] + C * (lambda[4, ] - lambda[2, ])
    optlambda <- lambda[3, ]

    tmpyimputewr <- images

    for (k in 1:cv.kfold) {
        yimpute2 <- images

        yimpute2[cv.index1[k, ], cv.index2[k, ]] <- imageimpute[cv.index1[k, ], cv.index2[k, ]]
        ywdth <- imwd(yimpute2, filter.number=filter.number)
        for (j in 1:ndim)
            ywdth[[slevel[j]]] <- threshld(ywdth[[slevel[j]]], optlambda[j], hard=FALSE)
        tmpyimputewr[cv.index1[k, ], cv.index2[k, ]] <- imwr(ywdth)[cv.index1[k, ], cv.index2[k, ]]
    }
    perr <- mean((tmpyimputewr - images)^2)

    j <- 0
    repeat{

        for (i in 1:ndim) {

            optlambda[i] <- lambda[2, i]
            for (k in 1:cv.kfold) {
                yimpute2 <- images
                yimpute2[cv.index1[k, ], cv.index2[k, ]] <- imageimpute[cv.index1[k, ], cv.index2[k, ]]
                ywdth <- imwd(yimpute2, filter.number=filter.number)
                for (m in 1:ndim)
                    ywdth[[slevel[m]]] <- threshld(ywdth[[slevel[m]]], optlambda[m], hard=FALSE)
                tmpyimputewr[cv.index1[k, ], cv.index2[k, ]] <- imwr(ywdth)[cv.index1[k, ], cv.index2[k, ]]
            }
            f2 <- mean((tmpyimputewr - images)^2)

            optlambda[i] <- lambda[3, i]
            for (k in 1:cv.kfold) {
                yimpute2 <- images
                yimpute2[cv.index1[k, ], cv.index2[k, ]] <- imageimpute[cv.index1[k, ], cv.index2[k, ]]
                ywdth <- imwd(yimpute2, filter.number=filter.number)
                for (m in 1:ndim)
                    ywdth[[slevel[m]]] <- threshld(ywdth[[slevel[m]]], optlambda[m], hard=FALSE)
                tmpyimputewr[cv.index1[k, ], cv.index2[k, ]] <- imwr(ywdth)[cv.index1[k, ], cv.index2[k, ]]
            }
            f3 <- mean((tmpyimputewr - images)^2)

            if(f3 < f2) {
                optlambda[i] <- lambda[3, i]
                optf <- f3
                lambda[1, i] <- lambda[2, i]
                lambda[2, i] <- lambda[3, i]
                lambda[3, i] <- R * lambda[2, i] + C * lambda[4, i]
            }
            else {
                optlambda[i] <- lambda[2, i]
                optf <- f2
                lambda[4, i] <- lambda[3, i]
                lambda[3, i] <- lambda[2, i]
                lambda[2, i] <- R * lambda[3, i] + C * lambda[1, i]
            }
            perr <- c(perr, optf)
        }

        stopping <- NULL
        for (i in 1:ndim)
            stopping <- c(stopping, abs(lambda[4, i] - lambda[1, i]) / (abs(lambda[2, i]) + abs(lambda[3, i])))
        if (all(stopping < cv.tol) || abs(perr[ndim*(j+1)+1]-perr[ndim*j+1]) < (cv.tol * sdimage) ||
            j > cv.maxiter) break
        j <- j + 1
    }

    for (j in 1:ndim)
        imagewd[[slevel[j]]] <- threshld(imagewd[[slevel[j]]], optlambda[j], hard=FALSE)

    imagecv <- imwr(imagewd)

    list(imagecv=imagecv, cvthresh=optlambda)
}

heav <- function(norx = 1024) {
    if (length(norx) == 1) {
        if(norx <= 0) stop("The number of data must be greater than 0.") 
        x <- seq(0, 1, length = norx + 1)[1:norx]           
    } else
        x <- norx
    meanf <- 4 * sin(4*pi*x) - sign(x-0.3)-sign(0.72-x)
    list(x = x, meanf = meanf, sdf=2.969959)
}

dopp <- function(norx = 1024) {
    if (length(norx) == 1) {
        if(norx <= 0) stop("The number of data must be greater than 0.") 
        x <- seq(0, 1, length = norx + 1)[1:norx]           
    } else
        x <- norx
    meanf <- (x*(1-x))^0.5 * sin(2*pi*1.05/(x+0.05))
    list(x = x, meanf = meanf, sdf=0.2889965)
}

ppoly <- function(norx = 1024) {
    if (length(norx) == 1) {
        if(norx <= 0) stop("The number of data must be greater than 0.") 
        x <- seq(0, 1, length = norx + 1)[1:norx]           
    } else
        x <- norx
    meanf <- rep(0, length(x))
    xsv <- (x <= 0.5)               # Left hand end
    meanf[xsv] <- -16 * x[xsv]^3 + 12 * x[xsv]^2
    xsv <- (x > 0.5) & (x <= 0.75)  # Middle section
    meanf[xsv] <- (x[xsv] * (16 * x[xsv]^2 - 40 * x[xsv] + 28))/3 - 1.5
    xsv <- x > 0.75                 # Right hand end
    meanf[xsv] <- (x[xsv] * (16 * x[xsv]^2 - 32 * x[xsv] + 16))/3
    list(x = x, meanf = meanf, sdf=0.3033111)
}

fg1 <- function(norx = 1024) {
    if (length(norx) == 1) {
        if(norx <= 0) stop("The number of data must be greater than 0.") 
        x <- seq(0, 1, length = norx + 1)[1:norx]           
    } else
        x <- norx
    meanf <- 0.25 * ((4 * x - 2.0) + 2 * exp(-16 * (4 * x - 2.0)^2))
    list(x = x, meanf = meanf, sdf=0.3159882)
}
 
