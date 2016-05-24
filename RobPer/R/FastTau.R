FastTau <- function(x, y, taucontrol=list(N=500, kk=2, tt=5, rr=2, approximate=0), beta_gamma) {
    # Fast-Tau algorithm for linear regression
    # This is a slightly changed version of the R-function FastTau published in
    # Salibian-Barrera, M., Willems, G. und Zamar, R. (2008): The fast-tau estimator for regression.
    # Journal of Computational and Graphical Statistics, 17 Nr. 3, 659-682.
    # For details concerning the changes see Thieler, Fried and Rathjens (2013).
	
    if(is.null(taucontrol$N)) N <- 500 else N <- taucontrol$N
    if(is.null(taucontrol$kk)) kk <- 2 else kk <- taucontrol$kk
    if(is.null(taucontrol$tt)) tt <- 5 else tt <- taucontrol$tt
    if(is.null(taucontrol$rr)) rr <- 22 else rr <- taucontrol$rr
    if(is.null(taucontrol$approximate)) approximate <- FALSE else approximate <- taucontrol$approximate
                            
    if (tt<1) stop("parameter tt should be at least 1")

    x <- as.matrix(x)
    x.red<- unique(x) # no double rows

    n <- nrow(x)
    p <- ncol(x)

    n.red <- nrow(x.red)
    if(n.red< n) {
        ff <- function(x){
            x<-as.factor(x)
            levels(x)<- (1: length(levels(x)))
            return(x)
        }
        index<-apply(apply(x,2,ff),1,paste, collapse="_")
        index.red<- unique(index)  # index for different regressors
    }
	
    isStepmodel <- all(apply(x!=0,1,sum)==1)&(p>1)

    c1 <- .4046
    b1 <- .5
    c2 <- 1.09
    b2 <- .1278

    RWLStol <- 1e-11

    bestbetas <- matrix(0, p, tt)
    bestscales <- 1e20 * rep(1, tt)
    besttauscales <- 1e20 * rep(1, tt)
    worsti <- 1
    rworst <- y

    rhoOpt <- function(x, cc) {
        tmp <- x^2 / 2 / (3.25*cc^2)
        tmp2 <- (1.792 - 0.972 * x^2 / cc^2 + 0.432 * x^4 / cc^4 - 0.052 * x^6 / cc^6 + 0.002 * x^8 / cc^8) / 3.25
        tmp[abs(x) > 2*cc] <- tmp2[abs(x) > 2*cc]
        tmp[abs(x) > 3*cc] <- 1
        tmp
    }

    Mscale <- function(u, b, c, initialsc) {
        if (initialsc==0) initialsc = median(abs(u))/.6745
        maxit <- 100
        sc <- initialsc
        i <- 0 
        eps <- 1e-10
        err <- 1
        while  (( i < maxit ) & (err > eps)) {
            sc2 <- sqrt( sc^2 * mean(rhoOpt(u/sc,c)) / b)
            err <- abs(sc2/sc - 1)
            sc <- sc2
            i <- i+1
        }
        return(sc)
    }

    stepsampler <- function(hnr, nichtnull) sample(rep(nichtnull[which(nichtnull[,2]==hnr),1],2),1)

    for (i in 1:N) {
        # find a p-subset in general position.
        if(isStepmodel) {# which means in this case: One Point from every step
            ranset <- apply(cbind(1:p), 1, stepsampler, nichtnull=which(x!=0, arr.ind=TRUE))
            xj <- x[ranset,]
            yj <- y[ranset] 
            bj <- as.matrix(qr.coef(qr(xj),yj))
        } else {
            singular <- 1; itertest <- 1
            while (singular==1 && itertest<100) {
                if(n.red<n) {
                    ranind <- sample(index.red, p, replace=FALSE) # choose p different regressor-label
                    zahlen <- which(index%in% ranind) # which indices belong to these?
                    i.red <- index[which(index%in% ranind)]  #which regressor-labels do they have
                    ranset<- tapply(zahlen, INDEX=i.red, sample,1) # only one per regressor-label
                } else {
                    ranset <- sample(n,p)
                }
                xj <- x[ranset,]
                yj <- y[ranset] 
                bj <- as.matrix(qr.coef(qr(xj),yj))
                singular <- any(!is.finite(bj))
                itertest <- itertest + 1
            }
            if (itertest==100) return(list(scale=NA))
        }
		
        # perform kk steps of IRLS on elemental start
        if (kk > 0) {
            tmp <- IWLSiteration(x, y, bj, 0, kk, RWLStol, b1, c1, c2)
            betarw <- tmp$betarw
            resrw <- y - x %*% betarw
            scalerw <- tmp$scalerw
        } else {
            betarw <- bj
            resrw <- y - x %*% betarw
            scalerw <- median(abs(resrw))/.6745
        }
		
        # long-term memory vector, for finding a special extra candidate at the end :
        if (i > 1) {
            LTMvec = LTMvec + abs(resrw)
        } else {
            LTMvec = abs(resrw)
        }
		
        # check whether new subsample yields one of best t tau-objective values
        tempres <- checkbest(y=y,x=x,approximate=approximate, rhoOpt=rhoOpt, resrw=resrw, bestscales=bestscales, besttauscales= besttauscales,
            c1=c1,c2=c2, b1=b1, bestbetas=bestbetas, Mscale= Mscale, scalerw=scalerw, rr=rr, rworst=rworst, worsti=worsti, betarw=betarw)
        bestscales    <-tempres$bestscales
        besttauscales <-tempres$besttauscales
        bestbetas     <-tempres$bestbetas
        rworst        <-tempres$rworst
        worsti        <-tempres$worsti
    }
  
    # consider an extra subsample, made up of badly fit observations

    IXLTM <- order(LTMvec, decreasing=T)
    singular <- 1 
    extrasize <- p
    while (singular==1) {
        xs <- x[IXLTM[1:extrasize],]
        ys <- y[IXLTM[1:extrasize]]
        bbeta <- as.matrix(qr.coef(qr(xs),ys))
        singular <- any(!is.finite(bbeta))
        extrasize <- extrasize + 1
    }

    # perform kk steps of IRLS on elemental start
    if (kk > 0) {
        tmp <- IWLSiteration(x, y, bbeta, 0, kk, RWLStol, b1, c1, c2)
        betarw <- tmp$betarw
        resrw <- y - x %*% betarw
        scalerw <- tmp$scalerw
    } else {
        betarw <- bbeta
        resrw <- y - x %*% betarw
        scalerw <- median(abs(resrw))/.6745
    }

    # check whether this candidate yields one of best t tau-objective values

    tempres <- checkbest(y=y,x=x,approximate=approximate, rhoOpt=rhoOpt, resrw=resrw, bestscales=bestscales, besttauscales= besttauscales,
        c1=c1,c2=c2, b1=b1, bestbetas=bestbetas, Mscale= Mscale, scalerw=scalerw, rr=rr, rworst=rworst, worsti=worsti, betarw=betarw)
    bestscales    <-tempres$bestscales
    besttauscales <-tempres$besttauscales
    bestbetas     <-tempres$bestbetas
    rworst        <-tempres$rworst
    worsti        <-tempres$worsti

    superbesttauscale <- 1e20

    # RWLS-iterate each of the best tt candidates until convergence, and retain the best result
    for (i in 1:tt) {
        tmp <- IWLSiteration(x, y, bestbetas[,i], bestscales[i], 500, RWLStol, b1, c1, c2)
        resrw <- y - x %*% tmp$betarw
        tauscalerw <- tmp$scalerw * sqrt(mean(rhoOpt(resrw/tmp$scalerw,c2)))
        if (tauscalerw < superbesttauscale) {
            superbesttauscale <- tauscalerw
            superbestbeta <- tmp$betarw
            superbestscale <- tmp$scalerw
        }
    }

    superbestscale <- Mscale(y - x%*%superbestbeta, b1, c1, superbestscale)
    superbesttauscale <- superbestscale * sqrt(mean(rhoOpt((y - x%*%superbestbeta)/superbestscale,c2)))

    # Olive's two extra candidates:

    # add LS candidate
    betaLS <- as.matrix(qr.coef(qr(x),y))
    resLS <- y - x %*% betaLS
    scaleLS <- median(abs(resLS))/.6745  
    scaletest1 <- mean(rhoOpt(resLS / superbestscale,c1)) < b1
    scaletest2 <- sum(rhoOpt(resLS / superbestscale,c2)) < sum(rhoOpt((y - x%*%superbestbeta)/superbestscale,c2))
    if (scaletest1 || scaletest2) {
        snew <- Mscale(resLS, b1, c1, scaleLS)
        taunew <- snew * sqrt(mean(rhoOpt(resLS/snew,c2)))
        if (taunew < superbesttauscale) {
            superbestscale <- snew
            superbestbeta <- betaLS
            superbesttauscale <- taunew
        }
    }

    if(missing(beta_gamma)){
        # add HB candidate
        IXmed <- order(abs(y - median(y)))
        xhalf <- x[IXmed[1:floor(n/2)],]
        yhalf <- y[IXmed[1:floor(n/2)]]
        bbeta <- as.matrix(qr.coef(qr(xhalf),yhalf))
        # + 10 C-steps
    } else {
        bbeta <- beta_gamma
    }

    tmp <- IWLSiteration(x, y, bbeta, 0, 10, RWLStol, b1, c1, c2)
    betaHB <- tmp$betarw
    resHB <- y - x %*% betaHB
    scaleHB <- tmp$scalerw
    scaletest1 <- mean(rhoOpt(resHB / superbestscale,c1)) < b1
    scaletest2 <- sum(rhoOpt(resHB / superbestscale,c2)) < sum(rhoOpt((y - x%*%superbestbeta)/superbestscale,c2))
    if (scaletest1 || scaletest2) {
        snew <- Mscale(resHB, b1, c1, scaleHB)
        taunew <- snew * sqrt(mean(rhoOpt(resHB/snew,c2)))
        if (taunew < superbesttauscale) {
            superbestbeta <- betaHB
            superbesttauscale <- taunew
        }
    }
	
    return(list( beta = superbestbeta, scale = superbesttauscale/sqrt(b2) ))
}
