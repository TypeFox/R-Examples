orthogram <- function (x, tre=NULL, orthobas = NULL, prox = NULL,
                        nrepet = 999, posinega = 0, tol = 1e-07, cdot = 1.5,
                        cfont.main = 1.5, lwd = 2, nclass,
                        high.scores = 0,alter=c("greater", "less", "two-sided")){

    ## some checks and preliminary assignements
    ## if(!require(ade4)) stop("The ade4 package is not installed.")

    nobs <- length(x)
    alter <- match.arg(alter)

    if(is.numeric(x)&is.vector(x)){
        type <- "numeric"
        ##  } else if(is.factor(x)){
        ##     type <- "factor"
        ##   } else if (inherits(x, "dudi")){
        ##     type <- "dudi"
    } else {
        ## stop("x must be a numeric vector, a factor or a dudi object")
        stop("x must be a numeric vector")
    }
    ##  if(type == "dudi") {
    ##     nobs <- nrow(x$tab)
    ##   } else {
    ##     nobs <- length(x)
    ##   }
    ##   if (!is.null(neig)) {
    ##     orthobas <- scores.neig(neig)
    ##   } else if (!is.null(phylog)) {
    ##     if (!inherits(phylog, "phylog")) stop ("'phylog' expected with class 'phylog'")
    ##     orthobas <- phylog$Bscores
    ##   }

    ## if (is.null(orthobas)){
    ##  stop ("'orthobas','neig','phylog' all NULL")
    ## }

    ## retrieve the orthobasis from a proximity matrix
    if(is.null(orthobas)){
        if(is.null(prox)) { # both orthobas and prox are not given -> default orthobasis
            ## check that tre is provided and valid
            if(is.null(tre)) stop("tre, orthobasis or prox must be provided")
            tre <- as(tre, "phylo4")
            if (is.character(checkval <- checkPhylo4(tre))) stop(checkval)
            orthobas <- treePart(tre, result="orthobasis")
        } else { # else orthobasis from the proxi matrix.
            orthobas <- orthobasis.phylo(prox=prox)
        }
    }

    if (!inherits(orthobas, "data.frame")) stop ("'orthobas' is not a data.frame")
    if (nrow(orthobas) != nobs) stop ("non convenient dimensions")
    if (ncol(orthobas) != (nobs-1)) stop (paste("'orthobas' has",ncol(orthobas),"columns, expected:",nobs-1))
    vecpro <- as.matrix(orthobas)
    npro <- ncol(vecpro)

    w <- t(vecpro/nobs)%*%vecpro
    if (any(abs(diag(w)-1)>tol)) {

        stop("'orthobas' is not orthonormal for uniform weighting")
    }
    diag(w) <- 0
    if ( any( abs(as.numeric(w))>tol) )
        stop("'orthobas' is not orthogonal for uniform weighting")
    if (nrepet < 99) nrepet <- 99
    if (posinega !=0) {
        if (posinega >= nobs-1) stop ("Non convenient value in 'posinega'")
        if (posinega <0) stop ("Non convenient value in 'posinega'")
    }
    if(type!="dudi"){
        if (any(is.na(x)))
            stop("missing value in 'x'")
    }
    if(type == "factor"){
        dudi1 <- dudi.acm(data.frame(x), scannf = FALSE, nf = min(nobs, nlevels(x)))
    }
    if(type == "dudi") {
        if (!all.equal(x$lw, rep(1/nobs, nobs)))
            stop("not implemented for non-uniform row weights")
        dudi1 <- redo.dudi(x, newnf = x$rank)
        if(any(colMeans(dudi1$li)>tol))
            stop("not implemented for non-centered analysis")
    }

    if(type == "numeric") {
        z <- x - mean(x)
        et <- sqrt(mean(z * z))
        if ( et <= tol*(max(z)-min(z))) stop ("No variance")
        z <- z/et
        w <- .C("VarianceDecompInOrthoBasis",
                param = as.integer(c(nobs,npro,nrepet,posinega)),
                observed = as.double(z),
                vecpro = as.double(vecpro),
                phylogram = double(npro),
                phylo95 = double(npro),
                sig025 = double(npro),
                sig975 = double(npro),
                R2Max = double(nrepet+1),
                SkR2k = double(nrepet+1),
                Dmax = double(nrepet+1),
                SCE = double(nrepet+1),
                ratio = double(nrepet+1),
                PACKAGE="adephylo"
                )
    } else {
        w <- .C("MVarianceDecompInOrthoBasis",
                param = as.integer(c(nobs,npro,nrepet,posinega)),
                observed = as.double(as.matrix(dudi1$li)),
                nvar = as.integer(ncol(dudi1$li)),
                inertot = as.double(sum(dudi1$eig)),
                vecpro = as.double(vecpro),
                phylogram = double(npro),
                phylo95 = double(npro),
                sig025 = double(npro),
                sig975 = double(npro),
                R2Max = double(nrepet+1),
                SkR2k = double(nrepet+1),
                Dmax = double(nrepet+1),
                SCE = double(nrepet+1),
                ratio = double(nrepet+1),
                PACKAGE="adephylo"
                )
    }
    ##return(w$phylogram)
    ## multiple graphical window (6 graphs)
    ## 1 pgram
    ## 2 cumulated pgram
    ## 3-6 Randomization tests

    def.par <- par(no.readonly = TRUE)
    on.exit(par(def.par))
    layout (matrix(c(1,1,2,2,1,1,2,2,3,4,5,6),4,3))
    par(mar = c(0.1, 0.1, 0.1, 0.1))
    par(usr = c(0,1,-0.05,1))


    ylim <- max(c(w$phylogram, w$phylo95))
    names(w$phylogram) <- as.character(1:npro)
    phylocum <- cumsum(w$phylogram)
    lwd0=2
    fun <- function (y, last=FALSE) {
        delta <- (mp[2]-mp[1])/3
        sel <- 1:(npro - 1)
        segments(mp[sel]-delta,y[sel],mp[sel]+delta, y[sel],lwd=lwd0)
        if(last) segments(mp[npro]-delta,y[npro],mp[npro]+delta, y[npro],lwd=lwd0)
    }
    sig50 <- (1:npro)/npro
    y0 <- phylocum - sig50
    h.obs <- max(y0)
    x0 <- min(which(y0 == h.obs))
    par(mar = c(3.1, 2.5, 2.1, 2.1))
    if(type == "numeric"){
        z0 <- apply(vecpro, 2, function(x) sum(z * x))
        mp <- barplot(w$phylogram, col = grey(1 - 0.3 * (sign(z0) > 0)), ylim = c(0, ylim * 1.05))
    } else {
        mp <- barplot(w$phylogram, ylim = c(0, ylim * 1.05))
    }
    scores.order <- (1:length(w$phylogram))[order(w$phylogram, decreasing=TRUE)[1:high.scores]]
    fun(w$phylo95,TRUE)
    abline(h = 1/npro)
    if (posinega!=0) {
        verti = (mp[posinega]+mp[posinega+1])/2
        abline (v=verti, col="red",lwd=1.5)
    }
    title(main = "Variance decomposition",font.main=1, cex.main=cfont.main)
    box()
    obs0 <- rep(0, npro)
    names(obs0) <- as.character(1:npro)
    barplot(obs0, ylim = c(-0.05, 1.05))
    abline(h=0,col="white")
    if (posinega!=0) {
        verti = (mp[posinega]+mp[posinega+1])/2
        abline (v=verti, col="red",lwd=1.5)
    }

    title(main = "Cumulative decomposition",font.main=1, cex.main=cfont.main)
    points(mp, phylocum, pch = 21, cex = cdot, type = "b")
    segments(mp[1], 1/npro, mp[npro], 1, lty = 1)
    fun(w$sig975)
    fun(w$sig025)
    arrows(mp[x0], sig50[x0], mp[x0], phylocum[x0], angle = 15, length = 0.15,
           lwd = 2)
    box()
    if (missing(nclass)) {
        nclass <- as.integer (nrepet/25)
        nclass <- min(c(nclass,40))
    }
    plot(as.randtest (w$R2Max[-1],w$R2Max[1],call=match.call()),main = "R2Max",nclass=nclass)
    if (posinega !=0) {
        plot(as.randtest (w$ratio[-1],w$ratio[1],call=match.call()),main = "Ratio",nclass=nclass)
    } else {
        plot(as.randtest (w$SkR2k[-1],w$SkR2k[1],call=match.call()),main = "SkR2k",nclass=nclass)
    }
    plot(as.randtest (w$Dmax[-1],w$Dmax[1],call=match.call()),main = "DMax",nclass=nclass)
    plot(as.randtest (w$SCE[-1],w$SCE[1],call=match.call()),main = "SCE",nclass=nclass)

    w$param <- w$observed <- w$vecpro <- NULL
    w$phylo95 <- w$sig025 <- w$sig975 <- NULL
    if (posinega==0) {
        w <- as.krandtest(obs=c(w$R2Max[1],w$SkR2k[1],w$Dmax[1],w$SCE[1]),sim=cbind(w$R2Max[-1],w$SkR2k[-1],w$Dmax[-1],w$SCE[-1]),names=c("R2Max","SkR2k","Dmax","SCE"),alter=alter,call=match.call())
    } else {
        w <- as.krandtest(obs=c(w$R2Max[1],w$SkR2k[1],w$Dmax[1],w$SCE[1],w$ratio[1]),sim=cbind(w$R2Max[-1],w$SkR2k[-1],w$Dmax[-1],w$SCE[-1],w$ratio[-1]),names=c("R2Max","SkR2k","Dmax","SCE","ratio"),alter=alter,call=match.call())
    }

    if (high.scores != 0)
        w$scores.order <- scores.order
    return(w)
} # end orthogram
