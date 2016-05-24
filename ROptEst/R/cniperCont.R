
.plotData <- function(
  ## helper function for cniper-type plots to plot in data
   data, # data to be plot in
   dots, # dots from the calling function
   origCl, # call from the calling function
   fun, # function to determine risk difference
   L2Fam, # L2Family
   IC # IC1 in cniperContPlot and eta in cniperPointPlot
){
               dotsP <- .makedotsP(dots)
               dotsP$col <- rep(eval(origCl$col.pts), length.out=n)
               dotsP$pch <- rep(eval(origCl$pch.pts), length.out=n)

               al <- eval(origCl$alpha.trsp)
               if(!is.na(al))
                   dotsP$col <- sapply(dotsP$col, addAlphTrsp2col, alpha=al)

               n <- if(!is.null(dim(data))) nrow(data) else length(data)
               if(!is.null(lab.pts))
                    lab.pts <- rep(origCl$lab.pts, length.out=n)

               sel <- .SelectOrderData(data, function(x)sapply(x,fun),
                                       eval(origCl$which.lbs),
                                       eval(origCl$which.Order))
               i.d <- sel$ind
               i0.d <- sel$ind1
               y.d <- sel$y
               x.d <- sel$data
               n <- length(i.d)

               resc.dat <- .rescalefct(x.d, function(x) sapply(x,fun),
                              eval(origCl$scaleX), origCl$scaleX.fct, origCl$scaleX.inv,
                              eval(origCl$scaleY), origCl$scaleY.fct,
                              dots$xlim, dots$ylim, dots)

               dotsP$x <- resc.dat$X
               dotsP$y <- resc.dat$Y

               trafo <- trafo(L2Fam@param)
               dims <- nrow(trafo)
               QF <- diag(dims)
               if(is(IC,"ContIC") & dims>1 )
                      {if (is(normtype(IC),"QFNorm"))
                           QF <- QuadForm(normtype(IC))}

               absInfoEval <- function(x,y) sapply(x, y@Map[[1]])
               IC.rv <- as(diag(dims) %*% IC@Curve, "EuclRandVariable")
               absy.f <- t(IC.rv) %*% QF %*% IC.rv
               absy <- absInfoEval(x.d, absy.f)

               if(is.null(origCl$cex.pts)) origCl$cex.pts <- par("cex")
               dotsP$cex <-  log(absy+1)*3*rep(origCl$cex.pts, length.out=n)

               dotsT <- dotsP
               dotsT$pch <- NULL
               dotsT$cex <- dotsP$cex/2
               dotsT$labels <- if(is.null(lab.pts)) i.d else lab.pts[i.d]
               do.call(points,dotsP)
               if(!is.null(origCl$with.lab))
                   if(origCl$with.lab)  do.call(text,dotsT)
               if(!is.null(origCl$return$order))
                   if(origCl$return.Order) return(i0.d)
        return(invisible(NULL))
        }


.getFunCnip <- function(IC1,IC2, risk, L2Fam, r, b20=NULL){

        riskfct <- getRiskFctBV(risk, biastype(risk))

       .getTrVar <- function(IC){
           R <- Risks(IC)[["trAsCov"]]
           if(is.null(R)) R <- getRiskIC(IC, risk = trAsCov(), L2Fam = L2Fam)
           if(length(R) > 1) R <- R$value
           return(R)
        }
        R1 <- .getTrVar (IC1)
        R2 <- .getTrVar (IC2)


        fun <- function(x){
            y1 <- evalIC(IC1,as.matrix(x,ncol=1))
            r1 <- riskfct(var=R1,bias=r*fct(normtype(risk))(y1))
            if(!is.null(b20))
               r2 <- riskfct(var=R1,bias=b20) else{
               y2 <- sapply(x,function(x0) evalIC(IC2,x0))
               r2 <- riskfct(var=R2,bias=r*fct(normtype(risk))(y2))
            }
            r1 - r2
        }

        return(fun)
}

cniperCont <- function(IC1, IC2, data = NULL, ...,
                           neighbor, risk, lower=getdistrOption("DistrResolution"),
                           upper=1-getdistrOption("DistrResolution"), n = 101,
                           scaleX = FALSE, scaleX.fct, scaleX.inv,
                           scaleY = FALSE, scaleY.fct = pnorm, scaleY.inv=qnorm,
                           scaleN = 9, x.ticks = NULL, y.ticks = NULL,
                           cex.pts = 1, col.pts = par("col"),
                           pch.pts = 19, jitter.fac = 1, with.lab = FALSE,
                           lab.pts = NULL, lab.font = NULL, alpha.trsp = NA,
                           which.lbs = NULL, which.Order  = NULL,
                           return.Order = FALSE){

        mc <- match.call(expand.dots = FALSE)
        dots <- as.list(mc$"...")
        if(!is(IC1,"IC")) stop ("IC1 must be of class 'IC'")
        if(!is(IC2,"IC")) stop ("IC2 must be of class 'IC'")
        if(!identical(IC1@CallL2Fam, IC2@CallL2Fam))
        stop("IC1 and IC2 must be defined on the same model")

        L2Fam <- eval(IC1@CallL2Fam)

        b20 <- NULL
        fCpl <- eval(dots$fromCniperPlot)
        if(!is.null(fCpl))
            if(fCpl) b20 <- neighbor@radius*Risks(IC2)$asBias$value
        dots$fromCniperPlot <- NULL
        
        fun <- .getFunCnip(IC1,IC2, risk, L2Fam, neighbor@radius, b20)

        if(missing(scaleX.fct)){
           scaleX.fct <- p(L2Fam)
           scaleX.inv <- q(L2Fam)
        }

        if(!is.null(as.list(mc)$lower)) lower <- p(L2Fam)(lower)
        if(!is.null(as.list(mc)$upper)) upper <- p(L2Fam)(upper)
        x <-  q(L2Fam)(seq(lower,upper,length=n))
        if(is(distribution(L2Fam), "DiscreteDistribution"))
           x <- seq(q(L2Fam)(lower),q(L2Fam)(upper),length=n)
        resc <- .rescalefct(x, fun, scaleX, scaleX.fct,
                     scaleX.inv, scaleY, scaleY.fct, dots$xlim, dots$ylim, dots)
        dots$x <- resc$X
        dots$y <- resc$Y


        dots$type <- "l"
        if(is.null(dots$main)) dots$main <- gettext("Cniper region plot")
        if(is.null(dots$xlab)) dots$xlab <- gettext("Dirac point")
        if(is.null(dots$ylab))
           dots$ylab <- gettext("Asymptotic Risk difference (IC1 - IC2)")

        colSet <- ltySet <- lwdSet <- FALSE
        if(!is.null(dots$col)) {colSet <- TRUE; colo <- eval(dots$col)}
        if(colSet) {
           colo <- rep(colo,length.out=2)
           dots$col <- colo[1]
        }
        if(!is.null(dots$lwd)) {lwdSet <- TRUE; lwdo <- eval(dots$lwd)}
        if(lwdSet) {
           lwdo <- rep(lwdo,length.out=2)
           dots$lwd <- lwdo[1]
        }
        if(!is.null(dots$lty)) {ltySet <- TRUE; ltyo <- eval(dots$lty)}
        if(ltySet && ((!is.numeric(ltyo) && length(ltyo)==1)||
                        is.numeric(ltyo))){
           ltyo <- list(ltyo,ltyo)
           dots$lty <- ltyo[[1]]
        }else{ if (ltySet && !is.numeric(ltyo) && length(ltyo)==2){
                   dots$lty <- ltyo[[1]]
            }
        }
        do.call(plot,dots)

        dots <- .makedotsLowLevel(dots)
        dots$x <- dots$y <- NULL
        if(colSet) dots$col <- colo[2]
        if(lwdSet) dots$lwd <- lwdo[2]
        if(ltySet) dots$lty <- ltyo[[2]]

        dots$h <- if(scaleY) scaleY.fct(0) else 0
        do.call(abline, dots)

        .plotRescaledAxis(scaleX, scaleX.fct, scaleX.inv, scaleY,scaleY.fct,
                          scaleY.inv, dots$xlim, dots$ylim, resc$X, ypts = 400,
                          n = scaleN, x.ticks = x.ticks, y.ticks = y.ticks)
        if(!is.null(data))
           return(.plotData(data, dots, mc, fun, L2Fam, IC1))
        invisible(NULL)
}

cniperPoint <- function(L2Fam, neighbor, risk= asMSE(),
                        lower=getdistrOption("DistrResolution"),
                        upper=1-getdistrOption("DistrResolution")){


        mc <- match.call(expand.dots = FALSE)

        if(!is.null(as.list(mc)$lower)) lower <- p(L2Fam)(lower)
        if(!is.null(as.list(mc)$upper)) upper <- p(L2Fam)(upper)
        lower <- q(L2Fam)(lower)
        upper <- q(L2Fam)(upper)

        robMod <- InfRobModel(center = L2Fam, neighbor = neighbor)

        psi <- optIC(model = L2Fam, risk = asCov())
        eta <- optIC(model = robMod, risk = risk)

        fun <- .getFunCnip(psi,eta, risk, L2Fam, neighbor@radius)

        res <- uniroot(fun, lower = lower, upper = upper)$root
        names(res) <- "cniper point"
        res
}

cniperPointPlot <- function(L2Fam, data=NULL, ..., neighbor, risk= asMSE(),
                        lower=getdistrOption("DistrResolution"),
                        upper=1-getdistrOption("DistrResolution"), n = 101,
                        withMaxRisk = TRUE,
                           scaleX = FALSE, scaleX.fct, scaleX.inv,
                           scaleY = FALSE, scaleY.fct = pnorm, scaleY.inv=qnorm,
                           scaleN = 9, x.ticks = NULL, y.ticks = NULL,
                           cex.pts = 1, col.pts = par("col"),
                           pch.pts = 19, jitter.fac = 1, with.lab = FALSE,
                           lab.pts = NULL, lab.font = NULL, alpha.trsp = NA,
                           which.lbs = NULL, which.Order  = NULL,
                           return.Order = FALSE){

        mc <- match.call(#call = sys.call(sys.parent(1)),
                       expand.dots = FALSE)
        mcl <- as.list(mc[-1])
        dots <- as.list(mc$"...")

        robMod <- InfRobModel(center = L2Fam, neighbor = neighbor)

        mcl$IC1 <- optIC(model = L2Fam, risk = asCov())
        mcl$IC2 <- optIC(model = robMod, risk = risk)
        mcl$L2Fam <- NULL
        if(is.null(dots$ylab))
           mcl$ylab <- gettext("Asymptotic Risk difference (classic - robust)")
        if(is.null(dots$main))
           mcl$main <- gettext("Cniper point plot")

        if(withMaxRisk) mcl$fromCniperPlot <- TRUE
        do.call(cniperCont, mcl)
}



