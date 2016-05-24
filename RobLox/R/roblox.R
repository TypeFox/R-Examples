###############################################################################
## computation of radius-minimax IC
## using predefined functions included in "sysdata.rda"
###############################################################################
.getlsInterval <- function(r, rlo, rup){
    if(r > 10){
        b <- 1.618128043
        const <- 1.263094656
        A2 <- b^2*(1+r^2)/(1+const)
        A1 <- const*A2
    }else{
        A1 <- .getA1.locsc(r)
        A2 <- .getA2.locsc(r)
        b <- .getb.locsc(r)
    }

    if(rlo == 0){
        efflo <- (A1 + A2 - b^2*r^2)/1.5
    }else{
        A1lo <- .getA1.locsc(rlo)
        A2lo <- .getA2.locsc(rlo)
        efflo <- (A1 + A2 - b^2*(r^2 - rlo^2))/(A1lo + A2lo)
    }

    if(rup > 10){
        bup <- 1.618128043
        const.up <- 1.263094656
        A2up <- bup^2*(1+rup^2)/(1+const.up)
        A1up <- const.up*A2up
        effup <- (A1 + A2 - b^2*(r^2 - rup^2))/(A1up + A2up)
    }else{
        A1up <- .getA1.locsc(rup)
        A2up <- .getA2.locsc(rup)
        effup <- (A1 + A2 - b^2*(r^2 - rup^2))/(A1up + A2up)
    }

    return(effup-efflo)
}
.getlInterval <- function(r, rlo, rup){
    if(r > 10){
        b <- sqrt(pi/2)
        A <- b^2*(1+r^2)
    }else{
        A <- .getA.loc(r)
        b <- .getb.loc(r)
    }

    if(rlo == 0){
        efflo <- A - b^2*r^2
    }else{
        Alo <- .getA.loc(rlo)
        efflo <- (A - b^2*(r^2 - rlo^2))/Alo
    }

    if(rup > 10){
        Aup <- pi/2*(1+rup^2)
        effup <- (A - b^2*(r^2 - rup^2))/Aup
    }else{
        Aup <- .getA.loc(rup)
        effup <- (A - b^2*(r^2 - rup^2))/Aup
    }

    return(effup-efflo)
}
.getsInterval <- function(r, rlo, rup){
    if(r > 10){
        b <- 1/(4*qnorm(0.75)*dnorm(qnorm(0.75)))
        A <- b^2*(1+r^2)
    }else{
        A <- .getA.sc(r)
        b <- .getb.sc(r)
    }

    if(rlo == 0){
        efflo <- (A - b^2*r^2)/0.5
    }else{
        Alo <- .getA.sc(rlo)
        efflo <- (A - b^2*(r^2 - rlo^2))/Alo
    }

    if(rup > 10){
        bup <- 1/(4*qnorm(0.75)*dnorm(qnorm(0.75)))
        Aup <- bup^2*(1+rup^2)
        effup <- (A - b^2*(r^2 - rup^2))/Aup
    }else{
        Aup <- .getA.sc(rup)
        effup <- (A - b^2*(r^2 - rup^2))/Aup
    }

    return(effup-efflo)
}


###############################################################################
## computation of k-step construction
###############################################################################
.onestep.loc <- function(x, initial.est, A, b, sd){
    u <- A*(x-initial.est)/sd^2
    IC <- mean(u*pmin(1, b/abs(u)), na.rm = TRUE)
    return(initial.est + IC)
}
.kstep.loc <- function(x, initial.est, A, b, sd, k){
    est <- initial.est
    for(i in 1:k){
        est <- .onestep.loc(x = x, initial.est = est, A = A, b = b, sd = sd)
    }
    return(est)
}
.onestep.sc <- function(x, initial.est, A, a, b, mean){
    v <- A*(((x-mean)/initial.est)^2-1)/initial.est - a
    IC <- mean(v*pmin(1, b/abs(v)), na.rm = TRUE)
    return(initial.est + IC)
}
.kstep.sc <- function(x, initial.est, A, a, b, mean, k){
    est <- .onestep.sc(x = x, initial.est = initial.est, A = A, a = a, b = b, mean = mean)
    if(k > 1){
        for(i in 2:k){
            A <- est^2*A/initial.est^2
            a <- est*a/initial.est
            b <- est*b/initial.est
            initial.est <- est
            est <- .onestep.sc(x = x, initial.est = est, A = A, a = a, b = b, mean = mean)
        }
    }
    A <- est^2*A/initial.est^2
    a <- est*a/initial.est
    b <- est*b/initial.est

    return(list(est = est, A = A, a = a, b = b))
}
.onestep.locsc <- function(x, initial.est, A1, A2, a, b){
    mean <- initial.est[1]
    sd <- initial.est[2]
    u <- A1*(x-mean)/sd^2
    v <- A2*(((x-mean)/sd)^2-1)/sd - a[2]
    w <- pmin(1, b/sqrt(u^2 + v^2))
    IC <- c(mean(u*w, na.rm = TRUE), mean(v*w, na.rm = TRUE))
    return(initial.est + IC)
}
.kstep.locsc <- function(x, initial.est, A1, A2, a, b, mean, k){
    est <- .onestep.locsc(x = x, initial.est = initial.est, A1 = A1, A2 = A2, a = a, b = b)
    if(k > 1){
        for(i in 2:k){
            A1 <- est[2]^2*A1/initial.est[2]^2
            A2 <- est[2]^2*A2/initial.est[2]^2
            a <- est[2]*a/initial.est[2]
            b <- est[2]*b/initial.est[2]
            initial.est <- est
            est <- .onestep.locsc(x = x, initial.est = est,
                                                A1 = A1, A2 = A2, a = a, b = b)
        }
    }
    A1 <- est[2]^2*A1/initial.est[2]^2
    A2 <- est[2]^2*A2/initial.est[2]^2
    a <- est[2]*a/initial.est[2]
    b <- est[2]*b/initial.est[2]
    a1 <- A1/est[2]^2
    a3 <- A2/est[2]^2
    a2 <- a[2]/est[2]/a3 + 1
    asVar <- est[2]^2*.ALrlsVar(b = b/est[2], a1 = a1, a2 = a2, a3 = a3)

    return(list(est = est, A1 = A1, A2 = A2, a = a, b = b, asvar = asVar))
}


###############################################################################
## optimally robust estimator for normal location and/or scale
###############################################################################
roblox <- function(x, mean, sd, eps, eps.lower, eps.upper, initial.est, k = 1L, 
                   fsCor = TRUE, returnIC = FALSE, mad0 = 1e-4, na.rm = TRUE){
    es.call <- match.call()
    if(missing(x))
        stop("'x' is missing with no default")
    if(!is.numeric(x)){
        if(is.data.frame(x))
            x <- data.matrix(x)
        else
            x <- as.matrix(x)
        if(!is.matrix(x))
            stop("'x' has to be a numeric vector resp. a matrix or data.frame
                  with one row resp. column/(numeric) variable")
        if(ncol(x) > 1 & nrow(x) > 1)
            stop("number of rows and columns/variables > 1. Please, do use 'rowRoblox'
                  resp. 'colRoblox'.")
    }

    completecases <- complete.cases(x)
    if(na.rm) x <- na.omit(x)

    if(length(x) <= 2){
        if(missing(mean) && missing(sd)){
            warning("Sample size <= 2! => Median and MAD are used for estimation.")
            robEst <- c(median(x, na.rm = TRUE), mad(x, na.rm = TRUE))
            names(robEst) <- c("mean", "sd")
            Info.matrix <- matrix(c("roblox", 
                                  paste("median and MAD")),
                                  ncol = 2, dimnames = list(NULL, c("method", "message")))
            return(new("ALEstimate", name = "Median and MAD", 
                       completecases = completecases,
                       estimate.call = es.call, estimate = robEst, 
                       samplesize = length(x), asvar = NULL,
                       asbias = NULL, pIC = NULL, Infos = Info.matrix))
        }
        if(missing(mean)){
            warning("Sample size <= 2! => Median is used for estimation.")
            robEst <- median(x, na.rm = TRUE)
            names(robEst) <- "mean"
            Info.matrix <- matrix(c("roblox", 
                                  paste("median")),
                                  ncol = 2, dimnames = list(NULL, c("method", "message")))
            return(new("ALEstimate", name = "Median", 
                       completecases = completecases,
                       estimate.call = es.call, estimate = robEst,
                       samplesize = length(x), asvar = NULL,
                       asbias = NULL, pIC = NULL, Infos = Info.matrix))
        }
        if(missing(sd)){
            warning("Sample size <= 2! => MAD is used for estimation.")
            if(length(mean) != 1)
                stop("mean has length != 1")
            robEst <- mad(x, center = mean, na.rm = TRUE)
            names(robEst) <- "sd"
            Info.matrix <- matrix(c("roblox", 
                                  paste("MAD")),
                                  ncol = 2, dimnames = list(NULL, c("method", "message")))
            return(new("ALEstimate", name = "MAD", 
                       completecases = completecases,
                       estimate.call = es.call, estimate = robEst,
                       samplesize = length(x), asvar = NULL,
                       asbias = NULL, pIC = NULL, Infos = Info.matrix))
        }
    }
    if(missing(eps) && missing(eps.lower) && missing(eps.upper)){
        eps.lower <- 0
        eps.upper <- 0.5
    }
    if(missing(eps)){
        if(!missing(eps.lower) && missing(eps.upper))
            eps.upper <- 0.5
        if(missing(eps.lower) && !missing(eps.upper))
            eps.lower <- 0
        if(length(eps.lower) != 1 || length(eps.upper) != 1)
            stop("'eps.lower' and 'eps.upper' have to be of length 1")
        if(!is.numeric(eps.lower) || !is.numeric(eps.upper) || eps.lower >= eps.upper) 
            stop("'eps.lower' < 'eps.upper' is not fulfilled")
        if((eps.lower < 0) || (eps.upper > 0.5))
            stop("'eps.lower' and 'eps.upper' have to be in [0, 0.5]")
    }else{
        if(length(eps) != 1)
            stop("'eps' has to be of length 1")
        if((eps < 0) || (eps > 0.5))
            stop("'eps' has to be in (0, 0.5]")
        if(eps == 0){
            if(missing(mean) && missing(sd)){
                warning("eps = 0! => Mean and sd are used for estimation.")
                n <- sum(!is.na(x))
                robEst <- c(mean(x, na.rm = TRUE), sqrt((n-1)/n)*sd(x, na.rm = TRUE))
                names(robEst) <- c("mean", "sd")
                Info.matrix <- matrix(c("roblox", 
                                      paste("mean and sd")),
                                      ncol = 2, dimnames = list(NULL, c("method", "message")))
                return(new("ALEstimate", name = "Mean and sd", 
                          completecases = completecases,
                          estimate.call = es.call, estimate = robEst,
                          samplesize = n, asvar = NULL,
                          asbias = NULL, pIC = NULL, Infos = Info.matrix))
            }
            if(missing(mean)){
                warning("eps = 0! => Mean is used for estimation.")
                robEst <- mean(x, na.rm = TRUE)
                names(robEst) <- "mean"
                Info.matrix <- matrix(c("roblox", 
                                      paste("mean")),
                                      ncol = 2, dimnames = list(NULL, c("method", "message")))
                return(new("ALEstimate", name = "Mean", 
                          completecases = completecases,
                          estimate.call = es.call, estimate = robEst,
                          samplesize = length(x), asvar = NULL,
                          asbias = NULL, pIC = NULL, Infos = Info.matrix))
            }
            if(missing(sd)){
                warning("eps = 0! => sd is used for estimation.")
                n <- sum(!is.na(x))
                robEst <- sqrt((n-1)/n)*sd(x, na.rm = TRUE)
                names(robEst) <- "sd"
                Info.matrix <- matrix(c("roblox", 
                                      paste("sd")),
                                      ncol = 2, dimnames = list(NULL, c("method", "message")))
                return(new("ALEstimate", name = "sd", 
                          completecases = completecases,
                          estimate.call = es.call, estimate = robEst,
                          samplesize = n, asvar = NULL,
                          asbias = NULL, pIC = NULL, Infos = Info.matrix))
            }
        }
    }
    if(!is.integer(k))
        k <- as.integer(k)
    if(k < 1){
        stop("'k' has to be some positive integer value")
    }
    if(length(k) != 1){
        stop("'k' has to be of length 1")
    }
    if(missing(mean) && missing(sd)){
        if(missing(initial.est)){
            mean <- median(x, na.rm = TRUE)
            sd <- mad(x, center = mean, na.rm = TRUE)
            if(sd == 0){
                warning("'mad(x, na.rm = TRUE) = 0' => cannot compute a valid initial estimate. 
                        To avoid division by zero 'mad0' is used. You could also specify 
                        a valid scale estimate via 'initial.est'. Note that you have to provide
                        a location and scale estimate.")
                sd <- mad0
            }
        }else{
            if(!is.numeric(initial.est) || length(initial.est) != 2)
                stop("'initial.est' needs to be a numeric vector of length 2 or missing")
            mean <- initial.est[1]
            sd <- initial.est[2]
            if(sd <= 0)
                stop("initial estimate for scale <= 0 which is no valid scale estimate")
        }
        mean.sd <- matrix(c(mean, sd),nrow=1,ncol=2)
        colnames(mean.sd) <- c("mean","sd")
        if(!missing(eps)){
            r <- sqrt(length(x))*eps
            if(fsCor) r <- finiteSampleCorrection(r = r, n = length(x), model = "locsc")
            if(r > 10){
                b <- sd*1.618128043
                const <- 1.263094656
                A2 <- b^2*(1+r^2)/(1+const)
                A1 <- const*A2
                a <- c(0, -0.6277527697*A2/sd)
                mse <- A1 + A2
            }else{
                A1 <- sd^2*.getA1.locsc(r)
                A2 <- sd^2*.getA2.locsc(r)
                a <- sd*c(0, .geta.locsc(r))
                b <- sd*.getb.locsc(r)
                mse <- A1 + A2
            }
            robEst <- .kstep.locsc(x = x, initial.est = c(mean, sd), A1 = A1, A2 = A2, a = a, b = b, k = k)
            names(robEst$est) <- c("mean", "sd")
            if(fsCor){
                Info.matrix <- matrix(c("roblox", 
                                        paste("finite-sample corrected optimally robust estimate for contamination 'eps' =", 
                                              round(eps, 3), "and 'asMSE'")),
                                      ncol = 2, dimnames = list(NULL, c("method", "message")))
            }else{
                Info.matrix <- matrix(c("roblox", 
                                        paste("optimally robust estimate for contamination 'eps' =", round(eps, 3),
                                              "and 'asMSE'")),
                                      ncol = 2, dimnames = list(NULL, c("method", "message")))
            }
            if(returnIC){
                w <- new("HampelWeight")
                clip(w) <- robEst$b
                cent(w) <- robEst$a/robEst$A2
                stand(w) <- diag(c(robEst$A1, robEst$A2))
                weight(w) <- getweight(w, neighbor = ContNeighborhood(radius = r), 
                                       biastype = symmetricBias(), 
                                       normW = NormType())
                mse <- robEst$A1 + robEst$A2
                modIC <- function(L2Fam, IC){
                    ICL2Fam <- eval(CallL2Fam(IC))
                    if(is(L2Fam, "L2LocationScaleFamily") && is(distribution(L2Fam), "Norm")){
                        sdneu <- main(L2Fam)[2]
                        sdalt <- main(ICL2Fam)[2]
                        r <- neighborRadius(IC)
                        w <- weight(IC)
                        clip(w) <- sdneu*clip(w)/sdalt
                        cent(w) <- sdalt*cent(w)/sdneu
                        stand(w) <- sdneu^2*stand(w)/sdalt^2
                        weight(w) <- getweight(w, neighbor = ContNeighborhood(radius = r), 
                                      biastype = biastype(IC), 
                                      normW = normtype(IC))
                        A <- sdneu^2*stand(IC)/sdalt^2
                        b <- sdneu*clip(IC)/sdalt
                        a <- sdneu*cent(IC)/sdalt
                        mse <- sum(diag(A))
                        a1 <- A[1, 1]/sdneu^2
                        a3 <- A[2, 2]/sdneu^2
                        a2 <- a[2]/sdneu/a3 + 1
                        asVar <- sdneu^2*.ALrlsVar(b = b/sdneu, a1 = a1, a2 = a2, a3 = a3)
                        res <- list(A = A, a = sdneu*cent(IC)/sdalt, b = b, d = NULL,
                                    risk = list(asMSE = mse, asBias = b, 
                                                trAsCov = mse - r^2*b^2,
                                                asCov = asVar), info = Infos(IC), w = w,
                                    normtype = normtype(IC), biastype = biastype(IC),
                                    modifyIC = modifyIC(IC))
                        IC <- generateIC(neighbor = ContNeighborhood(radius = r),
                                        L2Fam = L2Fam, res = res)
                        addInfo(IC) <- c("modifyIC", "The IC has been modified")
                        addInfo(IC) <- c("modifyIC", "The entries in 'Infos' may be wrong")
                        return(IC)
                    }else{
                        makeIC(L2Fam, IC)
                    }
                }
                L2Fam <- substitute(NormLocationScaleFamily(mean = m1, sd = s1), 
                                    list(m1 = robEst$est[1], s1 = robEst$est[2]))
                info <- c("roblox", "optimally robust IC for AL estimators and 'asMSE'")
                IC1 <- generateIC(neighbor = ContNeighborhood(radius = r), 
                                  L2Fam = eval(L2Fam), 
                                  res = list(A = diag(c(robEst$A1, robEst$A2)), a = robEst$a, 
                                      b = robEst$b, d = NULL, 
                                      risk = list(asMSE = mse, asBias = robEst$b, 
                                                  trAsCov = mse - r^2*robEst$b^2,
                                                  asCov = robEst$asVar), 
                                      info = info, w = w, biastype = symmetricBias(), 
                                      normtype = NormType(), modifyIC = modIC))
                Infos(IC1) <- Info.matrix
                return(new("kStepEstimate", name = "Optimally robust estimate", 
                           completecases = completecases,
                           estimate.call = es.call, estimate = robEst$est,
                           samplesize = length(x), asvar = robEst$asvar,
                           asbias = r*robEst$b, steps = k, pIC = IC1, Infos = Info.matrix,
                           start = mean.sd, startval = mean.sd, ustartval = mean.sd))
            }else
                return(new("kStepEstimate", name = "Optimally robust estimate", 
                           completecases = completecases,
                           estimate.call = es.call, estimate = robEst$est,
                           samplesize = length(x), asvar = robEst$asvar,
                           asbias = r*robEst$b, steps = k, pIC = NULL, Infos = Info.matrix,
                           start = mean.sd, startval = mean.sd, ustartval = mean.sd))
        }else{
            sqrtn <- sqrt(length(x))
            rlo <- sqrtn*eps.lower
            rup <- sqrtn*eps.upper
            if(rlo > 10){
                r <- (rlo + rup)/2
            }else{
                r <- uniroot(.getlsInterval, lower = rlo+1e-8, upper = rup, 
                             tol = .Machine$double.eps^0.25, rlo = rlo, rup = rup)$root
            }
            if(fsCor){
                r.as <- r
                r <- finiteSampleCorrection(r = r, n = length(x), model = "locsc")
            }
            if(r > 10){
                b <- sd*1.618128043
                const <- 1.263094656
                A2 <- b^2*(1+r^2)/(1+const)
                A1 <- const*A2
                a <- c(0, -0.6277527697*A2/sd)
                mse <- A1 + A2
            }else{
                A1 <- sd^2*.getA1.locsc(r)
                A2 <- sd^2*.getA2.locsc(r)
                a <- sd*c(0, .geta.locsc(r))
                b <- sd*.getb.locsc(r)
                mse <- A1 + A2
            }
            if(rlo == 0){
                ineff <- (A1 + A2 - b^2*r.as^2)/(1.5*sd^2)
            }else{
                if(rlo > 10){
                    ineff <- 1
                }else{
                    A1lo <- sd^2*.getA1.locsc(rlo)
                    A2lo <- sd^2*.getA2.locsc(rlo)
                    ineff <- (A1 + A2 - b^2*(r.as^2 - rlo^2))/(A1lo + A2lo)
                }
            }
            robEst <- .kstep.locsc(x = x, initial.est = c(mean, sd), A1 = A1, A2 = A2, a = a, b = b, k = k)
            names(robEst$est) <- c("mean", "sd")
            if(fsCor){
                Info.matrix <- matrix(c(rep("roblox", 3), 
                                      paste("finite-sample corrected radius-minimax estimate for contamination interval [", 
                                        round(eps.lower, 3), ", ", round(eps.upper, 3), "]", sep = ""),
                                      paste("least favorable (uncorrected) contamination: ", round(r.as/sqrtn, 3), sep = ""),
                                      paste("maximum asymptotic MSE-inefficiency: ", round(ineff, 3), sep = "")), 
                                      ncol = 2, dimnames = list(NULL, c("method", "message")))
            }else{
                Info.matrix <- matrix(c(rep("roblox", 3), 
                                      paste("radius-minimax estimate for contamination interval [", 
                                        round(eps.lower, 3), ", ", round(eps.upper, 3), "]", sep = ""),
                                      paste("least favorable contamination: ", round(r/sqrtn, 3), sep = ""),
                                      paste("maximum asymptotic MSE-inefficiency: ", round(ineff, 3), sep = "")), 
                                      ncol = 2, dimnames = list(NULL, c("method", "message")))
            }
            if(returnIC){
                w <- new("HampelWeight")
                clip(w) <- robEst$b
                cent(w) <- robEst$a/robEst$A2
                stand(w) <- diag(c(robEst$A1, robEst$A2))
                weight(w) <- getweight(w, neighbor = ContNeighborhood(radius = r), 
                                       biastype = symmetricBias(), 
                                       normW = NormType())
                mse <- robEst$A1 + robEst$A2
                modIC <- function(L2Fam, IC){
                    ICL2Fam <- eval(CallL2Fam(IC))
                    if(is(L2Fam, "L2LocationScaleFamily") && is(distribution(L2Fam), "Norm")){
                        sdneu <- main(L2Fam)[2]
                        sdalt <- main(ICL2Fam)[2]
                        r <- neighborRadius(IC)
                        w <- weight(IC)
                        clip(w) <- sdneu*clip(w)/sdalt
                        cent(w) <- sdalt*cent(w)/sdneu
                        stand(w) <- sdneu^2*stand(w)/sdalt^2
                        weight(w) <- getweight(w, neighbor = ContNeighborhood(radius = r), 
                                      biastype = biastype(IC), 
                                      normW = normtype(IC))
                        A <- sdneu^2*stand(IC)/sdalt^2
                        b <- sdneu*clip(IC)/sdalt
                        a <- sdneu*cent(IC)/sdalt
                        mse <- sum(diag(A))
                        a1 <- A[1, 1]/sdneu^2
                        a3 <- A[2, 2]/sdneu^2
                        a2 <- a[2]/sdneu/a3 + 1
                        asVar <- sdneu^2*.ALrlsVar(b = b/sdneu, a1 = a1, a2 = a2, a3 = a3)
                        res <- list(A = A, a = sdneu*cent(IC)/sdalt, b = b, d = NULL,
                                    risk = list(asMSE = mse, asBias = b, 
                                                trAsCov = mse - r^2*b^2,
                                                asCov = asVar), info = Infos(IC), w = w,
                                    normtype = normtype(IC), biastype = biastype(IC),
                                    modifyIC = modifyIC(IC))
                        IC <- generateIC(neighbor = ContNeighborhood(radius = r),
                                        L2Fam = L2Fam, res = res)
                        addInfo(IC) <- c("modifyIC", "The IC has been modified")
                        addInfo(IC) <- c("modifyIC", "The entries in 'Infos' may be wrong")
                        return(IC)
                    }else{
                        makeIC(L2Fam, IC)
                    }
                }
                L2Fam <- substitute(NormLocationScaleFamily(mean = m1, sd = s1), 
                                    list(m1 = robEst$est[1], s1 = robEst$est[2]))
                info <- c("roblox", "optimally robust IC for AL estimators and 'asMSE'")
                IC1 <- generateIC(neighbor = ContNeighborhood(radius = r), 
                                  L2Fam = eval(L2Fam), 
                                  res = list(A = diag(c(robEst$A1, robEst$A2)), a = robEst$a, 
                                      b = robEst$b, d = NULL, 
                                      risk = list(asMSE = mse, asBias = robEst$b, 
                                                  trAsCov = mse - r^2*robEst$b^2,
                                                  asCov = robEst$asvar), 
                                      info = info, w = w, biastype = symmetricBias(), 
                                      normtype = NormType(), modifyIC = modIC))
                Infos(IC1) <- Info.matrix
                return(new("kStepEstimate", name = "Optimally robust estimate", 
                           completecases = completecases,
                           estimate.call = es.call, estimate = robEst$est,
                           samplesize = length(x), asvar = robEst$asvar,
                           asbias = r*robEst$b, steps = k, pIC = IC1, Infos = Info.matrix,
                           start = mean.sd, startval = mean.sd, ustartval = mean.sd))
            }else
                return(new("kStepEstimate", name = "Optimally robust estimate", 
                           completecases = completecases,
                           estimate.call = es.call, estimate = robEst$est,
                           samplesize = length(x), asvar = robEst$asvar,
                           asbias = r*robEst$b, steps = k, pIC = NULL, Infos = Info.matrix,
                           start = mean.sd, startval = mean.sd, ustartval = mean.sd))
        }
    }else{
        if(missing(mean)){
            if(length(sd) != 1)
                stop("'sd' has length != 1")
            if(sd <= 0)
                stop("'sd' has to be positive")
            if(missing(initial.est)){
                mean <- median(x, na.rm = TRUE)
            }else{
                if(!is.numeric(initial.est) || length(initial.est) != 1)
                    stop("'initial.est' needs to be a numeric vector of length 1 or missing")
                mean <- initial.est
            }

            if(!missing(eps)){
                r <- sqrt(length(x))*eps
                if(fsCor) r <- finiteSampleCorrection(r = r, n = length(x), model = "loc")
                if(r > 10){
                    b <- sd*sqrt(pi/2)
                    A <- b^2*(1+r^2)
                }else{
                    A <- sd^2*.getA.loc(r)
                    b <- sd*.getb.loc(r)
                }
                robEst <- .kstep.loc(x = x, initial.est = mean, A = A, b = b, sd = sd, k = k)
                names(robEst) <- "mean"
                if(fsCor){
                    Info.matrix <- matrix(c("roblox", 
                                            paste("finite-sample corrected optimally robust estimate for contamination 'eps' =", round(eps, 3),
                                                  "and 'asMSE'")),
                                          ncol = 2, dimnames = list(NULL, c("method", "message")))
                }else{
                    Info.matrix <- matrix(c("roblox", 
                                            paste("optimally robust estimate for contamination 'eps' =", round(eps, 3),
                                                  "and 'asMSE'")),
                                          ncol = 2, dimnames = list(NULL, c("method", "message")))
                }
                if(returnIC){
                    w <- new("HampelWeight")
                    clip(w) <- b
                    cent(w) <- 0
                    stand(w) <- as.matrix(A)
                    weight(w) <- getweight(w, neighbor = ContNeighborhood(radius = r), 
                                           biastype = symmetricBias(), 
                                           normW = NormType())
                    modIC <- function(L2Fam, IC){
                        if(is(L2Fam, "L2LocationFamily") && is(distribution(L2Fam), "Norm")){
                            CallL2New <- call("NormLocationFamily", 
                                              mean = main(L2Fam))
                            CallL2Fam(IC) <- CallL2New
                            return(IC)
                        }else{
                            makeIC(L2Fam, IC)
                        }
                    }
                    L2Fam <- substitute(NormLocationFamily(mean = m1, sd = s1), 
                                        list(m1 = robEst, s1 = sd))
                    IC1 <- generateIC(neighbor = ContNeighborhood(radius = r), 
                                      L2Fam = eval(L2Fam), 
                                      res = list(A = as.matrix(A), a = 0, b = b, d = NULL, 
                                          risk = list(asMSE = A, asBias = b, asCov = A-r^2*b^2), 
                                          info = c("roblox", "optimally robust IC for AL estimators and 'asMSE'"),
                                          w = w, biastype = symmetricBias(), normtype = NormType(),
                                          modifyIC = modIC))
                    Infos(IC1) <- Info.matrix
                    return(new("kStepEstimate", name = "Optimally robust estimate",
                               completecases = completecases,
                               estimate.call = es.call, estimate = robEst,
                               samplesize = length(x), asvar = as.matrix(A-r^2*b^2),
                               asbias = r*b, steps = k, pIC = IC1, Infos = Info.matrix,
                           start = median, startval = matrix(mean,1,1), ustartval = matrix(mean,1,1)))
                }else
                    return(new("kStepEstimate", name = "Optimally robust estimate",
                               completecases = completecases,
                               estimate.call = es.call, estimate = robEst,
                               samplesize = length(x), asvar = as.matrix(A-r^2*b^2),
                               asbias = r*b, steps = k, pIC = NULL, Infos = Info.matrix,
                           start = median, startval = matrix(mean,1,1), ustartval = matrix(mean,1,1)))
            }else{
                sqrtn <- sqrt(length(x))
                rlo <- sqrtn*eps.lower
                rup <- sqrtn*eps.upper
                if(rlo > 10){ 
                    r <- (rlo+rup)/2
                }else{
                    r <- uniroot(.getlInterval, lower = rlo+1e-8, upper = rup, 
                                 tol = .Machine$double.eps^0.25, rlo = rlo, rup = rup)$root
                }
                if(fsCor){ 
                    r.as <- r
                    r <- finiteSampleCorrection(r = r, n = length(x), model = "loc")
                }
                if(r > 10){
                    b <- sd*sqrt(pi/2)
                    A <- b^2*(1+r^2)
                }else{
                    A <- sd^2*.getA.loc(r)
                    b <- sd*.getb.loc(r)
                }
                if(rlo == 0){
                    ineff <- (A - b^2*r^2)/sd^2
                }else{
                    if(rlo > 10){
                        ineff <- 1
                    }else{
                        Alo <- sd^2*.getA.loc(rlo)
                        ineff <- (A - b^2*(r^2 - rlo^2))/Alo
                    }
                }
                robEst <- .kstep.loc(x = x, initial.est = mean, A = A, b = b, sd = sd, k = k)
                names(robEst) <- "mean"
                if(fsCor){
                    Info.matrix <- matrix(c(rep("roblox", 3), 
                                          paste("finite-sample corrected radius-minimax estimate for contamination interval [", 
                                            round(eps.lower, 3), ", ", round(eps.upper, 3), "]", sep = ""),
                                          paste("least favorable (uncorrected) contamination: ", round(r.as/sqrtn, 3), sep = ""),
                                          paste("maximum asymptotic MSE-inefficiency: ", round(ineff, 3), sep = "")), 
                                          ncol = 2, dimnames = list(NULL, c("method", "message")))
                }else{
                    Info.matrix <- matrix(c(rep("roblox", 3), 
                                          paste("radius-minimax estimate for contamination interval [", 
                                            round(eps.lower, 3), ", ", round(eps.upper, 3), "]", sep = ""),
                                          paste("least favorable contamination: ", round(r/sqrtn, 3), sep = ""),
                                          paste("maximum MSE-inefficiency: ", round(ineff, 3), sep = "")), 
                                          ncol = 2, dimnames = list(NULL, c("method", "message")))
                }
                if(returnIC){
                    w <- new("HampelWeight")
                    clip(w) <- b
                    cent(w) <- 0
                    stand(w) <- as.matrix(A)
                    weight(w) <- getweight(w, neighbor = ContNeighborhood(radius = r), 
                                       biastype = symmetricBias(), 
                                       normW = NormType())
                    modIC <- function(L2Fam, IC){
                        if(is(L2Fam, "L2LocationFamily") && is(distribution(L2Fam), "Norm")){
                            CallL2New <- call("NormLocationFamily", 
                                              mean = main(L2Fam))
                            CallL2Fam(IC) <- CallL2New
                            return(IC)
                        }else{
                            makeIC(L2Fam, IC)
                        }
                    }
                    L2Fam <- substitute(NormLocationFamily(mean = m1, sd = s1), 
                                        list(m1 = robEst, s1 = sd))
                    IC1 <- generateIC(neighbor = ContNeighborhood(radius = r), 
                                      L2Fam = eval(L2Fam), 
                                      res = list(A = as.matrix(A), a = 0, b = b, d = NULL, 
                                          risk = list(asMSE = A, asBias = b, asCov = A-r^2*b^2), 
                                          info = c("roblox", "optimally robust IC for AL estimators and 'asMSE'"),
                                          w = w, biastype = symmetricBias(), normtype = NormType(),
                                          modifyIC = modIC))
                    Infos(IC1) <- Info.matrix
                    return(new("kStepEstimate", name = "Optimally robust estimate",
                               completecases = completecases,
                               estimate.call = es.call, estimate = robEst,
                               samplesize = length(x), asvar = as.matrix(A-r^2*b^2),
                               asbias = r*b, steps = k, pIC = IC1, Infos = Info.matrix,
                           start = median, startval = matrix(mean,1,1), ustartval = matrix(mean,1,1)))
                }else
                    return(new("kStepEstimate", name = "Optimally robust estimate",
                               completecases = completecases,
                               estimate.call = es.call, estimate = robEst,
                               samplesize = length(x), asvar = as.matrix(A-r^2*b^2),
                               asbias = r*b, steps = k, pIC = NULL, Infos = Info.matrix,
                           start = median, startval = matrix(mean,1,1), ustartval = matrix(mean,1,1)))
            }
        }
        if(missing(sd)){
            if(length(mean) != 1)
                stop("mean has length != 1")
            if(missing(initial.est)){ 
                sd <- mad(x, center = mean, na.rm = TRUE)
                if(sd == 0){
                    warning("'mad(x, na.rm = TRUE) = 0' => cannot compute a valid initial estimate. 
                            To avoid division by zero 'mad0' is used. You could also specify 
                            a valid scale estimate via 'initial.est'.")
                    sd <- mad0
                }
            }else{
                if(!is.numeric(initial.est) || length(initial.est) != 1)
                    stop("'initial.est' needs to be a numeric vector of length 1 or missing")
                sd <- initial.est
                if(initial.est <= 0)
                  stop("'initial.est <= 0'; i.e., is no valid scale estimate")
            }

            if(!missing(eps)){
                r <- sqrt(length(x))*eps
                if(fsCor) r <- finiteSampleCorrection(r = r, n = length(x), model = "sc")
                if(r > 10){
                    b <- sd/(4*qnorm(0.75)*dnorm(qnorm(0.75)))
                    A <- b^2*(1+r^2)
                    a <- (qnorm(0.75)^2 - 1)/sd*A
                }else{
                    A <- sd^2*.getA.sc(r)
                    a <- sd*.geta.sc(r)
                    b <- sd*.getb.sc(r)
                }
                robEst <- .kstep.sc(x = x, initial.est = sd, A = A, a = a, b = b, mean = mean, k = k)
                names(robEst$est) <- "sd"
                if(fsCor){
                    Info.matrix <- matrix(c("roblox", 
                                            paste("finite-sample corrected optimally robust estimate for contamination 'eps' =", round(eps, 3),
                                                  "and 'asMSE'")),
                                          ncol = 2, dimnames = list(NULL, c("method", "message")))
                }else{
                    Info.matrix <- matrix(c("roblox", 
                                            paste("optimally robust estimate for contamination 'eps' =", round(eps, 3),
                                                  "and 'asMSE'")),
                                          ncol = 2, dimnames = list(NULL, c("method", "message")))
                }
                if(returnIC){
                    w <- new("HampelWeight")
                    clip(w) <- robEst$b
                    cent(w) <- robEst$a/robEst$A
                    stand(w) <- as.matrix(robEst$A)
                    weight(w) <- getweight(w, neighbor = ContNeighborhood(radius = r), 
                                       biastype = symmetricBias(), 
                                       normW = NormType())
                    modIC <- function(L2Fam, IC){
                        ICL2Fam <- eval(CallL2Fam(IC))
                        if(is(L2Fam, "L2ScaleFamily") && is(distribution(L2Fam), "Norm")){
                            sdneu <- main(L2Fam)
                            sdalt <- main(ICL2Fam)
                            r <- neighborRadius(IC)
                            w <- weight(IC)
                            clip(w) <- sdneu*clip(w)/sdalt
                            cent(w) <- sdalt*cent(w)/sdneu
                            stand(w) <- sdneu^2*stand(w)/sdalt^2
                            weight(w) <- getweight(w, neighbor = ContNeighborhood(radius = r), 
                                          biastype = biastype(IC), 
                                          normW = normtype(IC))
                            A <- sdneu^2*stand(IC)/sdalt^2
                            b <- sdneu*clip(IC)/sdalt
                            res <- list(A = A, a = sdneu*cent(IC)/sdalt, b = b, d = NULL,
                                        risk = list(asMSE = A, asBias = b, asCov = A-r^2*b^2), 
                                        info = Infos(IC), w = w,
                                        normtype = normtype(IC), biastype = biastype(IC),
                                        modifyIC = modifyIC(IC))
                            IC <- generateIC(neighbor = ContNeighborhood(radius = r),
                                            L2Fam = L2Fam, res = res)
                            addInfo(IC) <- c("modifyIC", "The IC has been modified")
                            addInfo(IC) <- c("modifyIC", "The entries in 'Infos' may be wrong")
                            return(IC)
                        }else{
                            makeIC(L2Fam, IC)
                        }
                    }
                    L2Fam <- substitute(NormScaleFamily(mean = m1, sd = s1), 
                                        list(m1 = mean, s1 = robEst$est))
                    IC1 <- generateIC(neighbor = ContNeighborhood(radius = r), 
                                      L2Fam = eval(L2Fam), 
                                      res = list(A = as.matrix(robEst$A), a = robEst$a, b = robEst$b, d = NULL, 
                                          risk = list(asMSE = robEst$A, asBias = robEst$b, 
                                                      asCov = robEst$A-r^2*robEst$b^2), 
                                          info = c("roblox", "optimally robust IC for AL estimators and 'asMSE'"),
                                          w = w, biastype = symmetricBias(), normtype = NormType(),
                                          modifyIC = modIC))
                    Infos(IC1) <- Info.matrix
                    return(new("kStepEstimate", name = "Optimally robust estimate",
                               completecases = completecases,
                               estimate.call = es.call, estimate = robEst$est,
                               samplesize = length(x), asvar = as.matrix(robEst$A-r^2*robEst$b^2),
                               asbias = r*robEst$b, steps = k, pIC = IC1, Infos = Info.matrix,
                           start = mad, startval = matrix(sd,1,1), ustartval = matrix(sd,1,1)))
                }else
                    return(new("kStepEstimate", name = "Optimally robust estimate",
                               completecases = completecases,
                               estimate.call = es.call, estimate = robEst$est,
                               samplesize = length(x), asvar = as.matrix(robEst$A-r^2*robEst$b^2),
                               asbias = r*robEst$b, steps = k, pIC = NULL, Infos = Info.matrix,
                           start = mad, startval = matrix(sd,1,1), ustartval = matrix(sd,1,1)))
            }else{
                sqrtn <- sqrt(length(x))
                rlo <- sqrtn*eps.lower
                rup <- sqrtn*eps.upper
                if(rlo > 10){
                    r <- (rlo+rup)/2
                }else{
                    r <- uniroot(.getsInterval, lower = rlo+1e-8, upper = rup, 
                             tol = .Machine$double.eps^0.25, rlo = rlo, rup = rup)$root
                }
                if(fsCor){ 
                    r.as <- r
                    r <- finiteSampleCorrection(r = r, n = length(x), model = "sc")
                }
                if(r > 10){
                    b <- sd/(4*qnorm(0.75)*dnorm(qnorm(0.75)))
                    A <- b^2*(1+r^2)
                    a <- (qnorm(0.75)^2 - 1)/sd*A
                }else{
                    A <- sd^2*.getA.sc(r)
                    a <- sd*.geta.sc(r)
                    b <- sd*.getb.sc(r)
                }
                if(rlo == 0){
                    ineff <- (A - b^2*r^2)/(0.5*sd^2)
                }else{
                    if(rlo > 10){
                        ineff <- 1
                    }else{
                        Alo <- sd^2*.getA.sc(rlo)
                        ineff <- (A - b^2*(r^2 - rlo^2))/Alo
                    }
                }
                robEst <- .kstep.sc(x = x, initial.est = sd, A = A, a = a, b = b, mean = mean, k = k)
                names(robEst$est) <- "sd"
                if(fsCor){
                    Info.matrix <- matrix(c(rep("roblox", 3), 
                                          paste("finite-sample corrected radius-minimax estimate for contamination interval [", 
                                            round(eps.lower, 3), ", ", round(eps.upper, 3), "]", sep = ""),
                                          paste("least favorable (uncorrected) contamination: ", round(r.as/sqrtn, 3), sep = ""),
                                          paste("maximum asymptotic MSE-inefficiency: ", round(ineff, 3), sep = "")), 
                                          ncol = 2, dimnames = list(NULL, c("method", "message")))
                }else{
                    Info.matrix <- matrix(c(rep("roblox", 3), 
                                          paste("radius-minimax estimate for contamination interval [", 
                                            round(eps.lower, 3), ", ", round(eps.upper, 3), "]", sep = ""),
                                          paste("least favorable contamination: ", round(r/sqrtn, 3), sep = ""),
                                          paste("maximum MSE-inefficiency: ", round(ineff, 3), sep = "")), 
                                          ncol = 2, dimnames = list(NULL, c("method", "message")))
                }
                if(returnIC){
                    w <- new("HampelWeight")
                    clip(w) <- robEst$b
                    cent(w) <- robEst$a/robEst$A
                    stand(w) <- as.matrix(robEst$A)
                    weight(w) <- getweight(w, neighbor = ContNeighborhood(radius = r), 
                                       biastype = symmetricBias(), 
                                       normW = NormType())
                    modIC <- function(L2Fam, IC){
                        ICL2Fam <- eval(CallL2Fam(IC))
                        if(is(L2Fam, "L2ScaleFamily") && is(distribution(L2Fam), "Norm")){
                            sdneu <- main(L2Fam)
                            sdalt <- main(ICL2Fam)
                            r <- neighborRadius(IC)
                            w <- weight(IC)
                            clip(w) <- sdneu*clip(w)/sdalt
                            cent(w) <- sdalt*cent(w)/sdneu
                            stand(w) <- sdneu^2*stand(w)/sdalt^2
                            weight(w) <- getweight(w, neighbor = ContNeighborhood(radius = r), 
                                          biastype = biastype(IC), 
                                          normW = normtype(IC))
                            A <- sdneu^2*stand(IC)/sdalt^2
                            b <- sdneu*clip(IC)/sdalt
                            res <- list(A = A, a = sdneu*cent(IC)/sdalt, b = b, d = NULL,
                                        risk = list(asMSE = A, asBias = b, asCov = A-r^2*b^2), 
                                        info = Infos(IC), w = w,
                                        normtype = normtype(IC), biastype = biastype(IC),
                                        modifyIC = modifyIC(IC))
                            IC <- generateIC(neighbor = ContNeighborhood(radius = r),
                                            L2Fam = L2Fam, res = res)
                            addInfo(IC) <- c("modifyIC", "The IC has been modified")
                            addInfo(IC) <- c("modifyIC", "The entries in 'Infos' may be wrong")
                            return(IC)
                        }else{
                            makeIC(L2Fam, IC)
                        }
                    }
                    L2Fam <- substitute(NormScaleFamily(mean = m1, sd = s1), 
                                        list(m1 = mean, s1 = robEst$est))
                    IC1 <- generateIC(neighbor = ContNeighborhood(radius = r), 
                                      L2Fam = eval(L2Fam), 
                                      res = list(A = as.matrix(robEst$A), a = robEst$a, b = robEst$b, d = NULL, 
                                          risk = list(asMSE = robEst$A, asBias = robEst$b, 
                                                      asCov = robEst$A-r^2*robEst$b^2), 
                                          info = c("roblox", "optimally robust IC for AL estimators and 'asMSE'"),
                                          w = w, biastype = symmetricBias(), normtype = NormType(),
                                          modifyIC = modIC))
                    Infos(IC1) <- Info.matrix
                    return(new("kStepEstimate", name = "Optimally robust estimate",
                               completecases = completecases,
                               estimate.call = es.call, estimate = robEst$est,
                               samplesize = length(x), asvar = as.matrix(robEst$A-r^2*robEst$b^2),
                               asbias = r*robEst$b, steps = k, pIC = IC1, Infos = Info.matrix,
                           start = mad, startval = matrix(sd,1,1), ustartval = matrix(sd,1,1)))
                }else
                    return(new("kStepEstimate", name = "Optimally robust estimate",
                               completecases = completecases,
                               estimate.call = es.call, estimate = robEst$est,
                               samplesize = length(x), asvar = as.matrix(robEst$A-r^2*robEst$b^2),
                               asbias = r*robEst$b, steps = k, pIC = NULL, Infos = Info.matrix,
                           start = mad, startval = matrix(sd,1,1), ustartval = matrix(sd,1,1)))
            }
        }
    }
}
