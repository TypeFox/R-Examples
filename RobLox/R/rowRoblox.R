###############################################################################
## computation of k-step construction in case x is a matrix
###############################################################################
.onestep.loc.matrix <- function(x, initial.est, A, b, sd){
    u <- A*(x-initial.est)/sd^2
    ind <- b/abs(u) <= 1
    IC <- rowMeans(u*(ind*b/abs(u) + !ind), na.rm = TRUE)
    return(initial.est + IC)
}
.kstep.loc.matrix <- function(x, initial.est, A, b, sd, k){
    est <- initial.est
    for(i in 1:k){
        est <- .onestep.loc.matrix(x = x, initial.est = est, A = A, b = b, sd = sd)
    }
    return(est)
}
.onestep.sc.matrix <- function(x, initial.est, A, a, b, mean){
    v <- A*(((x-mean)/initial.est)^2-1)/initial.est - a
    ind <- b/abs(v) <= 1
    IC <- rowMeans(v*(ind*b/abs(v) + !ind), na.rm = TRUE)
    return(initial.est + IC)
}
.kstep.sc.matrix <- function(x, initial.est, A, a, b, mean, k){
    est <- .onestep.sc.matrix(x = x, initial.est = initial.est, A = A, a = a, b = b, mean = mean)
    if(k > 1){
        for(i in 2:k){
            A <- est^2*A/initial.est^2
            a <- est*a/initial.est
            b <- est*b/initial.est
            initial.est <- est
            est <- .onestep.sc.matrix(x = x, initial.est = est, A = A, a = a, b = b, mean = mean)
        }
    }
    A <- est^2*A/initial.est^2
    a <- est*a/initial.est
    b <- est*b/initial.est
    return(list(est = est, A = A, a = a, b = b))
}
.onestep.locsc.matrix <- function(x, initial.est, A1, A2, a, b){
    mean <- initial.est[,1]
    sd <- initial.est[,2]
    u <- A1*(x-mean)/sd^2
    v <- A2*(((x-mean)/sd)^2-1)/sd - a
    ind <- b/sqrt(u^2 + v^2) <= 1
    IC1 <- rowMeans(u*(ind*b/sqrt(u^2 + v^2) + !ind), na.rm = TRUE)
    IC2 <- rowMeans(v*(ind*b/sqrt(u^2 + v^2) + !ind), na.rm = TRUE)
    IC <- cbind(IC1, IC2)
    return(initial.est + IC)
}
.kstep.locsc.matrix <- function(x, initial.est, A1, A2, a, b, mean, k){
    est <- .onestep.locsc.matrix(x = x, initial.est = initial.est, A1 = A1, A2 = A2, a = a, b = b)
    if(k > 1){
        for(i in 2:k){
            A1 <- est[,2]^2*A1/initial.est[,2]^2
            A2 <- est[,2]^2*A2/initial.est[,2]^2
            a <- est[,2]*a/initial.est[,2]
            b <- est[,2]*b/initial.est[,2]
            initial.est <- est
            est <- .onestep.locsc.matrix(x = x, initial.est = est, A1 = A1, A2 = A2, a = a, b = b)
        }
    }
    A1 <- est[,2]^2*A1/initial.est[,2]^2
    A2 <- est[,2]^2*A2/initial.est[,2]^2
    a <- est[,2]*a/initial.est[,2]
    b <- est[,2]*b/initial.est[,2]

    return(list(est = est, A1 = A1, A2 = A2, a = a, b = b))
}


###############################################################################
## Evaluate roblox on rows of a matrix
###############################################################################
rowRoblox <- function(x, mean, sd, eps, eps.lower, eps.upper, initial.est, 
                      k = 1L, fsCor = TRUE, mad0 = 1e-4, na.rm = TRUE){
    es.call <- match.call()
    if(missing(x))
        stop("'x' is missing with no default")
    if(is.data.frame(x))
        x <- data.matrix(x)
    else
        x <- as.matrix(x)
    if(!is.matrix(x))
        stop("'x' has to be a matrix resp. convertable to a matrix by 'as.matrix'
              or 'data.matrix'")

    completecases <- complete.cases(x)
    if(na.rm) x <- na.omit(x)

    if(ncol(x) <= 2){
        if(missing(mean) && missing(sd)){
            warning("Sample size <= 2! => Median and MAD are used for estimation.")
            if(TRUE){#(require(Biobase))
                mean <- rowMedians(x, na.rm = TRUE)
                sd <- rowMedians(abs(x-mean), na.rm = TRUE)/qnorm(0.75)
            }else{
                warning("You can speed up the computations by installing Bioconductor package 'Biobase'!")
                mean <- apply(x, 1, median, na.rm = TRUE)
                sd <- apply(abs(x-mean), 1, median, na.rm = TRUE)/qnorm(0.75)
            }            
            robEst <- cbind(mean, sd)
            colnames(robEst) <- c("mean", "sd")
            Info.matrix <- matrix(c("rowRoblox", 
                                  paste("median and MAD")),
                                  ncol = 2, dimnames = list(NULL, c("method", "message")))
            return(new("ALEstimate", name = "Median and MAD", 
                       completecases = completecases,
                       estimate.call = es.call, estimate = robEst,
                       samplesize = ncol(x), asvar = NULL,
                       asbias = NULL, pIC = NULL, Infos = Info.matrix))
        }
        if(missing(mean)){
            warning("Sample size <= 2! => Median is used for estimation.")
            if(TRUE){#(require(Biobase))
                mean <- rowMedians(x, na.rm = TRUE)
            }else{
                warning("You can speed up the computations by installing Bioconductor package 'Biobase'!")
                mean <- apply(x, 1, median, na.rm = TRUE)
            }
            robEst <- matrix(mean)
            colnames(robEst) <- "mean"
            Info.matrix <- matrix(c("rowRoblox", 
                                  paste("median")),
                                  ncol = 2, dimnames = list(NULL, c("method", "message")))
            return(new("ALEstimate", name = "Median", 
                       completecases = completecases,
                       estimate.call = es.call, estimate = robEst,
                       samplesize = ncol(x), asvar = NULL,
                       asbias = NULL, pIC = NULL, Infos = Info.matrix))
        }
        if(missing(sd)){
            warning("Sample size <= 2! => MAD is used for estimation.")
            if(!is.numeric(mean) || (length(mean) != 1 && length(mean) != nrow(x)))
                stop("'mean' has wrong dimension")
            if(TRUE){#(require(Biobase))
                sd <- rowMedians(abs(x-mean), na.rm = TRUE)/qnorm(0.75)
            }else{
                warning("You can speed up the computations by installing Bioconductor package 'Biobase'!")
                sd <- apply(abs(x-mean), 1, median, na.rm = TRUE)/qnorm(0.75)
            }
            robEst <- matrix(sd)
            colnames(robEst) <- "sd"
            Info.matrix <- matrix(c("rowRoblox", 
                                  paste("MAD")),
                                  ncol = 2, dimnames = list(NULL, c("method", "message")))
            return(new("ALEstimate", name = "MAD", 
                       completecases = completecases,
                       estimate.call = es.call, estimate = robEst,
                       samplesize = ncol(x), asvar = NULL,
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
                mean <- rowMeans(x, na.rm = TRUE)
                n <- rowSums(!is.na(x))
                n[n < 1] <- NA
                sd <- rowSums((x - mean)^2, na.rm = TRUE)/n
                robEst <- cbind(mean, sd)

                colnames(robEst) <- c("mean", "sd")
                Info.matrix <- matrix(c("rowRoblox", 
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
                robEst <- matrix(rowMeans(x, na.rm = TRUE))
                colnames(robEst) <- "mean"
                Info.matrix <- matrix(c("rowRoblox", 
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
                if(!is.numeric(mean) || (length(mean) != 1 && length(mean) != nrow(x)))
                    stop("'mean' has wrong dimension")
                n <- rowSums(!is.na(x))
                n[n < 1] <- NA
                robEst <- matrix(rowSums((x - mean)^2, na.rm = TRUE)/n)
                colnames(robEst) <- "sd"
                Info.matrix <- matrix(c("rowRoblox", 
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
    if(k < 1){
        stop("'k' has to be some positive integer value")
    }
    if(length(k) != 1){
        stop("'k' has to be of length 1")
    }
    k <- as.integer(k)

    if(missing(mean) && missing(sd)){
        if(missing(initial.est)){
            if(TRUE){#(require(Biobase))
                mean <- rowMedians(x, na.rm = TRUE)
                sd <- rowMedians(abs(x-mean), na.rm = TRUE)/qnorm(0.75)
            }else{
                warning("You can speed up the computations by installing Bioconductor package 'Biobase'!")
                mean <- apply(x, 1, median, na.rm = TRUE)
                sd <- apply(abs(x-mean), 1, median, na.rm = TRUE)/qnorm(0.75)
            }
            if(any(sd == 0)){
                stop("'mad(x, na.rm = TRUE) = 0' for some samples => cannot compute a valid initial estimate.
                      To avoid division by zero 'mad0' is used. You could also specify 
                      a valid scale estimate via 'initial.est'. Note that you have to provide
                      a location and scale estimate.")
                sd[sd == 0] <- mad0
            }
            mean.sd <- cbind(mean=mean,sd=sd)
        }else{
            if(nrow(initial.est) != nrow(x) || ncol(initial.est) != 2)
                stop("'initial.est' has wrong dimension")
            mean <- initial.est[,1]
            sd <- initial.est[,2]
            if(any(sd <= 0))
                stop("initial estimate for scale <= 0 which is no valid scale estimate")
            mean.sd <- initial.est
        }

        if(!missing(eps)){
            r <- sqrt(ncol(x))*eps
            if(fsCor) r <- finiteSampleCorrection(r = r, n = ncol(x), model = "locsc")
            if(r > 10){
                b <- sd*1.618128043
                const <- 1.263094656
                A2 <- b^2*(1+r^2)/(1+const)
                A1 <- const*A2
                a <- -0.6277527697*A2/sd
                mse <- A1 + A2
            }else{
                A1 <- sd^2*.getA1.locsc(r)
                A2 <- sd^2*.getA2.locsc(r)
                a <- sd*.geta.locsc(r)
                b <- sd*.getb.locsc(r)
                mse <- A1 + A2
            }
            robEst <- .kstep.locsc.matrix(x = x, initial.est = cbind(mean, sd), 
                                          A1 = A1, A2 = A2, a = a, b = b, k = k)
            colnames(robEst$est) <- c("mean", "sd")
            Info.matrix <- matrix(c("roblox", 
                                    paste("optimally robust estimates for contamination 'eps' =", round(eps, 3),
                                          "and 'asMSE'")),
                                  ncol = 2, dimnames = list(NULL, c("method", "message")))
            return(new("kStepEstimate", name = "Optimally robust estimate",
                       completecases = completecases,
                       estimate.call = es.call, estimate = robEst$est,
                       samplesize = ncol(x), steps = k, 
                       pIC = NULL, Infos = Info.matrix,
                       start = mean.sd, startval = mean.sd, ustartval = mean.sd))
## we need a class like "list of estimates" to set asvar and asbias consistently ...
#            return(new("kStepEstimate", name = "Optimally robust estimate",
#                       estimate = robEst$est, samplesize = ncol(x), asvar = NULL, 
#                       asbias = r*robEst$b, steps = k, pIC = NULL, Infos = Info.matrix))
        }else{
            sqrtn <- sqrt(ncol(x))
            rlo <- sqrtn*eps.lower
            rup <- sqrtn*eps.upper
            if(rlo > 10){
                r <- (rlo + rup)/2
            }else{
                r <- uniroot(.getlsInterval, lower = rlo+1e-8, upper = rup, 
                             tol = .Machine$double.eps^0.25, rlo = rlo, rup = rup)$root
            }
            if(fsCor) r <- finiteSampleCorrection(r = r, n = ncol(x), model = "locsc")
            if(r > 10){
                b <- sd*1.618128043
                const <- 1.263094656
                A2 <- b^2*(1+r^2)/(1+const)
                A1 <- const*A2
                a <- -0.6277527697*A2/sd
                mse <- A1 + A2
            }else{
                A1 <- sd^2*.getA1.locsc(r)
                A2 <- sd^2*.getA2.locsc(r)
                a <- sd*.geta.locsc(r)
                b <- sd*.getb.locsc(r)
                mse <- A1 + A2
            }
            if(rlo == 0){
                ineff <- (A1 + A2 - b^2*r^2)/(1.5*sd^2)
            }else{
                if(rlo > 10){
                    ineff <- 1
                }else{
                    A1lo <- sd^2*.getA1.locsc(rlo)
                    A2lo <- sd^2*.getA2.locsc(rlo)
                    ineff <- (A1 + A2 - b^2*(r^2 - rlo^2))/(A1lo + A2lo)
                }
            }
            robEst <- .kstep.locsc.matrix(x = x, initial.est = cbind(mean, sd), 
                                          A1 = A1, A2 = A2, a = a, b = b, k = k)
            colnames(robEst$est) <- c("mean", "sd")
            Info.matrix <- matrix(c(rep("roblox", 3), 
                                  paste("radius-minimax estimates for contamination interval [", 
                                    round(eps.lower, 3), ", ", round(eps.upper, 3), "]", sep = ""),
                                  paste("least favorable contamination: ", round(r/sqrtn, 3), sep = ""),
                                  paste("maximum MSE-inefficiency: ", round(ineff[1], 3), sep = "")), 
                                  ncol = 2, dimnames = list(NULL, c("method", "message")))
            return(new("kStepEstimate", name = "Optimally robust estimate",
                       completecases = completecases,
                       estimate.call = es.call, estimate = robEst$est, #
                       samplesize = ncol(x), steps = k, 
                       pIC = NULL, Infos = Info.matrix,
                       start = mean.sd, startval = mean.sd, ustartval = mean.sd))
## we need a class like "list of estimates" to set asvar and asbias consistently ...
#            return(new("kStepEstimate", name = "Optimally robust estimate",
#                       estimate = robEst$est, samplesize = ncol(x), asvar = NULL, 
#                       asbias = r*robEst$b, steps = k, pIC = NULL, Infos = Info.matrix))
        }
    }else{
        if(missing(mean)){
            if(any(sd <= 0))
                stop("'sd' has to be positive")
            if(!is.numeric(sd) || (length(sd) != 1 && length(sd) != nrow(x)))
                stop("'sd' has wrong dimension")
            if(missing(initial.est)){
                if(TRUE){#(require(Biobase))
                    mean <- rowMedians(x, na.rm = TRUE)
                }else{
                    warning("You can speed up the computations by installing Bioconductor package 'Biobase'!")
                    mean <- apply(x, 1, median, na.rm = TRUE)
                }
            }else{
                if(!is.numeric(initial.est) || length(initial.est) != nrow(x))
                    stop("'initial.est' has wrong dimension")
                mean <- initial.est
            }
            if(length(sd) == 1)
                sd <- rep(sd, length(mean))

            mean.sd <- cbind(mean=mean,sd=sd)

            if(!missing(eps)){
                r <- sqrt(ncol(x))*eps
                if(fsCor) r <- finiteSampleCorrection(r = r, n = ncol(x), model = "loc")
                if(r > 10){
                    b <- sd*sqrt(pi/2)
                    A <- b^2*(1+r^2)
                }else{
                    A <- sd^2*.getA.loc(r)
                    b <- sd*.getb.loc(r)
                }
                robEst <- as.matrix(.kstep.loc.matrix(x = x, initial.est = mean, A = A, b = b, sd = sd, k = k))
                colnames(robEst) <- "mean"
                Info.matrix <- matrix(c("roblox", 
                                        paste("optimally robust estimates for contamination 'eps' =", round(eps, 3),
                                              "and 'asMSE'")),
                                      ncol = 2, dimnames = list(NULL, c("method", "message")))
                return(new("kStepEstimate", name = "Optimally robust estimate",
                       completecases = completecases,
                           estimate.call = es.call, estimate = robEst,
                           samplesize = ncol(x), steps = k, 
                           pIC = NULL, Infos = Info.matrix,
                       start = mean.sd[,1,drop=F], startval = mean.sd[,1,drop=F], ustartval = mean.sd))
## we need a class like "list of estimates" to set asvar and asbias consistently ...
#                return(new("kStepEstimate", name = "Optimally robust estimate",
#                           estimate = robEst$est, samplesize = ncol(x), asvar = as.matrix(A - r^2*b^2), 
#                           asbias = r*b, steps = k, pIC = NULL, Infos = Info.matrix))
            }else{
                sqrtn <- sqrt(ncol(x))
                rlo <- sqrtn*eps.lower
                rup <- sqrtn*eps.upper
                if(rlo > 10){ 
                    r <- (rlo+rup)/2
                }else{
                    r <- uniroot(.getlInterval, lower = rlo+1e-8, upper = rup, 
                                 tol = .Machine$double.eps^0.25, rlo = rlo, rup = rup)$root
                }
                if(fsCor) r <- finiteSampleCorrection(r = r, n = ncol(x), model = "loc")
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
                robEst <- as.matrix(.kstep.loc.matrix(x = x, initial.est = mean, A = A, b = b, sd = sd, k = k))
                colnames(robEst) <- "mean"
                Info.matrix <- matrix(c(rep("roblox", 3), 
                                      paste("radius-minimax estimates for contamination interval [", 
                                        round(eps.lower, 3), ", ", round(eps.upper, 3), "]", sep = ""),
                                      paste("least favorable contamination: ", round(r/sqrtn, 3), sep = ""),
                                      paste("maximum MSE-inefficiency: ", round(ineff[1], 3), sep = "")), 
                                      ncol = 2, dimnames = list(NULL, c("method", "message")))
                return(new("kStepEstimate", name = "Optimally robust estimate",
                       completecases = completecases,
                           estimate.call = es.call, estimate = robEst,
                           samplesize = ncol(x), steps = k, 
                           pIC = NULL, Infos = Info.matrix,
                       start = mean.sd[,1,drop=F], startval = mean.sd[,1,drop=F], ustartval = mean.sd))
## we need a class like "list of estimates" to set asvar and asbias consistently ...
#                return(new("kStepEstimate", name = "Optimally robust estimate",
#                           estimate = robEst$est, samplesize = ncol(x), asvar = as.matrix(A - r^2*b^2), 
#                           asbias = r*b, steps = k, pIC = NULL, Infos = Info.matrix))
            }
        }
        if(missing(sd)){
            if(!is.numeric(mean) || (length(mean) != 1 && length(mean) != nrow(x)))
                stop("'mean' has wrong dimension")
            if(missing(initial.est)){
                if(TRUE){#(require(Biobase))
                    sd <- rowMedians(abs(x-mean), na.rm = TRUE)/qnorm(0.75)
                }else{
                    warning("You can speed up the computations by installing Bioconductor package 'Biobase'!")
                    sd <- apply(abs(x-mean), 1, median, na.rm = TRUE)/qnorm(0.75)
                }
                if(any(sd == 0)){
                    stop("'mad(x, na.rm = TRUE) = 0' for some samples => cannot compute a valid initial estimate.
                          To avoid division by zero 'mad0' is used. You could also specify 
                          a valid scale estimate via 'initial.est'.")
                    sd[sd == 0] <- mad0
                }
            }else{
                if(!is.numeric(initial.est) || length(initial.est) != nrow(x))
                    stop("'initial.est' has wrong dimension")
                if(any(initial.est <= 0))
                    stop("'initial.est <= 0'; i.e., is no valid scale estimate")
                sd <- initial.est
            }

            mean.sd <- cbind(mean=mean,sd=sd)

            if(!missing(eps)){
                r <- sqrt(ncol(x))*eps
                if(fsCor) r <- finiteSampleCorrection(r = r, n = ncol(x), model = "sc")
                if(r > 10){
                    b <- sd/(4*qnorm(0.75)*dnorm(qnorm(0.75)))
                    A <- b^2*(1+r^2)
                    a <- (qnorm(0.75)^2 - 1)/sd*A
                }else{
                    A <- sd^2*.getA.sc(r)
                    a <- sd*.geta.sc(r)
                    b <- sd*.getb.sc(r)
                }
                robEst <- .kstep.sc.matrix(x = x, initial.est = sd, A = A, a = a, b = b, mean = mean, k = k)
                robEst$est <- as.matrix(robEst$est)
                colnames(robEst$est) <- "sd"
                Info.matrix <- matrix(c("roblox", 
                                        paste("optimally robust estimates for contamination 'eps' =", round(eps, 3),
                                              "and 'asMSE'")),
                                      ncol = 2, dimnames = list(NULL, c("method", "message")))
                return(new("kStepEstimate", name = "Optimally robust estimate",
                       completecases = completecases,
                           estimate.call = es.call, estimate = robEst$est,
                           samplesize = ncol(x), steps = k, 
                           pIC = NULL, Infos = Info.matrix,
                       start = mean.sd[,2,drop=F], startval = mean.sd[,2,drop=F], ustartval = mean.sd))
## we need a class like "list of estimates" to set asvar and asbias consistently ...
#                return(new("kStepEstimate", name = "Optimally robust estimate",
#                           estimate = robEst$est, samplesize = ncol(x), asvar = as.matrix(robEst$A - r^2*robEst$b^2), 
#                           asbias = r*robEst$b, steps = k, pIC = NULL, Infos = Info.matrix))
            }else{
                sqrtn <- sqrt(ncol(x))
                rlo <- sqrtn*eps.lower
                rup <- sqrtn*eps.upper
                if(rlo > 10){
                    r <- (rlo+rup)/2
                }else{
                    r <- uniroot(.getsInterval, lower = rlo+1e-8, upper = rup, 
                             tol = .Machine$double.eps^0.25, rlo = rlo, rup = rup)$root
                }
                if(fsCor) r <- finiteSampleCorrection(r = r, n = ncol(x), model = "sc")
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
                robEst <- .kstep.sc.matrix(x = x, initial.est = sd, A = A, a = a, b = b, mean = mean, k = k)
                robEst$est <- as.matrix(robEst$est)
                colnames(robEst$est) <- "sd"
                Info.matrix <- matrix(c(rep("roblox", 3), 
                                      paste("radius-minimax estimates for contamination interval [", 
                                        round(eps.lower, 3), ", ", round(eps.upper, 3), "]", sep = ""),
                                      paste("least favorable contamination: ", round(r/sqrtn, 3), sep = ""),
                                      paste("maximum MSE-inefficiency: ", round(ineff[1], 3), sep = "")), 
                                      ncol = 2, dimnames = list(NULL, c("method", "message")))
                return(new("kStepEstimate", name = "Optimally robust estimate",
                       completecases = completecases,
                           estimate.call = es.call, estimate = robEst$est,
                           samplesize = ncol(x), steps = k, 
                           pIC = NULL, Infos = Info.matrix,
                       start = mean.sd[,2,drop=F], startval = mean.sd[,2,drop=F],
                       ustartval = mean.sd))
## we need a class like "list of estimates" to set asvar and asbias consistently ...
#                return(new("kStepEstimate", name = "Optimally robust estimate",
#                           estimate = robEst$est, samplesize = ncol(x), asvar = as.matrix(robEst$A - r^2*robEst$b^2), 
#                           asbias = r*robEst$b, steps = k, pIC = NULL, Infos = Info.matrix))
            }
        }
    }
}
