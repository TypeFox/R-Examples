###############################################################################
## Use robloxbioc to compute optimally robust (rmx) estimator for rows of a 
## matrix
###############################################################################
setMethod("robloxbioc", signature(x = "matrix"),
    function(x, eps = NULL, eps.lower = 0, eps.upper = 0.05, steps = 3L, 
             fsCor = TRUE, mad0 = 1e-4){
        stopifnot(is.numeric(x))
        if(ncol(x) <= 2){
            mean <- rowMedians(x, na.rm = TRUE)
            sd <- rowMedians(abs(x-mean), na.rm = TRUE)/qnorm(0.75)
            robEst <- cbind(mean, sd)
            colnames(robEst) <- c("mean", "sd")
            return(robEst)
        }
        if(is.null(eps)){
            if(length(eps.lower) != 1 || length(eps.upper) != 1)
                stop("'eps.lower' and 'eps.upper' have to be of length 1")
            if(!is.numeric(eps.lower) || !is.numeric(eps.upper) || eps.lower >= eps.upper) 
                stop("'eps.lower' < 'eps.upper' is not fulfilled")
            if((eps.lower < 0) || (eps.upper > 0.5))
                stop("'eps.lower' and 'eps.upper' have to be in [0, 0.5]")
        }else{
            if(length(eps) != 1){
                warning("'eps' has to be of length 1 => only first element is used")
                eps <- eps[1]
            }
            if(!is.numeric(eps))
                stop("'eps' has to be a double in (0, 0.5]")
            if((eps < 0) || (eps > 0.5))
                stop("'eps' has to be in (0, 0.5]")
            if(eps == 0){
                warning("eps = 0! => Mean and sd are used for estimation.")
                mean <- rowMeans(x, na.rm = TRUE)
                n <- rowSums(!is.na(x))
                n[n < 1] <- NA
                sd <- rowSums((x - mean)^2, na.rm = TRUE)/n
                robEst <- cbind(mean, sd)
                colnames(robEst) <- c("mean", "sd")
                return(robEst)
            }
        }
        if(length(steps) != 1){
            warning("'steps' has to be of length 1 => only first element is used!")
            steps <- steps[1]
        }
        if(steps < 1)
            stop("'steps' has to be some positive integer value")
        steps <- as.integer(steps)
        if(steps > 10)
            warning("steps > 10 => numbers between 1 and 5 should be sufficient.")

        mean <- rowMedians(x, na.rm = TRUE)
        sd <- rowMedians(abs(x-mean), na.rm = TRUE)/qnorm(0.75)
        if(any(sd == 0)){
            warning("Some of the initial scale estimates were 0 => set to 'mad0'") 
            sd[sd == 0] <- mad0
        }

        if(!is.null(eps)){
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
                                          A1 = A1, A2 = A2, a = a, b = b, k = steps)$est
            colnames(robEst) <- c("mean", "sd")
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
                                          A1 = A1, A2 = A2, a = a, b = b, k = steps)$est
            colnames(robEst) <- c("mean", "sd")
        }
        return(robEst)
    })

