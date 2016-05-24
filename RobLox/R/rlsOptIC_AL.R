###############################################################################
## weight function
###############################################################################
.ALrlsGetw <- function(x, b, a1, a2, a3){
    hvkt <- sqrt(a3^2*x^4 + (a1^2 - 2*a2*a3^2)*x^2 + a2^2*a3^2)
    ind1 <- (hvkt < b)

    return(ind1 + (1-ind1)*b/hvkt)
}

###############################################################################
## computation of r
###############################################################################
.ALrlsGetr <- function(b, r, a1, a2, a3){
    integrandr <- function(x, b, a1, a2, a3){
        hvkt <- sqrt(a3^2*x^4 + (a1^2 - 2*a2*a3^2)*x^2 + a2^2*a3^2)/b - 1
        return((hvkt > 0)*hvkt*dnorm(x))
    }
    Int <- integrate(integrandr, lower = 0, upper = Inf, 
                rel.tol = .Machine$double.eps^0.5, a1 = a1, a2 = a2, 
                a3 = a3, b = b)$value

    return(r-sqrt(2*Int))
}


###############################################################################
## computation of a1, a2 and a3
###############################################################################
.ALrlsGeta1a2a3 <- function(b, a1, a2, a3){
    integrand1 <- function(x, b, a1, a2, a3){ 
        x^2*.ALrlsGetw(x, b, a1, a2, a3)*dnorm(x)
    }
    Int1 <- 2*integrate(integrand1, lower = 0, upper = Inf, 
                    rel.tol = .Machine$double.eps^0.5, b = b, a1 = a1, 
                    a2 = a2, a3 = a3)$value
    a1 <- 1/Int1

    integrand2 <- function(x, b, a1, a2, a3){
        .ALrlsGetw(x, b, a1, a2, a3)*dnorm(x)
    }
    Int2 <- 2*integrate(integrand2, lower = 0, upper = Inf, 
                    rel.tol = .Machine$double.eps^0.5, b = b, a1 = a1, 
                    a2 = a2, a3 = a3)$value
    a2 <- Int1/Int2

    integrand3 <- function(x, b, a1, a2, a3){
        (x^2 - a2)^2*.ALrlsGetw(x, b, a1, a2, a3)*dnorm(x)
    }
    Int3 <- 2*integrate(integrand3, lower = 0, upper = Inf, 
                    rel.tol = .Machine$double.eps^0.5, b = b, a1 = a1, 
                    a2 = a2, a3 = a3)$value
    a3 <- 1/Int3

    return(list(a1=a1, a2=a2, a3=a3))
}
.ALrlsVar <- function(b, a1, a2, a3){
    integrand1 <- function(x, b, a1, a2, a3){ 
        x^2*.ALrlsGetw(x, b, a1, a2, a3)^2*dnorm(x)
    }
    Int1 <- 2*integrate(integrand1, lower = 0, upper = Inf, 
                    rel.tol = .Machine$double.eps^0.5, b = b, a1 = a1, 
                    a2 = a2, a3 = a3)$value
    V1 <- a1^2*Int1

    integrand2 <- function(x, b, a1, a2, a3){
        (x^2 - a2)^2*.ALrlsGetw(x, b, a1, a2, a3)^2*dnorm(x)
    }
    Int2 <- 2*integrate(integrand2, lower = 0, upper = Inf, 
                    rel.tol = .Machine$double.eps^0.5, b = b, a1 = a1, 
                    a2 = a2, a3 = a3)$value
    V2 <- a3^2*Int2

    return(diag(c(V1, V2)))
}


###############################################################################
## optimal IC
###############################################################################
rlsOptIC.AL <- function(r, mean = 0, sd = 1, A.loc.start = 1, a.sc.start = 0, 
                        A.sc.start = 0.5, bUp = 1000, delta = 1e-6, itmax = 100, 
                        check = FALSE, computeIC = TRUE){
    a1 <- A.loc.start; a2 <- 1+a.sc.start; a3 <- A.sc.start
    b <- uniroot(.ALrlsGetr, lower = 1e-4, upper = bUp, 
            tol = .Machine$double.eps^0.5, r = r, a1 = a1, a2 = a2, 
            a3 = a3)$root

    iter <- 0
    repeat{
        iter <- iter + 1
        if(iter > itmax){
            stop("Algorithm did not converge!\n", 
                 "=> increase itmax or try different starting values",
                 "'A.loc.start', 'a.sc.start' and 'A.sc.start'\n")
        }
        a1.old <- a1; a2.old <- a2; a3.old <- a3; b.old <- b

        a1a2a3 <- .ALrlsGeta1a2a3(b = b, a1 = a1, a2 = a2, a3 = a3)
        a1 <- a1a2a3$a1
        a2 <- a1a2a3$a2
        a3 <- a1a2a3$a3

        b <- uniroot(.ALrlsGetr, lower = 1e-4, upper = bUp, 
                tol = .Machine$double.eps^0.5, r = r, a1 = a1, a2 = a2, 
                a3 = a3)$root
        if(max(abs(a1.old-a1), abs(a2.old-a2), abs(a3.old-a3), abs(b.old-b))<delta)
            break
    }

    if(check){
        integrand1 <- function(x, b, a1, a2, a3){
            x^2*.ALrlsGetw(x, b, a1, a2, a3)*dnorm(x)
        }
        Int1 <- 2*integrate(integrand1, lower = 0, upper = Inf, 
                        rel.tol = .Machine$double.eps^0.5, b = b, a1 = a1, 
                        a2 = a2, a3 = a3)$value
        ch1 <- a1*Int1

        integrand2 <- function(x, b, a1, a2, a3){
            (x^2 - a2)^2*.ALrlsGetw(x, b, a1, a2, a3)*dnorm(x)
        }
        Int2 <- 2*integrate(integrand2, lower = 0, upper = Inf, 
                        rel.tol = .Machine$double.eps^0.5, b = b, a1 = a1, 
                        a2 = a2, a3 = a3)$value
        ch2 <- a3*Int2

        integrand3 <- function(x, b, a1, a2, a3){
            (x^2 - a2)*.ALrlsGetw(x, b, a1, a2, a3)*dnorm(x)
        }
        Int3 <- 2*integrate(integrand3, lower=0, upper=Inf, 
                        rel.tol = .Machine$double.eps^0.5, b = b, a1 = a1, 
                        a2 = a2, a3 = a3)$value
        ch3 <- a3*Int3

        ch4 <- .ALrlsGetr(b = b, r = r, a1 = a1, a2 = a2, a3 = a3)

        cat("Fisher consistency of eta.loc:\t", ch1-1, "\n")
        cat("centering of eta.sc:\t", ch3, "\n")
        cat("Fisher consistency of eta.sc:\t", ch2-1, "\n")
        cat("MSE equation:\t", ch4, "\n")
    }

    asVar <- sd^2*.ALrlsVar(b = b, a1 = a1, a2 = a2, a3 = a3)
    A <- sd^2*diag(c(a1, a3))
    a <- sd*c(0, a3*(a2-1))
    b <- sd*b
    mse <- sd^2*(a1 + a3)


    if(computeIC){
        w <- new("HampelWeight")
        clip(w) <- b
        cent(w) <- c(0, a2-1)/sd
        stand(w) <- A
        weight(w) <- getweight(w, neighbor = ContNeighborhood(radius = r), 
                               biastype = symmetricBias(), 
                               normW = NormType())

        modIC <- function(L2Fam, IC){
            ICL2Fam <- eval(CallL2Fam(IC))
            if(is(L2Fam, "L2LocationScaleFamily") && is(distribution(L2Fam), "Norm")){
                sdneu <- main(L2Fam)[2]
                sdalt <- main(ICL2Fam)[2]
                w <- weight(IC)
                clip(w) <- sdneu*clip(w)/sdalt
                cent(w) <- sdalt*cent(w)/sdneu
                stand(w) <- sdneu^2*stand(w)/sdalt^2
                weight(w) <- getweight(w, neighbor = ContNeighborhood(radius = neighborRadius(IC)), 
                               biastype = biastype(IC), 
                               normW = normtype(IC))
                A <- sdneu^2*stand(IC)/sdalt^2
                b <- sdneu*clip(IC)/sdalt
                mse <- sum(diag(A))
                r <- neighborRadius(IC)
                a1 <- A[1, 1]/sdneu^2
                a3 <- A[2, 2]/sdneu^2
                a2 <- a[1]/sd/a3 + 1
                asVar <- sdneu^2*.ALrlsVar(b = b/sdneu, a1 = a1, a2 = a2, a3 = a3)
                res <- list(A = A, a = sdneu*cent(IC)/sdalt, b = b, d = NULL,
                            risk = list(asMSE = mse, asBias = b, trAsCov = mse - r^2*b^2,
                                        asCov = asVar), 
                            info = Infos(IC), w = w,
                            normtype = normtype(IC), biastype = biastype(IC),
                            modifyIC = modifyIC(IC))
                IC <- generateIC(neighbor = ContNeighborhood(radius = r),
                                 L2Fam = L2Fam, res = res)
                addInfo(IC) <- c("modifyIC", "The IC has been modified")
                addInfo(IC) <- c("modifyIC", "The entries in 'Infos' may be wrong")
                return(IC)
            }else{
                makeIC(IC, L2Fam)
            }
        }

        L2Fam <- substitute(NormLocationScaleFamily(mean = m1, sd = s1), 
                            list(m1 = mean, s1 = sd))
        return(generateIC(neighbor = ContNeighborhood(radius = r), 
                    L2Fam = eval(L2Fam), 
                    res = list(A = as.matrix(A), a = a, b = b, d = NULL, 
                               risk = list(asMSE = mse, asBias = b, trAsCov = mse - r^2*b^2,
                                           asCov = asVar), 
                               info = c("rlOptIC", "optimally robust IC for AL estimators and 'asMSE'"),
                               w = w, biastype = symmetricBias(), normtype = NormType(),
                               modifyIC = modIC)))
    }else{
        return(list(A = A, a = a, b = b))
    }
}
