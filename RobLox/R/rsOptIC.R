###############################################################################
## computation of centering constang
###############################################################################
.ALrsGetz <- function(z, b){
    b1 <- sqrt(pmax(z - b,0))
    b2 <- sqrt(b + z)

    return(b1*dnorm(b1) - b2*dnorm(b2) + (1 - z - b)*pnorm(b2) 
           - (1 - z + b)*pnorm(b1) + 1.5*b)
}

###############################################################################
## computation of clipping bound
###############################################################################
.ALrsGetr <- function(b, r, z){
    b1 <- sqrt(pmax(z - b,0))
    b2 <- sqrt(b + z)

    return(r^2*b - 2*((1 - z - b)*(1-pnorm(b2)) + b2*dnorm(b2) 
                      + (b - z + 1)*(1/2-pnorm(b1)) + b1*dnorm(b1)))
}


###############################################################################
## optimal IC
###############################################################################
rsOptIC <- function(r, mean = 0, sd = 1, bUp = 1000, delta = 1e-6, itmax = 100,
                    computeIC = TRUE){
    z <- 1
    b <- uniroot(.ALrsGetr, lower = 1e-4, upper = bUp, 
                 tol = .Machine$double.eps^0.5, r = r, z = z)$root

    iter <- 0
    repeat{
        iter <- iter + 1
        if(iter > itmax) 
            stop("Algorithm did not converge => increase 'itmax'!")

        z.old <- z; b.old <- b

        z <- uniroot(.ALrsGetz, lower = 0, upper = 1, tol = .Machine$double.eps^0.5, 
                     b = b)$root

        b <- uniroot(.ALrsGetr, lower = 1e-4, upper = bUp, 
                     tol = .Machine$double.eps^0.5, r = r, z = z)$root

        if(max(abs(z.old-z), abs(b.old-b))<delta)
            break
    }


    b1 <- sqrt(pmax(z-b,0))
    b2 <- sqrt(b + z)

    aa <- dnorm(b2)*(-b2^3 - 3*b2 + 2*z*b2) - dnorm(b1)*(-b1^3 - 3*b1 + 2*z*b1)
    aa <- aa + (pnorm(b2) - pnorm(b1))*(z^2 + 3 - 2*z) + b*(1-z)*(1.5 - pnorm(b2) - pnorm(b1))
    aa <- aa + b*(b2*dnorm(b2) + b1*dnorm(b1))

    A1 <- 1/(2*aa)
    b <- sd*A1*b
    a <- sd*A1*(z - 1)
    A <- sd^2*A1

    if(computeIC){
        w <- new("HampelWeight")
        clip(w) <- b
        cent(w) <- (z-1)/sd
        stand(w) <- as.matrix(A)
        weight(w) <- getweight(w, neighbor = ContNeighborhood(radius = r), 
                               biastype = symmetricBias(), 
                               normW = NormType())

        modIC <- function(L2Fam, IC){
            ICL2Fam <- eval(CallL2Fam(IC))
            if(is(L2Fam, "L2ScaleFamily") && is(distribution(L2Fam), "Norm")){
                sdneu <- main(L2Fam)
                sdalt <- main(ICL2Fam)
                w <- weight(IC)
                clip(w) <- sdneu*clip(w)/sdalt
                cent(w) <- sdalt*cent(w)/sdneu
                stand(w) <- sdneu^2*stand(w)/sdalt^2
                weight(w) <- getweight(w, neighbor = ContNeighborhood(radius = neighborRadius(IC)), 
                               biastype = biastype(IC), 
                               normW = normtype(IC))
                A <- sdneu^2*stand(IC)/sdalt^2
                b <- sdneu*clip(IC)/sdalt
                res <- list(A = as.matrix(A), a = sdneu*cent(IC)/sdalt, b = b, d = NULL,
                            risk = list(asMSE = A, asBias = b, asCov = A - r^2*b^2), 
                            info = Infos(IC), w = w,
                            normtype = normtype(IC), biastype = biastype(IC),
                            modifyIC = modifyIC(IC))
                IC <- generateIC(neighbor = ContNeighborhood(radius = neighborRadius(IC)),
                                 L2Fam = L2Fam, res = res)
                addInfo(IC) <- c("modifyIC", "The IC has been modified")
                addInfo(IC) <- c("modifyIC", "The entries in 'Infos' may be wrong")
                return(IC)
            }else{
                makeIC(IC, L2Fam)
            }
        }

        L2Fam <- substitute(NormScaleFamily(sd = s1, mean = m1), 
                            list(m1 = mean, s1 = sd))
        return(generateIC(neighbor = ContNeighborhood(radius = r), 
                    L2Fam = eval(L2Fam), 
                    res = list(A = as.matrix(A), a = a, b = b, d = NULL, 
                               risk = list(asMSE = A, asBias = b, asCov = A - r^2*b^2), 
                               info = c("rlOptIC", "optimally robust IC for AL estimators and 'asMSE'"),
                               w = w, biastype = symmetricBias(), normtype = NormType(),
                               modifyIC = modIC)))
    }else{
        return(list(A = A, a = a, b = b))
    }
}
