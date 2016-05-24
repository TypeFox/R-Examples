.LowerCaseMultivariate <- function(L2deriv, neighbor, biastype,
             normtype, Distr, Finfo, trafo, z.start = NULL,
             A.start = NULL, z.comp = NULL, A.comp = NULL, maxiter, tol,
             verbose = NULL){

        if(missing(verbose)|| is.null(verbose))
           verbose <- getRobAStBaseOption("all.verbose")

        w <- new("HampelWeight")

        if(is.null(z.start)) z.start <- numeric(ncol(trafo))
        if(is.null(A.start)) A.start <- trafo%*%solve(as.matrix(Finfo))
        if(is.null(A.comp)) 
           A.comp <- matrix(TRUE, nrow = nrow(trafo), ncol = ncol(trafo))
        if(is.null(z.comp)) 
           z.comp <- rep(TRUE, nrow(trafo))

        force(normtype)
        lA.comp <- sum(A.comp)
        
        abs.fct <- function(x, L2, stand, cent, normtype.0){
            X <- evalRandVar(L2, as.matrix(x))[,,1] - cent
            Y <- stand %*% X
            return(fct(normtype.0)(Y))
        }

        itermin <- 0
        bmin.fct <- function(param, L2deriv, Distr, trafo){
            itermin <<- itermin + 1
            p <- nrow(trafo)
            k <- ncol(trafo)
            A <- matrix(0, ncol = k, nrow = p)
            
            A[A.comp] <- param[1:lA.comp]
            A.max <- max(abs(A.comp))
            A <- A/A.max
            z <- numeric(k)
            z[z.comp] <- param[(lA.comp+1):length(param)]

#            if(is(normtype,"SelfNorm")) 
#               A <- A/max(A)
            
            w0 <- w
            cent(w0) <- z
            stand(w0) <- A
            weight(w0) <- minbiasweight(w0, neighbor = neighbor,
                                           biastype = biastype,
                                           normW = normtype)
            w <<- w0
            if (is(normtype,"SelfNorm")){
               normtype  <<- updateNorm(normtype = normtype, L2 = L2deriv,
                                        neighbor = neighbor, biastype = biastype,
                                        Distr = Distr, V.comp = A.comp,
                                        cent = z, stand = A,  w = w0)
               weight(w0) <- minbiasweight(w0, neighbor = neighbor,
                                           biastype = biastype,
                                           normW = normtype)

               w <<- w0
               }

            E1 <- E(object = Distr, fun = abs.fct, L2 = L2deriv, stand = A,
                     cent = z, normtype.0 = normtype, useApply = FALSE)
            stA <- if (is(normtype,"QFNorm"))
                       QuadForm(normtype)%*%A else A
#            erg <- E1/sum(diag(stA %*% t(trafo)))
#            print(list(A,stA,E1))
            erg <- E1/sum(diag(stA %*% t(trafo)))
            clip(w0) <- 1/erg
            w <<- w0
            if(verbose && itermin %% 15 == 1){
#            if(verbose && itermin %% 2 == 1){
                cat("trying to find lower case solution;\n")
               cat("current Lagrange Multiplier value:\n")
               print(list(A=A, z=z,erg=erg))
               }
  
            return(erg)
        }

        A.vec <- as.vector(A.start[A.comp])
        A.max <- max(abs(A.vec))
        A.vec <- A.vec/A.max
        p.vec <- c(A.vec, z.start[z.comp]/A.max)
        
        erg <- optim(p.vec, bmin.fct, method = "Nelder-Mead",
                    control = list(reltol = tol, maxit = 100*maxiter),
                    L2deriv = L2deriv, Distr = Distr, trafo = trafo)
        problem <- (erg$convergence > 0)
        A.max <- max(abs(stand(w)))
        stand(w) <- stand(w)/A.max
        weight(w) <- minbiasweight(w, neighbor = neighbor,
                                           biastype = biastype,
                                           normW = normtype)

        return(list(erg=erg, w=w, normtype = normtype, z.comp = z.comp, itermin = itermin,
                    problem = problem ))
    }


.LowerCaseMultivariateTV <- function(L2deriv, neighbor, biastype,
             normtype, Distr, Finfo, trafo,
             A.start = NULL,  maxiter, tol,
             verbose = NULL){

        if(missing(verbose)|| is.null(verbose))
           verbose <- getRobAStBaseOption("all.verbose")

        w <- new("BdStWeight")
        k <- ncol(trafo)

        if(is.null(A.start)) A.start <- trafo%*%solve(Finfo)

        pos.fct <- function(x, L2, stand){
            X <- evalRandVar(L2, as.matrix(x))[,,1]
            Y <- stand %*% X
            return(Y*(Y>0))
        }

        itermin <- 0

        bmin.fct <- function(param, L2deriv, Distr, trafo){
            itermin <<- itermin + 1
            p <- 1
            A <- matrix(param, ncol = k, nrow = 1)
         #   print(A)
            E1 <- E(object = Distr, fun = pos.fct, L2 = L2deriv, stand = A,
                    useApply = FALSE)
            erg <- E1/sum(diag(A %*% t(trafo)))
            return(erg)
        }

        erg <- optim(as.numeric(A.start), bmin.fct, method = "Nelder-Mead",
                    control = list(reltol = tol, maxit = 100*maxiter),
                    L2deriv = L2deriv, Distr = Distr, trafo = trafo)

        problem <- (erg$convergence > 0)
        A <- matrix(erg$par, ncol = k, nrow = 1)
        b <- 1/erg$value
        stand(w) <- A

        pr.fct <- function(x, L2, pr.sign=1){
                  X <- evalRandVar(L2, as.matrix(x)) [,,1]
                  Y <- as.numeric(A %*% X)
                  return(as.numeric(pr.sign*Y>0))
                  }
        p.p   <- E(object = Distr, fun = pr.fct, L2 = L2deriv,
                   useApply = FALSE, pr.sign =  1)
        m.p   <- E(object = Distr, fun = pr.fct, L2 = L2deriv,
                   useApply = FALSE, pr.sign = -1)


        a <- -b * p.p/(p.p+m.p)
        
        clip(w) <- c(0,b)+a
        weight(w) <- minbiasweight(w, neighbor = neighbor,
                                           biastype = biastype,
                                           normW = normtype)
        return(list(A=A,b=b, w=w, a=a, itermin = itermin, problem = problem))
    }

