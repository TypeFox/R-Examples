NMs.randomCRW <- function(ltraj, rangles=TRUE, rdist=TRUE, fixedStart=TRUE,
                          x0=NULL, rx=NULL, ry=NULL, treatment.func=NULL,
                          treatment.par=NULL, constraint.func=NULL,
                          constraint.par=NULL, nrep=999)
{
    ## First checks the arguments
    if (!inherits(ltraj, "ltraj"))
        stop("ltraj should be of class ltraj")
    if (is.null(x0)) {
        x0 <- lapply(ltraj, function(x) c(x$x[1],x$y[1]))
    }
    if (!is.list(x0)) {
        x0 <- lapply(ltraj, function(x) x0)
    }
    if (length(x0)!=length(ltraj))
        stop("x0 should be a list with the same length as ltraj")
    if (any(unlist(lapply(x0, length))!=2))
        stop("each x0 should be of length 2")

    if (is.null(rx)) {
        rx <- do.call("rbind",lapply(ltraj, function(x) range(na.omit(x$x))))
        rx <- c(min(rx[,1]), max(rx[,2]))
    }
    if (is.null(ry)) {
        ry <- do.call("rbind",lapply(ltraj, function(x) range(na.omit(x$y))))
        ry <- c(min(ry[,1]), max(ry[,2]))
    }
    if (!is.list(rx)) {
        rx <- lapply(ltraj, function(x) rx)
    }
    if (!is.list(ry)) {
        ry <- lapply(ltraj, function(x) ry)
    }

    if (is.null(treatment.func)) {
        treatment.func <- function(x, par) return(x)
    }
    if (is.null(treatment.par)) {
        treatment.par <- 0
    }
    if (is.null(constraint.func)) {
        constraint.func <- function(x, par) return(TRUE)
    }
    if (is.null(constraint.par)) {
        constraint.par <- 0
    }

    ## convert the functions as lists of functions (same length as ltraj)
    ## if not already a list
    if (!is.list(treatment.func)) {
        treatment.func <- lapply(1:length(ltraj), function(i) treatment.func)
        treatment.par <- lapply(1:length(ltraj), function(i) treatment.par)
    }

    if (!is.list(constraint.func)) {
        constraint.func <- lapply(1:length(ltraj), function(i) constraint.func)
        constraint.par <- lapply(1:length(ltraj), function(i) constraint.par)
    }
    if (length(treatment.func)!=length(ltraj))
        stop("treatment.func should have the same length as ltraj")
    if (length(treatment.par)!=length(ltraj))
        stop("treatment.par should have the same length as ltraj")
    if (length(constraint.func)!=length(ltraj))
        stop("constraint.func should have the same length as ltraj")
    if (length(constraint.par)!=length(ltraj))
        stop("constraint.par should have the same length as ltraj")

    ## checks the name of the arguments in the function
    lapply(treatment.func, function(x) {
        art <- names(as.list(args(x)))
        if ((!all(art[1:2]==c("x","par")))|(length(art)!=3))
            stop("The treatment function should have only two arguments named x and par")
    })
    lapply(constraint.func, function(x) {
        art <- names(as.list(args(x)))
        if ((!all(art[1:2]==c("x","par")))|(length(art)!=3))
            stop("The treatment function should have only two arguments named x and par")
    })

    ## checks if the treatment works at least on the argument
    res <- lapply(1:length(ltraj), function(i) {
        ii <- try(treatment.func[[i]](ltraj[[i]][,1:3], treatment.par[[i]]))
        if (inherits(ii, "try-error"))
            stop(paste("The treatment function fails for animal", i,
                       ":\nit should work at least for all observed trajectories"))
    })

    ## Checks whether the constraints are satisfied on the argument
    ## warning otherwise
    warningConstraint <- unlist(lapply(1:length(ltraj), function(i) {
        ii <- try(constraint.func[[i]](ltraj[[i]][,1:3], constraint.par[[i]]))
        if (inherits(ii, "try-error"))
            stop(paste("The constraint function fails for animal", i,
                       ":\nit should work at least for all observed trajectories"))
        if (!is.logical(ii)|(length(ii)>1))
            stop("the constraint function should return a logical value of length 1")
        return(ii)
    }))
    if (any(!warningConstraint))
        warning(paste("The current ltraj does not satisfy the constraint for bursts:\n",
                      paste(burst(ltraj)[!warningConstraint], collapse=",")))

    ## build the object needed for the simulations
    res <- lapply(1:length(ltraj), function(i) {
        env <- new.env()
        assign("treatment", treatment.func[[i]], envir=env)
        assign("constraint", constraint.func[[i]], envir=env)
        ## par1 contains:
        ## - the environment
        ## - rangles
        ## - rdist
        ## - the type of starting point
        ## - fixed starting point coordinates
        ## - for random starting point: rx
        ## - for random starting point: ry

        param <- list(par1=list(env=env, rangles=as.integer(rangles),
                      rdist=as.integer(rdist),
                      fixedStart=as.integer(fixedStart),
                      x0=x0[[i]], rx=rx[[i]], ry=ry[[i]]),
                      par2=treatment.par[[i]],
                      parcon=constraint.par[[i]],
                      treatment=treatment.func[[i]],
                      constraint=constraint.func[[i]])
        list(xyt=ltraj[[i]], nrep=as.integer(nrep), type=as.integer(1), par = param)
    })

    names(res) <- burst(ltraj)
    class(res) <- c("NMs", "randomCRW")
    attr(res, "warningConstraint") <- warningConstraint
    return(res)
}




NMs.randomShiftRotation <- function(ltraj, rshift=TRUE, rrot=TRUE,
                                    rx=NULL, ry=NULL, treatment.func=NULL,
                                    treatment.par=NULL, constraint.func=NULL,
                                    constraint.par=NULL, nrep=999)
{
    ## First checks the arguments
    if (!inherits(ltraj, "ltraj"))
        stop("ltraj should be of class ltraj")
    if ((!rshift)&(!rrot))
        stop(paste("What do you want to randomize?\n",
                   "if rshift=FALSE and rrot=FALSE, the simulations\n",
                   "will be identical to the actual trajectory"))

    if (is.null(rx)) {
        rx <- do.call("rbind",lapply(ltraj, function(x) range(na.omit(x$x))))
        rx <- c(min(rx[,1]), max(rx[,2]))
    }
    if (is.null(ry)) {
        ry <- do.call("rbind",lapply(ltraj, function(x) range(na.omit(x$y))))
        ry <- c(min(ry[,1]), max(ry[,2]))
    }
    if (!is.list(rx)) {
        rx <- lapply(ltraj, function(x) rx)
    }
    if (!is.list(ry)) {
        ry <- lapply(ltraj, function(x) ry)
    }

    if (is.null(treatment.func)) {
        treatment.func <- function(x, par) return(x)
    }
    if (is.null(treatment.par)) {
        treatment.par <- 0
    }
    if (is.null(constraint.func)) {
        constraint.func <- function(x, par) return(TRUE)
    }
    if (is.null(constraint.par)) {
        constraint.par <- 0
    }

    ## convert the functions as lists of functions (same length as ltraj)
    ## if not already a list
    if (!is.list(treatment.func)) {
        treatment.func <- lapply(1:length(ltraj), function(i) treatment.func)
        treatment.par <- lapply(1:length(ltraj), function(i) treatment.par)
    }

    if (!is.list(constraint.func)) {
        constraint.func <- lapply(1:length(ltraj), function(i) constraint.func)
        constraint.par <- lapply(1:length(ltraj), function(i) constraint.par)
    }
    if (length(treatment.func)!=length(ltraj))
        stop("treatment.func should have the same length as ltraj")
    if (length(treatment.par)!=length(ltraj))
        stop("treatment.par should have the same length as ltraj")
    if (length(constraint.func)!=length(ltraj))
        stop("constraint.func should have the same length as ltraj")
    if (length(constraint.par)!=length(ltraj))
        stop("constraint.par should have the same length as ltraj")

    ## checks the name of the arguments in the function
    lapply(treatment.func, function(x) {
        art <- names(as.list(args(x)))
        if ((!all(art[1:2]==c("x","par")))|(length(art)!=3))
            stop("The treatment function should have only two arguments named x and par")
    })
    lapply(constraint.func, function(x) {
        art <- names(as.list(args(x)))
        if ((!all(art[1:2]==c("x","par")))|(length(art)!=3))
            stop("The treatment function should have only two arguments named x and par")
    })

    ## checks if the treatment works at least on the argument
    res <- lapply(1:length(ltraj), function(i) {
        ii <- try(treatment.func[[i]](ltraj[[i]][,1:3], treatment.par[[i]]))
        if (inherits(ii, "try-error"))
            stop(paste("The treatment function fails for animal", i,
                       ":\nit should work at least for all observed trajectories"))
    })

    ## Checks whether the constraints are satisfied on the argument
    ## warning otherwise
    warningConstraint <- unlist(lapply(1:length(ltraj), function(i) {
        ii <- try(constraint.func[[i]](ltraj[[i]][,1:3], constraint.par[[i]]))
        if (inherits(ii, "try-error"))
            stop(paste("The constraint function fails for animal", i,
                       ":\nit should work at least for all observed trajectories"))
        if (!is.logical(ii)|(length(ii)>1))
            stop("the constraint function should return a logical value of length 1")
        return(ii)
    }))
    if (any(!warningConstraint))
        warning(paste("The current ltraj does not satisfy the constraint for bursts:\n",
                      paste(burst(ltraj)[!warningConstraint], collapse=",")))

    ## build the object needed for the simulations
    res <- lapply(1:length(ltraj), function(i) {
        env <- new.env()
        assign("treatment", treatment.func[[i]], envir=env)
        assign("constraint", constraint.func[[i]], envir=env)
        ## par1 contains:
        ## - the environment
        ## - rshift
        ## - rrot
        ## - the type of starting point
        ## - fixed starting point coordinates
        ## - for random starting point: rx
        ## - for random starting point: ry

        param <- list(par1=list(env=env, rshift=as.integer(rshift),
                      rrot=as.integer(rrot),
                      rx=rx[[i]], ry=ry[[i]]),
                      par2=treatment.par[[i]],
                      parcon=constraint.par[[i]],
                      treatment=treatment.func[[i]],
                      constraint=constraint.func[[i]])
        list(xyt=ltraj[[i]], nrep=as.integer(nrep), type=as.integer(0), par = param)
    })

    names(res) <- burst(ltraj)
    class(res) <- c("NMs", "randomShiftRotation")
    attr(res, "warningConstraint") <- warningConstraint
    return(res)
}



NMs.CRW <- function(N=1, nlocs=100, rho=0, h=1, x0=c(0,0), treatment.func=NULL,
                    treatment.par=NULL, constraint.func=NULL,
                    constraint.par=NULL, nrep=999)
{
    ## First checks the arguments
    if (!is.list(nlocs)) {
        if (length(nlocs)!=1)
            stop("nlocs should be of length 1 (same for all) or N (one per animal)")
        nlocs <- sapply(1:N, function(x) nlocs)
    }
    if (!is.list(rho)) {
        if (length(rho)!=1)
            stop("rho should be of length 1 (same for all) or N (one per animal)")
        rho <- sapply(1:N, function(x) rho)
    }
    if (!is.list(h)) {
        if (length(h)!=1)
            stop("h should be of length 1 (same for all) or N (one per animal)")
        h <- sapply(1:N, function(x) h)
    }
    if (!is.list(x0)) {
        if (length(x0)!=2)
            stop("h should be a vector of length 2 (same for all) or a list of length N (one per animal)")
        x0 <- lapply(1:N, function(x) x0)
    }

    if (is.null(treatment.func)) {
        treatment.func <- function(x, par) return(x)
    }
    if (is.null(treatment.par)) {
        treatment.par <- 0
    }
    if (is.null(constraint.func)) {
        constraint.func <- function(x, par) return(TRUE)
    }
    if (is.null(constraint.par)) {
        constraint.par <- 0
    }

    ## convert the functions as lists of functions (same length as N)
    ## if not already a list
    if (!is.list(treatment.func)) {
        treatment.func <- lapply(1:N, function(i) treatment.func)
        treatment.par <- lapply(1:N, function(i) treatment.par)
    }

    if (!is.list(constraint.func)) {
        constraint.func <- lapply(1:N, function(i) constraint.func)
        constraint.par <- lapply(1:N, function(i) constraint.par)
    }
    if (length(treatment.func)!=N)
        stop("treatment.func should have a length = N")
    if (length(treatment.par)!=N)
        stop("treatment.par should have a length = N")
    if (length(constraint.func)!=N)
        stop("constraint.func should have a length = N")
    if (length(constraint.par)!=N)
        stop("constraint.par should have a length = N")

    ## checks the name of the arguments in the function
    lapply(treatment.func, function(x) {
        art <- names(as.list(args(x)))
        if ((!all(art[1:2]==c("x","par")))|(length(art)!=3))
            stop("The treatment function should have only two arguments named x and par")
    })
    lapply(constraint.func, function(x) {
        art <- names(as.list(args(x)))
        if ((!all(art[1:2]==c("x","par")))|(length(art)!=3))
            stop("The treatment function should have only two arguments named x and par")
    })

    ## checks if the treatment works at least on the argument
    res <- lapply(1:N, function(i) {
        so <- simm.crw(1:nlocs[[i]], h[[i]], r=rho[[i]], x0=x0[[i]])[[1]][,1:3]
        ii <- try(treatment.func[[i]](so, treatment.par[[i]]))
        if (inherits(ii, "try-error"))
            stop(paste("The treatment function fails for animal", i,
                       ":\nit should work at least for all observed trajectories"))
    })

    ## Checks whether the constraints are satisfied on the argument
    ## warning otherwise
    warningConstraint <- FALSE

    ## build the object needed for the simulations
    res <- lapply(1:N, function(i) {
        env <- new.env()
        assign("treatment", treatment.func[[i]], envir=env)
        assign("constraint", constraint.func[[i]], envir=env)
        ## par1 contains:
        ## - the length of the trajectory
        ## - the concentration
        ## - the dispersion
        ## - the starting point

        param <- list(par1=list(env=env, n=as.integer(nlocs[[i]]),
                      rho=rho[[i]], h=h[[i]], x0=x0[[i]]),
                      par2=treatment.par[[i]],
                      parcon=constraint.par[[i]],
                      treatment=treatment.func[[i]],
                      constraint=constraint.func[[i]])
        list(xyt=0, nrep=as.integer(nrep), type=as.integer(2), par = param)
    })

    class(res) <- c("NMs", "CRW")
    attr(res, "warningConstraint") <- warningConstraint
    return(res)
}


NMs.randomCs <- function(ltraj, Cs=NULL, rDistCs=TRUE, rAngleCs=TRUE,
                         rCentroidAngle=TRUE, rCs=TRUE,
                         newCs=NULL, newDistances=NULL,
                         treatment.func=NULL, treatment.par=NULL,
                         constraint.func=NULL, constraint.par=NULL,
                         nrep=999)
{
    if (!inherits(ltraj, "ltraj"))
        stop("ltraj should be of class ltraj")
    if ((!is.null(Cs))&((!is.list(Cs))|(length(Cs)!=length(ltraj))))
        stop("Cs should be a list of length equal to ltraj")
    if (is.null(Cs))
        Cs <- lapply(1:length(ltraj), function(i) c(0,0))
    if (any(unlist(lapply(Cs, length))!=2))
        stop(paste("Cs should contain vectors of length 2\n",
                   "storing the x and y coordinates of capture sites"))
    if ((rDistCs+rAngleCs+rCentroidAngle)==0)
        stop("Nothing to randomize: all possibilities are set to FALSE")

    samCs <- Cs
    if (!is.null(newCs)) {
        if (!is.list(newCs)) {
            stop("When not NULL, newCs should be a list")
        }
        if (any(unlist(lapply(newCs, length))!=2))
            stop(paste("newCs should contain vectors of length 2\n",
                       "storing the x and y coordinates of capture sites"))
        samCs <- newCs
    }

    if (is.null(treatment.func)) {
        treatment.func <- function(x, par) return(x)
    }
    if (is.null(treatment.par)) {
        treatment.par <- 0
    }
    if (is.null(constraint.func)) {
        constraint.func <- function(x, par) return(TRUE)
    }
    if (is.null(constraint.par)) {
        constraint.par <- 0
    }

    ## convert the functions as lists of functions (same length as ltraj)
    ## if not already a list
    if (!is.list(treatment.func)) {
        treatment.func <- lapply(1:length(ltraj), function(i) treatment.func)
        treatment.par <- lapply(1:length(ltraj), function(i) treatment.par)
    }

    if (!is.list(constraint.func)) {
        constraint.func <- lapply(1:length(ltraj), function(i) constraint.func)
        constraint.par <- lapply(1:length(ltraj), function(i) constraint.par)
    }
    if (length(treatment.func)!=length(ltraj))
        stop("treatment.func should have the same length as ltraj")
    if (length(treatment.par)!=length(ltraj))
        stop("treatment.par should have the same length as ltraj")
    if (length(constraint.func)!=length(ltraj))
        stop("constraint.func should have the same length as ltraj")
    if (length(constraint.par)!=length(ltraj))
        stop("constraint.par should have the same length as ltraj")

    ## checks the name of the arguments in the function
    lapply(treatment.func, function(x) {
        art <- names(as.list(args(x)))
        if ((!all(art[1:2]==c("x","par")))|(length(art)!=3))
            stop("The treatment function should have only two arguments named x and par")
    })
    lapply(constraint.func, function(x) {
        art <- names(as.list(args(x)))
        if ((!all(art[1:2]==c("x","par")))|(length(art)!=3))
            stop("The treatment function should have only two arguments named x and par")
    })

    ## checks if the treatment works at least on the argument
    res <- lapply(1:length(ltraj), function(i) {
        ii <- try(treatment.func[[i]](ltraj[[i]][,1:3], treatment.par[[i]]))
        if (inherits(ii, "try-error"))
            stop(paste("The treatment function fails for animal", i,
                       ":\nit should work at least for all observed trajectories"))
    })

    ## Checks whether the constraints are satisfied on the argument
    ## warning otherwise
    warningConstraint <- unlist(lapply(1:length(ltraj), function(i) {
        ii <- try(constraint.func[[i]](ltraj[[i]][,1:3], constraint.par[[i]]))
        if (inherits(ii, "try-error"))
            stop(paste("The constraint function fails for animal", i,
                       ":\nit should work at least for all observed trajectories"))
        if (!is.logical(ii)|(length(ii)>1))
            stop("the constraint function should return a logical value of length 1")
        return(ii)
    }))
    if (any(!warningConstraint))
        warning(paste("The current ltraj does not satisfy the constraint for bursts:\n",
                      paste(burst(ltraj)[!warningConstraint], collapse=",")))

    ## vector of distances to the capture site
    if (is.null(newDistances)) {
        dis <- unlist(lapply(1:length(ltraj), function(i) {
            sqrt(sum((apply(ltraj[[i]][,1:2],2,mean)-Cs[[i]])^2))
        }))
    } else {
        dis <- newDistances
        if (!is.vector(newDistances)) {
            stop("when not NULL, newDistances should be a vector")
        }
    }

    ## build the object needed for the simulations
    res <- lapply(1:length(ltraj), function(i) {
        env <- new.env()
        assign("treatment", treatment.func[[i]], envir=env)
        assign("constraint", constraint.func[[i]], envir=env)
        ## par1 contains:
        ## - the environment
        ## - rDistCs
        ## - rAngleCs
        ## - rCentroidAngle
        ## - the capture site coordinates
        ## - the vector of distance
        ## - whether the capture site should be sampled
        ## - sampled capture sites
        param <- list(par1=list(env=env, rDistCs=as.integer(rDistCs),
                      rAngleCs=as.integer(rAngleCs),
                      rCentroidAngle=as.integer(rCentroidAngle), cs=Cs[[i]],
                      distances=dis, rCs=as.integer(rCs), samCs=samCs),
                      par2=treatment.par[[i]],
                      parcon=constraint.par[[i]],
                      treatment=treatment.func[[i]],
                      constraint=constraint.func[[i]])
        list(xyt=ltraj[[i]], nrep=as.integer(nrep), type=as.integer(3), par = param)
    })

    names(res) <- burst(ltraj)
    class(res) <- c("NMs", "randomCs")
    attr(res, "warningConstraint") <- warningConstraint
    return(res)

}



## We define a list of null models
NMs2NMm <- function(NMs, treatment.func=NULL,
                    treatment.par=NULL, constraint.func=NULL,
                    constraint.par=NULL, nrep=999)
{
    if (!inherits(NMs, "NMs"))
        stop("NMs should be of class NMs")
    if (length(NMs)==1)
        stop("Only one individual in the object NMS: no need to convert to NMm")

    ## Checks the existence of a treatment and constraint function
    if (is.null(treatment.func)) {
        treatment.func <- function(x, par) return(x)
    }
    if (is.null(treatment.par)) {
        treatment.par <- 0
    }
    if (is.null(constraint.func)) {
        constraint.func <- function(x, par) return(TRUE)
    }
    if (is.null(constraint.par)) {
        constraint.par <- 0
    }

    lixyt <- lapply(NMs, function(x) x$xyt)
    litype <- lapply(NMs, function(x) x$type)
    lipar <- lapply(NMs, function(x) x$par)

    ## Checks that the treatment function works correctly on the arguments
    ltrajdf <- data.frame(do.call("rbind", lixyt)[,1:3],
                          id=rep(names(NMs), unlist(lapply(lixyt, nrow))))
    ltraj <- as.ltraj(ltrajdf[,1:2], ltrajdf[,3], id=ltrajdf[,4])
    ii <- try(treatment.func(ltraj, treatment.par))
    if (inherits(ii, "try-error"))
        stop(paste("The treatment function failed",
                   ":\nit should work at least on the observed data set"))

    res <- list(lixyt=lixyt, litype=litype, lipar=lipar, nrep=nrep,
                constraint.func=constraint.func, constraint.par=constraint.par,
                treatment.func=treatment.func, treatment.par=treatment.par)
    class(res) <- c("NMm", "NMs2NMm")
    attr(res, "warningConstraint") <- attr(NMs, "warningConstraint")
    return(res)
}


print.NMm <- function(x, ...)
{
    if (!inherits(x, "NMm"))
        stop("x should be of class NMm")
    cat("*********************************\n")
    cat("Null Model object of type", class(x)[2]," (multiple)\n")
    cat(x[[1]]$nrep, " repetitions will be carried out\n")
    cat("Please consider the function testNM() for the simulations\n")
    if (any(!attr(x, "warningConstraint")))
        cat("The current object does not satisfy the constraint function\n")
}








print.NMs <- function(x, ...)
{
    if (!inherits(x, "NMs"))
        stop("x should be of class NMs")
    cat("*********************************\n")
    cat("Null Model object of type", class(x)[2]," (single)\n")
    cat(x[[1]]$nrep, " repetitions will be carried out\n")
    cat("Please consider the function testNM() for the simulations\n")
    if (any(!attr(x, "warningConstraint")))
        cat("The current object does not satisfy the constraint function\n")
}

testNM <- function(NM, count=TRUE)
{
    if ((!inherits(NM, "NMs"))&(!inherits(NM, "NMm")))
        stop("NM should be of class NMs or NMm")
    if (inherits(NM, "NMs")) {
        res <- lapply(NM, function(x) {
            x$par$treatment <- quote(treatment(x,par))
            x$par$constraint <- quote(constraint(x,par))
            .Call("simulmod", x$xyt, x$nrep, x$type,
                  x$par, as.integer(count), PACKAGE="adehabitatLT")
        })
    } else {
        lixyt <- NM$lixyt
        litype <- lapply(NM$litype, function(x) as.integer(x))
        lipar <- NM$lipar
        tr <- function(x, par) return(x)
        for (i in 1:length(lipar)) {
            con <- lipar[[i]]$constraint
            assign("treatment", tr, envir=lipar[[i]]$par1$env)
            assign("constraint", con, envir=lipar[[i]]$par1$env)
            lipar[[i]]$treatment <- quote(treatment(x,par))
            lipar[[i]]$constraint <- quote(constraint(x,par))
        }
        nrep <- NM$nrep
        treatment <- NM$treatment.func
        constraint <- NM$constraint.func
        tpar <- NM$treatment.par
        cpar <- NM$constraint.par

        env2 <- new.env()

        convltraj <- function(lixyt, na, nlo) {
            ltrajdf <- data.frame(do.call("rbind", lixyt)[,1:3],
                                  id=rep(na, nlo))
            as.ltraj(ltrajdf[,1:2], ltrajdf[,3], id=ltrajdf[,4])
        }
        na <- names(lixyt)
        nlo <- unlist(lapply(lixyt, nrow))
        assign("convltraj", convltraj, envir=env2)
        res <- .Call("simulmodmv", lixyt, as.integer(nrep), litype, lipar,
                     as.integer(count),
                     env2, quote(convltraj(lixyt, na, nlo)), na, nlo,
                     quote(treatment(x, par)), tpar, quote(constraint(x,par)),
                     cpar, PACKAGE="adehabitatLT")
    }
    return(res)
}
