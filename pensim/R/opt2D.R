opt2D <-
    function(nsim,L1range=c(0.001,100),L2range=c(0.001,100),dofirst="both",nprocessors=1,L1gridsize=10,L2gridsize=10,cl=NULL,...){
        opt.L1L2 <- function(lambdarange,...){
            minL1 <- min(lambdarange[1:2])
            maxL1 <- max(lambdarange[1:2])
            minL2 <- min(lambdarange[3:4])
            maxL2 <- max(lambdarange[3:4])
            opt1 <- try(optL1(lambda2=minL2,minlambda1=minL1,...)) #testcode
            opt2 <- try(optL2(lambda1=opt1$lambda,...)) #testcode
            if(class(opt1)!="try-error"&class(opt2)!="try-error"){
                L1=opt1$lambda
                L2=opt2$lambda
                cvall <- opt2$cvl
                is.standardized <- list(...)$standardize
                cc <- coefficients(opt2$fullfit,"all",standardize=is.standardized)
                output <- c(L1,L2,cvall,cc)
                names(output)[1:3] <- c("L1","L2","cvl")
                return(output)
            } else{
                return(rep(NA,3+ncol(list(...)$penalized)))
            }      
        }
        opt.L2L1 <- function(lambdarange,...){
            minL1 <- min(lambdarange[1:2])
            maxL1 <- max(lambdarange[1:2])
            minL2 <- min(lambdarange[3:4])
            maxL2 <- max(lambdarange[3:4])
            opt2 <- try(optL2(lambda1=minL1,minlambda2=minL2,...)) #testcode
            opt1 <- try(optL1(lambda2=opt2$lambda,...)) #testcode
            if(class(opt1)!="try-error"&class(opt2)!="try-error"){
                L2=opt2$lambda
                L1=opt1$lambda
                cvall <- opt1$cvl
                is.standardized <- list(...)$standardize
                cc <- coefficients(opt1$fullfit,"all",standardize=is.standardized)
                output <- c(L1,L2,cvall,cc)
                names(output)[1:3] <- c("L1","L2","cvl")
                return(output)
            } else{
                return(rep(NA,3+ncol(list(...)$penalized)))
            }
        }
        opt.both <- function(lambdarange,...){
            ##lambdarange[1:2] to be boundaries for lambda1,
            ##lambdarange[3:4] to be boundaries for lambda2
            ##lambdarange[5:6] to be the starting point for optimization,
            getfolds <- function(N,nfolds){
                evenly <- nfolds * N%/%nfolds
                extras <- N - evenly
                groupsize <- evenly/nfolds
                evenly.folds <- rep(1:nfolds,groupsize)
                evenly.folds <- sample(evenly.folds,size=length(evenly.folds))
                extra.folds <- sample(1:nfolds,size=extras)
                folds <- c(evenly.folds,extra.folds)
                return(folds)
            }
            myargs <- list(...)
            if(is.null(myargs$fold)){
                myargs$fold <- 1:nrow(myargs$penalized)
                warning("fold argument missing, assuming leave-one-out cross-validation")
            }
            if(all.equal(length(myargs$fold),1)){
                myfolds <- getfolds(N=nrow(myargs$penalized),nfolds=myargs$fold)
                myargs$fold <- myfolds
            }
            if(!all.equal(length(myargs$fold),nrow(myargs$penalized))){
                stop("fold must be a single integer or a string of length equal to the number of samples.")
            }
            if(is.null(myargs$standardize)){
                myargs$standardize <- FALSE
                warning("standardize not specified, so no standardization performed (standardize=FALSE for cvl() function")
            }
            cvl.simple <- function(x1,x2,...){
                trycvl <- try(cvl(lambda1=x1,lambda2=x2,...))
                if(class(trycvl)=="try-error"){
                    print(paste("warning: try-error, faking cvl =",-10000))
                    return(list(cvl=-10000))
                }else{
                    if(trycvl$cvl < (-1e12)){
                        print(paste("warning: cvl set to -1e12, was",trycvl$cvl))
                        trycvl$cvl <- (-1e12)
                    }
                    return(trycvl)
                }
            }
            fun.fit <- function(x,myargs){do.call(cvl.simple,c(x,myargs))}
            fun.cvlonly <- function(x,myargs){do.call(cvl.simple,c(x,myargs))$cvl}
            opt.out <- try(optim(par=lambdarange[5:6],
                                 fn=fun.cvlonly,
                                 method="L-BFGS-B",
                                 lower=c(lambdarange[1],lambdarange[3]),
                                 upper=c(lambdarange[2],lambdarange[4]),
                                 myargs=myargs,
                                 ##factr=1e12 for very fast convergence -  1e10 or 11 for very accurate
                                 control=list(fnscale=-1,maxit=50,factr=1e11,trace=0)))  #trace=6 for full info
            if(class(opt.out)=="try-error"){
                output <- rep(NA,nrow(myargs$penalized)+5)
                names(output) <- c("L1","L2","cvl","convergence","fncalls",rownames(myargs$penalized))
            }else{
                L1 <- opt.out$par[1]
                L2 <- opt.out$par[2]
                cvall <- opt.out$value
                convergence <- opt.out$convergence
                fncalls <- opt.out$counts[1]
                fit <- fun.fit(opt.out$par,myargs)$fullfit
                cc <- coefficients(fit,"all",standardize=myargs$standardize)
                output <- c(opt.out$par,opt.out$value,opt.out$convergence,opt.out$counts[1],cc)
                names(output)[1:5] <- c("L1","L2","cvl","convergence","fncalls")
            }
            return(output)
        }
        ##the next line seems like a complicated way just to copy L1range and L2range into nsim elements of the list "looplist", which will be used to replicate the simulations:
        looplist <- lapply(1:nsim,function(x,...){return(as.numeric(...))},c(L1range,L2range))
        ##now do the optimization, whichever method was chosen.
        if(identical(dofirst,"L1")){
            FUN <- opt.L1L2
        } else if(identical(dofirst,"L2")){
            FUN <- opt.L2L1
        } else if(identical(dofirst,"both")){
            FUN <- opt.both
            print("scanning for good starting points for 2-D optimization...")
            ##Use a clipped lambdarange for scanning, scanning with very small values of lambda1
            ##and lambda2 can be *very* slow.
            L1range.clipped <- L1range
            L2range.clipped <- L2range
            if(min(L1range)<1){
                L1range.clipped[which.min(L1range.clipped)] <- 1
            }
            if(min(L2range)<1){
                L2range.clipped[which.min(L2range.clipped)] <- 1
            }
            print("Scanning L1 and L2 for good start points for 2D optimization.")
            print("Note: not scanning values of L1 or L2 less than 1 to speed computation, but optimization can still converge on values less than 1 if the minimum range is less than 1.")
            print(paste("Scanning L1 between",paste(L1range.clipped,collapse="-")))
            print(paste("Scanning L2 between",paste(L2range.clipped,collapse="-")))
            cvlmatrix <- scan.l1l2(L1range=L1range.clipped,L2range=L2range.clipped,L1.ngrid=L1gridsize,L2.ngrid=L2gridsize,nprocessors=nprocessors,polydegree=1,...)$cvl
            ranking <- matrix(rank(-cvlmatrix,ties.method="random"),ncol=ncol(cvlmatrix))
            top.positions <- lapply(1:5,function(n) which(ranking==n,arr.ind=TRUE))
            top.positions <- lapply(top.positions,function(x) as.numeric(c(L1range,L2range,rownames(cvlmatrix)[x[1]],colnames(cvlmatrix)[x[2]])))
            startpositions <- t(sapply(top.positions,function(x) x[5:6]))
            colnames(startpositions) <- c("lambda1","lambda2")
            print("done scanning.  The following start positions will be used:")
            print(startpositions)
            ##randomly sample the top 5 positions to provide starting points for
            ##the nsim optimizations
            looplist <- sample(top.positions,nsim,replace=TRUE)
        }else stop("dofirst must be one of L1, L2, or both")
        clusterIsSet <- "cluster" %in% class(cl)
        if(nprocessors>1 | clusterIsSet){
            if(!clusterIsSet){
                nprocessors <- as.integer(round(nprocessors))
                cl <- makeCluster(nprocessors, type="PSOCK")
            }
            clusterSetRNGStream(cl, iseed=NULL)
            if(nprocessors > 1)
                print(paste("beginning simulations on", nprocessors, "processors..."))
            output <- parSapply(cl,looplist,function(x,...){FUN(x,...)},...)
            if(!clusterIsSet){
                stopCluster(cl)
            }
        } else{
            print("beginning simulations on one processor...")
            output <- sapply(looplist,function(x,...){FUN(x,...)},...)
        }
        print("finished simulations.")
        output <- t(output)
        return(output)
    }

