###################################################################
#                                                                 #
# Diffusion Kurtosis Imaging (DKI)                                #
#                                                                 # 
# Literature: Jensen, J. H. et al. (2005), MRM 53, 1432-1440      #
#             Jensen, J. H. et al. (2010), NMR Biomed 23, 698-710 #
#             Tabesh, A. et al. (2011), MRM 65, 823-836           #
#             Hui et al. (2008), NI 42, 122-134                   #
#                                                                 #
###################################################################

dkiTensor <- function(object,  ...) cat("No DKI tensor calculation defined for this class:", class( object), "\n")

setGeneric("dkiTensor", function(object, ...) standardGeneric("dkiTensor"))

setMethod("dkiTensor", "dtiData",
          function(object, method=c("CLLS-QP", "CLLS-H", "ULLS") , 
                   mc.cores=setCores(, reprt = FALSE),
                   verbose=FALSE) {
            
            if (verbose) cat("dkiTensor: entering function", format(Sys.time()), "\n")
            
            ## check method
            method <- match.arg(method)
            
            ## call history
            args <- c(object@call, sys.call(-1))
            
            ## check available number of cores for parallel computation
            if(is.null(mc.cores)) mc.cores <- 1
            mc.cores <- min(mc.cores, detectCores())
            
            ## define the design matrix of the estimation problem
            s0ind <- object@s0ind
            xxx <- dkiDesign(object@gradient[, - s0ind])
            
            ## container for estimation results
            ddim <- object@ddim  
            ntotal <- prod(ddim)
            D <- array(0, dim=c(6, ntotal))
            W <- array(0, dim=c(15, ntotal))
            bvalues <- object@bvalue[- s0ind]  
            meanbv <- mean(bvalues)
            bvalues <- bvalues/meanbv
            maxbv <- max(bvalues)
            
            ## check for outliers and
            ## we need the mean s0 image and a mask
            z <- sioutlier1(object@si,s0ind,object@level,mc.cores=mc.cores)
            ## z$si and z$s0 only contain voxel in the mask
            ## first dimension of matrix z$si corresponds to gradients 
            cat("sioutlier completed\n")
            mask <- z$mask
            nvox <- sum(mask)
            ttt <- -log1p(sweep(z$si[-s0ind,],2,as.vector(z$s0),"/")-1)
            #  suggestion by B. Ripley
            ttt[is.na(ttt)] <- 0
            ttt[(ttt == Inf)] <- 0
            ttt[(ttt == -Inf)] <- 0
            s0 <- array(0,ddim)
            s0[mask] <- z$s0
            maskind <- (1:ntotal)[mask] 
            cat("Data transformation completed ",format(Sys.time()),"\n")
            
            if (method == "CLLS-QP"||method == "NLR"||method == "QL") {
              
              #if (!require(quadprog)) return("dkiTensor: did not find package quadprog, please install for the CLLS-QP method")
              
              ## Tabesh Eq. [11, 14, 15]
              Tabesh_AD <- xxx[, 1:6]
              Tabesh_AK <- xxx[, 7:21]
              
              ## Tabesh Eq. [10]
              Tabesh_A <- cbind(sweep(Tabesh_AD, 1, - bvalues, "*"),
                                sweep(Tabesh_AK, 1, bvalues^2/6, "*"))
              ##
              ##   old variant gave condition number of 3268 for Tabesh_A instead of 35 
              ##                                     10684619 for Dmat instead of 1263
              ##   for Siawooshs MQ01286 Data
              ##
              ngrad <- object@ngrad - length(s0ind)
              
              ## Tabesh Eq. [13]
              Tabesh_C <- rbind( cbind(Tabesh_AD, matrix(0, ngrad, 15)), 
                                 cbind(matrix(0, ngrad, 6), Tabesh_AK), 
                                 cbind(3/maxbv*Tabesh_AD, - Tabesh_AK))
              
              Dmat <- t(Tabesh_A) %*% Tabesh_A
              Amat <- t(Tabesh_C)
              
              ## go through the dataset
              Dind <- c(1, 4, 5, 2, 6, 3)
              Dinvind <- c(1, 4, 6, 2, 3, 5)
              if(mc.cores==1){
                ## single core, just do the loop 
                for (i in 1:nvox){
                  mi <- maskind[i]
                  ## Tabesh Eq. [12]
                  Tabesh_B <- -ttt[,i]
                  ## TODO: PERFORM SOME TEST BEFORE IT!!
                  
                  ## QP solution
                  dvec <- as.vector(t(Tabesh_A) %*% Tabesh_B)
                  resQPsolution <- solve.QP(Dmat, dvec, Amat)$solution
                  D[, mi] <- meanbv * resQPsolution[Dind] # re-order DT estimate to comply with dtiTensor
                  
                  ## Tabesh Eq. [9]
                  W[, mi] <- resQPsolution[7:21] / mean(resQPsolution[1:3])^2 # kurtosis tensor
                }
              } else {
                ## many cores, just split
                param <- plmatrix(ttt,pdkiQP,TA=Tabesh_A,Dmat=Dmat,Amat=Amat)
                D[,mask] <-  meanbv * param[Dind,]
                W[,mask] <- param[7:21,]
              }
            }
#            if (method == "NLR" || method == "QL"){
#              if (!require(Rsolnp)) return("dkiTensor: did not find package Rsolnp, please install for the NLR method")
#              ttt <- sweep(z$si[-s0ind,],2,as.vector(z$s0),"/")
#              funopt <- function(param,siq,Tabesh_A,Tabesh_C,rho){
#                z <- Tabesh_C%*%param
#                sum((siq-exp(Tabesh_A%*%param))^2) + rho*sum(((z-1e-5)[z<1e-5])^2)
#              }
#              fun <- function(param,siq,Tabesh_A,Tabesh_C){
#                sum((siq-exp(Tabesh_A%*%param))^2)
#              }
#              ineqfun <- function(param,siq,Tabesh_A,Tabesh_C){
#                Tabesh_C%*%param
#              }
#              lb <- rep(0,3*ngrad)
#              ub <- rep(25,3*ngrad)
#              for(i in 1:nvox){
#                mi <- maskind[i]
#                param <- c(D[Dinvind, mi]/meanbv,W[, mi]*mean(D[c(1,4,6),mi]/meanbv)^2)
#                siq <- ttt[,i]
#                cat("i",i,"mi",mi,"\n","param",signif(param,3),"\n","fw",fun(param,siq,Tabesh_A,Tabesh_C))
#                lconstr <- TRUE
#                for(k in 1:length(rho)){
#                  if(lconstr){
#                    param <- optim(param, funopt, method="BFGS", siq=siq,Tabesh_A=Tabesh_A,Tabesh_C=Tabesh_C,rho=rho[k])$par
#                    lconstr <- min(ineqfun(param,siq,Tabesh_A,Tabesh_C)-lb) < 0
#                  }
#                }
#                cat("fw optim",fun(param,siq,Tabesh_A,Tabesh_C), "constr", constr <-min(ineqfun(param,siq,Tabesh_A,Tabesh_C)-lb),"\n","param",signif(param,3),"\n")
#                if(constr< 0){
#                  cat(format(Sys.time()),"\n")
#                  solnpres <- solnp(param, fun = fun, ineqfun = ineqfun, ineqLB = lb, ineqUB = ub, 
#                                    control = control, #list(tol=1e-6,delta=1e-5,rho=.01), 
#                                    siq=siq,Tabesh_A=Tabesh_A,Tabesh_C=Tabesh_C)
#                  param <- solnpres$pars
#                  cat("fw solnp",fun(param,siq,Tabesh_A,Tabesh_C), "constr",                   min(ineqfun(param,siq,Tabesh_A,Tabesh_C)-lb),"\n","param",signif(param,3),"\n")
#                  cat("diagnostics: vals",solnpres$values,"nfunc",solnpres$nfuneval,
#                      "time",solnpres$elapsed,"\n")
#                  cat(format(Sys.time()),"\n")
#                  
#                }
#                D[, mi] <- meanbv * param[Dind] # re-order DT estimate to comply with dtiTensor
#                W[, mi] <- param[7:21] / mean(param[1:3])^2 #  
#              }
#            }
#            if(method=="QL"){
#              CL <- sqrt(pi/2)*gamma(L+1/2)/gamma(L)/gamma(3/2)
#              funoptQL <- function(param, siq, sigma, L, CL, Tabesh_A, Tabesh_C, rho){
#                z <- Tabesh_C%*%param
#                fv <- exp(Tabesh_A%*%param)
#                sum((siq-exp(Tabesh_A%*%param))^2) + rho*sum(((z-1e-5)[z<1e-5])^2)
#              }
#              funQL <- function(param, siq, sigma, L, CL, Tabesh_A, Tabesh_C){
#                sum((siq-exp(Tabesh_A%*%param))^2)
#              }
#              ineqfunQL <- function(param, siq, sigma, L, CL, Tabesh_A, Tabesh_C){
#                Tabesh_C%*%param
#              }
#              lb <- rep(0,3*ngrad)
#              ub <- rep(25,3*ngrad)
#              
#            }
            if (method == "CLLS-H") {
              
              ## these are the distinct bvalues
              bv <- unique(bvalues)
              bv <- bv[order(bv)]
              if (length(bv) != 2) stop("Need exactly two shells for CLLS-H method, choose other method!")
              ## now we have bv[ 1] < bv[ 2] 
              
              indbv1 <- (1:length(bvalues))[bvalues == bv[ 1]] # smaller b-value
              indbv2 <- (1:length(bvalues))[bvalues == bv[ 2]] # larger b-value
              if ((ngrad <- length(indbv1)) != length(indbv2)) stop("Need same number of gradients vectors for both shells!")
              
              ## we must have two shell data and the gradient direction coincide for both shells
              gradient <- object@gradient[,-object@s0ind]
              gradtest <- abs(range(diag(t(gradient[, indbv1]) %*% gradient[, indbv2]) - 1))
              if (max(gradtest) > 1e-6) warning("dkiTensor: Gradient directions on the two shells are not identical (this might be a numerical problem).\n Proceed, but strange things may happen!")
              
              ## We need only a reduced design, as the gradient directions coincide
              xxx <- dkiDesign(gradient[, indbv1])
              Tabesh_AD <- xxx[, 1:6]
              Tabesh_AK <- xxx[, 7:21]
              PI_Tabesh_AD <- pseudoinverseSVD(Tabesh_AD)
              PI_Tabesh_AK <- pseudoinverseSVD(Tabesh_AK)
              
              ## some hard-coded constraints, see Tabesh Eq. [6] ff. 
              Kmin <- 0
              ## Tabesh Eq. [6]
              C <- 3
              ## TODO: Do we want this as user defined arguments?
              ## for all voxel
              ## TODO: Make this more efficient matrix operation!
              for(i in 1:nvox){
                mi <- maskind[i]
                ## Tabesh Eq. [20]
                D1 <- ttt[indbv1,i]/bv[1]
                D2 <- ttt[indbv2,i]/bv[2]
                ## Tabesh Eq. [18]
                Di <- (bv[2] * D1 - bv[1] * D2) / (bv[2] - bv[1])
                
                ## Tabesh Eq. [19]
                Ki <-  6 * (D1 - D2) / (bv[2] - bv[1]) / Di^2 
                
                ## now we apply some constraints Tabesh
                ## D1
                ind1 <- (Di <= 0)
                Di[ind1] <- 0
                ## D2
                ind2 <- (!ind1 & (D1 < 0))
                Di[ind2] <- 0
                ## D3
                ind3 <- (!ind1 & !ind2 & (Di > 0) & (Ki < Kmin))
                if (Kmin == 0) {
                  Di[ind3] = D1[ind3]
                } else {
                  x <- - Kmin * bv[1] / 3 
                  Di[ind3] = (sqrt(1 + 2 * x * D1[ind3]) - 1) / x
                }
                ## D4
                Kmax <- numeric(length(Di))
                dim(Kmax) <- dim(Di)
                Kmax[Di > 0] <- C / bv[2] / Di[Di > 0]
                ind4 <- (!ind1 & !ind2 & !ind3 & (Di > 0) & (Ki > Kmax))
                Di[ind4] <- D1[ind4] / (1 - C * bv[1] / 6 / bv[2])
                
                ## estimate D tensor! Tabesh Eq. [21]
                Dihat <- meanbv * PI_Tabesh_AD %*% Di
                D[c(1, 4, 6, 2, 3, 5), mi] <- meanbv * Dihat 
                
                ## re-estimate Di: Tabesh Eq. [22]
                DiR <- Tabesh_AD %*% Dihat
                
                ## constraints for the KT estimation step
                KiR = numeric(length(Ki))
                ## Tabesh Eq. [23]
                ind <- (DiR > 0)
                KiR[ind] <- 6 * (DiR[ind] - D2[ind]) / bv[2] / DiR[ind]^2
                
                ## K1 and K2
                Kmax <- numeric(length(DiR))
                dim(Kmax) <- dim(DiR)# DiR is a vector 
                Kmax[DiR > 0] <- C / bv[2] / DiR[DiR > 0]
                ind <- (KiR > Kmax)
                KiR[ind] <- Kmax[ind]
                ind <- (KiR < Kmin)
                KiR[ind] <- Kmin
                
                
                ## estimate KT Tabesh Eq. [24]
                W[, mi] <- PI_Tabesh_AK %*% (DiR^2 * KiR) / (mean(Dihat[1:3]))^2
                
              }
              #              if (verbose) close(pb)  
            }  
            if (method == "ULLS") {
              
              ## Tabesh Eq. [11, 14, 15]
              Tabesh_AD <- xxx[, 1:6]
              Tabesh_AK <- xxx[, 7:21]
              
              ## Tabesh Eq. [10]
              Tabesh_A <- cbind(sweep(Tabesh_AD, 1, bvalues, "*"),
                                sweep(Tabesh_AK, 1, bvalues^2/6, "*"))
              PI_Tabesh_A <- pseudoinverseSVD(Tabesh_A)
              #
              #   this can be written in compact form !!
              #
              ## go through the dataset
              for (i in 1:nvox){
                mi <- maskind[i]
                
                ## TODO: make this for all voxel
                Tabesh_B <- -ttt[,i] 
                
                ## TODO: correct assignment! Check matrix dimensions!
                Tabesh_X <- PI_Tabesh_A %*% Tabesh_B
                D[c(1, 4, 6, 2, 3, 5), mi] <- meanbv * Tabesh_X[1:6]
                W[, mi] <- Tabesh_X[7:21] / (mean(Tabesh_X[1:3]))^2
              }
              #              if (verbose) close(pb)  
              
            }
            
            dim(D) <- c(6,ddim)
            dim(W) <- c(15,ddim)
            
            if (verbose) cat("dkiTensor: finished estimation", format(Sys.time()), "\n")
            
            ## TODO: CHECK THESE!
            sigma2 <- array(0, dim=ddim)
            scorr <- array(0, dim=c(3, 3, 3))
            bw <- c(0, 0, 0)
            index <- 1
            scale <- 1
            ## do we need this?
            
            if (verbose) cat("dkiTensor: exiting function", format(Sys.time()), "\n")	 	  
            
            invisible(new("dkiTensor",
                          call  = args,
                          D     = D,
                          W     = W,
                          th0   = s0,
                          sigma = sigma2,
                          scorr = scorr, 
                          bw    = bw, 
                          mask  = mask,
                          hmax  = 1,
                          gradient = object@gradient,
                          bvalue = object@bvalue,
                          btb   = xxx,
                          ngrad = object@ngrad,
                          s0ind = object@s0ind,
                          replind = object@replind,
                          ddim  = object@ddim,
                          ddim0 = object@ddim0,
                          xind  = object@xind,
                          yind  = object@yind,
                          zind  = object@zind,
                          voxelext = object@voxelext,
                          level = object@level,
                          orientation = object@orientation,
                          rotation = object@rotation,
                          source = object@source,
                          outlier = index,
                          scale = scale,
                          method = method)
            )
          })


dkiIndices <- function(object, ...) cat( "No DKI indices calculation defined for this class:", class( object), "\n")

setGeneric( "dkiIndices", function( object, ...) standardGeneric( "dkiIndices"))

setMethod("dkiIndices", "dkiTensor",
          function(object, 
                   mc.cores = setCores(, reprt=FALSE),
                   verbose = FALSE) {
            
            if (verbose) cat("dkiTensor: entering function", format(Sys.time()), "\n")
            
            ## call history  
            args <- c(object@call, sys.call(-1))
            
            ## we need this for all the arrays
            ddim <- object@ddim
            nvox <- prod(ddim)
            nvox0 <- sum(object@mask)
            
            ## perform the DTI indices calculations
            D <- object@D
            W <- object@W
            z <- dtiind3D(object@D, object@mask, mc.cores=mc.cores, verbose=verbose)
            
            if (verbose) cat("dkiTensor: DTI indices calculated", format(Sys.time()), "\n")
            
            ## perform the DKI indices calculations
            ## DO WE NEED THIS?
            if (mc.cores > 1) {
              mc.cores.old <- setCores(, reprt = FALSE)
              setCores(mc.cores)
            }
            
            dim(D) <- c(6, nvox)
            dim(W) <- c(15, nvox)
            andir <- matrix(0, 9, nvox)
            lambda <- matrix(1e20, 3, nvox)
            
            t1 <- Sys.time()
            
            zz <- .Fortran("dti3DevAll",
                           as.double(D[ , object@mask]),
                           as.integer(nvox0),
                           andir=double(9*nvox0),
                           evalues=double(3*nvox0),
                           PACKAGE="dti")[c("andir", "evalues")]
            andir[, object@mask] <- zz$andir
            lambda[, object@mask] <- zz$evalues 
            
            t2 <- Sys.time()
            if (verbose) cat("dkiTensor: calculation took ", difftime(t2, t1), attr(difftime( t2, t1), "units"), " for", nvox0, "voxel\n")
            
            if (mc.cores > 1) setCores(mc.cores.old, reprt=FALSE)
            dim(andir) <- c(3, 3, nvox)
            
            xxx <- dkiDesign(object@gradient[, - object@s0ind])
            Tabesh_AD <- xxx[, 1:6]
            Tabesh_AK <- xxx[, 7:21]
            
            ## Mean Kurtosis definition following Hui et al. (2008), this is only an approximation
            
            ## Note: the DT entries in D are re-ordered compared to Tabesh_AD for backward comp.
            ## Maybe we dont need this!
            Dapp <- Tabesh_AD %*% D[c(1, 4, 6, 2, 3, 5), object@mask]
            MD2 <- apply(D[c(1, 4, 6), object@mask], 2, mean)^2
            Kapp <- sweep((Tabesh_AK %*% W[, object@mask]) / Dapp^2, 2, MD2, "*")
            ## remove pathological values!
            Kapp[Kapp < 0] <- 0
            Kapp[Dapp <= 0] <- 0
            ## define mean kurtosis
            mk <- array(0, ddim)
            mk[object@mask] <- apply(Kapp, 2, mean)
            
            ## START 
            ## Mean Kurtosis definition following Tabesh et al. (2011), this should be exact
            ## Note, these values already contain MD^2 as factor!
            ## Tabesh Eq. [26] needed for Tabesh Eq. [25]
            Wtilde1111 <- rotateKurtosis(andir[, , object@mask], W[, object@mask], 1, 1, 1, 1) 
            Wtilde2222 <- rotateKurtosis(andir[, , object@mask], W[, object@mask], 2, 2, 2, 2) 
            Wtilde3333 <- rotateKurtosis(andir[, , object@mask], W[, object@mask], 3, 3, 3, 3) 
            Wtilde1122 <- rotateKurtosis(andir[, , object@mask], W[, object@mask], 1, 1, 2, 2) 
            Wtilde1133 <- rotateKurtosis(andir[, , object@mask], W[, object@mask], 1, 1, 3, 3) 
            Wtilde2233 <- rotateKurtosis(andir[, , object@mask], W[, object@mask], 2, 2, 3, 3) 
            
            ## Tabesh Eq. [25]
            mk2 <- array(0, ddim)
            mk2[ object@mask] <- 
              kurtosisFunctionF1(lambda[1, object@mask], lambda[2, object@mask], lambda[3, object@mask]) * Wtilde1111 +
              kurtosisFunctionF1(lambda[2, object@mask], lambda[1, object@mask], lambda[3, object@mask]) * Wtilde2222 +
              kurtosisFunctionF1(lambda[3, object@mask], lambda[2, object@mask], lambda[1, object@mask]) * Wtilde3333 +
              kurtosisFunctionF2(lambda[1, object@mask], lambda[2, object@mask], lambda[3, object@mask]) * Wtilde2233 +
              kurtosisFunctionF2(lambda[2, object@mask], lambda[1, object@mask], lambda[3, object@mask]) * Wtilde1133 +
              kurtosisFunctionF2(lambda[3, object@mask], lambda[2, object@mask], lambda[1, object@mask]) * Wtilde1122
            ## END
            
            
            ## START
            ## Fractional kurtosis following Hui et al. (2008), this might be not useful!
            x1 <- dkiDesign(andir[, 1, ])
            x2 <- dkiDesign(andir[, 2, ])
            x3 <- dkiDesign(andir[, 3, ])
            
            ## cannot allocate memory for the following:
            ## w1111 <- diag( x1[ , 7:21] %*% D[ 7:21, ])
            w1111 <- (x1[, 7:21] * t(W)) %*% rep(1, 15)
            #            w1111 <- numeric( nvox)
            #            for ( i in 1:nvox) if (object@mask[ i]) w1111[ i] <- x1[ i, 7:21] %*% W[ , i]
            w2222 <- (x2[, 7:21] * t(W)) %*% rep(1, 15)
            #            w2222 <- numeric( nvox)
            #            for ( i in 1:nvox) if (object@mask[ i]) w2222[ i] <- x2[ i, 7:21] %*% W[ , i]
            w3333 <- (x3[, 7:21] * t(W)) %*% rep( 1, 15)
            #            w3333 <- numeric( nvox)
            #            for ( i in 1:nvox) if (object@mask[ i]) w3333[ i] <- x3[ i, 7:21] %*% W[ , i]
            
            k1 <- w1111 / lambda[1, ]^2
            k2 <- w2222 / lambda[2, ]^2
            k3 <- w3333 / lambda[3, ]^2
            
            kbar <- (k1 + k2 + k3) / 3
            kaxial <- k1
            kradial <- (k2 + k3) / 2
            fak <- sqrt(3/2 * ((k1-kbar)^2 + (k2-kbar)^2 + (k3-kbar)^2) / (k1^2 + k2^2 + k3^2))
            ## END
            
            
            ## we finally got k1, k2, k3, kaxial, kradial, mk, fak
            dim(k1) <- ddim
            dim(k2) <- ddim
            dim(k3) <- ddim
            dim(mk) <- ddim
            dim(mk2) <- ddim
            dim(kaxial) <- ddim
            dim(kradial) <- ddim
            dim(fak) <- ddim
            
            if (verbose) cat("dkiTensor: DKI indices calculated", format(Sys.time()), "\n")
            
            if (verbose) cat("dkiTensor: exiting function", format(Sys.time()), "\n")
            
            invisible(new("dkiIndices",
                          call = args,
                          fa = array(z$fa, ddim),
                          ga = array(z$ga, ddim),
                          md = array(z$md, ddim),
                          andir = array(z$andir, c(3, ddim)),
                          bary = array(z$bary, c(3, ddim)),
                          k1 = k1,
                          k2 = k2,
                          k3 = k3,
                          mk = mk,
                          mk2 = mk2,
                          kaxial = kaxial,
                          kradial = kradial,
                          fak = fak,
                          gradient = object@gradient,
                          bvalue = object@bvalue,
                          btb   = object@btb,
                          ngrad = object@ngrad,
                          s0ind = object@s0ind,
                          ddim  = ddim,
                          ddim0 = object@ddim0,
                          voxelext = object@voxelext,
                          orientation = object@orientation,
                          rotation = object@rotation,
                          xind  = object@xind,
                          yind  = object@yind,
                          zind  = object@zind,
                          method = object@method,
                          level = object@level,
                          source = object@source)
            )
          })

dkiDesign <- function(gradients) {
  
  ## start with some checks
  dgrad <- dim(gradients)
  
  ## do we have a matrix (check for class is not correct here, as array or dataframe would also do)
  if (length(dgrad) != 2) stop("dkiDesign: dimensionality of gradients is not 2!")
  
  ## we expect the columns of gradients to be the gradient vectors
  if (dgrad[1] != 3) stop("dkiDesign: not a valid gradient matrix, expect gradient vectors as columns!")
  
  ## now we have: - dgrad[2] is the number of gradients
  
  X <- matrix(0, dgrad[2], 21)
  
  X[, 1] <- gradients[1, ]^2 # D_11
  X[, 2] <- gradients[2, ]^2 # D_22
  X[, 3] <- gradients[3, ]^2 # D_33
  
  X[, 4] <- 2 * gradients[1, ] * gradients[2, ] # D_12
  X[, 5] <- 2 * gradients[1, ] * gradients[3, ] # D_13
  X[, 6] <- 2 * gradients[2, ] * gradients[3, ] # D_23
  
  X[, 7] <- gradients[1, ]^4 # V_1111
  X[, 8] <- gradients[2, ]^4 # V_2222
  X[, 9] <- gradients[3, ]^4 # V_3333
  
  X[, 10] <- 4 * gradients[1, ]^3 * gradients[2, ] # V_1112
  X[, 11] <- 4 * gradients[1, ]^3 * gradients[3, ] # V_1113
  X[, 12] <- 4 * gradients[2, ]^3 * gradients[1, ] # V_2221
  X[, 13] <- 4 * gradients[2, ]^3 * gradients[3, ] # V_2223
  X[, 14] <- 4 * gradients[3, ]^3 * gradients[1, ] # V_3331
  X[, 15] <- 4 * gradients[3, ]^3 * gradients[2, ] # V_3332
  
  X[, 16] <- 6 * gradients[1, ]^2 * gradients[2, ]^2 # V_1122
  X[, 17] <- 6 * gradients[1, ]^2 * gradients[3, ]^2 # V_1133
  X[, 18] <- 6 * gradients[2, ]^2 * gradients[3, ]^2 # V_2233
  
  X[, 19] <- 12 * gradients[1, ]^2 * gradients[2, ] * gradients[3, ] # V_1123
  X[, 20] <- 12 * gradients[2, ]^2 * gradients[1, ] * gradients[3, ] # V_2213
  X[, 21] <- 12 * gradients[3, ]^2 * gradients[1, ] * gradients[2, ] # V_3312
  
  X
}











