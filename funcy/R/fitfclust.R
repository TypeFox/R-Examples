#  Original  Copyright (C) 2003  Gareth M. James and Catherine A. Sugar
#  Modified Code: Copyright (C) 2011-2015 Christina Yassouridis
#  all functions for regular data were implemented by Yassouridis
#


"fitfclust" <-
    function(data, k=2, reg=reg, regTime=NULL, dimBase=dimBase, h=1, p=5, epsilon=0.01,
             maxiter=20, pert =0.01, hard=hard,
             seed=seed, init=init, nrep=nrep, fpcCtrl, baseType="splines"){
        if(baseType=="fourier")
            dimBase <- fourierWarn(dimBase)
        ##Check Input
        if(h>k)
            stop("Reduced dimension in mbcCtrl cannot be greater than k!")
        if(reg==1){
            if(is.null(regTime))
                grid <- 1:dim(data)[1]
            else
                grid <- regTime
            initFct <- initStepReg
            MstepFct <- fitfclustMstep
            EstepFct <- fitfclustEstep
        }else if(reg==0){
            grid <- unique(sort(data[,3]))
            initFct <- initStepIrreg
            MstepFct <- fitfclustMstepIrreg
            EstepFct <- fitfclustEstepIrreg
        }
        
        ##INITIAL CLUSTERING --------------------------------------------
        initFit <- initFct(data=data, pert=pert, grid=grid, h=h,
                           p=p, dimBase=dimBase, k=k, baseType,
                           seed=seed, init=init, nrep=nrep, fpcCtrl=fpcCtrl)

        ## Initialize the parameters ------------------------------------
        plotParams <- initFit$plotParams
        parameters <- initFit$parameters
        vars <- initFit$vars
        fullBase <- initFit$fullBase
        base <- initFit$base
        
        sigma.old <- 0
        sigma.new <- 1
        kk <- 1
        ## Main loop. Iterates between M and E steps and stops when
        ## sigma  has converged.
        while(abs(sigma.old - sigma.new)/sigma.new > epsilon & (kk <=
                                                                    maxiter)) {
            parameters <- MstepFct(parameters=parameters, data=data,
                                   vars=vars, base=fullBase,
                                   epsilon=epsilon, p=p, hard=hard)
            
            ##returns lambda0, Lambda, alpha, pi, Gamma, sigma
            vars <- EstepFct(parameters=parameters, data=data, vars=vars,
                             base=fullBase, hard=hard)
            ##returns gamma, piigivej, gprod, gcov 
            sigma.old <- sigma.new
            sigma.new <- parameters$sigma[1]
            kk <- kk + 1
        }

        ##added by Christina
        nc <- dim(vars$gamma)[1]
        loglik <- vars$loglik
        nrParams <- 1+dimBase^2+(k-1)+k*h+dimBase*h+k
        aic <- 2*nrParams-2*loglik
        bic <- -2*loglik+nrParams*log(nc)
        
        ## Enforce parameter constraint (7)
        cfit <- fitfclustconst(parameters, vars, base)
        list(data=data, plotParams=plotParams, parameters=cfit$parameters, vars=cfit$vars,
             fullBase=fullBase, base=base, grid=grid, nrIter=kk,
             aic=aic, bic=bic, loglik=loglik)
    }


##IRREGULAR DATA FUNCTIONS**************************************************************
"initStepIrreg" <- function(data, pert, grid, h, p,
                            dimBase, k, baseType, seed, init, nrep,
                            fpcCtrl){
    base <- fullBase <- NULL

    if(baseType=="eigenbasis"){
        mkdat <- formatFuncy(data=data, format="Format3")
        res <- with(mkdat,fpc(Yin=Yin, Tin=Tin, isobs=isobs,
                              dimBase=dimBase,fpcCtrl=fpcCtrl,
                              Tout=grid, reg=FALSE))
        plotParams <- res$plotParams
    }else{
        res <- makeCoeffs(data=data, reg=FALSE, dimBase=dimBase,
                          grid=grid, pert=pert,
                          baseType=baseType)
        plotParams <- NULL
    }
    coeffs <- res$coeffs
    base <- res$base
    fullBase <- res$fullBase

    ## Use k-means to get initial cluster memberships from coeffs.
    nc <- max(data[,1])
    if(k > 1) 
        class <- initClust(data=coeffs, k=k, init=init, seed=seed, nrep=nrep)$clusters
    else
        class <- rep(1, nc)

    ## Initial estimates for the posterior probs of cluster membership.
    piigivej <- matrix(0, nc, k)
    piigivej[col(piigivej) == class] <- 1 
    ## Calculate basis coefficients for cluster means.
    groupCoeffs <- matrix(0,k,dimBase)
    
    for (l in 1:k)
        {
            if(sum(class==l)>1)
                groupCoeffs[l,] <- apply(coeffs[class==l,], 2, mean)
            else
                groupCoeffs[l,] <- coeffs[class==l,] 
        }
    
    ## Calculate overall mean lambda0 
    lambda.zero <- apply(groupCoeffs, 2, mean)
    Lambda <- as.matrix(svd(scale(groupCoeffs, scale=F))$v[, 1:h])
    alpha <- scale(groupCoeffs, scale=F)%*%Lambda 
    gamma <- t(t(coeffs) - lambda.zero - (Lambda %*% t(alpha[class,  ])))
    gprod <- NULL
    for(i in 1:nc)
        gprod <- cbind(gprod, (gamma[i,  ]) %*% t(gamma[i,  ]))
    
    ##covariance matrices for the errors of the coefficients.
    gamma <- array(gamma[, rep(1:sum(dimBase), rep(k, sum(dimBase)))],
                   c(nc, k, sum(dimBase)))
    
    gcov <- matrix(0, sum(dimBase), nc * sum(dimBase))
    list(base=base, fullBase=fullBase, plotParams=plotParams, parameters =
             list(lambda.zero=lambda.zero, Lambda=Lambda, alpha =
                      alpha), vars=list(gamma=gamma, gprod=gprod, gcov=gcov, piigivej=piigivej))
}


"fitfclustMstepIrreg" <-
    function(parameters, data, vars, base, epsilon, p, hard){
        ## This function implements the M step of the EM algorithm. It
        ## computes the parameters that maximize the function.
        
        k <- dim(parameters$alpha)[1]
        alpha <- parameters$alpha
        Lambda <- parameters$Lambda
        gamma <- vars$gamma
        gcov <- vars$gcov
        curve <- data[,1]
        piigivej <- vars$piigivej
        nc <- dim(gamma)[1] 
        k<- dim(alpha)[1] 
        h <- dim(alpha)[2] 
        nTot <- length(curve) 
        dimBase <- dim(base)[2] 
        sigma.old <- 2
        sigma <- 1
        
        ##PI
        if(hard)
            parameters$pi <- rep(1/k, k) 
        else
            parameters$pi <- apply(vars$piigivej, 2, mean)
        
        ## Compute rank p estimate of Gamma
        ind <- matrix(rep(c(rep(c(1, rep(0, dimBase - 1)), nc), 0), dimBase)[1:(nc*dimBase^2)], nc * dimBase, dimBase)
        gsvd <- svd(vars$gprod %*% ind/nc)
        gsvd$d[ - (1:p)] <- 0 
        parameters$Gamma <- gsvd$u %*% diag(gsvd$d) %*% t(gsvd$u)
        
        ## This loop iteratively calculates alpha and then Lambda and stops
        ## when they have converged.
        loopnumber <- 1
        while((abs(sigma.old[1] - sigma[1])/sigma[1] > epsilon) & (loopnumber <10)){
            sigma.old <- sigma
            gamma.pi <- diag(base %*% t(apply(gamma * as.vector(piigivej), c(1, 3), sum)[curve,  ]))
            alpha.pi <- t(matrix(apply(alpha, 2, function(x, piigivej,k)
                {rep(1, k) %*% (piigivej * x)}
                                     , t(piigivej), k), nc, h)[curve,  ])

            lambda.zero <- solve(t(base) %*% base) %*% t(base) %*% (data[,2] - diag(base %*% Lambda %*% alpha.pi) - gamma.pi)
            
            x <- t(data[,2] - t(base %*% lambda.zero))
            for(i in 1.:k) {
                base.Lam <- base %*% Lambda
                base.Lam.pi <- base.Lam * piigivej[curve, i]
                if(sum(piigivej[, i]) > 10^(-4)){
                    alpha[i,  ] <- solve(t(base.Lam.pi) %*% base.Lam) %*% t(base.Lam.pi) %*% (x - diag(base %*% t(gamma[curve, i,  ])))
                }else print("Warning: empty cluster")
            }

            ## through each column of Lambda holding the other columns fixed.
            for(m in 1:h) {
                pi.alphasq <- apply(t(piigivej) * (alpha^2)[, m], 2,sum)[curve]
                pi.alpha <- apply(t(piigivej) * alpha[, m], 2, sum)[curve]
                base.Lambda <- t(base %*% Lambda)
                if(h != 1) {
                    temp <- NULL
                    for(i in 1:k) {
                        temp <- cbind(temp, as.vector(rep(1, h - 1) %*% matrix((base.Lambda * alpha[i,  ])[ - m,  ], h - 1,dim(base)[1])) * alpha[i, m])
                    }
                    otherLambda <- (as.vector(temp) * piigivej[curve,  ])%*%rep(1, k)
                }
                else otherLambda <- 0
                gamma.pi.alpha <- apply(gamma * as.vector(piigivej) *
                                            rep(alpha[, m], rep(nc, k)), c(1, 3), sum)[curve,  ]
                Lambda[, m] <- solve(t(base * pi.alphasq) %*% base) %*% t(base) %*%
                    (x * pi.alpha - otherLambda - (base *gamma.pi.alpha) %*% rep(1, sum(dimBase)))
            }
            
            ind <- matrix(rep(c(rep(F, dimBase), rep(T, nc * dimBase)),nc)
                          [1:(nc * nc * dimBase)], nc, nc * dimBase, byrow=T)[rep(1:nc, table(curve)),]
            mat1 <- matrix(rep(base, nc), nTot, nc * dimBase)
            mat3 <- t(mat1)
            mat3[t(ind)] <- 0
            ind2 <- matrix(rep(c(rep(F, dimBase), rep(T, nc * dimBase)),
                               nc)[1:(nc * nc * dimBase)], nc, nc * dimBase, byrow=T)[rep(1:nc, rep(dimBase, nc)),]
            mat2 <- matrix(rep(t(gcov), nc), nc * dimBase, nc * dimBase,byrow=T)
            mat2[ind2] <- 0
            econtrib <- 0
            for(i2 in 1:k) {
                vect1 <- x - base %*% Lambda %*% alpha[i2,  ] - (base *
                                                                     gamma[curve, i2,  ]) %*% rep(1, dimBase)
                econtrib <- econtrib + t(piigivej[curve,i2] * vect1) %*% vect1
            }
            sigma <- as.vector((econtrib + sum(diag(mat1 %*% mat2 %*% mat3)))/nTot)
            loopnumber <- loopnumber + 1
        }
        parameters$lambda.zero <- as.vector(lambda.zero)
        parameters$alpha <- alpha
        parameters$Lambda <- Lambda
        parameters$sigma <- sigma
        parameters
    }

"fitfclustEstepIrreg" <-
    function(parameters, data, vars, base, hard){
        ## This function performs the E step of the EM algorithm by
        ## calculating the expected values of gamma and gamma %*% t(gamma)
        ## given The currenT parameter estimates.
        par <- parameters
        nc <- dim(vars$gamma)[1]
        k <- dim(vars$gamma)[2]
        dimBase <- dim(vars$gamma)[3]
        Gamma <- par$Gamma
        Lambda.alpha <- par$lambda.zero + par$Lambda %*% t(par$alpha)
        vars$gprod <- vars$gcov <- NULL
        curveIndx <- data[,1]
        cost <- matrix(0, nrow=nc, ncol=k)
        
        for(j in 1:nc) {
            ## Calculate expected value of gamma.
            basej <- base[curveIndx==j,]
            if(is.null(dim(basej)))
                basej <- t(basej)
            nj <- sum(curveIndx == j)
            if(!nj==1)
                invvar <- diag(1/rep(par$sigma, nj))
            else
                invvar <- 1/par$sigma

            ## Variance of coefficients
            Cgamma <- Gamma - Gamma %*% t(basej) %*% solve(diag(nj) + basej %*%
                                                               Gamma %*% t(basej) %*% invvar) %*% invvar %*% basej %*% Gamma
            centx <- data[curveIndx==j,2]  - basej %*% Lambda.alpha
            vars$gamma[j,  ,  ] <- t(Cgamma %*% t(basej) %*% invvar %*% centx)
            
            ## Posterior probabilities 
            covx <- basej %*% par$Gamma %*% t(basej) + solve(invvar)
            d <- exp( - diag(t(centx) %*% solve(covx) %*% centx)/2) * par$pi
            cost[j,] <-exp( - diag(t(centx) %*% solve(covx) %*% centx)/2)
            
            vars$piigivej[j,  ] <- d/sum(d)
            if(hard) {
                m <- order( - d)[1]
                vars$piigivej[j,  ] <- 0
                vars$piigivej[j, m] <- 1
                
            }
            ## Calculate expected value of gamma and gamma %*% t(gamma)
            vars$gprod <- cbind(vars$gprod, t(matrix(vars$gamma[j,  ,  ],
                                                     k, dimBase)) %*% (matrix(vars$gamma[j,  ,  ], k, dimBase) * vars$
                                                                       piigivej[j,  ]) + Cgamma)
            vars$gcov <- cbind(vars$gcov, Cgamma)
        }
        vars$loglik <- sum(log(cost))
        vars
    }

##Additional Functions------------------------------------------------------
"fitfclust.predIrreg" <-
    function(fit, data=NULL, reweight=F){
        ## This function produces the alpha hats used to provide a low
        ## dimensional pictorial respresentation of each curve. It also
        ## produces a class prediction for each curve. It takes as
        ## input the fit from fldafit (for predictions on the original data)
        ## or the fit and a new data set (for predictions on new data).
        if (is.null(data))
            data <- fit$data
        fullBase <- fit$fullBase
        base <- fit$base
        par <- fit$parameters
        curve <- data[,1]#data$curveIndx
        grid <- sort(unique(data[,3]))
        time <- match(data[,3], grid)#data$timeIndx
        nc <- length(table(curve)) #number of curves
        h <- dim(par$alpha)[2]
        alpha.hat <- matrix(0, nc, h)
        k <- dim(fit$par$alpha)[1]
        distance <- matrix(0, nc, k)
        Calpha <- array(0, c(nc, h, h))
        
        for(i in 1:nc) {
            baseij <- base[time[curve == i],  ]
            if(is.null(dim(baseij)))
                baseij <- t(baseij)
            Y_ij <-data[curve==i,2] #data$yin[curve == i]
            n_ij<- length(Y_ij)
            
            Sigma <- par$sigma * diag(n_ij) + baseij %*% par$Gamma %*% t(baseij)
            ## Calculate covariance for each alpha hat.(11)
            InvCalpha <- t(par$Lambda) %*% t(baseij) %*% solve(Sigma) %*%
                baseij %*% par$Lambda
            Calpha[i,  ,  ] <- ginv(InvCalpha) ##Christina changed solve to ginv 
            ## Calculate each alpha hat(11).
            alpha.hat[i,  ] <- Calpha[i,  ,  ] %*% t(par$Lambda) %*% t(
                baseij) %*% solve(Sigma) %*% (Y_ij - baseij %*% par$lambda.zero)
            ## Calculate the matrix of distances, relative to the
            ## appropriate metric of each curve from each class centroid. 
            for (k in 1:k){
                y <- as.vector(alpha.hat[i,])-fit$par$alpha[k,] #fit$par$alpha=centroid
                distance[i,k] <- t(y)%*%InvCalpha %*%y
            }
        }
        
        ## Calculate final class predictions for each curve.
        class.pred <- rep(1, nc)
        log.pi <- log(fit$par$pi)
        if (!reweight)
            log.pi <- rep(0,k)
        probs <- t(exp(log.pi-t(distance)/2)) #(13)
        probs <- probs/apply(probs,1,sum) #probability that curve i belongs to class j
        ##class prediction for the curves (highest prob from above)
        m <- probs[,1] #prob that curve i belongs to clas 1
        if(k!= 1)
            for(l in 2:k) {
                test <- (probs[, l] > m)
                class.pred[test] <- l
                m[test] <- probs[test, l]
            }
        list(Calpha=Calpha, alpha.hat=alpha.hat, class.pred=class.pred,
             distance=distance, m=m, probs=probs)
    }

"fitfclust.curvepredIrreg" <-
    function(fit, data=NULL, index=NULL, tau=0.95, tau1=0.975){
        if (is.null(data))
            data <- fit$data
        if (is.null(index))
            index <- 1:max(data[,1])
        curveIndx <- data[,1]
        grid <- unique(sort(data[,3]))
        timeIndx <- match(data[,3],grid)
        tau2 <- tau/tau1
        sigma <- fit$par$sigma
        Gamma <- fit$par$Gamma
        Lambda <- fit$par$Lambda
        alpha <- fit$par$alpha
        lambda.zero <- as.vector(fit$par$lambda.zero)
        fullBase <- fit$fullBase
        base <- fit$base
        nc <- length(index)
        upci <-lowci <- uppi <- lowpi <- gpred <- matrix(0,nc,nrow(base))
        etapred <- matrix(0,nc,ncol(base))
        ind <- 1
        Lambda.alpha <- lambda.zero + Lambda %*% t(alpha)
        for (i in index){
            y <-data[curveIndx==i,2] 
            basei <-  fullBase[timeIndx[curveIndx == i],  ]
            ni <- dim(basei)[1]
            
            if(is.null(dim(basei))){
                ni <- 1
                basei <- t(basei)
                invvar <- 1/sigma
            } else
                invvar <- diag(1/rep(sigma, ni))
            covx <- basei %*% Gamma %*% t(basei) + solve(invvar)
            centx <- y- basei %*% Lambda.alpha
            d <- exp( - diag(t(centx) %*% solve(covx) %*% centx)/2) * fit$par$pi
            pi <- d/sum(d)
            k<- length(pi)
            mu <- lambda.zero + Lambda %*% t(alpha * pi) %*% rep(1, k)
            cov <- (Gamma - Gamma %*% t(basei) %*% solve(diag(ni) + basei %*% Gamma %*%
                                                             t(basei)/sigma) %*%
                                                                 basei %*% Gamma/sigma)/sigma
            
            etapred[ind,] <- mu + cov %*% t(basei) %*% (y - basei %*% mu)
            ord <- order( - pi) #sort from biggest to smallest
            numb <- sum(cumsum(pi[ord]) <= tau1) + 1
            v <- diag(base %*% (cov * sigma) %*% t(base))
            pse <- sqrt(v + sigma)
            se <- sqrt(v)
            lower.p <- upper.p <- lower.c <- upper.c <- matrix(0, nrow(base), numb)
            for(j in 1:numb) {
                mu <- lambda.zero + Lambda %*% alpha[ord[j],  ]
                mean <- base %*% (mu + cov %*% t(basei) %*% (y - basei %*% mu))
                upper.p[, j] <- mean + qnorm(tau2) * pse #mixture likelihood
                lower.p[, j] <- mean - qnorm(tau2) * pse
                upper.c[, j] <- mean + qnorm(tau2) * se #classification likelihood
                lower.c[, j] <- mean - qnorm(tau2) * se
            }
            upci[ind,] <- nummax(upper.c)$max
            lowci[ind,] <-  - nummax( - lower.c)$max
            uppi[ind,] <- nummax(upper.p)$max
            lowpi[ind,] <-  - nummax( - lower.p)$max
            gpred[ind,] <- as.vector(base %*%etapred[ind,])
            ind <- ind+1
        }
        
        meancurves <- base%*%Lambda.alpha
        list(etapred=etapred, gpred=gpred,  upci=upci,lowci=lowci,  uppi=uppi, lowpi=lowpi,index=index,grid=fit$grid,data=data,meancurves=meancurves)
    }


"fitfclust.plotcurvesIrreg" <-function(object=NULL, fit=NULL, index=NULL,
                                       ci=T, pi=T, clustermean=F){
    if (is.null(object))
        object <- fitfclust.curvepredIrreg(fit=fit)
    if (is.null(index))
        index <- 1:length(table(object$data[,1]))
    r <- min(ceiling(sqrt(length(index))),5)
    timeIndx <- match(object$data[,3],object$grid)
    curveIndx <- object$data[,1]
    par(mfrow=c(r,r))
    for (i in index){
        grid <- object$grid
        upci <- object$upci[i,]
        uppi <- object$uppi[i,]
        lowci <- object$lowci[i,]
        lowpi <- object$lowpi[i,]
        gpred <- object$gpred[i,]
        meancurves <- (object$mean)
        if (clustermean){
            yrange <- c(min(c(lowpi,meancurves)),max(c(uppi,meancurves)))
        }  else  yrange <- c(min(lowpi),max(uppi))
        plot(grid,grid,ylim=yrange,ylab="Predictions",xlab="Time",type='n',
             main=paste("Curve ",i))
        if (clustermean)
            for (k  in 1:ncol(meancurves))
                lines(grid,meancurves[,k],col=6,lty=2,lwd=2)
        if (ci){
            lines(grid,upci,col=3)
            lines(grid,lowci,col=3)}
        if (pi){
            lines(grid,uppi,col=4)
            lines(grid,lowpi,col=4)}
        lines(grid,gpred,col=2,lwd=2)
        lines(grid[timeIndx[curveIndx==i]],object$data[curveIndx==i,2],lwd=2)
        points(grid[timeIndx[curveIndx==i]],object$data[curveIndx==i,2],pch=19,cex=1.5)
    }
    par(mfrow=c(1,1))
}

############################################################
##Functions for regular data *******************************
"initStepReg" <- function(data, pert, grid, h, p,
                          dimBase, k, baseType, seed, init, nrep,
                          fpcCtrl){

    if(baseType=="eigenbasis"){
        mkdat <- formatFuncy(data, regTime=grid, format="Format3")
        res <- with(mkdat,fpc(Yin=Yin, Tin=Tin, isobs=isobs, dimBase=dimBase,
                              fpcCtrl=fpcCtrl, reg=TRUE))
        plotParams <- res$plotParams
    }else{
        res <- makeCoeffs(data=t(data), reg=TRUE, dimBase=dimBase,
                          grid=grid, pert=pert, baseType=baseType)
        plotParams <- NULL
    }

    coeffs <- res$coeffs
    base <- res$base
    fullBase <- res$fullBase
    nc <- dim(data)[2] 
    
    ## Use k-means to get initial cluster memberships from coeffs.
    if(k > 1) {
        class <- initClust(data=coeffs, k=k, init=init, seed=seed, nrep=nrep)$clusters 
    }  else
        class <- rep(1, nc)

    ## Initial estimates for the posterior probs of cluster membership.
    piigivej <- matrix(0, nc, k)
    piigivej[col(piigivej) == class] <- 1
    classmean <- matrix(0,k,dimBase)
    
    for (l in 1:k)
        {
            if(sum(class==l)>1)
                classmean[l,] <- colMeans(coeffs[class==l,])
            else
                classmean[l,] <- coeffs[class==l,] 
        }
    
    ## Calculate overall mean lambda0 
    lambda.zero <- colMeans(classmean) 
    Lambda <- as.matrix(svd(scale(classmean, scale=F))$v[, 1:h])
    alpha <- scale(classmean, scale=F)%*%Lambda
    gamma <- t(t(coeffs) - lambda.zero - (Lambda %*%
                                              t(alpha[class,  ])))
    
    gprod <- do.call(cbind,lapply(1:nc, function(x) gamma[x,]%*%t(gamma[x,])))
    gamma <- array(gamma[, rep(1:sum(dimBase), rep(k, sum(dimBase)))],
                   c(nc, k, sum(dimBase)))
    
    gcov <- matrix(0, sum(dimBase), nc * sum(dimBase))
    list(base=base, fullBase=fullBase, plotParams=plotParams, parameters=list(lambda.zero=lambda.zero,Lambda=Lambda, alpha=alpha), vars=list(gamma=gamma,gprod=gprod, gcov=gcov, piigivej=piigivej))
}


"fitfclustMstep" <-
    function(parameters, data, vars, base, epsilon, p, hard){
        ## This function implements the M step of the EM algorithm. It computes the parameters that maximize the function.
        k <- dim(parameters$alpha)[1]
        alpha <- parameters$alpha
        Lambda <- parameters$Lambda
        gamma <- vars$gamma
        gcov <- vars$gcov
        piigivej <- vars$piigivej
        nc <- dim(gamma)[1] #number of curves
        k <- dim(alpha)[1] #number of classes
        h <- dim(alpha)[2] #dim of projection
        nt <- dim(data)[1]#total number time coeffs of all curves together
        dimBase <- dim(base)[2] #number of basis functions for spline basis

        sigma.old <- 2
        sigma <- 1
        
        if(hard)
            parameters$pi <- rep(1/k, k) #each class gehts the same prob
        else
            parameters$pi <- apply(vars$piigivej, 2, mean)#percentage of curves in each class
        ## Compute rank p estimate of Gamma
        ind <- matrix(rep(c(rep(c(1, rep(0, dimBase - 1)), nc), 0), dimBase)[1:(nc*dimBase^2)], nc * dimBase, dimBase)
        
        ##Covariance is reconstructed after projection into p dimensions
        gsvd <- svd(vars$gprod %*% ind/nc)# takes the mean of the covariance
                                        # matrix entries in each dimension
                                        # and makes singular value
                                        # decomposition of the resulting
                                        # total covariance matrix
        gsvd$d[ - (1:p)] <- 0 #sets all exept first p eigenvalues to zero
        parameters$Gamma <- gsvd$u %*% diag(gsvd$d) %*% t(gsvd$u)
        ## This loop iteratively calculates alpha and then Lambda and stops
        ## when they have converged.
        loopnumber <- 1
        while((abs(sigma.old[1] - sigma[1])/sigma[1] > epsilon) & (loopnumber <10)){
            sigma.old <- sigma
            
            ## Calculate lambda.zero.
            gamma.pi <- (base %*% t(apply(gamma * as.vector(piigivej), c(1, 3), sum)))
            alpha.pi <- t(matrix(apply(alpha, 2, function(x, piigivej,k)
                {rep(1, k) %*% (piigivej * x)}
                                     , t(piigivej), k), nc, h))
            lambda.zero <- solve(t(base) %*% base) %*% t(base) %*% ((data) - (base %*% Lambda %*% alpha.pi) - gamma.pi)
            lambda.zero <- 1/nc*rowSums(lambda.zero)
            x <- (data) -as.vector(base %*% lambda.zero)
            
            for(i in 1.:k) {
                S.Lam <- base %*% Lambda
                S.Lam.pi <- S.Lam
                nr <- sum(piigivej[,i])
                
                if(sum(piigivej[, i]) > 10^(-4)){
                    alpha[i,  ] <- (solve((apply(t(S.Lam.pi)%o%(piigivej[,i]),1:2,sum))%*%S.Lam) %*% t(S.Lam.pi) %*% (x - (base %*% t(gamma[, i,  ]))))%*%(piigivej[,i])
                    
                }else print("Warning: empty cluster")
            }

            ## Calculate Lambda given alpha. This is done by iterating
            ## through each column of Lambda holding the other columns fixed. (31)
            for(m in 1:h) {
                pi.alphasq <- apply(t(piigivej) * (alpha^2)[, m], 2,sum)
                pi.alpha <- apply(t(piigivej) * alpha[, m], 2, sum)
                S.Lambda <- t(base %*% Lambda)
                if(h != 1) {
                    temp <- NULL
                    for(i in 1:k) {
                        temp <- cbind(temp, as.vector(rep(1, h - 1) %*% matrix((S.Lambda * alpha[i,  ])[ - m,  ], h - 1,dim(base)[1])) * alpha[i, m])
                    }
                    otherLambda <- apply(piigivej,1, function(x) colSums((x%o%temp)[,,]))
                }
                else otherLambda <- 0
                
            }
            ## Calculate sigma (32)
            ind <- matrix(rep(c(rep(F, dimBase), rep(T, nc * dimBase)),nc)
                          [1:(nc * nc * dimBase)], nc, nc * dimBase, byrow=T)[rep(1, nt),]
            mat1 <- matrix(rep(base, nc), nt, nc * dimBase)
            mat3 <- t(mat1)
            mat3[t(ind)] <- 0
            ind2 <- matrix(rep(c(rep(F, dimBase), rep(T, nc * dimBase)),
                               nc)[1:(nc * nc * dimBase)], nc, nc * dimBase, byrow=T)[rep(1:nc, rep(dimBase, nc)),]
            mat2 <- matrix(rep(t(gcov), nc), nc * dimBase, nc * dimBase,byrow=T)
            mat2[ind2] <- 0
            econtrib <- 0
            for(i2 in 1:k) {
                vect1 <- x - as.vector(base %*% Lambda %*% alpha[i2,  ]) - (base %*% t(gamma[, i2,  ])) #%*% rep(1, dimBase)
                econtrib <- econtrib + piigivej[,i2]%*%diag(t(vect1) %*% vect1)
            }
            sigma <- as.vector((econtrib + sum(diag(mat1 %*% mat2 %*% mat3))*nc)/(nt*nc))
            loopnumber <- loopnumber + 1
        }
        parameters$lambda.zero <- as.vector(lambda.zero)
        parameters$alpha <- alpha
        parameters$Lambda <- Lambda
        parameters$sigma <- sigma
        parameters
    }

"fitfclustEstep" <-
    function(parameters, data, vars, base, hard){
        ## This function performs the E step of the EM algorithm by
        ## calculating the expected values of gamma and gamma %*% t(gamma)
        ## given the current parameter estimates.
        par <- parameters
        nc <- dim(vars$gamma)[1]
        k <- dim(vars$gamma)[2]
        dimBase <- dim(vars$gamma)[3]
        Gamma <- par$Gamma
        Lambda.alpha <- par$lambda.zero + par$Lambda %*% t(par$alpha)
        vars$gprod <- vars$gcov <- NULL
        nt <- dim(data)[1]
        invvar <- diag(1/rep(par$sigma, nt))
        Cgamma <- Gamma - Gamma %*% t(base) %*% solve(diag(nt) + base %*%
                                                          Gamma %*% t(base) %*% invvar) %*% invvar %*% base %*% Gamma
        covx <- base %*% par$Gamma %*% t(base) + solve(invvar)
        cost <- matrix(0, nrow=nc, ncol=k)
        
        for(j in 1:nc) {
            ## Calculate expected value of gamma.
            centx <- data[,j] - base %*% Lambda.alpha
            vars$gamma[j,  ,  ] <- t(Cgamma %*% t(base) %*% invvar %*% centx)
            ## Calculate pi i given j.
            d <- exp( - diag(t(centx) %*% solve(covx) %*% centx)/2) *
                par$pi
            cost[j,] <-  exp( - diag(t(centx) %*% solve(covx) %*% centx)/2)
            if(sum(d)!=0)
                vars$piigivej[j,  ] <- d/sum(d)
            else
                vars$piigivej[j,  ] <- 0
            
            if(hard) {
                m <- order( - d)[1]
                vars$piigivej[j,  ] <- 0
                vars$piigivej[j, m] <- 1
            }
            
            ## Calculate expected value of gamma %*% t(gamma).
            vars$gprod <- cbind(vars$gprod, t(matrix(vars$gamma[j,  ,
                                                                ],k,
                                                     dimBase)) %*%
                                                         (matrix(vars$gamma[j,  ,  ], k, dimBase) * vars$piigivej[j,  ]) + Cgamma)
            vars$gcov <- cbind(vars$gcov, Cgamma)
        }
        vars$loglik <- sum(log(cost))
        vars
    }


"fitfclust.pred" <-
    function(fit, data=NULL, reweight=F){
                                        # This function produces the alpha hats used to provide a low
                                        # dimensional pictorial respresentation of each curve. It also
                                        # produces a class prediction for each curve. It takes as
                                        # input the fit from fldafit (for predictions on the original data)
                                        # or the fit and a new data set (for predictions on new data).
        if (is.null(data))
            data <- t(fit$data)
        fullBase <- fit$base
        par <- fit$parameters
                                        #curve <- data$curve
                                        #time <- data$timeindex
        nc <- dim(data)[1] #number of curves
        h <- dim(par$alpha)[2]
        alpha.hat <- matrix(0, nc, h)
        k <- dim(fit$par$alpha)[1]
        distance <- matrix(0, nc, k)
        Calpha <- array(0, c(nc, h, h))
        for(i in 1:nc) {
            baseij <- fullBase
            xij <- data[i,]
            nij <- length(xij)
            Sigma <- par$sigma * diag(nij) + baseij %*% par$Gamma %*% t(baseij)
                                        # Calculate covariance for each alpha hat.(11)
            InvCalpha <- t(par$Lambda) %*% t(baseij) %*% solve(Sigma) %*% baseij %*%
                par$Lambda
            Calpha[i,  ,  ] <- solve(InvCalpha)
                                        # Calculate each alpha hat(11).
            alpha.hat[i,  ] <- Calpha[i,  ,  ] %*% t(par$Lambda) %*% t(
                baseij) %*% solve(Sigma) %*% (xij - baseij %*% par$lambda.zero)
                                        # Calculate the matrix of distances, relative to the
                                        # appropriate metric of each curve from each class centroid. 
            for (l in 1:k){
                y <- as.vector(alpha.hat[i,])-fit$par$alpha[l,] #fit$par$alpha=centroid
                distance[i,l] <- t(y)%*%InvCalpha %*%y
            }
        }
                                        # Calculate final class predictions for each curve.
        class.pred <- rep(1, nc)
        log.pi <- log(fit$par$pi)
        if (!reweight)
            log.pi <- rep(0,k)
        probs <- t(exp(log.pi-t(distance)/2)) #(13)
        probs <- probs/apply(probs,1,sum) #probability that curve i belongs to class j
        ##class prediction for the curves (highest prob from above)
        m <- probs[,1] #prob that curve i belongs to clas 1
        if(k != 1)
            for(l in 2:k) {
                test <- (probs[, l] > m)
                class.pred[test] <- l
                m[test] <- probs[test, l]
            }
        list(Calpha=Calpha, alpha.hat=alpha.hat, class.pred=class.pred,
             distance=distance, m=m,probs=probs)
    }


"fitfclust.curvepred" <-
    function(fit, data=NULL, index=NULL, tau=0.95, tau1=0.975){
        if (is.null(data))
            data <- fit$data
        if (is.null(index))
            index <- 1:dim(data)[2]
        tau2 <- tau/tau1
        sigma <- fit$par$sigma
        Gamma <- fit$par$Gamma
        Lambda <- fit$par$Lambda
        alpha <- fit$par$alpha
        lambda.zero <- as.vector(fit$par$lambda.zero)
        base <- fit$base
        nc <- length(index)
        upci <-lowci <- uppi <- lowpi <- gpred <- matrix(0,nc,nrow(base))
        etapred <- matrix(0,nc,ncol(base))
        ind <- 1
        Lambda.alpha <- lambda.zero + Lambda %*% t(alpha)
        for (i in index){
            y <- data[,i]
            basei <- base
            ni <- dim(basei)[1]
            invvar <- diag(1/rep(sigma, ni))
            covx <- basei %*% Gamma %*% t(basei) + solve(invvar)
            
            centx <- data[,i] - basei %*% Lambda.alpha
            d <- exp( - diag(t(centx) %*% solve(covx) %*% centx)/2) *
                fit$par$pi

            if(all(d==0))
                pi <- rep(1,length(fit$par$pi))
            else
                pi <- d/sum(d)
            
            k <- length(pi)
            mu <- lambda.zero + Lambda %*% t(alpha * pi) %*% rep(1, k)
            cov <- (Gamma - Gamma %*% t(basei) %*% solve(diag(ni) + basei %*% Gamma %*%
                                                             t(basei)/sigma) %*% basei %*% Gamma/sigma)/sigma
            etapred[ind,] <- mu + cov %*% t(basei) %*% (y - basei %*% mu)
            ord <- order( - pi) #sort from biggest to smallest
            numb <- sum(cumsum(pi[ord]) <= tau1) + 1
            v <- diag(base %*% (cov * sigma) %*% t(base))
            pse <- sqrt(v + sigma)
            se <- sqrt(v)
            lower.p <- upper.p <- lower.c <- upper.c <- matrix(0, nrow(base), numb)
            for(j in 1:numb) {
                mu <- lambda.zero + Lambda %*% alpha[ord[j],  ]
                mean <- base %*% (mu + cov %*% t(basei) %*% (y - basei %*% mu))
                upper.p[, j] <- mean + qnorm(tau2) * pse #mixture likelihood
                lower.p[, j] <- mean - qnorm(tau2) * pse
                upper.c[, j] <- mean + qnorm(tau2) * se #classification likelihood
                lower.c[, j] <- mean - qnorm(tau2) * se
            }
            upci[ind,] <- nummax(upper.c)$max
            lowci[ind,] <-  - nummax( - lower.c)$max
            uppi[ind,] <- nummax(upper.p)$max
            lowpi[ind,] <-  - nummax( - lower.p)$max
            gpred[ind,] <- as.vector(base %*%etapred[ind,])
            ind <- ind+1
        }
        
        meancurves <- base%*%Lambda.alpha
        list(etapred=etapred, gpred=gpred,  upci=upci,lowci=lowci,  uppi=uppi, lowpi=lowpi,index=index,grid=fit$grid,data=data,meancurves=meancurves)
    }

"fitfclust.discrim" <-
    function(fit, absvalue=F){
        base <- fit$base
        sigma <- fit$par$sigma
        nt <- nrow(base)
        Gamma <- fit$par$Gamma
        Sigma <- base%*%Gamma%*%t(base)+sigma*diag(nt)
        Lambda <- fit$par$Lambda
        discrim <- solve(Sigma)%*%base%*%Lambda
        if (absvalue)
            discrim <- abs(discrim)
        n <- ncol(discrim)
        nrows <- ceiling(sqrt(n))
        par(mfrow=squareGrid(n))
        for (i in 1:n){
            plot(fit$grid,discrim[,i], ylab=paste("Discriminant Function ",i),
                 xlab="Time", type='n')
            lines(fit$grid, discrim[,i], lwd=3)
            abline(0,0)}
        par(mfrow=c(1,1))
    }

"fitfclust.plotcurves" <-
    function(object=NULL, fit=NULL, index=NULL, ci=T, pi=T, clustermean=F){
        
        if (is.null(object))
            object <- fitfclust.curvepred(fit)
        if (is.null(index))
            index <- 1:(dim(object$data)[2])
        r <- min(ceiling(sqrt(length(index))),5)
        par(mfrow=c(r,r))
        for (i in index){
            grid <- object$grid
            upci <- object$upci[i,]
            uppi <- object$uppi[i,]
            lowci <- object$lowci[i,]
            lowpi <- object$lowpi[i,]
            gpred <- object$gpred[i,]
            meancurves <- (object$mean)
            if (clustermean){
                yrange <- c(min(c(lowpi,meancurves)),max(c(uppi,meancurves)))
            }  else  yrange <- c(min(lowpi),max(uppi))
            plot(grid,grid,ylim=yrange,ylab="Predictions",xlab="Time",type='n',
                 main=paste("Curve ",i))
            if (clustermean)
                for (k  in 1:ncol(meancurves))
                    lines(grid,meancurves[,k],col=6,lty=2,lwd=2)
            if (ci){
                lines(grid,upci,col=3)
                lines(grid,lowci,col=3)}
            if (pi){
                lines(grid,uppi,col=4)
                lines(grid,lowpi,col=4)}
            lines(grid,gpred,col=2,lwd=2)
            lines(object$grid,object$data[,i],lwd=2)
            points(object$grid,object$data[,i],pch=19,cex=1.5)
        }
        par(mfrow=c(1,1))
    }

############################################################
##Functions independent of data*****************************
"fitfclustconst" <-
    function(parameters, vars, base){
                                        # This function enforces the constraint (7) from the paper on the
                                        # parameters. This means that the alphas can be interpreted as the
                                        # number of standard deviations that the groups are apart etc.
        par <- parameters
        A <- t(base) %*% solve(par$sigma * diag(dim(base)[1]) + base %*% par$Gamma %*%
                                   t(base)) %*% base # t(base)*inv(Sigma)*base
        svdA <- svd(A)
        sqrtA <- diag(sqrt(svdA$d)) %*% t(svdA$u)
        negsqrtA <- svdA$u %*% diag(1/sqrt(svdA$d))
        finalsvd <- svd(sqrtA %*% par$Lambda)
        par$Lambda <- negsqrtA %*% finalsvd$u
        if(dim(par$Lambda)[2] > 1)
            par$alpha <- t(diag(finalsvd$d) %*% t(finalsvd$v) %*% t(par$alpha))
        else par$alpha <- t(finalsvd$d * t(finalsvd$v) %*% t(par$alpha))
        meanalpha <- apply(par$alpha, 2, mean)
        par$alpha <- t(t(par$alpha) - meanalpha)
        par$lambda.zero <- par$lambda.zero + par$Lambda %*% meanalpha
        list(parameters=par, vars=vars)
    }


"nummax" <-
    function(X){
	ind <- rep(1, dim(X)[1])
	m <- X[, 1]
	if(dim(X)[2] > 1)
            for(i in 2:dim(X)[2]) {
                test <- X[, i] > m
                ind[test] <- i
                m[test] <- X[test, i]
            }
	list(ind=ind, max=m)
    }
