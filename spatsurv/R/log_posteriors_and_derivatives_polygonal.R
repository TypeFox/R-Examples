##' logPosterior_polygonal function
##'
##' A function to evaluate the log-posterior of a spatial parametric proportional hazards model. Not intended for general use.
##'
##' @param surv an object of class Surv 
##' @param X the design matrix, containing covariate information 
##' @param beta parameter beta 
##' @param omega parameter omega 
##' @param eta parameter eta 
##' @param gamma parameter gamma 
##' @param priors the priors, an object of class 'mcmcPriors' 
##' @param cov.model the spatial covariance model 
##' @param u vector of interpoint distances 
##' @param control a list containg various control parameters for the MCMC and post-processing routines    
##' @param gradient logical whether to evaluate the gradient
##' @param hessian logical whether to evaluate the Hessian 
##' @return evaluates the log-posterior and the gradient and hessian, if required.
##' @references 
##' \enumerate{
##'     \item Benjamin M. Taylor. Auxiliary Variable Markov Chain Monte Carlo for Spatial Survival and Geostatistical Models. Benjamin M. Taylor. Submitted. \url{http://arxiv.org/abs/1501.01665}
##' }
##' @export

logPosterior_polygonal <- function(surv,X,beta,omega,eta,gamma,priors,cov.model,u,control,gradient=FALSE,hessian=FALSE){

    censoringtype <- control$censoringtype
    censored <- control$censored
    notcensored <- control$notcensored
    Ctest <- control$Ctest
    Utest <- control$Utest
    rightcensored <- control$rightcensored
    notcensored <- control$notcensored
    leftcensored <- control$leftcensored
    intervalcensored <- control$intervalcensored
    Rtest <- control$Rtest
    Utest <- control$Utest
    Ltest <- control$Ltest
    Itest <- control$Itest

    if(hessian){
        gradient <- TRUE
    }

    priorcontrib <- -(1/2)*sum(gamma^2) + do.call(priors$call,args=list(beta=beta,omega=omega,eta=eta,priors=priors))
    if(gradient){
        deriv <- do.call(priors$derivative,args=list(beta=beta,omega=omega,eta=eta,priors=priors))$deriv1 # first derivaties of priors
        deriv <- c(deriv,-gamma) # tag on gamma
        deriv[(length(beta)+length(omega)+1):(length(beta)+length(omega)+length(eta))] <- 0 # random walk for eta ...
    } 
    
    omegaorig <- omega # recall we are working with omega on the transformed scale
    omega <- control$omegaitrans(omega) # this is omega on the correct scale
    
    haz <- setupHazard(dist=control$dist,pars=omega,grad=gradient,hess=hessian)
    
    
    etapars <- cov.model$itrans(eta) 
    sigma <- matrix(EvalCov(cov.model,u=u,parameters=etapars),control$n,control$n)
    cholsigma <- t(chol(sigma))
    MU <- -etapars[control$sigmaidx]^2/2
    Y <- MU + cholsigma%*%gamma
    
    Xbeta <- X%*%beta
    XbetaplusY <- Xbeta + Y[control$idx]
    expXbetaplusY <- exp(XbetaplusY)

    # setup function J=exp(X%*%beta + Y)*H_0(t)
    if(censoringtype=="left" | censoringtype=="right"){
        H <- haz$H(surv[,"time"])
        h <- haz$h(surv[,"time"])
        J <- expXbetaplusY*H
        
        J <- as.vector(J)
    }
    else{ # else interval censored
        H1 <- haz$H(surv[,"time1"])
        H2 <- haz$H(surv[,"time2"])
        J1 <- expXbetaplusY*H1
        J2 <- expXbetaplusY*H2
        
        J1 <- as.vector(J1)
        J2 <- as.vector(J2)
    }    
    
    # setup function S=exp(-J(t))
    if(censoringtype=="right"|censoringtype=="left"){ 
        S <- exp(-J)
    }    
    else{
        S1 <- exp(-J1)
        S2 <- exp(-J2)
    }
    
    
     
    if(censoringtype=="right"){
        h <- haz$h(surv[,"time"])
        
        loglik <-  (if(Utest){sum(XbetaplusY[notcensored] + log(h)[notcensored] - J[notcensored])}else{0}) + 
                    (if(Ctest){sum(-J[censored])}else{0})

        logpost <- loglik + priorcontrib                    
    }
    else if(censoringtype=="left"){
        h <- haz$h(surv[,"time"])
        
        loglik <-  (if(Utest){sum(XbetaplusY[notcensored] + log(h)[notcensored] - J[notcensored])}else{0}) + 
                    (if(Ctest){sum(log(1-S[censored]))}else{0})
        
        logpost <- loglik + priorcontrib
    }
    else{ #censoringtype=="interval"
        h <- haz$h(surv[,"time1"])
        
        loglik <-  (if(Utest){sum(XbetaplusY[notcensored] + log(h)[notcensored] - J1[notcensored])}else{0}) + 
                    (if(Rtest){sum(-J1[rightcensored])}else{0}) + 
                    (if(Ltest){sum(log(1-S1[leftcensored]))}else{0}) + 
                    (if(Itest){sum(log(S1[intervalcensored]-S2[intervalcensored]))}else{0})
        
        logpost <- loglik + priorcontrib
    }

    #browser()

    
    if(gradient){
    
        if(censoringtype=="left" | censoringtype=="right"){
            dJ_dbeta <- J*X
            dJ_domega <- matrix(as.vector(expXbetaplusY)*haz$gradH(surv[,"time"]),ncol=length(omega))
            #dJ_dgamma <- J * cholsigma
        }
        else{ # interval censoring
            dJ_dbeta1 <- J1*X
            dJ_domega1 <- matrix(as.vector(expXbetaplusY)*haz$gradH(surv[,"time1"]),ncol=length(omega))
            #dJ_dgamma1 <- J1 * cholsigma
            
            dJ_dbeta2 <- J2*X
            dJ_domega2 <- matrix(as.vector(expXbetaplusY)*haz$gradH(surv[,"time2"]),ncol=length(omega))
            #dJ_dgamma2 <- J2 * cholsigma
        }               
            
        if(censoringtype=="right"){
            dP_dbeta <- (if(Utest){colSums(X[notcensored,,drop=FALSE] - dJ_dbeta[notcensored,,drop=FALSE])}else{0}) + 
                        (if(Ctest){colSums(-dJ_dbeta[censored,,drop=FALSE])}else{0})
            dh_domega <- matrix(haz$gradh(surv[,"time"]),ncol=length(omega))
            dP_domega <-    (if(Utest){colSums(dh_domega[notcensored,,drop=FALSE] / h[notcensored] - dJ_domega[notcensored,])}else{0}) + 
                            (if(Ctest){colSums(-dJ_domega[censored,,drop=FALSE])}else{0})
            dP_domega <- control$omegajacobian(omegaorig)*dP_domega # this puts the derivative back on the correct scale dL/dpsi = dL/dtheta * dtheta/dpsi, e.g. psi=log(theta)

            bitsnbobs <- rep(0,control$n)
            browser()
            bitsnbobs[control$uqidx] <- bitsnbobs[control$uqidx] + sapply(control$uqidx,function(i){control$sumidxinotcensored[[i]]-sum(J[control$idxi[[i]]])})            
            dP_dgamma <- cholsigma%*%bitsnbobs #as.vector(Re((1/(control$Mext*control$Next))*fft(invrootQeigs*fft(bitsnbobs,inverse=TRUE))))

            #dP_dgamma <-    (if(Utest){colSums(cholsigma[notcensored,,drop=FALSE] - dJ_dgamma[notcensored,,drop=FALSE])}else{0}) + 
            #                (if(Ctest){colSums(-dJ_dgamma[censored,,drop=FALSE])}else{0})
            
            grad <- deriv + c(dP_dbeta,dP_domega,rep(0,length(eta)),dP_dgamma)
            
        }
        else if(censoringtype=="left"){
            dP_dbeta <- (if(Utest){colSums(X[notcensored,,drop=FALSE] - dJ_dbeta[notcensored,,drop=FALSE])}else{0}) + 
                        (if(Ctest){colSums((S[censored]/(1-S[censored]))*dJ_dbeta[censored,,drop=FALSE])}else{0})
            dh_domega <- matrix(haz$gradh(surv[,"time"]),ncol=length(omega))
            dP_domega <-    (if(Utest){colSums(dh_domega[notcensored,,drop=FALSE] / h[notcensored] - dJ_domega[notcensored,,drop=FALSE])}else{0}) + 
                            (if(Ctest){colSums((S[censored]/(1-S[censored]))*dJ_domega[censored,,drop=FALSE])}else{0})
            dP_domega <- control$omegajacobian(omegaorig)*dP_domega # this puts the derivative back on the correct scale dL/dpsi = dL/dtheta * dtheta/dpsi, e.g. psi=log(theta)
            
            bitsnbobs <- rep(0,control$n)
            bitsnbobs[control$uqidx] <- bitsnbobs[control$uqidx] + 
                                        (if(Utest){sapply(control$uqidx,function(i){sum(1-J[control$idxinotcensored[[i]]])})}else{0}) +
                                        (if(Ctest){sapply(control$uqidx,function(i){sum((S[control$idxicensored[[i]]]/(1-S[control$idxicensored[[i]]]))*J[control$idxicensored[[i]]])})}else{0})            
            dP_dgamma <- cholsigma%*%bitsnbobs # as.vector(Re((1/(control$Mext*control$Next))*fft(invrootQeigs*fft(bitsnbobs,inverse=TRUE))))

            #dP_dgamma <-    (if(Utest){colSums(cholsigma[notcensored,,drop=FALSE] - dJ_dgamma[notcensored,,drop=FALSE])}else{0}) + 
            #                (if(Ctest){colSums((S[censored]/(1-S[censored]))*dJ_dgamma[censored,,drop=FALSE])}else{0})
            
            grad <- deriv + c(dP_dbeta,dP_domega,rep(0,length(eta)),dP_dgamma)
        }
        else{ #censoringtype=="interval" 
            dP_dbeta <- (if(Utest){colSums(X[notcensored,,drop=FALSE] - dJ_dbeta1[notcensored,,drop=FALSE])}else{0}) + 
                        (if(Rtest){colSums(-dJ_dbeta1[rightcensored,,drop=FALSE])}else{0}) + 
                        (if(Ltest){colSums((S1[leftcensored]/(1-S1[leftcensored]))*dJ_dbeta1[leftcensored,,drop=FALSE])}else{0}) + 
                        (if(Itest){colSums((1/(S1[intervalcensored]-S2[intervalcensored]))*(dJ_dbeta2[intervalcensored,]*S2[intervalcensored]-dJ_dbeta1[intervalcensored,,drop=FALSE]*S1[intervalcensored]))}else{0})
            dh_domega <- matrix(haz$gradh(surv[,"time1"]),ncol=length(omega))
            dP_domega <-    (if(Utest){colSums(dh_domega[notcensored,,drop=FALSE] / h[notcensored] - dJ_domega1[notcensored,,drop=FALSE])}else{0}) + 
                            (if(Rtest){colSums(-dJ_domega1[rightcensored,,drop=FALSE])}else{0}) + 
                            (if(Ltest){colSums((S1[leftcensored]/(1-S1[leftcensored]))*dJ_domega1[leftcensored,,drop=FALSE])}else{0}) + 
                            (if(Itest){colSums((1/(S1[intervalcensored]-S2[intervalcensored]))*(dJ_domega2[intervalcensored,]*S2[intervalcensored]-dJ_domega1[intervalcensored,,drop=FALSE]*S1[intervalcensored]))}else{0})
            dP_domega <- control$omegajacobian(omegaorig)*dP_domega # this puts the derivative back on the correct scale dL/dpsi = dL/dtheta * dtheta/dpsi, e.g. psi=log(theta)

            bitsnbobs <- rep(0,control$n)
            bitsnbobs[control$uqidx] <- bitsnbobs[control$uqidx] + 
                                        (if(Utest){sapply(control$uqidx,function(i){sum(1-J1[control$idxinotcensored[[i]]])})}else{0}) -
                                        (if(Rtest){sapply(control$uqidx,function(i){sum(-J1[control$idxirightcensored[[i]]])})}else{0}) +  
                                        (if(Ltest){sapply(control$uqidx,function(i){sum((S1[control$idxileftcensored[[i]]]/(1-S1[control$idxileftcensored[[i]]]))*J1[control$idxileftcensored[[i]]])})}else{0}) +                                                   
                                        (if(Itest){sapply(control$uqidx,function(i){sum(((1/(S1-S2))*(J2*S2-J1*S1))[control$idxiintervalcensored[[i]]])})}else{0})
        

            dP_dgamma <- cholsigma%*%bitsnbobs #as.vector(Re((1/(control$Mext*control$Next))*fft(invrootQeigs*fft(bitsnbobs,inverse=TRUE))))  

            # dP_dgamma <-    (if(Utest){colSums(cholsigma[notcensored,,drop=FALSE] - dJ_dgamma1[notcensored,])}else{0}) + 
            #                 (if(Rtest){colSums(-dJ_dgamma1[rightcensored,,drop=FALSE])}else{0}) + 
            #                 (if(Ltest){colSums((S1[leftcensored]/(1-S1[leftcensored]))*dJ_dgamma1[leftcensored,,drop=FALSE])}else{0}) + 
            #                 (if(Itest){colSums((1/(S1[intervalcensored]-S2[intervalcensored]))*(dJ_dgamma2[intervalcensored,,drop=FALSE]*S2[intervalcensored]-dJ_dgamma1[intervalcensored,,drop=FALSE]*S1[intervalcensored]))}else{0})
            
            grad <- deriv + c(dP_dbeta,dP_domega,rep(0,length(eta)),dP_dgamma)
        }
    
    }
    
    
    if(hessian){
    
        cross_gradh <- lapply(1:nrow(X),function(i){outer(dh_domega[i,],dh_domega[i,])}) 
    
        if(censoringtype=="left" | censoringtype=="right"){
            d2h_domega1_domega2 <- haz$hessh(surv[,"time"])                    
                    
            d2J_dbetak_dbetaj <- lapply(1:nrow(X),function(i){outer(X[i,],X[i,])*J[i]})
            d2J_domega_dbeta <- lapply(1:length(omega),function(i){dJ_domega[,i]*X})
            d2J_domega1_domega2 <- mapply("*", haz$hessH(surv[,"time"]), expXbetaplusY, SIMPLIFY = FALSE)
            #d2J_dgamma2 <- J*sigma

            dS_dbeta <- -S*dJ_dbeta
            dS_domega <- -S*dJ_domega
            #dS_dgamma <- -S*dJ_dgamma             
            
            d2S_dbetak_dbetaj <- mapply('*',mapply('-',lapply(1:nrow(X),function(i){outer(dJ_dbeta[i,],dJ_dbeta[i,])}),d2J_dbetak_dbetaj,SIMPLIFY=FALSE),S,SIMPLIFY=FALSE)
            if(length(beta)>1){
                d2S_domega_dbeta <- mapply('*',lapply(1:nrow(X),function(i){outer(dJ_domega[i,],dJ_dbeta[i,])-t(sapply(d2J_domega_dbeta,function(x){x[i,]}))}),S,SIMPLIFY=FALSE)
            }
            else{
                d2S_domega_dbeta <- mapply('*',lapply(1:nrow(X),function(i){outer(dJ_domega[i,],dJ_dbeta[i,])-t(t(sapply(d2J_domega_dbeta,function(x){x[i,]})))}),S,SIMPLIFY=FALSE)
            }

            d2S_domega1_domega2 <- mapply('*',mapply('-',lapply(1:nrow(X),function(i){outer(dJ_domega[i,],dJ_domega[i,])}),d2J_domega1_domega2,SIMPLIFY=FALSE),S,SIMPLIFY=FALSE)
            #d2S_dgamma2 <- S*(J^2*sigma - d2J_dgamma2) #S*(dJ_dgamma^2 - d2J_dgamma2)

            cross_dS_dbeta <- lapply(1:nrow(X),function(i){outer(dS_dbeta[i,],dS_dbeta[i,])})
            cross_dS_domega_dbeta <- lapply(1:nrow(X),function(i){outer(dS_domega[i,],dS_dbeta[i,])})
            cross_dS_domega <- lapply(1:nrow(X),function(i){outer(dS_domega[i,],dS_domega[i,])})
            #cross_dS_dgamma <- J^2*S^2*sigma #dS_dgamma^2 # not considering cross product terms here
        } 
        else{ # censoringtype=="interval"
            d2h_domega1_domega2 <- haz$hessh(surv[,"time1"])        
        
            d2J_dbetak_dbetaj_1 <- lapply(1:nrow(X),function(i){outer(X[i,],X[i,])*J1[i]})
            d2J_domega_dbeta_1 <- lapply(1:length(omega),function(i){dJ_domega1[,i]*X})
            d2J_domega1_domega2_1 <- mapply("*", haz$hessH(surv[,"time1"]), expXbetaplusY, SIMPLIFY = FALSE)
            #d2J_dgamma2_1 <- J1*sigma
            
            d2J_dbetak_dbetaj_2 <- lapply(1:nrow(X),function(i){outer(X[i,],X[i,])*J2[i]})
            d2J_domega_dbeta_2 <- lapply(1:length(omega),function(i){dJ_domega2[,i]*X})
            d2J_domega1_domega2_2 <- mapply("*", haz$hessH(surv[,"time2"]), expXbetaplusY, SIMPLIFY = FALSE)
            #d2J_dgamma2_2 <- J2*sigma
            
            
            
            dS_dbeta_1 <- -S1*dJ_dbeta1
            dS_domega_1 <- -S1*dJ_domega1
            #dS_dgamma_1 <- -S1*dJ_dgamma1             
            
            d2S_dbetak_dbetaj_1 <- mapply('*',mapply('-',lapply(1:nrow(X),function(i){outer(dJ_dbeta1[i,],dJ_dbeta1[i,])}),d2J_dbetak_dbetaj_1,SIMPLIFY=FALSE),S1,SIMPLIFY=FALSE)
            if(length(beta)>1){
                d2S_domega_dbeta_1 <- mapply('*',lapply(1:nrow(X),function(i){outer(dJ_domega1[i,],dJ_dbeta1[i,])-t(sapply(d2J_domega_dbeta_1,function(x){x[i,]}))}),S1,SIMPLIFY=FALSE)
            }
            else{
                d2S_domega_dbeta_1 <- mapply('*',lapply(1:nrow(X),function(i){outer(dJ_domega1[i,],dJ_dbeta1[i,])-t(t(sapply(d2J_domega_dbeta_1,function(x){x[i,]})))}),S1,SIMPLIFY=FALSE)   
            }
            d2S_domega1_domega2_1 <- mapply('*',mapply('-',lapply(1:nrow(X),function(i){outer(dJ_domega1[i,],dJ_domega1[i,])}),d2J_domega1_domega2_1,SIMPLIFY=FALSE),S1,SIMPLIFY=FALSE)
            #d2S_dgamma2_1 <- S1*(J1^2*sigma - d2J_dgamma2_1) #S1*(dJ_dgamma1^2 - d2J_dgamma2_1)
            
            
            
            dS_dbeta_2 <- -S2*dJ_dbeta2
            dS_domega_2 <- -S2*dJ_domega2
            #S_dgamma_2 <- -S2*dJ_dgamma2             
            
            d2S_dbetak_dbetaj_2 <- mapply('*',mapply('-',lapply(1:nrow(X),function(i){outer(dJ_dbeta2[i,],dJ_dbeta2[i,])}),d2J_dbetak_dbetaj_2,SIMPLIFY=FALSE),S2,SIMPLIFY=FALSE)
            if(length(beta)>1){
                d2S_domega_dbeta_2 <- mapply('*',lapply(1:nrow(X),function(i){outer(dJ_domega2[i,],dJ_dbeta2[i,])-t(sapply(d2J_domega_dbeta_2,function(x){x[i,]}))}),S2,SIMPLIFY=FALSE)
            }
            else{
                d2S_domega_dbeta_2 <- mapply('*',lapply(1:nrow(X),function(i){outer(dJ_domega2[i,],dJ_dbeta2[i,])-t(t(sapply(d2J_domega_dbeta_2,function(x){x[i,]})))}),S2,SIMPLIFY=FALSE)   
            }
            d2S_domega1_domega2_2 <- mapply('*',mapply('-',lapply(1:nrow(X),function(i){outer(dJ_domega2[i,],dJ_domega2[i,])}),d2J_domega1_domega2_2,SIMPLIFY=FALSE),S2,SIMPLIFY=FALSE)
            #d2S_dgamma2_2 <- S2*(J2^2*sigma - d2J_dgamma2_2)
            
    
            
            cross_dS_dbeta_1 <- lapply(1:nrow(X),function(i){outer(dS_dbeta_1[i,],dS_dbeta_1[i,])})
            cross_dS_domega_dbeta_1 <- lapply(1:nrow(X),function(i){outer(dS_domega_1[i,],dS_dbeta_1[i,])})
            cross_dS_domega_1 <- lapply(1:nrow(X),function(i){outer(dS_domega_1[i,],dS_domega_1[i,])})
            #cross_dS_dgamma_1 <- J1^2*S1^2*sigma #dS_dgamma_1^2 # not considering cross product terms here
            
            cross_dS_dbeta_2 <- lapply(1:nrow(X),function(i){outer(dS_dbeta_2[i,],dS_dbeta_2[i,])})
            cross_dS_domega_dbeta_2 <- lapply(1:nrow(X),function(i){outer(dS_domega_2[i,],dS_dbeta_2[i,])})
            cross_dS_domega_2 <- lapply(1:nrow(X),function(i){outer(dS_domega_2[i,],dS_domega_2[i,])})
            #cross_dS_dgamma_2 <- J2^2*S2^2*sigma #dS_dgamma_2^2 # not considering cross product terms here
            
            cross_dS_dbeta_12 <- lapply(1:nrow(X),function(i){outer(dS_dbeta_1[i,],dS_dbeta_2[i,])})
            cross_dS_domega_dbeta_12 <- lapply(1:nrow(X),function(i){outer(dS_domega_1[i,],dS_dbeta_2[i,])})
            cross_dS_domega_12 <- lapply(1:nrow(X),function(i){outer(dS_domega_1[i,],dS_domega_2[i,])})
            #cross_dS_dgamma_12 <- J1*J2*S1*S2*sigma #dS_dgamma_1*dS_dgamma_2
            
            cross_dS_dbeta_21 <- lapply(1:nrow(X),function(i){outer(dS_dbeta_2[i,],dS_dbeta_1[i,])})
            cross_dS_domega_dbeta_21 <- lapply(1:nrow(X),function(i){outer(dS_domega_2[i,],dS_dbeta_1[i,])})
            cross_dS_domega_21 <- lapply(1:nrow(X),function(i){outer(dS_domega_2[i,],dS_domega_1[i,])})
            #cross_dS_dgamma_21 <- J2*J1*S2*S1*sigma #dS_dgamma_2*dS_dgamma_1
        }
        
        
        if(censoringtype=="right"){
            hess_beta <- (if(Utest){-Reduce('+',d2J_dbetak_dbetaj[notcensored])}else{0}) - 
                                    (if(Ctest){Reduce('+',d2J_dbetak_dbetaj[censored])}else{0})
            
            hess_omega <- (if(Utest){Reduce('+',mapply('*',d2h_domega1_domega2[notcensored],1/h[notcensored],SIMPLIFY=FALSE)) - Reduce('+',mapply('*',cross_gradh[notcensored],1/h[notcensored]^2,SIMPLIFY=FALSE)) - Reduce('+',d2J_domega1_domega2[notcensored])}else{0}) -
                                    (if(Ctest){Reduce('+',d2J_domega1_domega2[censored])}else{0})                                             
            hess_omega <- diag((dP_domega/control$omegajacobian(omegaorig))*sapply(1:length(omega),function(i){control$omegahessian[[i]](omegaorig[i])}),length(omega)) + hess_omega * outer(control$omegajacobian(omegaorig),control$omegajacobian(omegaorig))
            # note, have dP_domega/control$omegajacobian(omegaorig) to get onto correct scale since further above in the computation of dP_domega we have control$omegajacobian(omegaorig)*dP_domega                                      

            hess_omega_beta <- (if(Utest){-t(sapply(d2J_domega_dbeta,function(x){colSums(x[notcensored,,drop=FALSE])}))}else{0}) - 
                                    (if(Ctest){t(sapply(d2J_domega_dbeta,function(x){colSums(x[censored,,drop=FALSE])}))}else{0})
            hess_omega_beta <- control$omegajacobian(omegaorig)*hess_omega_beta                        
            
            bitsnbobs <- rep(0,control$n)
            bitsnbobs[control$uqidx] <- bitsnbobs[control$uqidx] + sapply(control$uqidx,function(i){sum(-J[control$idxi[[i]]])})            
            hess_gamma <- sigma%*%bitsnbobs #as.vector(Re((1/(control$Mext*control$Next))*fft(covbase*fft(bitsnbobs,inverse=TRUE))))

            #hess_gamma <- (if(Utest){-colSums(d2J_dgamma2[notcensored,,drop=FALSE])}else{0}) - 
            #                        (if(Ctest){colSums(d2J_dgamma2[censored,,drop=FALSE])}else{0})
        }
        else if(censoringtype=="left"){
            hess_beta <- (if(Utest){-Reduce('+',d2J_dbetak_dbetaj[notcensored])}else{0}) - 
                                    (if(Ctest){Reduce('+',mapply('*',d2S_dbetak_dbetaj[censored],1/(1-S[censored]),SIMPLIFY=FALSE)) - Reduce('+',mapply('*',cross_dS_dbeta[censored],1/(1-S[censored])^2,SIMPLIFY=FALSE))}else{0})            
            
            hess_omega <- (if(Utest){Reduce('+',mapply('*',d2h_domega1_domega2[notcensored],1/h[notcensored],SIMPLIFY=FALSE)) - Reduce('+',mapply('*',cross_gradh[notcensored],1/h[notcensored]^2,SIMPLIFY=FALSE)) - Reduce('+',d2J_domega1_domega2[notcensored])}else{0}) -
                                    (if(Ctest){Reduce('+',mapply('*',d2S_domega1_domega2[censored],1/(1-S[censored]),SIMPLIFY=FALSE)) - Reduce('+',mapply('*',cross_dS_domega[censored],1/(1-S[censored])^2,SIMPLIFY=FALSE))}else{0})
            hess_omega <- diag((dP_domega/control$omegajacobian(omegaorig))*sapply(1:length(omega),function(i){control$omegahessian[[i]](omegaorig[i])}),length(omega)) + hess_omega * outer(control$omegajacobian(omegaorig),control$omegajacobian(omegaorig))
            # note, have dP_domega/control$omegajacobian(omegaorig) to get onto correct scale since further above in the computation of dP_domega we have control$omegajacobian(omegaorig)*dP_domega                                    
            
            hess_omega_beta <- (if(Utest){-t(sapply(d2J_domega_dbeta,function(x){colSums(x[notcensored,,drop=FALSE])}))}else{0}) - 
                                    (if(Ctest){Reduce('+',mapply('*',d2S_domega_dbeta[censored],1/(1-S[censored]),SIMPLIFY=FALSE)) - Reduce('+',mapply('*',cross_dS_domega_dbeta[censored],1/(1-S[censored])^2,SIMPLIFY=FALSE))}else{0})
            hess_omega_beta <- control$omegajacobian(omegaorig)*hess_omega_beta                                     
            
            bitsnbobs <- rep(0,control$n)
            bitsnbobs[control$uqidx] <- bitsnbobs[control$uqidx] + 
                                    (if(Utest){sapply(control$uqidx,function(i){sum(-J[control$idxinotcensored[[i]]])})}else{0}) -
                                    (if(Ctest){sapply(control$uqidx,function(i){sum((S*(J^2-J)/(1-S)+(J*S/(1-S))^2)[control$idxicensored[[i]]])})}else{0})           
            hess_gamma <- sigma%*%bitsnbobs #as.vector(Re((1/(control$Mext*control$Next))*fft(covbase*fft(bitsnbobs,inverse=TRUE))))

            #hess_gamma <- (if(Utest){-colSums(d2J_dgamma2[notcensored,,drop=FALSE])}else{0}) -
            #                        (if(Ctest){colSums((1/(1-S[censored]))*d2S_dgamma2[censored,,drop=FALSE]) - colSums((1/(1-S[censored])^2)*cross_dS_dgamma[censored,,drop=FALSE])}else{0})
        }
        else{ # censoringtype=="interval"        
            hess_beta <- (if(Utest){-Reduce('+',d2J_dbetak_dbetaj_1[notcensored])}else{0}) - 
                                    (if(Rtest){Reduce('+',d2J_dbetak_dbetaj_1[rightcensored])}else{0}) - 
                                    (if(Ltest){Reduce('+',mapply('*',d2S_dbetak_dbetaj_1[leftcensored],1/(1-S1[leftcensored]),SIMPLIFY=FALSE)) - Reduce('+',mapply('*',cross_dS_dbeta_1[leftcensored],1/(1-S1[leftcensored])^2,SIMPLIFY=FALSE))}else{0}) +
                                    (if(Itest){Reduce('+',mapply('*',mapply('-',d2S_dbetak_dbetaj_1[intervalcensored],d2S_dbetak_dbetaj_2[intervalcensored],SIMPLIFY=FALSE),1/(S1[intervalcensored]-S2[intervalcensored]),SIMPLIFY=FALSE)) - Reduce('+',mapply('*',mapply(function(a,b,c,d){a-b-c+d},cross_dS_dbeta_1[intervalcensored],cross_dS_dbeta_12[intervalcensored],cross_dS_dbeta_21[intervalcensored],cross_dS_dbeta_2[intervalcensored],SIMPLIFY=FALSE),1/(S1[intervalcensored]-S2[intervalcensored])^2,SIMPLIFY=FALSE))}else{0})
            
            
            hess_omega <- (if(Utest){Reduce('+',mapply('*',d2h_domega1_domega2[notcensored],1/h[notcensored],SIMPLIFY=FALSE)) - Reduce('+',mapply('*',cross_gradh[notcensored],1/h[notcensored]^2,SIMPLIFY=FALSE)) - Reduce('+',d2J_domega1_domega2_1[notcensored])}else{0}) - 
                                    (if(Rtest){Reduce('+',d2J_domega1_domega2_1[rightcensored])}else{0}) - 
                                    (if(Ltest){Reduce('+',mapply('*',d2S_domega1_domega2_1[leftcensored],1/(1-S1[leftcensored]),SIMPLIFY=FALSE)) - Reduce('+',mapply('*',cross_dS_domega_1[leftcensored],1/(1-S1[leftcensored])^2,SIMPLIFY=FALSE))}else{0}) +
                                    (if(Itest){Reduce('+',mapply('*',mapply('-',d2S_domega1_domega2_1[intervalcensored],d2S_domega1_domega2_2[intervalcensored],SIMPLIFY=FALSE),1/(S1[intervalcensored]-S2[intervalcensored]),SIMPLIFY=FALSE)) - Reduce('+',mapply('*',mapply(function(a,b,c,d){a-b-c+d},cross_dS_domega_1[intervalcensored],cross_dS_domega_12[intervalcensored],cross_dS_domega_21[intervalcensored],cross_dS_domega_2[intervalcensored],SIMPLIFY=FALSE),1/(S1[intervalcensored]-S2[intervalcensored])^2,SIMPLIFY=FALSE))}else{0}) 
            hess_omega <- diag((dP_domega/control$omegajacobian(omegaorig))*sapply(1:length(omega),function(i){control$omegahessian[[i]](omegaorig[i])}),length(omega)) + hess_omega * outer(control$omegajacobian(omegaorig),control$omegajacobian(omegaorig))
            # note, have dP_domega/control$omegajacobian(omegaorig) to get onto correct scale since further above in the computation of dP_domega we have control$omegajacobian(omegaorig)*dP_domega
            
            hess_omega_beta <- (if(Utest){-t(sapply(d2J_domega_dbeta_1,function(x){colSums(x[notcensored,,drop=FALSE])}))}else{0}) - 
                                    (if(Rtest){t(sapply(d2J_domega_dbeta_1,function(x){colSums(x[rightcensored,,drop=FALSE])}))}else{0}) - 
                                    (if(Ltest){Reduce('+',mapply('*',d2S_domega_dbeta_1[leftcensored],1/(1-S1[leftcensored]),SIMPLIFY=FALSE)) - Reduce('+',mapply('*',cross_dS_domega_dbeta_1[leftcensored],1/(1-S1[leftcensored])^2,SIMPLIFY=FALSE))}else{0}) +
                                    (if(Itest){Reduce('+',mapply('*',mapply('-',d2S_domega_dbeta_1[intervalcensored],d2S_domega_dbeta_2[intervalcensored],SIMPLIFY=FALSE),1/(S1[intervalcensored]-S2[intervalcensored]),SIMPLIFY=FALSE)) - Reduce('+',mapply('*',mapply(function(a,b,c,d){a-b-c+d},cross_dS_domega_dbeta_1[intervalcensored],cross_dS_domega_dbeta_12[intervalcensored],cross_dS_domega_dbeta_21[intervalcensored],cross_dS_domega_dbeta_2[intervalcensored],SIMPLIFY=FALSE),1/(S1[intervalcensored]-S2[intervalcensored])^2,SIMPLIFY=FALSE))}else{0})
            hess_omega_beta <- control$omegajacobian(omegaorig)*hess_omega_beta 
            
            bitsnbobs <- rep(0,control$n)
            
            bitsnbobs[control$uqidx] <- bitsnbobs[control$uqidx] + 
                                    (if(Utest){sapply(control$uqidx,function(i){sum(-J1[control$idxinotcensored[[i]]])})}else{0}) +
                                    (if(Rtest){sapply(control$uqidx,function(i){sum(-J1[control$idxirightcensored[[i]]])})}else{0}) -
                                    (if(Ltest){sapply(control$uqidx,function(i){sum((S1*(J1^2-J1)/(1-S1)+(J1*S1/(1-S1))^2)[control$idxileftcensored[[i]]])})}else{0}) +
                                    (if(Itest){sapply(control$uqidx,function(i){sum(((1/(S1-S2))*(S1*(J1^2-J1)-S2*(J2^2-J2))-(1/(S1-S2)^2)*(J1^2*S1^2-2*J1*J2*S1*S2+J2^2*S2^2))[control$idxiintervalcensored[[i]]])})}else{0})          
            hess_gamma <- sigma%*%bitsnbobs # as.vector(Re((1/(control$Mext*control$Next))*fft(covbase*fft(bitsnbobs,inverse=TRUE))))

            # hess_gamma <- (if(Utest){-colSums(d2J_dgamma2_1[notcensored,,drop=FALSE])}else{0}) - 
            #                         (if(Rtest){colSums(d2J_dgamma2_1[rightcensored,,drop=FALSE])}else{0}) -
            #                         (if(Ltest){colSums((1/(1-S1[leftcensored]))*d2S_dgamma2_1[leftcensored,,drop=FALSE]) - colSums((1/(1-S1[leftcensored])^2)*cross_dS_dgamma_1[leftcensored,,drop=FALSE])}else{0}) +
            #                         (if(Itest){colSums((1/(S1[intervalcensored]-S2[intervalcensored]))*(d2S_dgamma2_1-d2S_dgamma2_2)[intervalcensored,,drop=FALSE]) - colSums((1/(S1[intervalcensored]-S2[intervalcensored])^2)*(cross_dS_dgamma_1-cross_dS_dgamma_12-cross_dS_dgamma_21+cross_dS_dgamma_2)[intervalcensored,,drop=FALSE])}else{0})
        }
        
        
        deriv2 <- do.call(priors$derivative,args=list(beta=beta,omega=omegaorig,eta=eta,priors=priors))$deriv2
        
        #tag on contributions from the prior ...
        hess_beta <- hess_beta + deriv2[1:length(beta),1:length(beta)]
        hess_omega <- hess_omega + deriv2[(length(beta)+1):(length(beta)+length(omega)),(length(beta)+1):(length(beta)+length(omega))]
        hess_gamma <- hess_gamma - 1 # -1 comes from the prior
    }  
      
    
    if(!gradient & !hessian){
        return(logpost)
    } 
    else{
        retlist <- list()
        retlist$logpost <- logpost
        retlist$loglik <- loglik
        retlist$Y <- Y
        if(gradient){
            retlist$grad <- grad
        }
        if(hessian){
            retlist$hess_beta <- hess_beta
            retlist$hess_omega <- hess_omega
            retlist$hess_omega_beta <- hess_omega_beta
            retlist$hess_gamma <- hess_gamma
        }
        return(retlist)
    }
}

