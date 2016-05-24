#
#Copyright (c) 2012 Madison Giacofci, Sophie Lambert-Lacroix, Guillemette Marot, Franck Picard
#
#

setGeneric(name = "getFCMM", def = function(.Object,CCO){standardGeneric("getFCMM")})

setMethod(
    f          = "getFCMM",
    signature  = "CClustData",
    definition =
        function(.Object,CCO) {

            if (CCO@nbclust==1){
                stop("[if CCO@nbclust==1 use getFMM instead]")
            }
            if (CCO@Gamma2.structure=="none"){
                stop("[if CCO@Gamma2.structure==\"none\" use getFCM instead]")
            }
            
            #cat("[ initialization with",CCO@burn,"burns using",CCO@init,"]\n") 
            llik      = c()
            loglikbm1 = -Inf
            CCR       = new("CClustRes",CCO,.Object)                        
            CCRu      = new("CClustRes",CCO,.Object)
            CCRum1    = new("CClustRes",CCO,.Object)
            
            if (CCO@init=="SEM"){
                b = 0
                while ( (b < CCO@burn)& (min(CCR@prop)>1e-10) & ( CCR@varE>1e-10 ) ){
                    b   = b+1
                    CCR = SstepFCMM(CCR,CCO,.Object)
                    CCR = MstepFCMM(CCR,CCO,.Object)
                    CCR = EstepFCMM(CCR,CCO,.Object)

                    if ((CCR@loglik> loglikbm1)) {
                        CCRum1    = CCRu
                        CCRu      = CCR
                        loglikbm1 = CCR@loglik
                    }
                    llik = c(llik,CCR@loglik)
                }
                if (min(CCR@prop)<1e-10){stop("[SEM initializations resulted in empty clusters]\n")}
                if (CCR@varE<1e-10){stop("[SEM initializations resulted in null variance]\n")}
                
            } else if (CCO@init=="rEM"){
                dv = 0
                for (nbseeds in c(1:CCO@burn)){
                    CCR  = new("CClustRes",CCO,.Object)
                    CCR  = SstepFCMM(CCR,CCO,.Object)
                    i    = 0
                    while ( (i < 10) & (min(CCR@prop)>1e-10) & ( CCR@varE>1e-10 ) ){
                        i   = i+1
                        CCR = MstepFCMM(CCR,CCO,.Object)
                        CCR = EstepFCMM(CCR,CCO,.Object)
                    }
                    if (min(CCR@prop)<1e-10){CCR@loglik=-Inf; dv=dv+1}
                    if (CCR@varE <1e-10){CCR@loglik=-Inf; dv=dv+1}
                    if (CCR@loglik> loglikbm1){
                        CCRum1    = CCRu
                        CCRu      = CCR
                        loglikbm1 = CCR@loglik
                    }
                    llik = c(llik,CCR@loglik)
                }
                if (dv>=CCO@burn){stop("[All rEM initializations resulted in empty clusters or null variance]\n")}
            }

            #cat("[EM algorithm]\n")            
            CCR           = CCRu            
            theta.h       = gettheta(CCRu,CCO)
            theta.hm1     = gettheta(CCRum1,CCO)
            theta.dot.hm2 = gettheta(CCRu,CCO)
            epsilon       = Inf
            
            while ( (epsilon>CCO@eps) & (min(CCR@prop)>1e-6) & ( CCR@varE>1e-10 )){
                
                CCR           = MstepFCMM(CCR,CCO,.Object)
                CCR           = EstepFCMM(CCR,CCO,.Object)
                
                theta.hp1     = gettheta(CCR,CCO)              
                theta.dot.hm1 = theta.h + invnorm( invnorm(theta.hm1-theta.h) + invnorm( theta.hp1- theta.h) )
                epsilon       = sum( (theta.dot.hm1-theta.dot.hm2)^2 / (theta.dot.hm1)^2 )
                theta.hm1     = theta.h
                theta.h       = theta.hp1
                theta.dot.hm2 = theta.dot.hm1
                llik          = c(llik,CCR@loglik)

            }

            if (min(CCR@prop)<1e-10){cat("[EM stopped because of empty clusters]")}
            if (CCR@varE <1e-10){cat("[EM stopped because of null variance]")}
            if (CCO@loglikplot){plot(llik,type="b");abline(v=CCO@burn)}
            
            return(CCR)
        })

invnorm <- function(x){
    y = x/sum(x^2)
    y[is.na(y)] = 0
    return(y)
}

setGeneric( name = "MstepFCMM"     ,def =
               function(.Object,CCO,CCD){standardGeneric("MstepFCMM")})

setMethod(
    f          = "MstepFCMM",
    signature  = "CClustRes",
    definition =
        function(.Object,CCO,CCD) {

            ## avoid label switching, sort groups by levels
            tmp           = order(unlist(.Object@Alpha))
            .Object@Beta  = .Object@Beta[tmp]
            .Object@Alpha = .Object@Alpha[tmp]
            .Object@prop  = .Object@prop[tmp]
            .Object@Tau   = .Object@Tau[,tmp]
            .Object@Blup  = .Object@Blup[tmp]
            
            n             = length(CCD@wavecoef)
            M             = length(CCD@wavecoef[[1]]$D) + 1
            Wd            = t(sapply(CCD@wavecoef,FUN=function(x){x$D},USE.NAMES=TRUE))
            Wc            = t(t(sapply(CCD@wavecoef,FUN=function(x){x$C},USE.NAMES=TRUE)))
            Hcoef         = 2^-(CCD@jklevels[,1]*.Object@eta)

            clustfreq     = apply(.Object@Tau,2,sum)
            .Object@prop  = as.vector(clustfreq) / n 
            .Object@Beta  = lapply(1:CCO@nbclust,FUN = function(ell){
                apply( .Object@Tau[,ell]* (Wd-.Object@Blup[[ell]]$Theta),2,sum)/clustfreq[ell]
            })
            .Object@Alpha = lapply(1:CCO@nbclust,FUN = function(ell){
                sum(.Object@Tau[,ell]*(Wc-.Object@Blup[[ell]]$Nu))/clustfreq[ell]
            })

            if( (CCO@Gamma2.structure=="group.scale.location") | (CCO@Gamma2.structure=="group") ){

                .Object@Gamma2.Theta = .Object@Gamma2.Theta[tmp]
                .Object@Gamma2.Nu    = .Object@Gamma2.Nu[tmp]
                
                lambda.Theta = lapply(1:CCO@nbclust,FUN = function(ell){.Object@varE /.Object@Gamma2.Theta[[ell]]})
                lambda.Nu    = lapply(1:CCO@nbclust,FUN = function(ell){.Object@varE /.Object@Gamma2.Nu[[ell]]})
                
                varE.D =  sapply(1:CCO@nbclust,FUN = function(ell){
                    RSS = apply(Wd-.Object@Blup[[ell]]$Theta,1,FUN = function(x){x-.Object@Beta[[ell]]})^2 + .Object@varE/(1+lambda.Theta[[ell]]/Hcoef) 
                    sum(.Object@Tau[,ell]*t(RSS))
                },simplify=FALSE,USE.NAMES=TRUE)
                
                varE.C =  sapply(1:CCO@nbclust,FUN = function(ell){
                    RSS = (Wc-.Object@Blup[[ell]]$Nu-.Object@Alpha[[ell]])^2 + .Object@varE/(1+lambda.Nu[[ell]])
                    sum(.Object@Tau[,ell]*RSS)
                },simplify=FALSE,USE.NAMES=TRUE)
                
                .Object@varE = (Reduce("+",varE.D) + Reduce("+",varE.C)) / (n*M)
                
                .Object@Gamma2.Nu = sapply(1:CCO@nbclust,FUN = function(ell){
                    Z = .Object@Blup[[ell]]$Nu^2 + .Object@varE /(1+lambda.Nu[[ell]])
                    sum(.Object@Tau[,ell]*Z) /clustfreq[ell]
                },simplify=FALSE,USE.NAMES=TRUE)
                
                Gamma2.Theta = sapply(1:CCO@nbclust,FUN = function(ell){
                    Z = (t(.Object@Blup[[ell]]$Theta^2) + .Object@varE /(1+lambda.Theta[[ell]]/Hcoef) ) / Hcoef
                    return(t(Z)*.Object@Tau[,ell])
                },simplify=FALSE,USE.NAMES=TRUE)

                if( (CCO@Gamma2.structure=="group.scale.location") ){
                    
                    .Object@Gamma2.Theta = sapply(1:CCO@nbclust,FUN = function(ell){
                        apply(Gamma2.Theta[[ell]],2,sum)/clustfreq[ell]
                    },simplify=FALSE,USE.NAMES=TRUE)
                    
                    
                } else if (CCO@Gamma2.structure=="group"){
                    
                    .Object@Gamma2.Theta = sapply(1:CCO@nbclust,FUN = function(ell){
                        sum(Gamma2.Theta[[ell]])/(clustfreq[ell]*(M-1))
                    },simplify=FALSE,USE.NAMES=TRUE)
                    
                }
            }
            
            if ((CCO@Gamma2.structure=="scale.location") | (CCO@Gamma2.structure=="constant")) {
                
                lambda.Theta = .Object@varE/.Object@Gamma2.Theta
                lambda.Nu    = .Object@varE/.Object@Gamma2.Nu                             
                
                varE.D =  sapply(1:CCO@nbclust,FUN = function(ell){
                    RSS = apply(Wd-.Object@Blup[[ell]]$Theta,1,FUN = function(x){x-.Object@Beta[[ell]]})^2 + .Object@varE/(1+lambda.Theta/Hcoef) 
                    sum(.Object@Tau[,ell]*t(RSS))
                },simplify=FALSE,USE.NAMES=TRUE)
                
                varE.C =  sapply(1:CCO@nbclust,FUN = function(ell){
                    RSS = (Wc-.Object@Blup[[ell]]$Nu-.Object@Alpha[[ell]])^2 + .Object@varE/(1+lambda.Nu)
                    sum(.Object@Tau[,ell]*RSS)
                },simplify=FALSE,USE.NAMES=TRUE)
                
                .Object@varE = (Reduce("+",varE.D) + Reduce("+",varE.C)) / (n*M)
                
                Gamma2.Nu = sapply(1:CCO@nbclust,FUN = function(ell){
                    sum(.Object@Tau[,ell]*.Object@Blup[[ell]]$Nu^2 + .Object@varE /(1+lambda.Nu))
                },simplify=FALSE,USE.NAMES=TRUE)
                .Object@Gamma2.Nu    = Reduce("+",Gamma2.Nu)/n

                Gamma2.Theta = sapply(1:CCO@nbclust,FUN = function(ell){
                    Z = (t(.Object@Blup[[ell]]$Theta^2) + .Object@varE /(1+lambda.Theta/Hcoef) ) / Hcoef
                    return(t(Z)*.Object@Tau[,ell])
                },simplify=FALSE,USE.NAMES=TRUE)
                
                if (CCO@Gamma2.structure=="scale.location"){
                    .Object@Gamma2.Theta = apply(Reduce("+",Gamma2.Theta),2,sum)/n
                }
                if (CCO@Gamma2.structure=="constant"){                
                    .Object@Gamma2.Theta = sum(Reduce("+",Gamma2.Theta))/(n*(M-1))
                }
            }
            .Object@eta = getetaFCMM(.Object,CCD,CCO)
            return(.Object) 
        })


setGeneric( name = "getllkFCMM" ,def = function(.Object,CCD,CCO,eta){standardGeneric("getllkFCMM")})
setMethod(
    f          = "getllkFCMM",
    signature  = "CClustRes",
    definition =
        function(.Object,CCD,CCO,eta) {
            
            Hcoef = 2^-(CCD@jklevels$scale*eta)
            
            if (CCO@nbclust==1){
                vardens.C           = .Object@Gamma2.Nu+.Object@varE
                vardens.D           = .Object@Gamma2.Theta*Hcoef+.Object@varE    
                logdens   = sapply(CCD@wavecoef, FUN = function(x){
                    sum(dnorm( c(x$C,x$D),c(.Object@Alpha,.Object@Beta) , sd=sqrt(c(vardens.C,vardens.D)),log=TRUE))
                })
                llk= sum(logdens)
            } else{        
                if( (CCO@Gamma2.structure=="group.scale.location") | (CCO@Gamma2.structure=="group") ){
                    vardens.C   = sapply(1:CCO@nbclust,FUN=function(ell){.Object@Gamma2.Nu[[ell]]+.Object@varE})
                    vardens.D   = sapply(1:CCO@nbclust,FUN=function(ell){.Object@Gamma2.Theta[[ell]]*Hcoef+.Object@varE})
                }    
                if ((CCO@Gamma2.structure=="scale.location") | (CCO@Gamma2.structure=="constant")) {
                    vardens.C   = .Object@Gamma2.Nu+.Object@varE
                    vardens.D   = .Object@Gamma2.Theta*Hcoef+.Object@varE
                }    
                logdens   = sapply(CCD@wavecoef, FUN = function(x){
                    apply( dnorm( c(x$C,x$D),sapply(1:CCO@nbclust,FUN = function(ell){c(.Object@Alpha[[ell]],.Object@Beta[[ell]])}) , sd=sqrt(c(vardens.C,vardens.D)),log=TRUE),2,sum)
                })
                
                tau1    = log(.Object@prop)+logdens
                tau.max = apply(tau1,2,max)
                tau1    = exp(t(tau1)-tau.max)
                llk     = sum(log( apply(tau1,1,sum)) + tau.max)
            }  
            return(llk)
        })


setGeneric( name = "getetaFCMM" ,def = function(.Object,CCD,CCO){standardGeneric("getetaFCMM")})
setMethod(
    f          = "getetaFCMM",
    signature  = "CClustRes",
    definition =
        function(.Object,CCD,CCO) {
            x1  = 1+1e-4          
            x3  = 6
            r   = (sqrt(5)-1)/2
            x2  = x1+(1-r)*(x3-x1)
            x4  = x1+r*(x3-x1)
            Jx2 = -getllkFCMM(.Object,CCD,CCO,x2)
            Jx4 = -getllkFCMM(.Object,CCD,CCO,x4)
            while (abs(x4-x2)>10^-4){
                if (Jx2<Jx4){
                    x3        = x4
                    x4        = x2
                    Jx4       = Jx2
                    x2        = x1+(1-r)*(x3-x1)
                    Jx2       = -getllkFCMM(.Object,CCD,CCO,x2)
                } else {
                    if (Jx2>Jx4){
                        x1         = x2
                        x2         = x4
                        Jx2        = Jx4
                        x4         = x1+r*(x3-x1)
                        Jx4       = -getllkFCMM(.Object,CCD,CCO,x4)
                    } else {
                        x1        = x2
                        x3        = x4
                        x2        = x1+(1-r)*(x3-x1)
                        x4        = x1+r*(x3-x1)        
                        Jx2       = -getllkFCMM(.Object,CCD,CCO,x2)
                        Jx4       = -getllkFCMM(.Object,CCD,CCO,x4)
                    }
                    
                }
            }
            etamax=max(round((x2+x4)/2,4),1+1e-4)
            etamax
        })


setGeneric( name = "EstepFCMM"     ,def = function(.Object,CCO,CCD){standardGeneric("EstepFCMM")})

setMethod(
    f          = "EstepFCMM",
    signature  = "CClustRes",
    definition =
        function(.Object,CCO,CCD){
            
            
            Hcoef = 2^-(CCD@jklevels$scale*.Object@eta)
            
            if( (CCO@Gamma2.structure=="group.scale.location") | (CCO@Gamma2.structure=="group") ){
                
                lambda.Theta = lapply(1:CCO@nbclust,FUN = function(ell){.Object@varE/.Object@Gamma2.Theta[[ell]]})
                lambda.Nu    = lapply(1:CCO@nbclust,FUN = function(ell){.Object@varE/.Object@Gamma2.Nu[[ell]]})
                
                .Object@Blup  = sapply(1:CCO@nbclust, FUN=function(ell){
                    Theta = t(sapply(CCD@wavecoef,FUN = function(x){
                        (x$D-.Object@Beta[[ell]])/(1+lambda.Theta[[ell]]/Hcoef)
                    }))
                    Nu    = t(t(sapply(CCD@wavecoef,FUN = function(x){
                        (x$C-.Object@Alpha[[ell]])/(1+lambda.Nu[[ell]])
                    })))
                    return(list(Theta=Theta,Nu=Nu))
                },simplify=FALSE,USE.NAMES=TRUE)
                
                vardens.C           = sapply(1:CCO@nbclust,FUN=function(ell){.Object@Gamma2.Nu[[ell]]+.Object@varE})                
                vardens.D           = sapply(1:CCO@nbclust,FUN=function(ell){.Object@Gamma2.Theta[[ell]]*Hcoef+.Object@varE})
                .Object@varBlup.Nu    = lapply(1:CCO@nbclust,FUN=function(ell){.Object@varE/(1+lambda.Nu[[ell]])})              
                .Object@varBlup.Theta = lapply(1:CCO@nbclust,FUN=function(ell){.Object@varE/(1+lambda.Theta[[ell]]/Hcoef)})
                if (CCO@Gamma2.structure=="group"){.Object@varBlup.Theta = lapply(.Object@varBlup.Theta,FUN=function(y){mean(y)})}        	              
                
                
            }
            
            
            if ((CCO@Gamma2.structure=="scale.location") | (CCO@Gamma2.structure=="constant")) {
                
                lambda.Theta = .Object@varE/.Object@Gamma2.Theta
                lambda.Nu    = .Object@varE/.Object@Gamma2.Nu
                
                .Object@Blup  = sapply(1:CCO@nbclust, FUN=function(ell){
                    Theta = t(sapply(CCD@wavecoef,FUN = function(x){
                        (x$D-.Object@Beta[[ell]])/(1+lambda.Theta/Hcoef)
                    }))
                    Nu    = t(t(sapply(CCD@wavecoef,FUN = function(x){
                        (x$C-.Object@Alpha[[ell]])/(1+lambda.Nu)
                    })))
                    return(list(Theta=Theta,Nu=Nu))
                },simplify=FALSE,USE.NAMES=TRUE)
                
                vardens.C           = .Object@Gamma2.Nu+.Object@varE
                vardens.D           = .Object@Gamma2.Theta*Hcoef+.Object@varE      
                .Object@varBlup.Nu    = .Object@varE/(1+lambda.Nu)
                .Object@varBlup.Theta = .Object@varE/(1+lambda.Theta/Hcoef)
                if (CCO@Gamma2.structure=="constant"){
                    .Object@varBlup.Theta = mean(.Object@varBlup.Theta)
                }
                
            }
            
            logdens   = sapply(CCD@wavecoef, FUN = function(x){
                apply( dnorm( c(x$C,x$D),sapply(1:CCO@nbclust,FUN = function(ell){c(.Object@Alpha[[ell]],.Object@Beta[[ell]])}) , sd=sqrt(c(vardens.C,vardens.D)),log=TRUE),2,sum)
            })
            
            tau1           = log(.Object@prop)+logdens 
            tau.max        = apply(tau1,2,max)
            tau1           = exp(t(tau1)-tau.max)
            .Object@loglik = sum(log( apply(tau1,1,sum)) + tau.max)         
            .Object@Tau    = t(apply(tau1,1,FUN=function(x) x/sum(x)))
            
            return(.Object)
        }
)


setGeneric( name = "SstepFCMM"     ,def = function(.Object,CCO,CCD){standardGeneric("SstepFCMM")})

setMethod(
    f          = "SstepFCMM",
    signature  = "CClustRes",
    definition =
        function(.Object,CCO,CCD) {
            n              = length(CCD@wavecoef)            
            M              = length(CCD@wavecoef[[1]]$D) + 1
            clustfreq      = apply(.Object@Tau,2,sum)
            .Object@prop   = as.vector(clustfreq) / n 
            Wd             = t(sapply(CCD@wavecoef,  FUN=function(x){x$D},USE.NAMES=TRUE))
            Wc             = t(t(sapply(CCD@wavecoef,FUN=function(x){x$C},USE.NAMES=TRUE)))
            Hcoef          = 2^-(CCD@jklevels[,1]*.Object@eta)
            eps            = 1e-10

            set.seed(CCO@seed)
            E     = t(apply(.Object@Tau,1,FUN=function(x) rmultinom(1,1,x)))  
            Nell  = apply(E,2,sum)
            count = 0
            while( (min(Nell)<eps) & (count<100)){
                E     = t(apply(.Object@Tau,1,FUN=function(x) rmultinom(1,1,x)))
                Nell  = apply(E,2,sum)
                count = count + 1
            }

            if (count==100){
                count = 0
                while( (min(Nell)<eps) & (count<10)){
                    count = count + 1
                    E     = t(apply(matrix(1/CCO@nbclust,ncol=CCO@nbclust,nrow=n),1,FUN=function(x) rmultinom(1,1,x)))   #n*L
                    Nell  = apply(E,2,sum)
                }                
            }
            
            .Object@Tau = E            
            if( (CCO@Gamma2.structure=="group.scale.location") | (CCO@Gamma2.structure=="group") ){
                .Object@Blup = sapply(1:CCO@nbclust,
                    FUN=function(ell){
                        Theta = t(apply(.Object@Blup[[ell]]$Theta,1,FUN =
                                            function(x){rnorm(M-1,x,sqrt(Hcoef*.Object@varBlup.Theta[[ell]]))}))
                        Nu    = t(t(apply(.Object@Blup[[ell]]$Nu,1,FUN = function(x){rnorm(1,x,sqrt(.Object@varBlup.Nu[[ell]]))})))
                        return(list(Theta=Theta,Nu=Nu))
                    },simplify=FALSE,USE.NAMES=TRUE)
            }
            if ((CCO@Gamma2.structure=="scale.location") |
                (CCO@Gamma2.structure=="constant")) {
                .Object@Blup = sapply(1:CCO@nbclust, FUN=function(ell){
                    Theta = t(apply(.Object@Blup[[ell]]$Theta,1,FUN =
                                        function(x){rnorm(M-1,x,sqrt(Hcoef*.Object@varBlup.Theta))}))
                    Nu    = t(t(apply(.Object@Blup[[ell]]$Nu,1,FUN = function(x){rnorm(1,x,sqrt(.Object@varBlup.Nu))})))
                    return(list(Theta=Theta,Nu=Nu))
                },simplify=FALSE,USE.NAMES=TRUE)
            }            
            return(.Object)
        }
)




setGeneric( name = "gettheta"      ,def =
               function(.Object,CCO){standardGeneric("gettheta")})

setMethod(
    f          = "gettheta",
    signature  = "CClustRes",
    definition =
        function(.Object,CCO) {
            
            if ((CCO@nbclust>1) & (CCO@Gamma2.structure!="none")){              
                theta = c(.Object@prop,unlist(.Object@Beta),unlist(.Object@Alpha),.Object@varE,unlist(.Object@Gamma2.Theta),unlist(.Object@Gamma2.Nu),.Object@eta)              
            }
            if ((CCO@nbclust==1) & (CCO@Gamma2.structure!="none")){
                theta = c(unlist(.Object@Beta),unlist(.Object@Alpha),.Object@varE,unlist(.Object@Gamma2.Theta),unlist(.Object@Gamma2.Nu),.Object@eta)              
            }
            if ((CCO@nbclust>1) & (CCO@Gamma2.structure=="none")){
                theta = c(.Object@prop,unlist(.Object@Beta),unlist(.Object@Alpha),.Object@varE,.Object@eta)               
            }
            return(theta)
        }
)


setGeneric( name = "getAICBIC"        ,def =
               function(.Object,CCD){standardGeneric("getAICBIC")})

setMethod(
    f          = "getAICBIC",
    signature  = "CClustRes",
    definition =
        function(.Object,CCD) {

            n = length(CCD@wavecoef)            
            M = length(CCD@wavecoef[[1]]$D) + 1            

            if (is.null(.Object@Tau)){
                nbclust = 1
            } else {
                nbclust = dim(.Object@Tau)[2]
            }
            
            if (is.null(.Object@Gamma2.Theta)){
                dimL = (M+1)*nbclust
            } else {
                if (is.list(.Object@Gamma2.Theta)){
                    if (length(.Object@Gamma2.Theta[[1]])==1){
                        Gamma2.structure = "group"
                    } else {
                        Gamma2.structure = "group.scale.location"                
                    }
                } else {
                    if (length(.Object@Gamma2.Theta)==1){
                        Gamma2.structure = "constant"
                    } else {
                        Gamma2.structure = "scale.location"                
                    }
                }
                
                if( (Gamma2.structure=="group.scale.location") ){
                    dimL = ( (M+1)*nbclust + M*nbclust )                 
                } else if (Gamma2.structure=="group"){
                    dimL = ( (M+1)*nbclust + 2*nbclust )                        
                } else if (Gamma2.structure=="scale.location"){
                    dimL = ( (M+1)*nbclust + M )                        
                } else if (Gamma2.structure=="constant"){
                    dimL = ( (M+1)*nbclust + 2 )                        
                }
            }            
            
            AIC <-  2*dimL-2*.Object@loglik
            BIC <- -2*.Object@loglik +dimL*log(n)  
            return(list(BIC=BIC, AIC=AIC)) 
        })


setGeneric( name = "getwr.mu"      ,def = function(.Object,CCO,CCD){standardGeneric("getwr.mu")})

setMethod(
    f          = "getwr.mu",
    signature  = "CClustRes",
    definition =
        function(.Object,CCO,CCD) {
            
            M     = length(CCD@wavecoef[[1]]$D) + 1           
            Jmin  = floor(log2(CCD@lengthsignal))
            needs = ceiling(CCD@lengthsignal/2^Jmin)*2^Jmin-CCD@lengthsignal
            M2J   = CCD@lengthsignal+needs
            
            if (CCO@nbclust>1){
                
                mu = lapply(1:CCO@nbclust,FUN = function(ell){
                    Wmu.ell = wd(rep(0,M2J),filter.number=CCD@filter.number,family="DaubExPhase")                                
                    Wmu.ell = putC.wd(Wmu.ell,level=0,.Object@Alpha[[ell]])
                    for (j in unique((CCD@jklevels$scale))){
                        new.coef = rep(0,2^j)
                        pos4fun  = CCD@jklevels[CCD@jklevels$scale==j,2]+1
                        pos4coef = which(CCD@jklevels$scale==j)
                        new.coef[pos4fun] = (.Object@Beta[[ell]][pos4coef])
                        Wmu.ell  = putD.wd(Wmu.ell,level=j,new.coef)
                    }
                    return(wr(Wmu.ell,filter.number=CCD@filter.number,family="DaubExPhase")[1:CCD@lengthsignal])
                })
            } else {
                
                Wmu = wd(rep(0,M2J),filter.number=CCD@filter.number,family="DaubExPhase")
                Wmu = putC.wd(Wmu,level=0,.Object@Alpha)
                for (j in unique((CCD@jklevels$scale))){
                    new.coef = rep(0,2^j)
                    pos4fun  = CCD@jklevels[CCD@jklevels$scale==j,2]+1
                    pos4coef = which(CCD@jklevels$scale==j)
                    new.coef[pos4fun] = (.Object@Beta[pos4coef])
                    Wmu  = putD.wd(Wmu,level=j,new.coef)
                }
                mu = wr(Wmu,filter.number=CCD@filter.number,family="DaubExPhase")[1:CCD@lengthsignal]
            }                          
            return(mu)
        })
