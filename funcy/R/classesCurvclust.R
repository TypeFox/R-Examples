#
#Copyright (c) 2012 Madison Giacofci, Sophie Lambert-Lacroix, Guillemette Marot, Franck Picard
#
#

setClassUnion("charOrNULL", c("character","NULL"))
setClassUnion("listOrNULL", c("list","NULL"))
setClassUnion("listOrNumericOrNULL", c("list","NULL","numeric"))
setClassUnion("matrixOrNULL", c("matrix","NULL"))
setClassUnion("numericOrNULL", c("numeric","NULL"))

setClass(Class = "CClustO",representation(
             nbclust            = "numeric",
             Gamma2.structure   = "charOrNULL",
             burn               = "numeric",
             eps                = "numeric",
             init               = "character",
             loglikplot         = "logical",
             seed               = "numeric"
),
         prototype(nbclust = 2,
                   Gamma2.structure = "constant",
                   burn             = 100,
                   eps              = 0.1,
                   init             = "SEM",
                   loglikplot       = FALSE
                   ),         
         validity=function(object){
             stopifnot(object@burn>=0)
             stopifnot(object@eps>=0)
             stopifnot(object@Gamma2.structure %in% c("none","constant", "group", "scale.location", "group.scale.location"))
             stopifnot(object@init %in% c("rEM","SEM"))
             return(TRUE)         
         })


setMethod(f="[","CClustO",function(x,i,j,drop){
    switch(EXPR=i,
           "nbclust"         ={return(x@nbclust)},
           "Gamma2.structure"={return(x@Gamma2.structure)},
           "burn"            ={return(x@burn)},
           "eps"             ={return(x@eps)},
           "loglikplot"      ={return(x@loglikplot)},
           "init"            ={return(x@init)},
           "seed"            ={return(x@seed)},
           stop("This argument does not exist!"))})


setReplaceMethod(f="[","CClustO",function(x,i,j,value){
    switch(EXPR=i,
           "nbclust"         ={x@nbclust<-value},
           "Gamma2.structure"={x@Gamma2.structure<-value},
           "burn"            ={x@burn<-value},
           "eps"             ={x@eps<-value},
           "loglikplot"      ={x@loglikplot<-value},
           "init"            ={x@init<-value},
           "seed"            ={x@seed<-value},
           stop("This argument does not exist!"))
    validObject(x)
    return(x)})


setClass(Class = "CClustRes",
         representation(
             prop          = "numericOrNULL",
             Beta          = "listOrNumericOrNULL",
             Alpha         = "listOrNumericOrNULL",
             Blup          = "listOrNULL",
             varBlup.Nu    = "listOrNumericOrNULL",
             varBlup.Theta = "listOrNumericOrNULL",
             varE          = "numeric",  
             Gamma2.Theta  = "listOrNumericOrNULL",
             Gamma2.Nu     = "listOrNumericOrNULL",
             Tau           = "matrixOrNULL",
             eta           = "numericOrNULL",
             loglik        = "numeric"
         ))


setMethod(
    f          = "initialize",
    signature  = "CClustRes",
    definition =          
        function(.Object, CCO, CCD) {

            n              = length(CCD@wavecoef)            
            M              = length(CCD@wavecoef[[1]]$D) + 1
            .Object@varE   = (mean(CCD@seMAD))^2
            .Object@loglik = 0
            .Object@eta    = 2
            gamma2.init    = (mean(CCD@seMAD))^2
            
            if (CCO@nbclust > 1){
                
                .Object@Tau    = matrix(1/CCO@nbclust,ncol=CCO@nbclust,nrow=n)
                .Object@Beta   = lapply(1:CCO@nbclust,FUN = function(ell){rep(0,M-1)})
                .Object@Alpha  = lapply(1:CCO@nbclust,FUN = function(ell){0})
                .Object@prop   = rep(1/CCO@nbclust,CCO@nbclust)
                
                if (CCO@Gamma2.structure=="group.scale.location"){
                    .Object@Blup          = lapply(1:CCO@nbclust,FUN = function(ell){list(Theta = matrix(0,ncol=M-1,nrow=n), Nu=matrix(0,ncol=1,nrow=n))})
                    .Object@Gamma2.Theta  = lapply(1:CCO@nbclust,FUN = function(ell){rep(gamma2.init,M-1)})
                    .Object@Gamma2.Nu     = lapply(1:CCO@nbclust,FUN = function(ell){gamma2.init})
                    .Object@varBlup.Theta = lapply(1:CCO@nbclust,FUN = function(ell){rep(gamma2.init,M-1)})
                    .Object@varBlup.Nu    = lapply(1:CCO@nbclust,FUN = function(ell){gamma2.init})
                    .Object@eta           = 1+1e-4
                }
                
                if (CCO@Gamma2.structure=="group"){
                    .Object@Blup          = lapply(1:CCO@nbclust,FUN = function(ell){list(Theta = matrix(0,ncol=M-1,nrow=n), Nu=matrix(0,ncol=1,nrow=n))})
                    .Object@Gamma2.Theta  = lapply(1:CCO@nbclust,FUN = function(ell){rep(gamma2.init,1)})
                    .Object@Gamma2.Nu     = lapply(1:CCO@nbclust,FUN = function(ell){rep(gamma2.init,1)})
                    .Object@varBlup.Theta = lapply(1:CCO@nbclust,FUN = function(ell){rep(gamma2.init,1)})
                    .Object@varBlup.Nu    = lapply(1:CCO@nbclust,FUN = function(ell){rep(gamma2.init,1)})
                    .Object@eta           = 1+1e-4
                }
                
                if (CCO@Gamma2.structure=="scale.location"){
                    .Object@Blup          = lapply(1:CCO@nbclust,FUN = function(ell){list(Theta = matrix(0,ncol=M-1,nrow=n), Nu=matrix(0,ncol=1,nrow=n))})
                    .Object@Gamma2.Theta  = rep(gamma2.init,M-1)
                    .Object@Gamma2.Nu     = gamma2.init
                    .Object@varBlup.Theta = rep(gamma2.init,M-1)
                    .Object@varBlup.Nu    = gamma2.init
                    .Object@eta           = 1+1e-4
                }            
                
                if (CCO@Gamma2.structure=="constant"){
                    .Object@Blup          = lapply(1:CCO@nbclust,FUN = function(ell){list(Theta = matrix(0,ncol=M-1,nrow=n), Nu=matrix(0,ncol=1,nrow=n))})
                    .Object@Gamma2.Theta  = gamma2.init
                    .Object@Gamma2.Nu     = gamma2.init
                    .Object@varBlup.Theta = gamma2.init
                    .Object@varBlup.Nu    = gamma2.init
                    .Object@eta           = 1+1e-4
                }
                if (CCO@Gamma2.structure=="none"){
                    .Object@Blup          = NULL
                    .Object@Gamma2.Theta  = NULL
                    .Object@Gamma2.Nu     = NULL
                    .Object@varBlup.Theta = NULL
                    .Object@varBlup.Nu    = NULL
                    .Object@eta           = NULL
                }
                
            } else {
                
                .Object@Tau    = NULL
                .Object@Beta   = rep(0,M-1)
                .Object@Alpha  = 0
                .Object@prop   = NULL
                
                if ( (CCO@Gamma2.structure=="scale.location") ){
                    .Object@Blup          = list(Theta = matrix(0,ncol=M-1,nrow=n), Nu=matrix(0,ncol=1,nrow=n))
                    .Object@Gamma2.Theta  = rep(gamma2.init,M-1)
                    .Object@Gamma2.Nu     = gamma2.init
                    .Object@varBlup.Theta = rep(gamma2.init,M-1)
                    .Object@varBlup.Nu    = gamma2.init
                    .Object@eta           = 1+1e-4
                }
                
                if ( (CCO@Gamma2.structure=="constant") ){
                    .Object@Blup          = list(Theta = matrix(0,ncol=M-1,nrow=n), Nu=matrix(0,ncol=1,nrow=n))
                    .Object@Gamma2.Theta  = gamma2.init
                    .Object@Gamma2.Nu     = gamma2.init
                    .Object@varBlup.Theta = gamma2.init
                    .Object@varBlup.Nu    = gamma2.init
                    .Object@eta           = 1+1e-4
                }
                if (CCO@Gamma2.structure=="none"){
                    .Object@Blup          = NULL
                    .Object@Gamma2.Theta  = NULL
                    .Object@Gamma2.Nu     = NULL
                    .Object@varBlup.Theta = NULL
                    .Object@varBlup.Nu    = NULL
                    .Object@eta           = NULL
                }
            }
            return(.Object)
        }
)


setMethod(f="[","CClustRes",function(x,i,j,drop){
    switch(EXPR=i,
           "prop"         ={return(x@prop)},
           "Beta"         ={return(x@Beta)},
           "Alpha"        ={return(x@Alpha)},
           "Blup"         ={return(x@Blup)},
           "varBlup.Nu"   ={return(x@varBlup.Nu)},
           "varBlup.Theta"={return(x@varBlup.Theta)},
           "varE"         ={return(x@varE)},
           "Gamma2.Theta" ={return(x@Gamma2.Theta)},
           "Gamma2.Nu"    ={return(x@Gamma2.Nu)},
           "Tau"          ={return(x@Tau)},
           "eta"          ={return(x@eta)},
           "loglik"       ={return(x@loglik)},
           stop("This argument does not exist!"))})


setReplaceMethod(f="[","CClustRes",function(x,i,j,value){
    switch(EXPR=i,
           "prop"         ={x@prop<-value},
           "Beta"         ={x@Beta<-value},
           "Alpha"        ={x@Alpha<-value},
           "Blup"         ={x@Blup<-value},
           "varBlup.Nu"   ={x@varBlup.Nu<-value},
           "varBlup.Theta"={x@varBlup.Theta<-value},
           "varE"         ={x@varE<-value},
           "Gamma2.Theta" ={x@Gamma2.Theta<-value},
           "Gamma2.Nu"    ={x@Gamma2.Nu<-value},
           "Tau"          ={x@Tau<-value},
           "eta"          ={x@eta<-value},
           "loglik"       ={x@loglik<-value},
           stop("This argument does not exist!"))
    validObject(x)
    return(x)}) 


setMethod(
    f = "summary",
    signature = "CClustRes",
    definition = function(object){			
        if (is.null(object@Tau)){
            nbclust = 1
        } else {
            nbclust = dim(object@Tau)[2]
        }
        cat("Number of clusters ",nbclust,"\n")
        if (is.null(object@Tau)){
            csize = 1
        } else {
            csize = object@prop
            cat("cluster size ",csize,"\n")
            cat("labels ","\n")
            print(apply(object@Tau,1,which.max))
        }

    }          
)


setClass(Class = "CClustData", representation  (
             wavecoef      = "list",
             jklevels      = "data.frame",
             seMAD         = "numeric",
             lengthsignal  = "numeric",
             filter.number = "numeric"
)
         )

setMethod(
    f          = "initialize",
    signature  = "CClustData",
    definition =          
        function(.Object,Y,filter.number) {

            if(missing(filter.number)){filter.number=8}

            n  = length(Y)
            if(is.null(names(Y))){
                names(Y) = paste("ID",1:n,sep="")
            }

            mi = sapply(Y,FUN = function(y){length(y)})
            if (length(unique(mi))>1){
                stop("All curves should have the same number of recorded positions")
            }
            if(any(is.na(unlist(Y)))){stop("please remove missing values")}

            Minit = length(Y[[1]])
            .Object@lengthsignal = length(Y[[1]])
            
            Jmin  = floor(log2(Minit))
            needs = ceiling(Minit/2^Jmin)*2^Jmin-Minit

            if (needs>0){
                Y = sapply(Y,FUN = function(y){y = c(y,rep(0,needs)); y},simplify=FALSE,USE.NAMES=TRUE)
            }
            
            M      = length(Y[[1]])
            J      = log2(M)-1
            
            wavdec = sapply(Y,FUN=function(y){wd(y,filter.number=filter.number,family="DaubExPhase")},simplify=F,USE.NAMES=TRUE)
            D    = sapply(wavdec, FUN=function(x){x$D},simplify=F,USE.NAMES=TRUE)
            C    = unlist(sapply(wavdec, FUN=function(x){accessC(x,level=0)},simplify=F,USE.NAMES=TRUE))
            
            jcoef              = rep(J:0,2^(J:0))
            kcoef              = unlist(lapply(J:0,FUN=function(x) 0:(2^x-1)))
            jklevels           = data.frame(scale=jcoef,location=kcoef)

            seMAD = sapply(Y, FUN = function(y){
                u         = floor(log2(length(y)))
                ywd       = wd(y,filter.number=filter.number,family="DaubExPhase")
                finecoef  = accessD(ywd,lev=(u-1))
                sigma     = mad(finecoef) 
            })

            .Object@seMAD = seMAD
            
            if(needs>0){
                pos2keep  = apply(Reduce("cbind",D),1,FUN=function(x){sum(x==0)})!=n
                jklevels  = jklevels[pos2keep,]
            } else {
                pos2keep = rep(TRUE,M-1)
            }

            .Object@wavecoef = sapply(names(D),FUN = function(ID){
                list(C=C[[ID]],D=D[[ID]][pos2keep])
            },simplify=F,USE.NAMES=TRUE)

            .Object@jklevels      = jklevels            
            .Object@filter.number = filter.number
            
            return(.Object)
        }
)

setGeneric(name = "getUnionCoef", def=function(.Object){standardGeneric("getUnionCoef")})
setMethod(
    f          = "getUnionCoef",
    signature  = "CClustData",
    definition =
        function(.Object) {

            M  = dim(.Object@jklevels)[1]+1
            n  = length(names(.Object@wavecoef))
            se = mean(.Object@seMAD)            
            zeros = sapply(.Object@wavecoef,FUN = function(x){
                coef = pmax(unlist(x)-se*sqrt(2*log(M-1)),0)
                return(coef)
            })
            union    = (apply(zeros,1,sum)!=0)
            union[1] = TRUE
            newcoef = lapply(.Object@wavecoef,FUN = function(x){
                z = c(x$C,x$D)[union]
                C = z[1]
                D = z[-1]
                list(C=C,D=D)
                list(C=C,D=D)
            })            
            .Object@wavecoef = newcoef
            .Object@jklevels = .Object@jklevels[union[-1],]
            return(.Object)
        }          
)



