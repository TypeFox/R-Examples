exact2x2CI<-function(x,tsmethod="minlike",  conf.level=.95,tol=.00001,orRange=c(10^-10,10^10)){
    alpha<-1-conf.level
    # can make tol for uniroot functions less than tol if you want, edit urtol
    # but it is not helpful for the current algorithm
    urtol<-tol     
    m<-sum(x[1,])
    n<-sum(x[2,])
    k<-sum(x[,1])
    x<-x[1,1]
    lo<-max(0,k-n)
    hi<-min(m,k)
    support<- lo:hi
    ns<-length(support)
    logf0 <- dhyper(support, m, n, k, log = TRUE)
    ## use dnhyper as in fisher.test
    dnhyper <- function(OR) {
        if (OR==0){
            d<-c(1,rep(0,ns-1))
        } else if (OR==Inf){
            d<-c(rep(0,ns-1),1)
        } else {
            d <- logf0 + log(OR) * support
            d <- exp(d - max(d))
            d <-d/sum(d)
        }
        d
    }
    pnhyper<- function(x,OR,lower.tail=TRUE){
        nOR<-length(OR)
        out<-rep(NA,nOR)
        if (lower.tail){ X<- support<=x 
        } else X<- support>=x

        for (i in 1:nOR){  
            f<-dnhyper(OR[i])
            out[i]<-sum(f[X])
        }
        out
    }
    intercept<-function(xlo,xhi,ORRange=orRange,TSmethod=tsmethod){
        ## old intercept function did not always work correctly with 
        ## minlike method. The problem was when the density function 
        ## was very very small for one of the ends of the orRange
        ## then it would give f[Xlo]==0 and f[Xhi]==0 and the uniroot 
        ## would stop at the limit of the range.

        ## After thinking about it, we do not even need to use uniroot
        ## we just need to solve the following for beta:
        ## logf0[Xlo] + xlo*log(beta) = logf0[Xhi]+xhi*log(beta)
        ## 
        if (TSmethod=="minlike"){
            Xlo<- support==xlo
            Xhi<- support==xhi
            beta<- exp( (logf0[Xlo] - logf0[Xhi])/(xhi-xlo) )

        } else if (TSmethod=="blaker"){
            Xlo<- support<=xlo
            Xhi<- support>=xhi
            rootfunc<-function(beta){
                nb<-length(beta)
                out<-rep(NA,nb)
                for (i in 1:nb){
                    f<-dnhyper(beta[i])
                    out[i]<- sum(f[Xlo]) - sum(f[Xhi]) 
                }
                out
            }
            beta<-uniroot(rootfunc,ORRange,tol=urtol )$root
        }
        beta
    }

    Bnds<-function(xlo,xhi,ORRange,ndiv=1){
        orlo<-min(ORRange)
        orhi<-max(ORRange)
        OR<- orlo + (orhi-orlo)*((0:ndiv)/ndiv)
        F<- pnhyper(xlo,OR,lower.tail=TRUE)
        Fbar<-pnhyper(xhi,OR,lower.tail=FALSE)
        estimate<- F+Fbar
        L<- F[-1] + Fbar[-(ndiv+1)]
        U<- F[-(ndiv+1)] + Fbar[-1]
        list(or=OR,estimate=estimate,bndlo=L,bndhi=U)
    }
    refine<-function(xlo,xhi,ORRange,NDIV=100,maxiter=50,limit="upper"){
        getCLbnds<-function(b){
            nb<-length(b$bndhi)
            HI<-max(b$bndhi)
            LO<-min(b$bndlo)
            if (HI<=alpha){
                CLbnds<-NULL
                continue<-TRUE
            } else if (LO>alpha){
                if (limit=="upper"){
                    CLbnds<-c(max(b$or)-tol/2,max(b$or)+tol/2)
                } else {
                    CLbnds<-c(min(b$or)-tol/2,min(b$or)+tol/2)
                }
                continue<-FALSE
            } else {
                CLbnds<-c(NA,NA)
                if (limit=="upper"){
                    if (any(b$bndlo>alpha)){
                        CLbnds[1]<-max( b$or[2:(nb+1)][b$bndlo>alpha] )
                    } else { 
                        CLbnds[1]<- min(b$or)
                    }
                    CLbnds[2]<- max( b$or[2:(nb+1)][b$bndhi>alpha])
                } else {
                # limit=lower
                    if (any(b$bndlo>alpha)){
                        CLbnds[2]<-min( b$or[1:nb][b$bndlo>alpha] )
                    } else {
                        CLbnds[2]<-max(b$or)
                    }
                    CLbnds[1]<-min( b$or[1:nb][b$bndhi>alpha] )
                }
                continue<-TRUE
            }
            out<-list(CLbnds=CLbnds,continue=continue)
            out
        }
        b<-Bnds(xlo,xhi,ORRange,ndiv=1)
        clb<-getCLbnds(b)     
        if (!is.null(clb$CLbnds) & clb$continue){
            ORRANGE<-clb$CLbnds
            for (i in 1:maxiter){
                b<-Bnds(xlo,xhi,ORRANGE,ndiv=NDIV) 
                clb<-getCLbnds(b)
                if (!clb$continue | (clb$continue & is.null(clb$CLbnds))) break()
                ORRANGE<-clb$CLbnds
                if (ORRANGE[2]-ORRANGE[1]>tol){ 
                    NDIV<-2*NDIV
                    if (i==maxiter){
                        warning("Could not estimate confidence interval to within tol level, see conf.limit.prec attr of conf.int")
                    }
                } else if (ORRANGE[2]-ORRANGE[1]<=tol){
                    clb$continue<-FALSE
                    break()
                }
            }
        }
        clb
    }
    CINT<-c(NA,NA)
    if (x==hi){
        CINT[2]<-Inf
        upper.prec<-c(Inf,Inf)
    } 
    if (x==lo){
        CINT[1]<-0
        lower.prec<-c(0,0)
    }
    xless<-lo:(x-1)
    if (is.na(CINT[2])){
        xgreater<-hi:(x+1)
        ngreater<- length(xgreater)
        ints<-rep(NA,ngreater)
        #bndlo<-bndhi<-pend1<-pend2<-orend1<-orend2<-rep(NA,ngreater)
        for (i in 1:ngreater){
            ints[i]<-intercept(x,xgreater[i])
            F<-pnhyper(x,ints[i],lower.tail=TRUE)
            if (i==1){
                if (F>alpha){
                    rootfunc<-function(or){
                        nor<-length(or)
                        out<-rep(NA,nor)
                        for (i in 1:nor){
                            out[i]<- alpha - pnhyper(x,or[i],lower.tail=TRUE)
                        }
                        out
                    }
                    if (rootfunc(orRange[2])<0) stop("very large odds ratio, modify orRange in exact2x2CI")
                    CINT[2]<-uniroot(rootfunc,c(ints[i],orRange[2]),tol=urtol)$root
                    upper.prec<-c(CINT[2]-urtol/2,CINT[2]+urtol/2)
                    break()
                } else if (F==alpha){
                    CINT[2]<-ints[i]
                    upper.prec<-c(CINT[2]-urtol/2,CINT[2]+urtol/2)
                    break()
                }
            } else if (i>1){
                rout<- refine(x,xgreater[i-1],c(ints[i],ints[i-1]),limit="upper") 
              if (!rout$continue){ 
                    CINT[2]<-rout$CLbnds[2]
                    upper.prec<-rout$CLbnds
                    break()
                } else if (i==ngreater){
                    CINT[2]<-ints[i]
                    upper.prec<-c(ints[i]-urtol/2,ints[i]+urtol/2)
                }
            }
        }
    }
    if (is.na(CINT[1])){
        xless<-lo:(x-1)
        nless<- length(xless)
        ints<-rep(NA,nless)
        for (i in 1:nless){
            ints[i]<-intercept(xless[i],x)
            Fbar<-pnhyper(x,ints[i],lower.tail=FALSE)
            if (i==1){
                if (Fbar>alpha){
                    rootfunc<-function(or){
                        nor<-length(or)
                        out<-rep(NA,nor)
                        for (i in 1:nor){
                            out[i]<- alpha - pnhyper(x,or[i],lower.tail=FALSE)
                        }
                        out
                    }
                    if (rootfunc(orRange[1])<0) stop("very small odds ratio, modify orRange in exact2x2CI")
                    CINT[1]<-uniroot(rootfunc,c(orRange[1],ints[i]),tol=urtol)$root
                    lower.prec<-c(CINT[1]-urtol/2,CINT[1]+urtol/2)
                    break()
                } else if (F==alpha){
                    CINT[1]<-ints[i]
                    lower.prec<-c(CINT[1]-urtol/2,CINT[1]+urtol/2)
                    break()
                }
            } else if (i>1){
                rout<- refine(xless[i-1],x,c(ints[i-1],ints[i]),limit="lower") 
                if (!rout$continue){ 
                    CINT[1]<-rout$CLbnds[1]
                    lower.prec<-rout$CLbnds
                    break()
                } else if (i==nless){
                    CINT[1]<-ints[i]
                    lower.prec<-c(ints[i]-urtol/2,ints[i]+urtol/2)
                }
            }
        }
    }
    attr(CINT,"conf.level")<-conf.level
    attr(CINT,"conf.limit.prec")<-list(estimate=CINT,lower=lower.prec,upper=upper.prec)
    # round to the digit above tol level
    CINT<-round(CINT,floor(-log10(tol))-1)
    CINT
}

## example that shows problem with version 1.0-1.1
##x<-matrix(c(17,126,64,769),2,2)
##exact2x2CI(x,tsmethod="minlike")
##exact2x2CI(x,tsmethod="minlike",orRange=c(10^-3,10^3))