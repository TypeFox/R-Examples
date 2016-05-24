exactbinomCI<-function(x, n, tsmethod="minlike",conf.level=.95,tol=.00001,
    pRange=c(1e-10,1-1e-10)){

    ## NOTE: no midp option yet...needs to be coded

    ## code is a modification of exact2x2CI

    alpha<-1-conf.level
    # can make tol for uniroot functions less than tol if you want, edit urtol
    # but it is not helpful for the current algorithm
    urtol<-tol   
    lo<-0
    hi<-n  
    support<- 0:n    
    ns<-length(support)
    logf0 <- dbinom(support, n, .5, log=TRUE) - n*log(.5)
    ## use dnhyper as in fisher.test
    d<-function(p){
        if (p==0){
            out<-c(1,rep(0,n))
       } else if (p==1){
            out<-c(rep(0,n),1)
        } else if (p<1 & p>0){
           out<-exp(logf0 + support*log(p) + (n-support)*log(1-p))
        }
        out
    }

    ## the intercept function finds the places where there is a jump 
    ## in the evidence (i.e., p-value) function

    intercept<-function(xlo,xhi,TSmethod=tsmethod){
        if (TSmethod=="minlike"){
            ## root is the parameter  
            ## where dbinom(xlo,n,root)=dbinom(xhi,n,root)
            ## we can solve this algebraically
            root<-1/(1+exp((lchoose(n,xhi) - 
                            lchoose(n,xlo))/(xhi-xlo)))
        } else if (TSmethod=="blaker"){
            ## root is the parameter where the tails are equal, i.e.,
            ## where pbinom(xlo,n,root)=pbinom(xhi-1,n,root,lower.tail=FALSE)
            ## we solve using uniroot
            rootfunc<-function(beta){
                nb<-length(beta)
                out<-rep(NA,nb)
                for (i in 1:nb){
                    out[i]<- pbinom(xlo,n,beta[i]) - pbinom(xhi-1,n,beta[i],lower.tail=FALSE)
                }
                out   
            }
            root<-uniroot(rootfunc,pRange,tol=urtol )$root
        } else {
            stop("method must equal 'minlike' or 'blaker' ")
        }
        root
    }

    # get upper and lower bounds for p-value within range=pRange
    # divide it into ndiv equal pieces and get the bounds within 
    # each piece
 
    Bnds<-function(xlo,xhi,pRange,ndiv=1){
        plo<-min(pRange)
        phi<-max(pRange)
        P<- plo + (phi-plo)*((0:ndiv)/ndiv)
        F<- pbinom(xlo,n,P,lower.tail=TRUE)
        Fbar<-pbinom(xhi-1,n,P,lower.tail=FALSE)
        estimate<- F+Fbar
        L<- F[-1] + Fbar[-(ndiv+1)]
        U<- F[-(ndiv+1)] + Fbar[-1]
        list(p=P,estimate=estimate,bndlo=L,bndhi=U)
    }

   refine<-function(xlo,xhi,PRange,NDIV=100,maxiter=50,limit="upper"){
        getCLbnds<-function(b){
            nb<-length(b$bndhi)
            HI<-max(b$bndhi)
            LO<-min(b$bndlo)
            if (HI<=alpha){
                CLbnds<-NULL
                continue<-TRUE
            } else if (LO>alpha){
                if (limit=="upper"){
                    CLbnds<-c(max(b$p)-tol/2,max(b$p)+tol/2)
                } else {
                    CLbnds<-c(min(b$p)-tol/2,min(b$p)+tol/2)
                }
                continue<-FALSE
            } else {
                CLbnds<-c(NA,NA)
                if (limit=="upper"){
                    if (any(b$bndlo>alpha)){
                        CLbnds[1]<-max( b$p[2:(nb+1)][b$bndlo>alpha] )
                    } else { 
                        CLbnds[1]<- min(b$p)
                    }
                    CLbnds[2]<- max( b$p[2:(nb+1)][b$bndhi>alpha])
                } else {
                # limit=lower
                    if (any(b$bndlo>alpha)){
                        CLbnds[2]<-min( b$p[1:nb][b$bndlo>alpha] )
                    } else {
                        CLbnds[2]<-max(b$p)
                    }
                    CLbnds[1]<-min( b$p[1:nb][b$bndhi>alpha] )
                }
                continue<-TRUE
            }
            out<-list(CLbnds=CLbnds,continue=continue)
            out
        } # end of getCLbnds
        b<-Bnds(xlo,xhi,PRange,ndiv=1)
        clb<-getCLbnds(b)     
        if (!is.null(clb$CLbnds) & clb$continue){
            PRANGE<-clb$CLbnds
            for (i in 1:maxiter){
                b<-Bnds(xlo,xhi,PRANGE,ndiv=NDIV) 
                clb<-getCLbnds(b)
                if (!clb$continue | (clb$continue & is.null(clb$CLbnds))) break()
                PRANGE<-clb$CLbnds
                if (PRANGE[2]-PRANGE[1]>tol){ 
                    NDIV<-2*NDIV
                    if (i==maxiter){
                        warning("Could not estimate confidence interval to within tol level, see conf.limit.prec attr of conf.int")
                    }
                } else if (PRANGE[2]-PRANGE[1]<=tol){
                    clb$continue<-FALSE
                    break()
                }
            }
        }
        clb
    } # end refine
    CINT<-c(NA,NA)
    if (x==n){
        CINT[2]<-1
        upper.prec<-c(1,1)
    } 
    if (x==0){
        CINT[1]<-0
        lower.prec<-c(0,0)
    }
    if (is.na(CINT[2])){
        xgreater<-n:(x+1)
        ngreater<- length(xgreater)
        ints<-rep(NA,ngreater)
        #bndlo<-bndhi<-pend1<-pend2<-pend1<-pend2<-rep(NA,ngreater)
        for (i in 1:ngreater){
            ints[i]<-intercept(x,xgreater[i])
            F<-pbinom(x,n,ints[i],lower.tail=TRUE)
            if (i==1){
                if (F>alpha){
                    rootfunc<-function(p){
                        np<-length(p)
                        out<-rep(NA,np)
                        for (i in 1:np){
                            #out[i]<- alpha - pbinom(x,n,p[i],lower.tail=TRUE)
                            out[i]<-alpha - exactbinomPvals(x,n,p[i],tsmethod=tsmethod)$pvals
                        }
                        out
                    }
                    if (rootfunc(pRange[2])<0) stop("very large odds ratio, modify pRange in exact2x2CI")
                    CINT[2]<-uniroot(rootfunc,c(ints[i],pRange[2]),tol=urtol)$root
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
        xless<-0:(x-1)
        nless<- length(xless)
        ints<-rep(NA,nless)
        for (i in 1:nless){
            ints[i]<-intercept(xless[i],x)
            ### fixed below mistake may 6, 2010
            #Fbar<-pbinom(x,n,ints[i],lower.tail=FALSE)
            Fbar<-pbinom(x-1,n,ints[i],lower.tail=FALSE)
            if (i==1){
                if (Fbar>alpha){
                    rootfunc<-function(p){
                        np<-length(p)
                        out<-rep(NA,np)
                        for (i in 1:np){
                            #out[i]<- alpha - pbinom(x,n,p[i],lower.tail=FALSE)
                            out[i]<- alpha - exactbinomPvals(x,n,p[i],tsmethod=tsmethod)$pvals
                         }
                        out
                    }
                    if (rootfunc(pRange[1])<0) stop("very small odds ratio, modify pRange in exact2x2CI")
                    CINT[1]<-uniroot(rootfunc,c(pRange[1],ints[i]),tol=urtol)$root
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
    attr(CINT,"conf.limit.prec")<-list(lower=lower.prec,upper=upper.prec)
    # round to the digit above tol level
    CINT<-round(CINT,floor(-log10(tol))-1)
    CINT
}

#p<-(1:1000)/1001
#pval<-exactbinomPvals(3,10,p,tsmethod="blaker")
#plot(p,pval$pvals,xlim=c(.61,.64),ylim=c(.04,.06))
#lines(c(0,1),c(.05,.05))
#lines(c(.0873,.0873),c(0,1))
#lines(c(.6066,.6066),c(0,1))
#exactbinomCI(3,10,tsmethod="blaker")
#exactbinomCI(18,99,tsmethod="minlike")
