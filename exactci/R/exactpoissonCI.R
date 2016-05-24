exactpoissonCI<-function(x, tsmethod="minlike",conf.level=.95,tol=.00001,
    pRange=c(1e-10,1-1e-10)){

    ## code is a modification of exact2x2CI

    alpha<-1-conf.level
    # can make tol for uniroot functions less than tol if you want, edit urtol
    # but it is not helpful for the current algorithm
    urtol<-tol   

       ## since there is no upper limit for Poisson we need to find n=xmax
        ## such that at some large value of the parameter, say 
        ## the 1-alpha/100 upper confidence limit,  
        ## the tail is much less than alpha, say alpha/100
        ## this is arbitrary, there is probably room for improvement
        n<- qpois(1-alpha/100,qgamma(1-alpha/100,x))   
 

    ## the intercept function finds the places where there is a jump 
    ## in the evidence (i.e., p-value) function
    intercept<-function(xlo,xhi,TSmethod=tsmethod){
        if (TSmethod=="minlike"){
            ## root is the parameter  
            ## where dpois(xlo,root)=dpois(xhi,root)
            ## we can solve this algebraically
            root<-exp( (lfactorial(xhi) - lfactorial(xlo))/(xhi-xlo) )
        } else if (TSmethod=="blaker"){
            ## root is the parameter where the tails are equal, i.e.,
            ## where ppois(xlo,root)=ppois(xhi-1,root,lower.tail=FALSE)
            ## we solve using uniroot
            rootfunc<-function(beta){
                nb<-length(beta)
                out<-rep(NA,nb)
                for (i in 1:nb){
                    out[i]<- ppois(xlo,beta[i]) - ppois(xhi-1,beta[i],lower.tail=FALSE)
                }
                out   
            }
            ## note that for intergers,i, ppois(i-1,x,lower.tail=FALSE)=pgamma(x,a)
            ## so we use qgamma to get range from pRange
            root<-uniroot(rootfunc,c(qgamma(pRange[1],xlo),qgamma(pRange[2],xhi)),tol=urtol )$root

        } else {
            stop("method must equal 'minlike' or 'blaker' ")
        }
        root
    }

    # get upper and lower bounds for p-value within range=pRange
    # divide it into ndiv equal pieces and get the bounds within 
    # each piece
    rRange<-c(qgamma(min(pRange),0),qgamma(max(pRange),n))

    Bnds<-function(xlo,xhi,rRange,ndiv=1){
        plo<-min(rRange)
        phi<-max(rRange)
        P<- plo + (phi-plo)*((0:ndiv)/ndiv)
        F<- ppois(xlo,P,lower.tail=TRUE)
        Fbar<-ppois(xhi-1,P,lower.tail=FALSE)
        estimate<- F+Fbar
        L<- F[-1] + Fbar[-(ndiv+1)]
        U<- F[-(ndiv+1)] + Fbar[-1]
        list(p=P,estimate=estimate,bndlo=L,bndhi=U)
    }

   refine<-function(xlo,xhi,RRange,NDIV=100,maxiter=50,limit="upper"){
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
        b<-Bnds(xlo,xhi,RRange,ndiv=1)
        clb<-getCLbnds(b)     
        if (!is.null(clb$CLbnds) & clb$continue){
            RRANGE<-clb$CLbnds
            for (i in 1:maxiter){
                b<-Bnds(xlo,xhi,RRANGE,ndiv=NDIV) 
                clb<-getCLbnds(b)
                if (!clb$continue | (clb$continue & is.null(clb$CLbnds))) break()
                RRANGE<-clb$CLbnds
                if (RRANGE[2]-RRANGE[1]>tol){ 
                    NDIV<-2*NDIV
                    if (i==maxiter){
                        warning("Could not estimate confidence interval to within tol level, see conf.limit.prec attr of conf.int")
                    }
                } else if (RRANGE[2]-RRANGE[1]<=tol){
                    clb$continue<-FALSE
                    break()
                }
            }
        }
        clb
    } # end refine
    CINT<-c(NA,NA)
    if (x==0){
        CINT[1]<-0
        lower.prec<-c(0,0)
    }
    if (is.na(CINT[2])){
        ## since there is no upper limit for Poisson we need to find n=xmax
        ## such that at some large value of the parameter, say 
        ## the 1-alpha/100 upper confidence limit,  
        ## the tail is much less than alpha, say alpha/100
        ## this is arbitrary, there is probably room for improvement
        ## Sept 21,2012: fixed error, gives n=0 if x=0
        if (x==0){ gammaParm<- 1
        } else gammaParm<-x
        n<- qpois(1-alpha/100,qgamma(1-alpha/100,gammaParm))   
        xgreater<-n:(x+1)
        ngreater<- length(xgreater)
        ints<-rep(NA,ngreater)
        #bndlo<-bndhi<-pend1<-pend2<-pend1<-pend2<-rep(NA,ngreater)
        for (i in 1:ngreater){
            ints[i]<-intercept(x,xgreater[i])
            F<-ppois(x,ints[i],lower.tail=TRUE)
            if (i==1){
                if (F>=alpha){
                    stop("original F greater than alpha, rewrite n<- code") 
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
            Fbar<-ppois(x-1,ints[i],lower.tail=FALSE)
            if (i==1){
                if (Fbar>alpha){
                    rootfunc<-function(p){
                        alpha - exactpoissonPvals(x,p,tsmethod=tsmethod)$pvals
                     }

                    if (rootfunc(rRange[1])<0) stop("very small rate, modify pRange")
                    CINT[1]<-uniroot(rootfunc,c(rRange[1],ints[i]),tol=urtol)$root
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

#exactpoissonCI(5,tsmethod="minlike")
#exactpoissonCI(0,tsmethod="blaker")
