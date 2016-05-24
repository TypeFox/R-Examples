`icfit` <-
function(L,...){
    UseMethod("icfit")
}

`icfit.formula` <-
function (formula, data,...) 
{
    ## Most of this function is copied or slightly modified from survfit
    ## Copied starting from here:
    call <- match.call()
    if ((mode(call[[2]]) == "call" && call[[2]][[1]] == as.name("Surv")) || 
        inherits(formula, "Surv")) {
        formula <- eval(parse(text = paste(deparse(call[[2]]), 
            1, sep = "~")))
        environment(formula) <- parent.frame()
    }
    ## change survfit code to UseMethod "icfit" 
    if (!inherits(formula, "formula")) 
        temp <- UseMethod("icfit")
    else {
        m <- match.call(expand.dots = FALSE)
        m$... <- NULL
        Terms <- terms(formula, "strata")
        ord <- attr(Terms, "order")
        if (length(ord) & any(ord != 1)) 
            stop("Interaction terms are not valid for this function")
        m$formula <- Terms
        m[[1]] <- as.name("model.frame")
        m <- eval(m, parent.frame())
        n <- nrow(m)
        ## left-hand-side of formula is "response"
        Y <- model.extract(m, "response")
        ## change survfit code next few lines
        ## allow response to be numeric vector, treated as known time to event
        if (!is.Surv(Y)){
            if (is.numeric(Y) & is.vector(Y)) Y<-Surv(Y,rep(1,length(Y))) 
            else stop("Response must be a survival object or numeric vector")
        }
        casewt <- model.extract(m, "weights")
        if (is.null(casewt)) 
            casewt <- rep(1, n)
        if (!is.null(attr(Terms, "offset"))) 
            warning("Offset term ignored")
        ll <- attr(Terms, "term.labels")
        if (length(ll) == 0) 
            X <- factor(rep(1, n))
        else X <- strata(m[ll])
        ### end of copied code from survfit.

        ### do a separate fit for each level of the factor
        group<-levels(X)
        nstrata<-length(group)
        sbind<-function(x,y){
            if (is.vector(x) & is.vector(y)){
                 if (class(x)=="list" & class(y)=="list"){ 
                     out<-list(x,y)
                 } else out<-c(x,y)
            } else if (is.matrix(x) & is.matrix(y)) out<-cbind(x,y) 
            if (!is.null(attr(x,"LRin")) & !is.null(attr(y,"LRin"))){
                xLRin<-attr(x,"LRin")
                yLRin<-attr(y,"LRin")
                attr(out,"LRin")<-cbind(xLRin,yLRin)
            } 
            return(out)
        }
        ## change original left-hand-side of formula to list with L and R vectors representing left and right endpoints
        Y<- SurvLR(Y)
        for (i in 1:nstrata){
             tempout<-icfit.default(Y$L[X==group[i]],Y$R[X==group[i]],Lin=Y$Lin[X==group[i]],Rin=Y$Rin[X==group[i]],...)
             tempout$strata<-length(tempout$pf)
             names(tempout$strata)<-group[i]
             if (i==1){ icout<-tempout
             } else{
                 icout$A<-NULL
                 tempout$A<-NULL
                 icout<-mapply(sbind,icout,tempout)
             }
        }
        class(icout) <- c("icfit")
        if (!is.null(attr(m, "na.action"))) 
            icout$na.action <- attr(m, "na.action")
    }
    icout$call <- call
    icout
}




`Aintmap`<-function(L,R,Lin=NULL,Rin=NULL){
    n<-length(L)
    if (is.null(Lin) & is.null(Rin)){
        Lin<-rep(FALSE,n)
        Rin<-rep(TRUE,n)
        Lin[L==R]<-TRUE
        Rin[R==Inf]<-FALSE
    } else if (length(Lin)==1 & length(Rin)==1 & is.logical(Lin) & is.logical(Rin) ){
        Lin<-rep(Lin,n)
        Rin<-rep(Rin,n)
    } else if (length(Lin)!=n | length(Rin)!=n | !all(is.logical(Lin)) | !all(is.logical(Rin)) ){
        stop("Lin and Rin should be either NULL, logical length 1 or length same as L,R")
    } 
    if(n != length(R))
        stop("length of L and R must be the same")
    # calculate a small number, eps, to differentiate between e.g.,  [L,R] and (L,R]
    # we will treat, (L,R] as [L+eps,R], and [L,R) as [L,R-eps] 
    # since eps is only 
    # used in ranking we do not need to make it super small
    # just smaller than the smallest difference
    LRvalues<-sort(unique(c(0,L,R,Inf)))
    eps<- min(diff(LRvalues))/2
    Le<-L
    Re<-R
    Le[!Lin]<-L[!Lin]+eps
    Re[!Rin]<-R[!Rin]-eps
    # let s be the vector of ordered L and R values with 
    # R values later when there are ties
    # then intmap are values s[i] and s[i+1] where s[i] is 
    # associated with L and s[i+1] is associated with R
    oLR<-order(c(Le,Re+eps/2) )
    # find the Turnbull intervals, or innermost intervals
    # this is the same as the primary reduction of 
    ### Aragon and Eberly (1992) J of Computational and Graphical
    ###     Statistics 1:129-140
    # label L=1 and R=2
    Leq1.Req2<-c(rep(1,n),rep(2,n))
    # order and see if an R is followed by an L
    # take difference of Leq1.Req2 after putting them in 
    # order, then if the difference is 1 then the R=2 is followed by L=1 
    flag<- c(0,diff( Leq1.Req2[oLR] ))
    R.right.of.L<- (1:(2*n))[flag==1]
    intmapR<- c(L,R)[oLR][R.right.of.L]
    intmapL<- c(L,R)[oLR][R.right.of.L - 1]
    intmapRin<- c(Lin,Rin)[oLR][R.right.of.L]
    intmapLin<- c(Lin,Rin)[oLR][R.right.of.L - 1]
    intmap<-matrix(c(intmapL,intmapR),byrow=TRUE,nrow=2)
    attr(intmap,"LRin")<-matrix(c(intmapLin,intmapRin),byrow=TRUE,nrow=2)
    k<-dim(intmap)[[2]]
    Lbracket<-rep("(",k)
    Lbracket[intmapLin]<-"["
    Rbracket<-rep(")",k)
    Rbracket[intmapRin]<-"]"
    intname<-paste(Lbracket,intmapL,",",intmapR,Rbracket,sep="")
    A<-matrix(0,n,k,dimnames=list(1:n,intname))
    intmapLe<-intmapL
    intmapLe[!intmapLin]<-intmapL[!intmapLin]+eps
    intmapRe<-intmapR
    intmapRe[!intmapRin]<-intmapR[!intmapRin]-eps
    for (i in 1:n){
        tempint<- Le[i]<=intmapRe & Re[i]>=intmapLe
        A[i,tempint]<-1
    }

   # previous versions (<=0.9-9.1) did primary reduction twice,
   # once as described in Turnbull (see above) and once as 
   # described in Aragon and Eberly (1992, J of Computational and Graphical
   #    Statistics 1:129-140) 
   # both do same thing, so we do not need to do it twice

    ## fix error when intmap=(0,Inf) and k=1, previously A was column matrix of 0, should be a column matrix of 1
    if (k==1 & intmap[1,1]==0 & intmap[2,1]==Inf) A[A==0]<-1  

    out<-list(A=A,intmap=intmap)
    out
}



icfitBootCI<-function(L,R,conf.level=.95,B=100,timeEpsilon=10^-8,seed=19439101,messages=TRUE,...){
    if (B<10) stop("B must be at least 10")
    if (!is.null(seed)) set.seed(seed)
    if (messages){
        message("Confidence intervals use modified bootstrap, can be very time consuming.")
        message("See icfitControl help, argument B, for changing the number of bootstrap replicates.")
        utils::flush.console()
    }
    fit0<-icfit(L,R)
    ## just get the bootstrap value at a little bit before and a little bit after each time
    ## use timeEpsilon to define the little bit
    eps<-timeEpsilon
    times<-unique(c(0,as.vector(fit0$intmap)))
    times<-c(0,times-eps,times+eps)
    times<-sort(times[times>=0])
    n<-length(L)

    nt<-length(times)
    LOWER<-UPPER<-matrix(NA,B,nt)

    #fit0<-icfit(L,R)


    ### time it so that have some idea of how long it will take
    t0<-proc.time()
    for (i in 1:10){
        I<-sample(1:n,replace=TRUE)
        fiti<-icfit(L[I],R[I],...)
        LOWER[i,]<-getsurv(times,fiti,nonUMLE.method="right")[[1]]$S
        UPPER[i,]<-getsurv(times,fiti,nonUMLE.method="left")[[1]]$S
    }
    t1<-proc.time()
    if (messages){
        message("Estimated time of one iteration (one stratum only): ",round((t1-t0)[1]/10,2)," sec")
        message("Estimated time for all ",B," iterations  (one stratum only): ",
             round((t1-t0)[1]*B/10,1)," sec")
        utils::flush.console()
    }


    for (i in 11:B){
        I<-sample(1:n,replace=TRUE)
        fiti<-icfit(L[I],R[I],...)
        LOWER[i,]<-getsurv(times,fiti,nonUMLE.method="right")[[1]]$S
        UPPER[i,]<-getsurv(times,fiti,nonUMLE.method="left")[[1]]$S
    }

    percci<-function(Ti,conf.level=.95){
        ### get percentile bootstrap confidence intervals
        ### see Efron and Tibshirani, p. 160 bottom
        alpha<- (1-conf.level)/2
        B<-length(Ti)
        k<-floor((B+1)*alpha)
        if (k==0){
             warning("increase number of bootstrap samples")
             ci<-c(-Inf,Inf)
        } else {
            oTi<-Ti[order(Ti)]
            ci<-oTi[c(k,B+1-k)]
        }
        ci
    }
    calclower<-function(x,CL=conf.level){ percci(x,conf.level=CL)[1] }
    calcupper<-function(x,CL=conf.level){ percci(x,conf.level=CL)[2] }
    lower<-apply(LOWER,2,calclower)
    upper<-apply(UPPER,2,calcupper)

    maxlower<-binom.test(n,n,conf.level=conf.level)$conf.int[1]
    lower[lower>maxlower]<-maxlower
    lower[lower<0]<-0

    minupper<-binom.test(0,n,conf.level=conf.level)$conf.int[2]
    upper[upper<minupper]<-minupper
    upper[upper>1]<-1

    list(time=times, lower=lower, upper= upper, confMethod="modboot", conf.level=conf.level)
}



`icfit.default` <-
function(L, R, initfit = NULL, control=icfitControl(), Lin=NULL, Rin=NULL, conf.int=FALSE,...)
{
    out<-icfitCalc(L, R, initfit, control, Lin, Rin,...)
    if (conf.int){
        if (control$confMethod!="modboot") stop("only modified bootstrap confidence method available")
        out$CI<-icfitBootCI(L,R,conf.level=control$conf.level, B=control$B, timeEpsilon= control$timeEpsilon, seed=control$seed, 
            messages=control$timeMessage,...) 
    }
    out
}


`icfitCalc` <-
function(L, R, initfit = NULL, control=icfitControl(), Lin=NULL, Rin=NULL,...)
{
    epsilon<-control$epsilon
    maxit<-control$maxit
    AI<-Aintmap(L,R,Lin,Rin)
    A<-AI$A
    if (any(apply(A,1,sum)==0)) stop("A row all zeros. Appears that there are some R<L")
    n<-dim(A)[[1]]
    k<-dim(A)[[2]]
    intmap<-AI$intmap
    if (k==1){
        pf<-1
        ### fix error 1/24/2011: needed to change name from numit to count in emout list 
        emout<-list(error=0,count=0,converge=TRUE,message="normal convergence")
        anypzero<-FALSE
    } else {
        ### come up with the initial estimates from the initfit option
        ### may be (1) NULL, (2) character vector giving name of function, 
        ### or (3) icfit or similar object
        ### If of type (1) or (2) we first convert it to type (3) 
        if(is.null(initfit)) {
            pbar <- apply(A/apply(A, 1, sum), 2, mean)
            initfit<-list(pf=pbar,intmap=intmap)
        ## the following else section is for initfit functions
        } else if (is.character(initfit) & length(initfit)==1){
            ## because some initfit functions will input A and some will 
            ## input L,R,Lin, and Rin, we input all 5 variables
            ## but any initfit function need not use all 5
            ## Get options for initfit function from control
            initfitOpts<-control$initfitOpts

            ## since initfit functions may not know how to interpret Lin=NULL and Rin=NULL create the values
            if (is.null(Lin) & is.null(Rin)){
                Lin<-rep(FALSE,length(L))
                Rin<-rep(TRUE,length(L))
                Lin[L==R]<-TRUE
                Rin[R==Inf]<-FALSE
            }
            ## use try function in case initfit function fails
            if (is.null(initfitOpts)){
                initfit<-try( do.call(initfit,args=list(L=L,R=R,Lin=Lin,Rin=Rin,A=A)) )
            } else {
                initfit<-try( do.call(initfit,args=c(list(L=L,R=R,Lin=Lin,Rin=Rin,A=A),initfitOpts)) )
            }
            if (class(initfit)=="try-error"){
                warning("initfit was a character, treated as a function name, and when called gave an error so will not be used")
                pbar <- apply(A/apply(A, 1, sum), 2, mean)
            } else {
                if (is.null(initfit$pf)) stop("initfit treated as function and did not produce list with pf element")
                ## if the initfit function outputs an intmap, check that it matches 
                ## what we have already calculated
                ## do not check attributes of intmap to allow functions that do not output that
                initintmap<-initfit$intmap
                pbar<-initfit$pf
            } 
            if (is.null(initfit$intmap)){ initfit<-list(pf=pbar,intmap=intmap)
            } else initfit<-list(pf=pbar,intmap=initfit$intmap)
        }

        ### Check the initfit:
        ## if the initfit has a different intmap but it has some elements that match 
        ## (as would happen if the initfit function deleted values from the intmap with 0 mass)
        ## the current intmap then we can still try this initfit, 
        ## we use pbar proportional to the initfit$pf values of those intmaps values that match 
        if (is.null(initfit$pf) | is.null(initfit$intmap)) stop("initfit should be either a character function name or a list with elements pf and intmap")
        nkeep<- dim(intmap)[[2]]
        pbar<-rep(0,nkeep)
        for (i in 1:nkeep){
            index<-initfit$intmap[1,]==intmap[1,i] & initfit$intmap[2,]==intmap[2,i]
            if (any(index)){
                if (length(initfit$pf[index])>1) stop("initfit has non-unique intmap columns")
                pbar[i]<- initfit$pf[index]
            }               
        }
        if (sum(pbar)==0) stop("initfit has no matching intmap elements with L and R")
        pbar<-pbar/sum(pbar)

        ## em is the em-algorithm, at any iteration if the estimate of the probability 
        ## mass at any point is lower than the lower.bound, then that probability 
        ## mass is set to zero. Then the Kuhn-Tucker conditions are checked, if they are 
        ## not met then a small mass is added back to those values set to zero. This is the 
        ## polishing methoded of 
        ## Gentleman and Geyer (1994, Biometrika, 618-623)

        ## start count keeps a running total of the number of iterations, since 
        ## em may be called more than once with a different lower.bound 
        em<-function(A,pbar,lower.bound=.01,startcount=1){
            converge<-FALSE 
            message<-"normal convergence"
            A.pbar<-as.vector(A %*% pbar)
            if (any(A.pbar==0)) stop("initfit$pf does not have nonzero values associated with needed intmap spaces")
            J1<-matrix(1,n,1)
            Jn<-matrix(1/n,n,1)
            for (i in startcount:maxit){
                tA.div.A.pbar <- t(A/A.pbar)
                newpbar<- drop((tA.div.A.pbar * pbar) %*% Jn)
                d<-drop(tA.div.A.pbar %*% J1)
                # below is an older slower version 
                #A.div.A.pbar <- A/A.pbar
                #newpbar<- apply(t(A.div.A.pbar)*pbar,1,mean)
                #d <- apply(A.div.A.pbar, 2, sum)
                u <-  - d + n
                u[newpbar > 0] <- 0
                error <- max(d + u - n)
                if (error<epsilon){
                    pbar<-newpbar
                    converge<-TRUE
                    if(any(u < 0)){
                        message<-"Kuhn-Tucker conditions not met, self-consistent estimator not MLE"
                        #pbar[u<0]<-min(pbar[pbar>0])
                        pbar[u<0]<-lower.bound
                        pbar<-pbar/sum(pbar)
                    }  
                    break()
                }
                newpbar[newpbar<lower.bound]<-0
                A.pbar<-as.vector(A %*% newpbar)
                if (any(A.pbar==0)){ 
                    message<-"lower bound too high"
                    break()}
                pbar<-newpbar
                if (i==maxit) message<-"maxit reached"
            }
            out<-list(A=A,pbar=pbar,count=i,message=message,converge=converge,error=error)
            out
        }

        ## try different values to "polish" the estimates, where if 
        ## the jump in the distribution (mass at any iterval) is less than the lower.bound then 
        ## set it equal to zero. If it turns out that a jump should not 
        ## have been set to zero, the em function checks that and spits out 
        ## the last pbar before the jumps were set to zero 
        ## (see Gentleman and Geyer, 1994)
        lower.bounds<- 10^(0:ceiling(log10(epsilon)))
        lower.bounds<-lower.bounds[lower.bounds<=max(1/n,epsilon)]
        emout<-list(pbar=pbar,count=0)
        for (lb in lower.bounds){
            emout<-em(A,pbar=emout$pbar,lower.bound=lb,startcount=emout$count+1)
            if (emout$message=="normal convergence") break()
            if (emout$count==maxit) {
                emout$message<-"problem with convergence, increase maxit"
                break()
            }
            #print(emout$message)
        }
        keep<- !emout$pbar==0
        anypzero<-FALSE
        if (!all(keep)){
            LRin<-attr(intmap,"LRin")[,keep]
            intmap<-intmap[,keep]
            attr(intmap,"LRin")<-LRin
            A<-A[,keep]
            anypzero<-TRUE
        }
        pf<-emout$pbar[keep]
    } # end else for k>1
    strata<-length(pf)
    ## if the A matrix corresponding to the non-zero pf values is full rank, then the NPMLE is mixture unique
    ## see Gentleman and Geyer, 1994, Biometrika 618- or Gentleman and Vandal, 2002 Can J Stat, 557- 
    ## DO NOT NEED THIS for univariate interval censored data because it will always be mixture unique
    ## munique<-qr(A)$rank==dim(A)[2]
    # if there is only one strata, title describes output:  NPMLE
    names(strata)<-"NPMLE"
    out <- list(A=A, strata=strata, error = emout$error, numit = emout$count, pf = pf, intmap = 
        intmap, converge= emout$converge, message= emout$message, anypzero=anypzero)
    class(out)<-c("icfit","list")
    return(out)
}

