###
# do we want to use rlogit linear comb to initialize? it may or may not be a good thing

# s=1; lambda=1; kernel="linear"; para=NULL; start.method="rlogit"; eta0.init=NULL; eta.diff.init=NULL; maxit=10; tol=1e-5; verbose=FALSE; minQuad.tol=NULL; ret.vcov=FALSE; minQuad.control = control.minQuad();mem.efficient=TRUE;init.alpha.from.previous=TRUE #default
rauc=function (formula, dat, s=1, lambda=1, kernel="linear", para=NULL, 
    start.method="rlogit", eta0.init=NULL, beta.init=NULL, eta.diff.init=NULL, # initial values
    maxit=50, tol=1e-5, # stopping criteria
    minQuad.control = control.minQuad(), 
    init.alpha.from.previous=TRUE, mem.efficient = TRUE, # performance gain in speed and space
    ret.vcov=FALSE, garbage.collection = TRUE, verbose=FALSE, ...)
{
        
    time0=Sys.time()
    
    tmp=model.frame(formula, dat)
    y1=tmp[,1]
    X1=model.matrix(formula, tmp[y1==1,])[,-1,drop=FALSE]
    X2=model.matrix(formula, tmp[y1==0,])[,-1,drop=FALSE]
    n1=nrow(X1)
    n2=nrow(X2)
    n=n1+n2
    p=ncol(X1)
    
    X.diff = get.X.diff(X1, X2)
    
    kernel=substr(kernel,1,1)
    if (kernel!="l" & is.null(para)) stop("kernel parameter not set") 
    
    # We stack case on top of control. To access j in the control, we add n1
    X=rbind(X1,X2)
    K=getK(X,kernel,para)/lambda
    if (!is.double(K)) K <- as.double(K) #this check is free but we must make
    if (!is.integer(n1)) n1 <- as.integer(n1)
    if (!is.integer(n2)) n2 <- as.integer(n2)
    
    Q <- NULL
    if(!mem.efficient){
        Q <-try(matrix(0.0,n1*n2, n1*n2)) 
        if(inherits(Q,'try-error')){ #failing to allocate large matrix 'Q'
            Q <- NULL; gc(); mem.efficient <- TRUE
        }else{
            aux=.C("get_Q", K=K, n1=n1, n2=n2, Q=Q, DUP = TRUE, NAOK = FALSE,PACKAGE = "aucm") # Q is influenced by lambda through K
        }
    }
    
    # set eta.init
    if (!is.null(eta.diff.init)) {
        
        if (start.method=="svm") {
#            # use svm fit to initialize
#            s=1
#            fit0 = svml (formula, alldat, kernel="r", fitted=TRUE, cost=1) 
#            eta.init= logit(fit0$fitted.values$posterior[,1])
#            n1=sum(alldat$y); n2=nrow(alldat)-n1
#            eta.diff.init=eta.init[rep(1:n1,each=n2)] - eta.init[rep(1:n2,n1)]
#            s.f=15 # this based on having saturation reach 70-80%
#            mean(eta.diff.init*s.f>.5*s | eta.diff.init*s.f< -.5*s)
#            eta.diff.init=eta.diff.init*s.f
            stop("not supported yet")
        } else {
            eta.init=eta.diff.init
        }        
        
    } else if (!is.null(eta0.init)) {
    
        eta.i.init = eta0.init[y1==1]
        eta.j.init = eta0.init[y1==0]
        eta.init=eta.i.init[rep(1:n1,each=n2)] - eta.j.init[rep(1:n2,n1)]
    
    } else {
        
        if (kernel=="l"){
            if (is.null(beta.init)){
                # use beta.init to create eta.init
                if (start.method=="rlogit") {
                    fit.rlogit=rlogit(formula, dat)
                    if (fit.rlogit$convergence) {
                        beta.init=coef(fit.rlogit)[-1]
                    } else {
                        beta.init=rep(1, ncol(X1))
                    }                
                } else if (start.method=="1") {
                    beta.init=rep(1, ncol(X1))
                } else if (start.method=="0") {
                    beta.init=rep(0, ncol(X1))
                } else stop("start.method not supported: "%+%start.method)
            }                             
            if (verbose) myprint(beta.init)
            eta.init=X.diff%*%beta.init       
        } else {
            eta.init=rep(0,nrow(X.diff))
        }
    
    }
    val.init=mean(ramp.f(eta.init,s))
    

    arg <- list(...)   
    # args=list() # useful when running this functine line by lin
    n1n2 <- n1*n2
    alpha0 <- arg$alpha #alpha0 initializes minQuad at each iteration 
    if(!length(alpha0))alpha0 <- 0.5
    alpha0 <- rep(0.5,length.out = n1n2)
    Qx <- double(n1n2);
    b <- double(n1n2);
    eta.new <- double(n1n2)
    alpha.pred <- double(n1n2)
    dh2.deta <- double(n1n2)
    
    iter=0
    converge=0
    minQuad.nonconvergence.list=c()    
    minQuad.elapsed <- double(maxit)
    penalized.losses <- losses <- NULL
    
    while (TRUE){
    
        iter=iter+1; if ( iter>maxit ) {converge=1; break;}
        
        # if dh2.deta are not changed from before, no change will happen
        dh2.deta=ifelse(eta.init < -.5*s, -1/s, 0.0)

        b <- NULL
        if(!mem.efficient){
            b=s^2*(1 - 2/s* Q %*% dh2.deta) 
        }else{
            # the following needs to be changed to .Call
            Qx=.C("get_Qx",as.double(.Machine$double.eps), K=K, n1=n1, n2=n2, x=dh2.deta, y = double(n1n2), DUP = TRUE, NAOK = FALSE,PACKAGE = "aucm")$y
            b=s^2*(1 - 2/s*Qx)
        }
        
        # Krisz: for debugging 
        if(length(arg$return.data) && (arg$return.data == -1)) return(list(Q=Q,K=K,b=-0.5*b,n1=n1,n2=n2,Qx=Qx))
        if(length(arg$SAVE) && arg$SAVE)dput(b,"b_iter_"%+%iter%+%".RDput")
    
##      minQuad: optimize f(a) = 0.5 * a'Qa + b'a, 0 <= a <= C. Note the different parametrization from a'Qa-b'a        
        begin <- Sys.time()
        fit = minQuad(H=if(mem.efficient) K else Q, b = -0.5*b, C = 1.0, n1 = n1, n2 = n2,alpha = alpha0,control = minQuad.control,mem.efficient = mem.efficient)
        end <- Sys.time()
        minQuad.elapsed[iter] <- as.double(difftime(end,begin,units = 'secs'))
        if(init.alpha.from.previous) alpha0 <- fit$alpha

        if(length(arg$SAVE) && arg$SAVE)dput(alpha0,"a_iter_"%+%iter%+%".RDput")

        if (verbose) cat(iter)
        if (fit$convergence!=0){
            minQuad.nonconvergence.list = c(minQuad.nonconvergence.list, iter)
        }

        # compute eta
        alpha.pred=dh2.deta + 1/s*fit$alpha # always using the newly fitted alpha and init eta
        if(!mem.efficient){
            eta.new = Q %*% alpha.pred
        }else{
            # the following needs to be changed to .Call
            eta.new=.C("get_Qx", as.double(.Machine$double.eps), K=K, n1=n1, n2=n2, x=alpha.pred,y = double(n1n2), DUP = TRUE, NAOK = FALSE,PACKAGE = "aucm")$y
        }
        
        if (verbose) cat (", saturated ",signif(mean(eta.new>.5*s | eta.new< -.5*s),3) )
        
        val.new=mean(ramp.f(eta.new,s))
        pen = drop(crossprod(alpha.pred, eta.new))/n1/n2/2
        penalized.losses=c(penalized.losses, val.new+pen)
        losses = c(losses, val.new)
        
        if (verbose) cat (", loss ",val.new, ", pen ",pen)
        
        # check stopping criteria and udpate eta.init
        if (kernel=="l"){
            beta.new =  drop(crossprod (X.diff, alpha.pred)/lambda)
            delta=max(abs((beta.new-beta.init)/beta.init))
            if(verbose) cat (", beta[1] ",signif(beta.new[1],4),", beta[-1]/beta[1] ",signif(beta.new[-1]/beta.new[1],4)," delta ",delta)
            if(delta < tol) {
                converge=0
                break;
            } else {
                beta.init=beta.new
                eta.init=eta.new 
            }
        } else {
            val.new=val.new + pen
            delta=val.new - val.init
            if(verbose) cat (" delta ",delta)            
            if (abs(delta) < tol) {
                converge=0
                break;
            } else {
                val.init=val.new
                eta.init=eta.new 
            }
        }
        if (verbose) cat (" eps ",fit$epsilon,"\n")
        
    }  # end while loop
    
    if (verbose & length(minQuad.nonconvergence.list)>0) print("minQuad did not converge in these iterations: "%+%concatList(minQuad.nonconvergence.list,","))    
    if (fit$convergence!=0) {
        print("minQuad did not converge in the last iteration")
        if (!verbose) print("minQuad did not converge in these iterations: "%+%concatList(minQuad.nonconvergence.list,","))        
    }
    if (verbose) cat ("\n")
    if (verbose) print("total time used: "%+%(Sys.time()-time0)%+%" "%+%attr((Sys.time()-time0),"units") )
    
    
    
#    #get linear combination
#    #R-code below has extra bits in 'extended double' to correct for accumulation of error in 
#    #matrix multiplication. As our sample size 'n' <= ~500 is expected, we do not worry about 
#    #this and use 'C' code that is more memory efficient and faster - takes advantage of 
#    #alpha's being zero.  
#    Q.pred = getQ(K,n1=n1,n2=n2,call.C=T,do.pred=TRUE)/lambda
#    linear.combination=Q.pred %*% alpha.pred    

     #fast.auc(Q.pred %*% (dh2.deta + alpha),)
    if (length(alpha.pred)!=n1*n2) stop("length(alpha.pred)!=n1*n2") # for whatever reason (e.g. the aborted pAUC development) if this happens, it is bad
    #.C("get_Q_pred_x",.Machine$double.eps, K=K, n1=n1, n2=n2,npred = n,x=alpha.pred,Qx = lin.comb, DUP = FALSE, NAOK = FALSE,PACKAGE = "aucm")
    lin.comb = .Call("get_Q_pred_x", K=K, n1=n1, n2=n2, x=alpha.pred, .Machine$double.eps)
    linear.combination = lin.comb / lambda # DEEP COPY
        
    result=list(
        convergence=converge,
        converged=ifelse(converge==0,TRUE,FALSE),
        iterations=iter,
        time.elapsed = difftime(Sys.time(),time0,units = 'secs'),
        #time.elapsed.minQuad = minQuad.elapsed,
        
        coefficients=NA,
        ratios=NA,
        alpha.pred=drop(alpha.pred),
        linear.combination=drop(linear.combination),
        last.minQuad.fit=fit,
        train.auc=fast.auc(drop(linear.combination), c(rep(1,n1),rep(0,n2))),
        saturation=mean(eta.new>.5*s | eta.new< -.5*s),
        losses=losses,
        penalized.losses=penalized.losses,
        
        formula=formula,
        X=X,
        y=c(rep(1,n1),rep(0,n2)),
        n.case=n1,
        n.control=n2,
        lambda=lambda,
        kernel=kernel,
        para=para
    )
    class(result) = c("rauc","auc",class(result))
    
    if (kernel=="l") {
        result$coefficients = drop(t(t(alpha.pred) %*% X.diff)/lambda)
        result$ratios = result$coefficients/result$coefficients[1]
    }
    
    if (result$convergence!=0) {
        print("rauc did not converge")
#        result$coefficients=NULL
#        result$alpha.pred=NULL
#        result$linear.combination=NULL
    }
    if(length(arg$return.data) && arg$return.data) result <- c(result,list(Q=Q,K=K,b=-0.5*b,n1=n1,n2=n2,Qx=Qx))

    if(garbage.collection)gc()
 
    if (ret.vcov) result$vcov = get.vcov.sauc ("normal", beta.hat=result$coefficients/result$coefficients[1], h=0.4/result$coefficients[1], X1, X2)     
    return(result)
    
}


# ramp function or h
ramp.f=function(eta,s,loss=TRUE){
    if (s==0) {
        out=ifelse(eta < 0, 1, 0)
    } else {
        out=ifelse(eta > -.5*s & eta < .5*s, .5 - 1/s*eta, 0)
        out[eta <= -.5*s] = 1
    }
    if (loss) drop(out) else drop(1-out)
}
# test ramp.f
# eta=seq(-5,5,length=100); plot(eta, ramp.f(eta,s=1))



## deprecated
## minQuad.tol=NULL; verbose=FALSE; maxit=1e2; reltol=1e-3; fix.beta1=FALSE # default arguments
#rauc.linear=function (formula, dat, lambda, s=1, beta.init=NULL, fix.beta1=FALSE, t0=NULL, t1=NULL, minQuad.tol=NULL, verbose=FALSE, maxit=1e2, reltol=1e-3){
#    
#    time0=Sys.time()
#    
#    tmp=model.frame(formula, dat)
#    x1=model.matrix(formula, tmp[tmp[,1]==1,])[,-1]
#    x2=model.matrix(formula, tmp[tmp[,1]==0,])[,-1]
#    n1=nrow(x1)
#    n2=nrow(x2)
#    X=rbind(x1,x2)
#    
#    X.diff=get.X.diff(x1, x2)
#    if(fix.beta1){
#        z=X.diff[,1,drop=FALSE]
#        X.diff=X.diff[,-1,drop=FALSE]
#    }
#    
#    # form Q, note that (X.diff%*%t(X.diff))/lambda would throw an error: cannot allocate vector of size 759.1 Mb
#    Q = X.diff %*% diag(1/lambda, ncol(X.diff)) %*%t(X.diff)
#    
#    pAUC=ifelse(!is.null(t0) | !is.null(t1), TRUE, FALSE)
#    if (pAUC) {
#        if (verbose) print("pAUC")
#        if (is.null(t0)) t0=0
#        if (is.null(t1)) t1=1            
#    } else {
#        A=rep(TRUE, n1*n2)         
#    }        
#    
#    iter=0
#    converge=0
#    minQuad.time=0
#    counter1=-1
#    rauc.old=Inf
#    while (TRUE){
#    
#        iter=iter+1
#        if (iter>maxit){
#            converge=1
#            break;
#        }
#        if (verbose) cat("iter ",iter)
#        
#        period=1 # changing period will change how often restricted set gets updated relative to updating model coefficients
#        counter1=(counter1+1) %% period 
#        if (pAUC & counter1==0) {
#            #if(verbose) print("update restricted set")
#            # compute restricted set
#            if (is.null(beta.init)) eta.j.init=rep(0,nrow(x2)) else eta.j.init=drop(x2 %*% beta.init)
#            # this is faster than ind1, a little different from ind1 due to different definition of quantile
#            ind2=order(eta.j.init)[ceiling((1-t1)*n2):floor((1-t0)*n2)] # note that if t1 is 0, the begining of the index is 0, it is ok in R, but be careful in C
#            #if (verbose) myprint(length(ind2))
#            #if (verbose) myprint(ind2)
##           # an alternative implemenation for ind2, slower    
##            q0=quantile(eta.j.init,(1-t0))
##            q1=quantile(eta.j.init,(1-t1))    
##            ind1=eta.j.init<=q0 & eta.j.init>=q1
##            unname(which(ind1))
#            A=rep(FALSE,n2)
#            A[ind2]=TRUE            
#            A= rep(A, n1)    
#            X.diff.sub=X.diff[A,]        
#            Q.sub=Q[A,A]               
#        }
#        
#        if (!fix.beta1) {
#            # beta1 is not fixed to be 1
#            if(!pAUC){
#                if (!is.null(beta.init)) eta.init=X.diff %*% beta.init else eta.init=rep(0, nrow(X.diff))
#                dh2.deta=ifelse(eta.init < -.5/s, -s, 0)
#                b=s^(-2)*( 1 - 2*s* Q %*% dh2.deta )
#            } else {
#                if (!is.null(beta.init)) eta.init=X.diff.sub %*% beta.init else eta.init=rep(0, nrow(X.diff.sub))
#                dh2.deta=ifelse(eta.init < -.5/s, -s, 0)
#                b=s^(-2)*( 1 - 2*s* Q.sub %*% dh2.deta )
#            }
#        } else {
#            # beta1 is fixed to be 1
#            eta.init=z + X.diff %*% beta.init
#            dh2.deta=ifelse(eta.init < -.5/s, -s, 0)
#            b=s^(-2)*( (1-2*s*z) - 2*s* Q %*% dh2.deta )
#        }
#        
#        if (pAUC) {
#            # this can only be computed after a new working set is computed from eta.init
#            rauc.new = mean(ramp.f(eta.init,s,loss=TRUE))
#            if (rauc.new > rauc.old) {
#                converge=0
#                print("stop due to increasing rauc loss")
#                break;
#            } else {
#                myprint(rauc.new, digits=10)
#                rauc.old = rauc.new
#            }
#        }                   
#    
#        start.time=Sys.time()
#        # doing [A,A] is very expensive
#        if(!pAUC) {
#            fit=minQuad(H=Q,b=b, C=1,control = control.minQuad(tol = minQuad.tol))
#        } else {
#            fit=minQuad(H=Q.sub,b=b, C=1,control = control.minQuad(tol = minQuad.tol))        
#        }
#        if (fit$convergence!=0) {
#            print("minQuad does not converge")
#            warning("minQuad does not converge")
#        }
#        minQuad.time=minQuad.time+Sys.time()-start.time
#        
#        if (!pAUC) {
#            beta.new = t(t(dh2.deta + s*fit$alpha) %*% X.diff)/lambda 
#        }else {
#            beta.new = t(t(dh2.deta + s*fit$alpha) %*% X.diff.sub)/lambda 
#        }
#        
#        if(verbose) {
#            if (fix.beta1) {
#                myprint(c(beta.new))
#            } else {
#                #myprint(c(beta.new[1], beta.new[-1]/beta.new[1]))
#                cat (", beta[1] ",beta.new[1],", beta[-1]/beta[1] ",beta.new[-1]/beta.new[1])
#            }
#            cat (", saturated ",mean(eta.init>.5/s | eta.init< -.5/s))
#            #myprint(mean(eta.init>.5/s | eta.init< -.5/s))
#            cat ("\n")
#        }
#        
#        if (is.null(beta.init)) beta.init=rep(1,length(beta.new))
#        if(max(abs((beta.new-beta.init)/beta.init))<reltol) {
#            converge=0
#            break;
#        } else {            
#            beta.init=beta.new
#        }
#            
#    }    
#    
#    if(verbose) print("total time used: "%+%(Sys.time()-time0)%+%" "%+%attr((Sys.time()-time0),"units")%+%", minQuad time used: "%+%minQuad.time%+%" "%+%attr(minQuad.time,"units"))
#    
#    result=list(
#        converge=converge,
#        iteration=iter,     
#        coefficients=c(beta.new),
#        last.minQuad.fit=fit,
#        lambda=lambda,
#        X=X,
#        y=c(rep(1,n1),rep(0,n2)),
#        n.case=n1,
#        n.control=n2,
#        formula=formula
#    )
#    class(result) = c("rauc","auc",class(result))
#   
#    if (result$converge!=0) result$coefficients=NA
#    
#    return(result)
#    
#}


# some debugging code to use within rauc
#        # check out trouble cases
#        fit1=minQuad(Q,b,C=1,tol=minQuad.tol, maxit=9972, trace=2)
#        fit2=minQuad(Q,b,C=1,tol=minQuad.tol, maxit=1, trace=2, alpha=fit$alpha)
#        fit3=C_minQuad(Q,b,C=1,maxit=1, trace=2, alpha=fit$alpha)
#
#        fit4=minQuad(Q,b,C=1,trace=2,q=1)#solveQuad
#        fit5=C_solveQuad(Q,b,C=1)#very slow
#
#        c=Q%*%fit$alpha
#        der=2*c-b
#        aux=rep(-Inf,nrow(der))
#        aux[fit$alpha<C]=-der[fit$alpha<C]
#        i=which.max(aux)
#        i; aux[i]; fit$alpha[i]
#        aux=rep(Inf,nrow(der))
#        aux[fit$alpha>0]=-der[fit$alpha>0]
#        j=which.min(aux)
#        j; aux[j]; fit$alpha[j]
#        
#        B=c(4702,8397)
#        Q[B,B]
#        solve(Q[B,B])
#        
#        # end
        

###function penalty gives the matrix P
###lamda and a are the penalty parameters
### beta is the estimated parameter
###type could be "L2" or SCAD 
###penalty=function(lambda,a,beta,type)
###{
###dim=length(beta)
###result=matrix(0,dim,dim)
###if (type=="L2")
### {
###  diag(result)=1
###  } else if (type=="SCAD")
###     {diag(result)=lambda*(ifelse(abs(beta)<=lambda,1,0)+ifelse(a*lambda-abs(beta)>0,a*lambda-abs(beta),0)*ifelse(abs(beta)>lambda,1,0)/((a-1)*lambda))*ifelse(abs(beta)<0.0001,0,1)/abs(beta)
###      }
###return(result)
###}


# sv.rauc needs to be passed in for now
plot.rauc=function(x, dat.train, fit, intercept, nsvs, xlab=expression(x[i1]), ylab=expression(x[i2]), main=expression(x[i]), cex=.8, cex.sv=.8, pch.sv=19, ...) {
    plot(x, dat.train, col=ifelse(dat.train$y==1,2,1), cex=cex, xlab=xlab, ylab=ylab, main=main)
    points(x, dat.train[nsvs,], col=4, cex=cex.sv, pch=pch.sv)
    abline(intercept, -coef(fit)[1]/coef(fit)[2]) # separating plane from rauc
    mylegend(x=7,col=c(2,1,4),legend=c("case","non-case","non-support vector"),pch=c(1,1,pch.sv),cex=cex)    
}




vcov.rauc=function(object, ...) {
    object$vcov
}

nsv <- function(object) UseMethod("nsv") 

# this method only works if all cases are on top of controls
nsv.rauc=function(object) {
    n1=sum(object$y)
    n2=length(object$y)-n1    
    ind=cbind(rep(1:n1,each=n2), rep(n1+1:n2,n1))
    
    #alpha=object$last.minQuad.fit$alpha
    alpha=object$alpha.pred
    col.=sapply(alpha, function (a) ifelse (a==0, 1, ifelse (a==1, 2, 4)))
    nsv.rauc=setdiff(1:(n1+n2), unique(c(ind[col.==2 | col.==4,]))) # observations that do not affect alpha completely
    nsv.rauc
}
