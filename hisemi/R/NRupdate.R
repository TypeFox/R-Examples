NRupdate=function(f, starts, gradient, hessian, ..., ridge0=1e-6, tolerance=sqrt(.Machine$double.eps), 
    iter.max=1500, halving.max=Inf, relative=FALSE, return.hessian=FALSE,debugging=FALSE)
{
    nparms=length(starts)
    theta=starts
    fval=f(starts,...)
    having.grad=!is.null(attr(fval,'gradient'))
    having.hess=!is.null(attr(fval,'hessian'))
    old.grad=NA

    quit.flag=FALSE
    iter.count=0L
    repeat{ ## NR loop
        grad=if(having.grad) attr(fval,'gradient') else gradient(theta)
        hess=if(having.hess) attr(fval,'hessian')  else hessian(theta)
        if(debugging)cat("\tmax.grad=",max(abs(grad)),fill=TRUE)

############    ridging
        ridge.fact=0
        ridge.count=0L
        repeat{ ## ridging loop

####################### method 1            
#            meval=min.eigval(hess)
#            if(meval<tolerance){
#                ridge.fact=if(ridge.fact==0) ridge0 else 9*ridge.fact
##                ridge.fact=if(ridge.fact==0) abs(meval)+tolerance else 9*ridge.fact
#                diag(hess)=diag(hess)+ridge.fact
#            }else{
#                change.amount=drop(solve(hess,grad))
#                break
#            }
####################### method 1 end

####################### method 2  (reasonably good for small problems)
#            try.err=try({
#                hess.chol=suppressWarnings(chol(hess))
#            },silent=TRUE)                                   
#            if(class(try.err)=='try-error'){
#                ridge.fact=if(ridge.fact==0) ridge0 else 9*ridge.fact
#                diag(hess)=diag(hess)+ridge.fact
#            }else{
#                change.amount=drop(solve(hess.chol,solve(t(hess.chol),grad)))
#                break
#            }
####################### method 2 end

####################### method 3  (not good when hess is not definite)
#            try.err=try({
#                change.amount=drop(solve(hess, grad))
#            },silent=TRUE)                                 
#            if(class(try.err)=='try-error'){
##############   convergence check with non P.D. Hessian
#                rel.conv.H=abs(drop(crossprod(grad)))/abs(fval)
#                if(rel.conv.H<tolerance) {quit.flag=TRUE; break}
#
#                ridge.fact=if(ridge.fact==0) ridge0 else 9*ridge.fact
#                diag(hess)=diag(hess)+ridge.fact
#            }else{
#                break
#            }
####################### method 3 end

####################### method 4
#            ridge.fact=-4
#            change.amount=modchl.solve(hess,grad)
#            break
####################### method 4 end

####################### method 4.1
#            try.err=try({
#                hess.chol=suppressWarnings(chol(hess))
#            },silent=TRUE)                                   
#            if(class(try.err)=='try-error'){
#                ridge.fact=-4.1
#                change.amount=modchl.solve(hess,grad)
#            }else{
#                change.amount=drop(solve(hess.chol,solve(t(hess.chol),grad)))
#            }
#            break
####################### method 4.1 end

####################### method 4.2
            try.err=try({
                hess.chol=suppressWarnings(chol(hess))
            },silent=TRUE)                                   
            if(class(try.err)=='try-error'){
                ridge.fact=-4.2
                change.amount=drop(solve(nearPD(as.matrix(hess),posd.tol=1e-4,eig.tol=1e-4)$mat, grad))
            }else{
                change.amount=drop(solve(hess.chol,solve(t(hess.chol),grad)))
            }
            break
####################### method 4.2 end

####################### method 5 (steepest descent, not using hessian information) 
#            try.err=try({
#                hess.chol=suppressWarnings(chol(hess))
#            },silent=TRUE)                                   
#            if(class(try.err)=='try-error'){
#                ridge.fact=-5
#                change.amount=drop(grad)/512
#            }else{
#                change.amount=drop(solve(hess.chol,solve(t(hess.chol),grad)))
#            }
#            break
####################### method 5 end

####################### method 6 (BFGS if hessian is not p.d.)
#            try.err=try({
#                hess.chol=suppressWarnings(chol(hess))
#            },silent=TRUE)                                   
#            if(class(try.err)=='try-error'){
#                ridge.fact=-6
#                if(is.na(old.grad[1])){B.inv=diag(1,nparms)/512
#                }else{
#                    y=grad-old.grad
#                    sty=drop(crossprod(s,y))
#                    By=B.inv%*%y
#                    Bys=tcrossprod(By,s)
#                    B.inv=B.inv+tcrossprod(s)*(sty+sum(y*By))/(sty*sty)-(Bys+t(Bys))/sty
#                }
#                change.amount=drop(B.inv%*%grad)
#            }else{
#                B.inv=as.matrix(solve(hess.chol,solve(t(hess.chol),diag(1,nparms))))
#                change.amount=drop(B.inv%*%grad)
#            }
#            break
####################### method 6 end


            ridge.count=ridge.count+1L
        }

#        rel.conv.H=abs(drop(grad%*%change.amount))/if(relative)abs(fval)else 1
        rel.conv.H=max(abs(grad))
        if(rel.conv.H<tolerance || iter.count>=iter.max) {break}

############    update with step-halving
#        theta.new=theta-change.amount
#        halving.count=0L
#        repeat{
#            fval.new=f(theta.new,...)
#            if(fval.new<=fval +1e-4*crossprod(grad,theta.new-theta) || halving.count>=halving.max) break
#            theta.new=0.5*(theta.new+theta)
#            halving.count=halving.count+1L
#        }

############    update with golden-section
        theta.new=theta-change.amount
        fval.new=f(theta.new,...)
        if(fval.new> fval +1e-4*crossprod(grad,theta.new-theta)){
            ff=function(a,...)f(theta-a*change.amount,...)
            opt.a=optimize(ff,c(1e-3, 1),tol=1e-2)
            fval.new=opt.a$objective
            halving.count=-opt.a$minimum
            theta.new=theta-opt.a$minimum*change.amount
        }
        halving.count=-1L
############    for use with method 6 (BFGS) only
        if(ridge.fact==-6){
            s=theta.new-theta
            old.grad=grad
        }
############    end (BFGS)
        fval=fval.new
        theta=theta.new
        if(having.grad && is.null(attr(fval,'gradient'))) fval=f(theta.new,...)
        iter.count=iter.count+1L

        if(debugging)cat("iter=",iter.count, "\tf=", fval, "\tridge.fact=",ridge.fact, "\tridge.count=", ridge.count, "\thalving.count",halving.count)
    }

    grad=if(having.grad) attr(fval,'gradient') else gradient(theta)
    hess=if(having.hess) attr(fval,'hessian')  else hessian(theta)
    ans=theta
    attr(ans, 'objective')=fval
    attr(ans, 'gradient')=grad
    attr(ans, 'iter')=iter.count
    if(return.hessian) attr(ans, 'hessian')=hess
    ans
}
