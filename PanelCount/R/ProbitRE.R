LL_ProbitRE = function(par,z,w,group,H,verbose=1){
    if(length(par) != ncol(w)+1)  stop("Number of parameters incorrect")
    alpha = par[1:ncol(w)]
    delta = par[length(par)]

    N = length(group)-1
    obs = length(z)

    rule = gauss.quad(H, "hermite")
    Omega = rule$weights
    v = rule$nodes

    wa = as.vector(w %*% alpha)
    d = 2*z - 1

    Li = rep(0, N)
    for(k in 1:H){
        q = wa + sqrt(2)*delta*v[k]
        Phi = pnorm(d*q)
        Phi_prod = as.vector(groupProd(Phi, group)) # @@@
        Li = Li + Omega[k]*Phi_prod
    }
    LL = sum(log(Li/sqrt(pi)))
    if(verbose>=1){
        writeLines(paste("==== Iteration ", tmp.env$iter, ": LL=",round(LL,digits=5)," =====", sep=""))
        print(round(par,digits=3))
    }
    tmp.env$iter = tmp.env$iter + 1
    if(is.na(LL) || !is.finite(LL)){
        if(verbose>=2) writeLines("NA or infinite likelihood, will try others")
        LL = -1e300
    }
    return (LL)
}


# 2. Gradient function
Gradient_ProbitRE = function(par,z,w,group,H,variance=F,verbose=1){
    if(length(par) != ncol(w)+1)  stop("Number of parameters incorrect")
    alpha = par[1:ncol(w)]
    delta = par[length(par)]

    N = length(group)-1
    obs = length(z)

    rule = gauss.quad(H, "hermite")
    Omega = rule$weights
    v = rule$nodes

    wa = as.vector(w %*% alpha)
    d = 2*z - 1
    w_ext = cbind(w,delta=0)

    Li = rep(0, N)
    sumK_alp = matrix(0, N, ncol(w_ext))
    for(k in 1:H){
        dq = d * (wa + sqrt(2)*delta*v[k])
        phi = dnorm(dq)
        Phi = pnorm(dq)
        Phi_prod = as.vector(groupProd(Phi, group)) # @@@
        Li = Li + Omega[k]*Phi_prod

        sumT_alp = matVecProdSum(w_ext, sqrt(2)*v[k], d*phi/Phi, group) # @@@
        sumK_alp = sumK_alp + matVecProd(sumT_alp, Omega[k]*Phi_prod)
    }
    Li = pmax(Li, 1e-100) # in case that some Li=0
    dLogLi = matVecProd(sumK_alp, 1/Li) # sqrt(pi) cancells out
    colnames(dLogLi) = names(par)
    gradient = colSums(dLogLi)
    if(variance){
        if(verbose>=1){
            writeLines("----Converged, the gradient at optimum:")
            print(round(gradient,digits=3))
        }
        var = solve(crossprod(dLogLi))
        return (list(g=gradient, var=var, I=crossprod(dLogLi)))
    }
    if(verbose>=2){
        writeLines("----Gradient:")
        print(round(gradient,digits=3))
    }
    if(any(is.na(gradient) | !is.finite(gradient))){
        if(verbose>=2) writeLines("NA or infinite gradient, reset to all -1")
        gradient = rep(-1, length(gradient))
    }
    return (gradient)
}

#' A Probit Model with Random Effects
#' @description Estimate a Probit model with random effects
#' @param formula Formula of the model
#' @param data Input data, a data frame
#' @param id A vector that represents the identity of individuals, numeric or character
#' @param delta Variance of random effects on the individual level for ProbitRE
#' @param max_delta Largest allowed initial delta
#' @param method Searching algorithm, don't change default unless you know what you are doing
#' @param lower Lower bound for estiamtes
#' @param upper Upper bound for estimates
#' @param H A vector of length 2, specifying the number of points for inner and outer Quadratures
#' @param accu L-BFGS-B only, 1e12 for low accuracy; 1e7 for moderate accuracy; 10.0 for extremely high accuracy. See optim
#' @param reltol Relative convergence tolerance. default typically 1e-8
#' @param verbose Level of output during estimation. Lowest is 0.
#' @param tol_gtHg tolerance on gtHg, not informative for L-BFGS-B
#' @return A list containing the results of the estimated model
#' @examples
#' \donttest{
#' data(rt)
#' est = ProbitRE(isRetweet~fans+tweets+as.factor(tweet.id),
#'                     id=rt$user.id, data=rt)
#' }
#' @export
#' @family PanelCount
ProbitRE = function(formula, id, data=NULL, delta=1, max_delta=3, method='BFGS',lower=NULL,upper=NULL,H=20,accu=10,reltol=1e-8,verbose=0,tol_gtHg=Inf){
    # 1.1 Sort data based on id
    ord = order(id)
    data = data[ord,]
    id = id[ord]
    group = c(0,cumsum(table(as.integer(factor(id)))))
    # 1.2 Parse w and z
    mf = model.frame(formula, data=data, na.action=NULL, drop.unused.levels=T)
    z = model.response(mf, "numeric")
    w = model.matrix(attr(mf, "terms"), data=mf)
    if(length(id)!=nrow(w)) stop("Error: id and w length don't match, potentially due to removal of data points with missing values!")
    # 1.3 use probit estimates as initial values for parameters
    probit = summary(glm(formula, family=binomial(link='probit'), data=data))
    par = c(probit$coefficients[,1], delta=delta)
    # 2. Estimation
    tmp.env <<- new.env()
    tmp.env$iter = 1
    begin = Sys.time()
    if(method=="L-BFGS-B"){
        res = optim(par=par, fn=LL_ProbitRE, gr=Gradient_ProbitRE, method="L-BFGS-B", control=list(factr=accu,fnscale=-1), lower=lower, upper=upper, z=z, w=w, group=group, H=H,verbose=verbose)
    } else {
        res = optim(par=par, fn=LL_ProbitRE, gr=Gradient_ProbitRE, method="BFGS", control=list(reltol=reltol,fnscale=-1), z=z, w=w, group=group, H=H, verbose=verbose)
    }
    # 3. Likelihood, standard error, and p values
    res$LL = LL_ProbitRE(res$par,z,w,group,H,verbose=verbose-1)
    res$AIC = -2*res$LL + 2 * length(res$par)
    res$BIC = -2*res$LL + log(length(z)) * length(res$par)

    gvar = Gradient_ProbitRE(res$par,z,w,group,H,variance=T,verbose=verbose-1)
    res$var = gvar$var
    res$g = gvar$g
    # The precision of gtHg is very sensitive to H^-1 (var)
    res$gtHg = matrix(res$g, nrow=1) %*% res$var %*% matrix(res$g, ncol=1)
    # 3.1 Check convergence
    if(any(is.na(res$g) | !is.finite(res$g)) | res$gtHg>tol_gtHg){
        new_delta = runif(1,0,max_delta)
        writeLines(paste("====Failed to converge gtHg =", res$gtHg, ", trying delta =", new_delta,", current gradient shown below===="))
        print(res$g)
        return(ProbitRE(formula, id, data=data, delta=new_delta, max_delta=max_delta, method=method,lower=lower,upper=upper,H=H,accu=accu,reltol=reltol,verbose=verbose,tol_gtHg=tol_gtHg))
    }
    res$se = sqrt(diag(res$var))
    res$z = res$par/res$se
    res$p = 1 - pchisq(res$z^2, 1)
    res$estimates = cbind(estimates=res$par,se=res$se,z=res$z,p=res$p)
    writeLines(paste("\n *** Estimation of ProbitRE model finished, LL=",res$LL," ***", sep=""))
    print(res$estimates)
    print(Sys.time()-begin)
    # 4. Convergence criterions
    res$scgrad=tryCatch(solve(chol(gvar$I),gvar$g), error=function(e)e)
    if('numeric' %in% class(res$scgrad)[1]) {
        max_grad = round(max(abs(res$scgrad)), digits=3)
        max_grad_var = names(which.max(abs(res$scgrad)))
        if(verbose>=1)
            writeLines(paste0('Max absolute scaled gradient: ', max_grad, ' on ', paste(max_grad_var, collapse=', ')))
    } else {
        writeLines('Error occured while computing scaled gradient, details below:')
        print(res$scgrad)
    }
    writeLines(paste0('Convergence criterion gtHg: ', round(res$gtHg, digits=6)))
    rm(tmp.env)
    return (res)
}
