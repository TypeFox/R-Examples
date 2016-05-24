LL_PLN_RE = function(par,y,x,group,H,verbose=1){
    if(length(par) != ncol(x)+2)  stop("Number of parameters incorrect")
    beta = par[1:ncol(x)]
    sigma = par[length(par)-1]
    gamma = par[length(par)]

    N = length(group)-1
    obs = length(y)

    rule = gauss.quad(H, "hermite")
    Omega = rule$weights
    v = rule$nodes

    xb = as.vector(x %*% beta)

    Li = rep(0, N)
    for(h in 1:H){
        sumK = rep(0, obs)
        for(k in 1:H){
            lamh = exp(xb + sqrt(2)*sigma*v[h] + sqrt(2)*gamma*v[k])
            PhiP = Omega[k] / sqrt(pi) * dpois(y, lamh)
            sumK = sumK + PhiP
        }
        prodT = groupProd(sumK, group)
        Li = Li + Omega[h] * prodT
    }
    Li = pmax(Li, 1e-100) # in case that some Li=0
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
Gradient_PLN_RE = function(par,y,x,group,H,variance=F,verbose=1){
    if(length(par) != ncol(x)+2)  stop("Number of parameters incorrect")
    beta = par[1:ncol(x)]
    sigma = par[length(par)-1]
    gamma = par[length(par)]

    N = length(group)-1
    obs = length(y)

    rule = gauss.quad(H, "hermite")
    Omega = rule$weights
    v = rule$nodes

    xb = as.vector(x %*% beta)
    x_ext = cbind(x,sigma=0,gamma=0)

    Li = rep(0, N)
    sumH = matrix(0, N, ncol(x_ext))
    for(h in 1:H){
        sumK = rep(0, obs)
        sumK_ext = matrix(0, obs, ncol(x_ext))
        for(k in 1:H){
            lamh = exp(xb + sqrt(2)*sigma*v[h] + sqrt(2)*gamma*v[k])
            PhiP = Omega[k] / sqrt(pi) * dpois(y, lamh)
            sumK = sumK + PhiP
            sumK_ext = sumK_ext + matVecProdSum(x_ext, c(sqrt(2)*v[h],sqrt(2)*v[k]), PhiP*(y-lamh), numeric(0))
        }
        prodT = groupProd(sumK, group)
        Li = Li + Omega[h] * prodT
        # if sumK=0, the corresponding rows in sumK_ext are likely to be zero as well
        sumT = matVecProdSum(sumK_ext, numeric(0), 1/pmax(sumK,1e-6), group)
        sumH = sumH + matVecProd(sumT, Omega[h] * prodT)
    }
    Li = pmax(Li, 1e-100) # in case that some Li=0
    dLogLi = matVecProd(sumH, 1/Li)
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

#' A Poisson Lognormal Model with Random Effects
#' @description Estimate a Poisson Lognormal model with random effects in panel counting data. This model accounts for heterogeneity on the individual level, and heterogeneity on the <individual, time> level.
#' @param formula Formula of the model
#' @param data Input data, a data frame
#' @param id A vector that represents the identity of individuals, numeric or character
#' @param method Searching algorithm, don't change default unless you know what you are doing
#' @param par Starting values for estimates
#' @param sigma Variance of random effects on the individual level for PLN_RE
#' @param gamma Variance of random effects on the <individual, time> level for PLN_RE
#' @param max_sigma Largest allowed initial sigma
#' @param max_gamma Largest allowed initial gamma
#' @param lower Lower bound for estiamtes
#' @param upper Upper bound for estimates
#' @param H A vector of length 2, specifying the number of points for inner and outer Quadratures
#' @param psnH Number of Quadrature points for Poisson RE model
#' @param accu L-BFGS-B only, 1e12 for low accuracy; 1e7 for moderate accuracy; 10.0 for extremely high accuracy. See optim
#' @param reltol Relative convergence tolerance. default typically 1e-8
#' @param verbose Level of output during estimation. Lowest is 0.
#' @param tol_gtHg tolerance on gtHg, not informative for L-BFGS-B
#' @return A list containing the results of the estimated model
#' @examples
#' \donttest{
#' data(rt)
#' est = PLN_RE(num.words~fans+tweets+as.factor(tweet.id),
#'               id=rt$user.id[rt$isRetweet==1],
#'               data=rt[rt$isRetweet==1,])
#' }
#' @export
#' @family PanelCount
#' @references 1. Jing Peng and Christophe Van den Bulte. Participation vs. Effectiveness of Paid Endorsers in Social Advertising Campaigns: A Field Experiment. Working Paper.
#' @references 2. Jing Peng and Christophe Van den Bulte. How to Better Target and Incent Paid Endorsers in Social Advertising Campaigns: A Field Experiment. In Proceedings of the 2015 International Conference on Information Systems.
PLN_RE = function(formula, id, data=NULL, par=NULL, gamma=1, max_gamma=5, sigma=1, max_sigma=3, method='BFGS',lower=NULL,upper=NULL,H=20, psnH=20,accu=10,reltol=1e-8,verbose=0,tol_gtHg=Inf){
    # 1.1 Sort data based on id
    ord = order(id)
    data = data[ord,]
    id = id[ord]
    group = c(0,cumsum(table(as.integer(factor(id)))))
    # 1.2 Parse x and y
    mf = model.frame(formula, data=data, na.action=NULL, drop.unused.levels=T)
    y = model.response(mf, "numeric")
    x = model.matrix(attr(mf, "terms"), data=mf)
    if(length(id)!=nrow(x)) stop("Error: id and x length don't match, potentially due to removal of data points with missing values!")
    # 1.3 Initial values
    if(is.null(par)){
        # use Poisson estimates as initial values for parameters
        writeLines('========Initializing estimates with Possion RE model===========')
        pre = PoissonRE(formula, id=id, data=data, sigma=sigma, max_sigma=max_sigma, method=method,lower=lower,upper=upper,H=psnH,accu=accu,reltol=reltol,verbose=verbose-1)
        par = c(pre$estimates[,1], gamma=gamma)
    }
    # 2. Estimation
    tmp.env <<- new.env()
    tmp.env$iter = 1
    begin = Sys.time()
    if(method=="L-BFGS-B"){
        res = optim(par=par, fn=LL_PLN_RE, gr=Gradient_PLN_RE, method="L-BFGS-B", control=list(factr=accu,fnscale=-1), lower=lower, upper=upper, y=y, x=x, group=group, H=H, verbose=verbose)
    } else {
        res = optim(par=par, fn=LL_PLN_RE, gr=Gradient_PLN_RE, method="BFGS", control=list(reltol=reltol,fnscale=-1), y=y, x=x, group=group, H=H, verbose=verbose)
    }
    # 3. Likelihood, standard error, and p values
    res$LL = LL_PLN_RE(res$par,y,x,group,H,verbose=verbose-1)
    res$AIC = -2*res$LL + 2 * length(res$par)
    res$BIC = -2*res$LL + log(length(y)) * length(res$par)
    gvar = Gradient_PLN_RE(res$par,y,x,group,H,variance=T,verbose=verbose-1)
    res$var = gvar$var
    res$g = gvar$g
    res$gtHg = matrix(res$g, nrow=1) %*% res$var %*% matrix(res$g, ncol=1)
    # 3.1 Check convergence
    if(any(is.na(res$g) | !is.finite(res$g)) | res$gtHg>tol_gtHg){
        par['gamma'] = runif(1,0,max_gamma)
        writeLines(paste("====Failed to converge gtHg =", res$gtHg, ", trying gamma =", par['gamma'],", current gradient shown below===="))
        print(res$g)
        return (PLN_RE(formula, id, data=data, par=par, max_gamma=max_gamma, sigma=sigma, max_sigma=max_sigma, method=method,lower=lower,upper=upper,H=H, psnH=psnH,accu=accu,reltol=reltol,verbose=verbose,tol_gtHg=tol_gtHg))
    }
    res$se = sqrt(diag(res$var))
    res$z = res$par/res$se
    res$p = 1 - pchisq(res$z^2, 1)
    res$estimates = cbind(estimates=res$par,se=res$se,z=res$z,p=res$p)
    writeLines(paste("\n *** Estimation of PLN_RE model finished, LL=",res$LL," ***", sep=""))
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