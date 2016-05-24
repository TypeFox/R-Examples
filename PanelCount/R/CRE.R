LL_CRE = function(par,y,z,x,w,group,rule,verbose=1){
    if(length(par) != ncol(x)+ncol(w)+3 && length(par) != ncol(x)+ncol(w)+2)
        stop("Number of parameters incorrect")
    alpha = par[1:ncol(w)]
    beta = par[ncol(w)+1:ncol(x)]
    sigma = par["sigma"]
    delta = par["delta"]
    rho = ifelse("rho" %in% names(par), par["rho"], 0)

    # Omega & v (inner), Weight & r (outer)
    Omega = rule$Omega
    v = rule$v
    Weight = rule$Weight
    r = rule$r

    N = length(group)-1
    obs = length(y)

    wa = as.vector(w %*% alpha)
    xb = as.vector(x %*% beta)
    d = 2*z - 1

    Li = rep(0, N)
    for(h in 1:length(r)){
        lamh = exp(xb + sqrt(2)*sigma*r[h])
        ix = (z==1)
        Pz = rep(1, obs)
        Pz[ix] = dpois(y[ix], lamh[ix])
        Pz_prod = as.vector(groupProd(Pz, group))
        sumk = rep(0, N)
        for(k in 1:length(v)){
            q = wa + sqrt(2)*delta*rho*r[h] + sqrt(2*(1-rho^2))*delta*v[k]
            Phi = pnorm(d*q)
            Phi_prod = as.vector(groupProd(Phi, group)) # @@@
            sumk = sumk + Omega[k]*Phi_prod
        }
        Li = Li + Weight[h] * Pz_prod * sumk
    }
    Li = pmax(Li, 1e-100) # in case that some Li=0
    LL = sum(log(Li/pi))
    if(verbose>=1){
        writeLines(paste("==== Iteration ", tmp.env$iter, ": LL=",round(LL,digits=5)," =====", sep=""))
        print(round(par,digits=3))
    }
    tmp.env$iter = tmp.env$iter + 1
    if(is.na(LL) || !is.finite(LL)){
        if(verbose>=2) writeLines("NA or infinite likelihood, will try others")
        LL = -1e300
    }
    if(tmp.env$iter==1) tmp.env$initLL = LL
    return (LL)
}


# 2. Gradient function
Gradient_CRE = function(par,y,z,x,w,group,rule,variance=F,verbose=1){
    if(length(par) != ncol(x)+ncol(w)+3 && length(par) != ncol(x)+ncol(w)+2)
        stop("Number of parameters incorrect")
    alpha = par[1:ncol(w)]
    beta = par[ncol(w)+1:ncol(x)]
    sigma = par["sigma"]
    delta = par["delta"]
    rho = ifelse("rho" %in% names(par), par["rho"], 0)

    N = length(group)-1
    obs = length(y)

    # Omega & v (inner), Weight & r (outer)
    Omega = rule$Omega
    v = rule$v
    Weight = rule$Weight
    r = rule$r

    wa = as.vector(w %*% alpha)
    xb = as.vector(x %*% beta)
    d = 2*z - 1

    w_ext = cbind(w,delta=0,rho=0) # preallocate memory for delta and rho
    x_ext = cbind(x,sigma=0)
    # renames variables for uniqueness
    colnames(w_ext) = c(paste0("w",1:ncol(w)),"delta","rho")
    colnames(x_ext) = c(paste0("x",1:ncol(x)),"sigma")

    Li = rep(0, N)
    sumH_alp = matrix(0, N, ncol(w_ext))
    sumH_beta = matrix(0, N, ncol(x_ext))
    colnames(sumH_alp) = colnames(w_ext)
    colnames(sumH_beta) = colnames(x_ext)

    # saved coef
    a1 = sqrt(2)*rho
    a2 = sqrt(2*(1-rho^2))
    a3 = sqrt(2)*delta
    a4 = sqrt(2/(1-rho^2))*rho*delta

    b1 = sqrt(2)*delta*rho
    b2 = sqrt(2*(1-rho^2))*delta

    for(h in 1:length(r)){
        lamh = exp(xb + sqrt(2)*sigma*r[h])
        ix = (z==1)
        Pz = rep(1, obs)
        Pz[ix] = dpois(y[ix], lamh[ix])
        Pz_prod = as.vector(groupProd(Pz, group))

        zyl = rep(0, obs)
        zyl[ix] = y[ix] - lamh[ix]

        # meta data for beta and sigma
        sumK_beta = rep(0, N)
        sumT_beta = matVecProdSum(x_ext, sqrt(2)*r[h], zyl, group)

        # meta data for alpha, delta, and rho
        sumK_alp = matrix(0, N, ncol(w_ext))
        for(k in 1:length(v)){
            dq = d * (b1*r[h] + b2*v[k] + wa)
            phi = dnorm(dq)
            Phi = pnorm(dq)
            Phi_prod = as.vector(groupProd(Phi, group)) # @@@
            sumK_beta = sumK_beta + Omega[k]*Phi_prod

            # delta and rho
            sumT_alp = matVecProdSum(w_ext,c(a1*r[h] + a2*v[k], a3*r[h] - a4*v[k]), d*phi/Phi, group) # @@@
            sumK_alp = sumK_alp + matVecProd(sumT_alp, Omega[k]*Phi_prod)
        }
        Li = Li + Weight[h] * Pz_prod * sumK_beta
        sumH_alp = sumH_alp + matVecProd(sumK_alp, Weight[h] * Pz_prod)
        sumH_beta = sumH_beta + matVecProd(sumT_beta, Weight[h]*Pz_prod*sumK_beta)
    }
    Li = pmax(Li, 1e-100) # in case that some Li=0
    # whether rho is included depends on length of par, pi cancells out
    dLogLi = matVecProd(cbind(sumH_alp, sumH_beta), 1/Li)
    colnames(dLogLi) = c(colnames(sumH_alp), colnames(sumH_beta))
    dLogLi = dLogLi[, c(paste0("w",1:ncol(w)),paste0("x",1:ncol(x)),names(par[(ncol(x)+ncol(w)+1):length(par)]))]
    colnames(dLogLi) = names(par)
    gradient = colSums(dLogLi)
    if(variance){
        if(verbose>=1){
            writeLines("----Converged, the gradient at optimum:")
            print(round(gradient,digits=3))
        }
        LL = sum(log(Li/pi))
        var = solve(crossprod(dLogLi))
        return (list(LL=LL, g=gradient, var=var, I=crossprod(dLogLi)))
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

# 3. Partial Effects
Partial_CRE = function(res,w,xnames,rule,intercept=F){
    wnames = colnames(w)
    par = res$estimates[, 1]
    alpha = par[1:length(wnames)]
    beta = par[length(wnames)+1:length(xnames)]
    sigma = par["sigma"]
    delta = par["delta"]
    rho = ifelse("rho" %in% names(par), par["rho"], 0)

    com = xnames[xnames %in% wnames & xnames!="(Intercept)"]
    if(intercept){ # no effect on other variables
        xunq = c("(Intercept)", xnames[!xnames %in% wnames])
        wunq = c("(Intercept)", wnames[!wnames %in% xnames])
    } else {
        xunq = xnames[!xnames %in% wnames]
        wunq = wnames[!wnames %in% xnames]
    }

    Omega = rule$Omega
    v = rule$v

    # ignoring beta in theta for now, as it's zero
    wa = as.vector(w %*% alpha)
    w_ext = cbind(w,delta=0,sigma=0,rho=0)

    obs = nrow(w)
    A = rep(0, obs)
    B = rep(0, obs)
    dA = matrix(0, obs, ncol(w_ext))
    dB = matrix(0, obs, ncol(w_ext))
    for(h in 1:length(v)){
        q = sqrt(2)*delta*v[h] + rho*sigma*delta + wa
        phi = dnorm(q)
        Phi = pnorm(q)
        A = A + Omega[h] * phi
        B = B + Omega[h] * Phi

        ext = c(sqrt(2)*v[h]+rho*sigma, rho*delta, sigma*delta)
        phi_dq = matVecProdSum(w_ext, ext, Omega[h]*phi, numeric(0))
        dA = dA - matVecProd(phi_dq, q)
        dB = dB + phi_dq
    }

    # Partial effects: average over input obs
    ratio = A/B
    Gw = mean(ratio) * alpha
    Gx = beta
    G = Gw[com]+Gx[com]
    if(length(wunq)>0) G = c(G, Gw[wunq])
    if(length(xunq)>0) G = c(G, Gx[xunq])
    names(G) = c(com, wunq, xunq)
    # Jw
    Jw_1 = mean(ratio) * cbind(diag(length(alpha)), matrix(0,length(alpha),3))
    nmr = dA - matVecProd(dB, A/B)
    frac = matVecProd(nmr, 1/B)
    Jw_2 = outer(alpha, colMeans(frac))
    Jw = Jw_1 + Jw_2
    # Now insert zeros for beta
    Jw = cbind(Jw[, 1:(ncol(Jw)-3)], matrix(0,length(alpha),length(beta)), Jw[, tail(1:ncol(Jw), 3)])
    rownames(Jw) = wnames
    Jx = cbind(matrix(0,length(beta),length(alpha)), diag(length(beta)), matrix(0,length(beta),3))
    rownames(Jx) = xnames

    J = Jx[com, ] + Jw[com, ]
    if(length(wunq)>0) J = rbind(J, Jw[wunq,])
    if(length(xunq)>0) J = rbind(J, Jx[xunq,])
    rownames(J) = c(com, wunq, xunq)
    colnames(J) = c(wnames,xnames,tail(colnames(w_ext), 3))

    # rearrange to allign with covariance matrix
    J = cbind(J[,1:(ncol(J)-3)], J[, tail(names(par), 3)])
    se = sqrt(diag(J %*% res$var %*% t(J)))
    z = G/se
    p = 1 - pchisq(z^2, 1)
    pe = cbind(estimates=G,se=se,z=z,p=p)
    return (pe)
}

#' A Model with Correlated Random Effects in Poisson and Probit Equations
#' @description Estimate a model in panel counting data, in which the selection equation is a Probit model with random effects on individuals, and the outcome equation is a Poisson model with random effects on the same individuals. The random effects on the same individual are correlated across two equations.
#' @param sel_form Formula for selection equation, a probit model with random effects
#' @param out_form Formula for outcome equation, a Poisson model with random effects
#' @param data Input data, a data frame
#' @param id A vector that represents the identity of individuals, numeric or character
#' @param par Starting values for estimates
#' @param par_files Loading initial values from saved ProbitRE and PoissonRE estimates
#' @param delta Variance of random effects in Probit model
#' @param sigma Variance of random effects in Poisson model
#' @param max_delta Largest allowed initial delta
#' @param max_sigma Largest allowed initial sigma
#' @param rho Correlation between random effects in Probit and Poisson models
#' @param lower Lower bound for estiamtes
#' @param upper Upper bound for estimates
#' @param method Searching algorithm, don't change default unless you know what you are doing
#' @param H A vector of length 2, specifying the number of points for inner and outer Quadratures
#' @param psnH Number of Quadrature points for Poisson RE model
#' @param prbH Number of Quddrature points for Probit RE model
#' @param accu L-BFGS-B only, 1e12 for low accuracy; 1e7 for moderate accuracy; 10.0 for extremely high accuracy. See optim
#' @param reltol Relative convergence tolerance. default typically 1e-8
#' @param verbose Level of output during estimation. Lowest is 0.
#' @param tol_gtHg tolerance on gtHg, not informative for L-BFGS-B
#' @return A list containing the results of the estimated model
#' @examples
#' \donttest{
#' data(rt)
#' # Note: estimation may take 2~3 minutes
#' est = CRE(isRetweet~fans+tweets+as.factor(tweet.id),
#'                    num.words~fans+tweets+as.factor(tweet.id),
#'                    id=rt$user.id, data=rt)
#' }
#' @export
#' @family PanelCount
#' @references 1. Jing Peng and Christophe Van den Bulte. Participation vs. Effectiveness of Paid Endorsers in Social Advertising Campaigns: A Field Experiment. Working Paper.
#' @references 2. Jing Peng and Christophe Van den Bulte. How to Better Target and Incent Paid Endorsers in Social Advertising Campaigns: A Field Experiment. In Proceedings of the 2015 International Conference on Information Systems.
CRE = function(sel_form, out_form, id, data=NULL, par=NULL, par_files=NULL, delta=1, max_delta=3, sigma=1, max_sigma=3, rho=0, lower=c(rho=-1), upper=c(rho=1), method='L-BFGS-B',H=c(10,10),psnH=20,prbH=20,accu=1e4,reltol=1e-8,verbose=0,tol_gtHg=Inf){
    # 1.1 Sort data based on id
    ord = order(id)
    data = data[ord,]
    id = id[ord]
    group = c(0,cumsum(table(as.integer(factor(id)))))
    # 1.2 Quadrature rules
    rule1 = gauss.quad(H[1], "hermite")
    rule2 = gauss.quad(H[2], "hermite")
    rule = list(Weight=rule1$weights, r=rule1$nodes, Omega=rule2$weights, v=rule2$nodes)
    # 1.3 parse w,z: need to make sure if data is reordered
    mf = model.frame(sel_form, data=data, na.action=NULL, drop.unused.levels=T)
    z = model.response(mf, "numeric")
    w = model.matrix(attr(mf, "terms"), data=mf)
    # 1.4 parse x,y
    mf2 = model.frame(out_form, data=data, na.action=NULL, drop.unused.levels=T)
    y = model.response(mf2, "numeric")
    x = model.matrix(attr(mf2, "terms"), data=mf2)
    if(length(id)!=nrow(w) || length(id)!=nrow(x)) stop("Error: id and w or x length don't match, potentially due to removal of data points with missing values!")
    # 1.5 Initialize parameters
    if(is.null(par)){
        if(length(par_files)==2){
            writeLines("========Loading initial parameters from estimated models=======")
            load(par_files$sel, .GlobalEnv)
            sel_est = res$estimates[,1]
            load(par_files$out, .GlobalEnv)
            out_est = res$estimates[,1]
            par = c(sel_est[-length(sel_est)], out_est, sel_est[length(sel_est)], rho=rho)
        } else {
            writeLines("========Initializing outcome equation parameters===========")
            psn_re = PoissonRE(out_form, id=id[!is.na(y)], data=data[!is.na(y),], sigma=sigma, max_sigma=max_sigma, H=psnH, method='BFGS',reltol=reltol,verbose=verbose-1)
            writeLines("========Initializing selection equation parameters=========")
            probit = ProbitRE(sel_form, id=id, data=data, delta=delta, max_delta=max_delta, method='BFGS',H=prbH,reltol=reltol,verbose=verbose-1)
            sel_est = probit$estimates[,1]
            out_est = psn_re$estimates[,1]
            par = c(sel_est[-length(sel_est)], out_est, sel_est[length(sel_est)], rho=rho)
        }
    }
    # 2. Estimation
    tmp.env <<- new.env()
    tmp.env$iter = 1
    tmp.env$initLL = -Inf
    tmp.env$begin = Sys.time()
    if(method=="L-BFGS-B"){
        lb = rep(-Inf, length(par))
        names(lb) = names(par)
        ub = rep(Inf, length(par))
        names(ub) = names(par)
        lb[match(names(lower), names(lb))] = lower
        ub[match(names(upper), names(ub))] = upper
        if(verbose>=2){
            writeLines("Lower and upper bounds, and initial values:")
            print(lb)
            print(ub)
            print(par)
        }
        res = optim(par=par, fn=LL_CRE, gr=Gradient_CRE, method="L-BFGS-B", control=list(factr=accu,fnscale=-1), lower=lb, upper=ub, y=y, z=z, x=x, w=w, group=group,rule=rule,verbose=verbose)
    } else {
        res = optim(par=par, fn=LL_CRE, gr=Gradient_CRE, method="BFGS", control=list(reltol=reltol,fnscale=-1), y=y, z=z, x=x, w=w, group=group,rule=rule,verbose=verbose)
    }
    # 3. Likelihood, standard error, and p values
    gvar = Gradient_CRE(res$par,y,z,x,w,group,rule,variance=T,verbose=verbose-1)
    res$LL = gvar$LL
    res$AIC = -2*res$LL + 2 * length(res$par)
    res$BIC = -2*res$LL + log(length(z)) * length(res$par)
    res$var = gvar$var
    res$g = gvar$g
    res$gtHg = matrix(res$g, nrow=1) %*% res$var %*% matrix(res$g, ncol=1)
    # 3.1 Check convergence
    if(any(is.na(res$g) | !is.finite(res$g)) | res$gtHg>tol_gtHg){
        par['rho'] = runif(1,-1,1)
        writeLines(paste("====Failed to converge gtHg =", res$gtHg, ", trying rho =", par['rho'], ", current gradient shown below===="))
        print(res$g)
        return (CRE(sel_form, out_form, id, data=data, par=par, par_files=par_files, delta=delta, max_delta=max_delta, sigma=sigma, max_sigma=max_sigma, lower=lower, upper=upper, method=method,H=H,psnH=psnH,prbH=prbH,accu=accu,reltol=reltol,verbose=verbose,tol_gtHg=tol_gtHg))
    }
    res$se = sqrt(diag(res$var))
    res$z = res$par/res$se
    res$p = 1 - pchisq(res$z^2, 1)

    # 3.2 Force hetergeneity terms to positive
    if('delta' %in% names(res$par) && res$par['delta']<0){
        res$par['delta'] = -res$par['delta']
        if('rho' %in% names(res$par)) res$par['rho'] = -res$par['rho']
    }
    if('sigma' %in% names(res$par) && res$par['sigma']<0){
        res$par['sigma'] = -res$par['sigma']
        if('rho' %in% names(res$par)) res$par['rho'] = -res$par['rho']
    }
    res$estimates = cbind(estimates=res$par,se=res$se,z=res$z,p=res$p)
    # 4. Partila Effects
    res$partial = Partial_CRE(res,w,colnames(x),rule)
    wavg = t(colMeans(w)) # convert vector to matrix with a single row (names kept)
    res$partialAvgObs = Partial_CRE(res,wavg,colnames(x),rule)
    # 5. Meta data
    res$iter = tmp.env$iter
    res$input = list(sel_form=sel_form, out_form=out_form, par=par, par_files=par_files, delta=delta, max_delta=max_delta, sigma=sigma, max_sigma=max_sigma, rho=rho, lower=lower, upper=upper, method=method,H=H,psnH=psnH,prbH=prbH,accu=accu,reltol=reltol,verbose=verbose,tol_gtHg=tol_gtHg)
    writeLines(paste0("\n *** Estimation of CRE model finished, LL=",res$LL," ***"))
    print(res$estimates)
    print(Sys.time()-tmp.env$begin)
    # 6. Convergence criterions
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
