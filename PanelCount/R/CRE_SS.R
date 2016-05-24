LL_CRE_SS = function(par,y,z,x,w,group,rule,verbose=1){
    if(length(par) < ncol(x)+ncol(w)+3 || length(par) > ncol(x)+ncol(w)+5)
        stop("Number of parameters incorrect")
    alpha = par[1:ncol(w)]
    beta = par[ncol(w)+1:ncol(x)]
    # names of variables in x and w should not use sigma, ..., tau
    if(any(!c("sigma","gamma","delta") %in% names(par))){
        print(names(par))
        stop("sigma, gamma, or delta not provided")
    }
    sigma = par["sigma"]
    gamma = par["gamma"]
    delta = par["delta"]
    rho = ifelse("rho" %in% names(par), par["rho"], 0)
    tau = ifelse("tau" %in% names(par), par["tau"], 0)

    # Omega & v (inner), Weight & r (outer)
    Omega = rule$Omega
    v = rule$v
    Weight = rule$Weight
    r = rule$r

    N = length(group)-1
    obs = length(y)

    wa = as.vector(w %*% alpha)
    xb = as.vector(x %*% beta)
    d = (2*z - 1)/sqrt(1-tau^2)
    ix = (z==1)

    Li = rep(0, N)
    for(h in 1:length(r)){
        sumK = rep(0, N)
        for(k in 1:length(v)){
            sumM = rep(0, obs)
            for(m in 1:length(v)){
                lam = exp(xb[ix] + sqrt(2)*sigma*r[h] + sqrt(2)*gamma*v[m])
                Pz = rep(1, obs)
                Pz[ix] = dpois(y[ix], lam)

                q = d * (sqrt(2)*rho*delta*r[h] + sqrt(2*(1-rho^2))*delta*v[k] + sqrt(2)*tau*v[m] + wa)
                Phi = pnorm(q)
                sumM = sumM + Omega[m]/sqrt(pi) * Pz * Phi
            }
            prodT = groupProd(sumM, group)
            sumK = sumK + Omega[k] * prodT
        }
        Li = Li + Weight[h] * sumK
    }
    Li = pmax(Li, 1e-100) # in case that some Li=0
    LL = sum(log(Li/pi))
    if(verbose>=1){
        writeLines(paste0("==== Iteration ", tmp.env$iter, ": LL=",round(LL,digits=5)," ====="))
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
Gradient_CRE_SS = function(par,y,z,x,w,group,rule,variance=F,verbose=1){
    if(length(par) < ncol(x)+ncol(w)+3 || length(par) > ncol(x)+ncol(w)+5)
        stop("Number of parameters incorrect")
    alpha = par[1:ncol(w)]
    beta = par[ncol(w)+1:ncol(x)]
    # names of variables in x and w should not use sigma, ..., tau
    sigma = par["sigma"]
    gamma = par["gamma"]
    delta = par["delta"]
    rho = ifelse("rho" %in% names(par), par["rho"], 0)
    tau = ifelse("tau" %in% names(par), par["tau"], 0)

    # Omega & v (inner), Weight & r (outer)
    Omega = rule$Omega
    v = rule$v
    Weight = rule$Weight
    r = rule$r

    N = length(group)-1
    obs = length(y)

    wa = as.vector(w %*% alpha)
    xb = as.vector(x %*% beta)
    d = (2*z - 1)/sqrt(1-tau^2)
    ix = (z==1)
    watau = wa * tau / (1-tau^2)

    w_ext = cbind(w,delta=0,rho=0,tau=0) # preallocate memory for delta and rho
    x_ext = cbind(x,sigma=0,gamma=0)
    # renames variables for uniqueness
    colnames(w_ext) = c(paste0("w",1:ncol(w)),"delta","rho","tau")
    colnames(x_ext) = c(paste0("x",1:ncol(x)),"sigma","gamma")

    Li = rep(0, N)
    sumH_alp = matrix(0, N, ncol(w_ext))
    sumH_beta = matrix(0, N, ncol(x_ext))
    colnames(sumH_alp) = colnames(w_ext)
    colnames(sumH_beta) = colnames(x_ext)
    for(h in 1:length(r)){
        sumK = rep(0, N)
        sumK_alp = matrix(0, N, ncol(w_ext))
        sumK_beta = matrix(0, N, ncol(x_ext))
        for(k in 1:length(v)){
            sumM = rep(0, obs)
            sumM_alp = matrix(0, obs, ncol(w_ext))
            sumM_beta = matrix(0, obs, ncol(x_ext))
            for(m in 1:length(v)){
                lam = exp(xb[ix] + sqrt(2)*sigma*r[h] + sqrt(2)*gamma*v[m])
                Pz = rep(1, obs)
                Pz[ix] = dpois(y[ix], lam)
                zyl = rep(0, obs)
                zyl[ix] = y[ix] - lam

                q = d * (sqrt(2)*rho*delta*r[h] + sqrt(2*(1-rho^2))*delta*v[k] + sqrt(2)*tau*v[m] + wa)
                phi = dnorm(q)
                Phi = pnorm(q)
                OPP = Omega[m]/sqrt(pi) * Pz * Phi

                sumM = sumM + OPP
                vecW = c(sqrt(2)*rho*r[h]+sqrt(2*(1-rho^2))*v[k], sqrt(2)*delta*r[h]-sqrt(2/(1-rho^2))*rho*delta*v[k])
                vecW2 = (sqrt(2)*v[m] + sqrt(2)*rho*delta*tau*r[h] + sqrt(2*(1-rho^2))*delta*tau*v[k])/(1-tau^2) + watau
                vecX = c(sqrt(2)*r[h], sqrt(2)*v[m])
                sumM_alp = sumM_alp + matVecProdSumExt(w_ext, vecW, vecW2, Omega[m]/sqrt(pi) * Pz * phi * d, numeric(0))
                sumM_beta = sumM_beta + matVecProdSum(x_ext, vecX, zyl * OPP, numeric(0))
            }
            prodT = groupProd(sumM, group)
            sumK = sumK + Omega[k] * prodT
            # If sumM=0, the corresponding rows in sumM_alp/sumM_beta are likely to be zero
            sumT_alp = matVecProdSum(sumM_alp, numeric(0), 1/pmax(sumM,1e-6), group)
            sumT_beta = matVecProdSum(sumM_beta, numeric(0), 1/pmax(sumM,1e-6), group)
            sumK_alp = sumK_alp + matVecProd(sumT_alp, Omega[k]*prodT)
            sumK_beta = sumK_beta + matVecProd(sumT_beta, Omega[k]*prodT)
        }
        Li = Li + Weight[h] * sumK
        sumH_alp = sumH_alp + Weight[h] * sumK_alp
        sumH_beta = sumH_beta + Weight[h] * sumK_beta
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
Partial_CRE_SS = function(res,w,xnames,rule,intercept=F){
    wnames = colnames(w)
    par = res$estimates[, 1]
    alpha = par[1:length(wnames)]
    beta = par[length(wnames)+1:length(xnames)]
    delta = par["delta"]
    sigma = par["sigma"]
    gamma = par["gamma"]
    rho = ifelse("rho" %in% names(par), par["rho"], 0)
    tau = ifelse("tau" %in% names(par), par["tau"], 0)

    com = xnames[xnames %in% wnames & xnames!="(Intercept)"]
    if(intercept){ # no effect on other variables
        xunq = c("(Intercept)", xnames[!xnames %in% wnames])
        wunq = c("(Intercept)", wnames[!wnames %in% xnames])
    } else {
        xunq = xnames[!xnames %in% wnames]
        wunq = wnames[!wnames %in% xnames]
    }

    # Omega & v (inner), Weight & r (outer)
    Omega = rule$Omega
    v = rule$v
    Weight = rule$Weight
    r = rule$r

    # Relative partial effects
    oneMT = 1-tau^2
    sqrtOneMT = sqrt(oneMT)
    wa = as.vector(w %*% alpha)
    qwa = wa / sqrtOneMT
    dqw = w /sqrtOneMT
    dqwa = 2*tau / (oneMT * sqrtOneMT) * wa
    w_ext = cbind(dqw,delta=0,sigma=0,gamma=0,rho=0,tau=0)

    obs = nrow(w)
    A = rep(0, obs)
    B = rep(0, obs)
    dA = matrix(0, obs, length(alpha)+5)
    dB = matrix(0, obs, length(alpha)+5)
    for(h in 1:length(r)){
        Ah = rep(0, obs)
        Bh = rep(0, obs)
        dAh = matrix(0, obs, length(alpha)+5)
        dBh = matrix(0, obs, length(alpha)+5)
        for(k in 1:length(v)){
            q = (sqrt(2)*delta*r[h] + sqrt(2)*tau*v[k] + rho*sigma*delta + tau*gamma) / sqrtOneMT + qwa
            phi = dnorm(q)
            Phi = pnorm(q)
            Ah = Ah + Omega[k] * phi
            Bh = Bh + Omega[k] * Phi

            # for variance
            ext1 = c(sqrt(2)*r[h]+rho*sigma, rho*delta, tau, sigma*delta) / sqrtOneMT
            ext2 = ( (gamma+sqrt(2)*v[k])*(1+tau*tau) + 2*tau*(sqrt(2)*delta*r[h] + rho*sigma*delta) ) / (oneMT * sqrtOneMT) + dqwa
            phi_dq = matVecProdSumExt(w_ext, ext1, ext2, phi * Omega[k], numeric(0))
            dAh = dAh - matVecProd(phi_dq, q)
            dBh = dBh + phi_dq
        }
        A = A + Weight[h] * Ah
        B = B + Weight[h] * Bh
        dA = dA + Weight[h] * dAh
        dB = dB + Weight[h] * dBh
    }
    # Partial effects
    ratio = A / (sqrtOneMT * B)
    Gw = mean(ratio) * alpha
    Gx = beta
    G = Gw[com]+Gx[com]
    if(length(wunq)>0) G = c(G, Gw[wunq])
    if(length(xunq)>0) G = c(G, Gx[xunq])
    names(G) = c(com, wunq, xunq)
    # Jw
    Jw_1 = mean(ratio) * cbind(diag(length(alpha)), matrix(0,length(alpha),4), 2*tau/oneMT*alpha)
    nmr = dA - matVecProd(dB, A/B)
    frac = matVecProd(nmr, 1/(sqrtOneMT*B))
    Jw_2 = outer(alpha, colMeans(frac))
    Jw = Jw_1 + Jw_2
    # Now insert zeros for beta
    Jw = cbind(Jw[, 1:(ncol(Jw)-5)], matrix(0,length(alpha),length(beta)), Jw[, tail(1:ncol(Jw), 5)])
    rownames(Jw) = wnames
    Jx = cbind(matrix(0,length(beta),length(alpha)), diag(length(beta)), matrix(0,length(beta),5))
    rownames(Jx) = xnames

    J = Jx[com, ] + Jw[com, ]
    if(length(wunq)>0) J = rbind(J, Jw[wunq,])
    if(length(xunq)>0) J = rbind(J, Jx[xunq,])
    rownames(J) = c(com, wunq, xunq)
    colnames(J) = c(wnames,xnames,"delta","sigma","gamma","rho","tau")

    # rearrange to allign with covariance matrix, remove missing structural parameters
    plain_par = length(wnames) + length(xnames)
    J = cbind(J[,1:plain_par], J[, tail(names(par), length(par)-plain_par)])
    se = sqrt(diag(J %*% res$var %*% t(J)))
    z = G/se
    p = 1 - pchisq(z^2, 1)
    pe = cbind(estimates=G,se=se,z=z,p=p)
    return (pe)
}

# Prediction
predict_CRE_SS = function(par,y,z,x,w,group,rule){
    if(length(par) < ncol(x)+ncol(w)+3 || length(par) > ncol(x)+ncol(w)+5)
        stop("Number of parameters incorrect")
    alpha = par[1:ncol(w)]
    beta = par[ncol(w)+1:ncol(x)]
    # names of variables in x and w should not use sigma, ..., tau
    if(any(!c("sigma","gamma","delta") %in% names(par))){
        print(names(par))
        stop("sigma, gamma, or delta not provided")
    }
    sigma = par["sigma"]
    gamma = par["gamma"]
    delta = par["delta"]
    rho = ifelse("rho" %in% names(par), par["rho"], 0)
    tau = ifelse("tau" %in% names(par), par["tau"], 0)

    # Omega & v (inner), Weight & r (outer)
    Omega = rule$Omega
    v = rule$v
    Weight = rule$Weight
    r = rule$r

    N = length(group)-1
    obs = length(y)

    wa = as.vector(w %*% alpha)
    xb = as.vector(x %*% beta)

    prob = rep(0, obs)
    outcome = rep(0, obs)
    for(h in 1:length(r)){
        prob = prob + Weight[h] * pnorm(wa + sqrt(2)*delta*r[h])
        outcome_h = rep(0, obs)
        for(k in 1:length(v)){
            outcome_h = outcome_h + Omega[k] * exp(sqrt(2)*sigma*r[h]+sqrt(2)*gamma*v[k]+xb)
        }
        outcome = outcome + Weight[h] * outcome_h
    }
    return (list(prob=prob/sqrt(pi), outcome=outcome/pi))
}

#' A Sample Selection Model with Correlated Random Effects
#' @description Estimate a sample selection model in panel counting data, in which the selection equation is a Probit model with random effects on individuals, and the outcome equation is a Poisson Lognormal model with random effects on the same individuals. The random effects on the same individual and the error terms on the same <individual, time> dyad are both correlated across two equations.
#' @param sel_form Formula for selection equation, a probit model with random effects
#' @param out_form Formula for outcome equation, a Poisson model with random effects
#' @param data Input data, a data frame
#' @param id A vector that represents the identity of individuals, numeric or character
#' @param par Starting values for estimates
#' @param killed_par correlation parameters to swtich off
#' @param par_files Loading initial values from saved ProbitRE and PoissonRE estimates
#' @param delta Variance of random effects on the individual level for ProbitRE
#' @param sigma Variance of random effects on the individual level for PLN_RE
#' @param gamma Variance of random effects on the <individual, time> level for PLN_RE
#' @param rho Correlation between random effects on the individual level
#' @param max_delta Largest allowed initial delta
#' @param max_sigma Largest allowed initial sigma
#' @param max_gamma Largest allowed initial gamma
#' @param tau Correlation between error terms on the <individual, time> level
#' @param lower Lower bound for estiamtes
#' @param upper Upper bound for estimates
#' @param method Searching algorithm, don't change default unless you know what you are doing
#' @param H A vector of length 2, specifying the number of points for inner and outer Quadratures
#' @param psnH Number of Quadrature points for Poisson RE model
#' @param prbH Number of Quddrature points for Probit RE model
#' @param plnreH Number of Quddrature points for PLN_RE model
#' @param accu L-BFGS-B only, 1e12 for low accuracy; 1e7 for moderate accuracy; 10.0 for extremely high accuracy. See optim
#' @param reltol Relative convergence tolerance. default typically 1e-8
#' @param verbose Level of output during estimation. Lowest is 0.
#' @param tol_gtHg tolerance on gtHg, not informative for L-BFGS-B
#' @return A list containing the results of the estimated model
#' @examples
#' \donttest{
#' data(rt)
#' # Note: estimation may take up 10~15 minutes
#' est = CRE_SS(isRetweet~fans+tweets+as.factor(tweet.id),
#'                        num.words~fans+tweets+as.factor(tweet.id),
#'                        id=rt$user.id, data=rt)
#' }
#' @export
#' @family PanelCount
#' @references 1. Jing Peng and Christophe Van den Bulte. Participation vs. Effectiveness of Paid Endorsers in Social Advertising Campaigns: A Field Experiment. Working Paper.
#' @references 2. Jing Peng and Christophe Van den Bulte. How to Better Target and Incent Paid Endorsers in Social Advertising Campaigns: A Field Experiment. In Proceedings of the 2015 International Conference on Information Systems.
CRE_SS = function(sel_form, out_form, id, data=NULL, par=NULL, killed_par=NULL, par_files=NULL, delta=1, sigma=1, gamma=1,max_delta=3,max_sigma=3,max_gamma=5, rho=0, tau=0, lower=c(rho=-1, tau=-1), upper=c(rho=1, tau=1), method='L-BFGS-B',H=c(10,10),psnH=20,prbH=20,plnreH=20,accu=1e4,reltol=sqrt(.Machine$double.eps),verbose=0,tol_gtHg=Inf){
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
            writeLines("======Loading initial parameters from estimated models========")
            load(par_files$sel, .GlobalEnv)
            sel_est = res$estimates[,1]
            load(par_files$out, .GlobalEnv)
            out_est = res$estimates[,1]
            par = c(sel_est[-length(sel_est)], out_est, sel_est[length(sel_est)], rho=rho,tau=tau)
        } else {
            writeLines("=========Initializing outcome equation parameters===========")
            pln_re = PLN_RE(out_form, id=id[!is.na(y)], data=data[!is.na(y),], gamma=gamma, max_gamma=max_gamma, sigma=sigma, max_sigma=max_sigma, psnH=psnH, method='BFGS',H=plnreH,reltol=reltol,verbose=verbose-1)
            writeLines("=========Initializing selection equation parameters=========")
            probit = ProbitRE(sel_form, id=id, data=data, delta=delta, max_delta=max_delta, method='BFGS',H=prbH,reltol=reltol,verbose=verbose-1)
            sel_est = probit$estimates[,1]
            out_est = pln_re$estimates[,1]
            par = c(sel_est[-length(sel_est)], out_est, sel_est[length(sel_est)], rho=rho,tau=tau)
        }
    }
    par = par[!names(par) %in% killed_par]
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
        lower = lower[names(lower) %in% names(par)]
        upper = upper[names(upper) %in% names(par)]
        lb[match(names(lower), names(lb))] = lower
        ub[match(names(upper), names(ub))] = upper
        if(verbose>=2){
            writeLines("Lower and upper bounds, and initial values:")
            print(lb)
            print(ub)
            print(par)
        }
        res = optim(par=par, fn=LL_CRE_SS, gr=Gradient_CRE_SS, method="L-BFGS-B", control=list(factr=accu,fnscale=-1), lower=lb, upper=ub, y=y, z=z, x=x, w=w, group=group,rule=rule,verbose=verbose)
    } else {
        res = optim(par=par, fn=LL_CRE_SS, gr=Gradient_CRE_SS, method="BFGS", control=list(reltol=reltol,fnscale=-1), y=y, z=z, x=x, w=w, group=group,rule=rule,verbose=verbose)
    }

    # 3. Likelihood, standard error, and p values
    gvar = Gradient_CRE_SS(res$par,y,z,x,w,group,rule,T,verbose=verbose-1)
    res$LL = gvar$LL
    res$AIC = -2*res$LL + 2 * length(res$par)
    res$BIC = -2*res$LL + log(length(z)) * length(res$par)
    res$var = gvar$var
    res$g = gvar$g
    res$gtHg = matrix(res$g, nrow=1) %*% res$var %*% matrix(res$g, ncol=1)
    # 3.1 Check convergence
    if(any(is.na(res$g) | !is.finite(res$g)) | res$gtHg>tol_gtHg){
        if('rho' %in% names(par)) par['rho'] = runif(1,-1,1)
        if('tau' %in% names(par)) par['tau'] = runif(1,-1,1)
        writeLines(paste("====Failed to converge gtHg =", res$gtHg, ", trying rho =", par['rho'], ", tau =", par['tau'],", current gradient shown below===="))
        print(res$g)
        return (CRE_SS(sel_form, out_form, id, data=data, par=par, par_files=par_files, delta=delta, sigma=sigma, gamma=gamma,max_delta=max_delta,max_sigma=max_sigma,max_gamma=max_gamma, lower=lower, upper=upper, method=method,H=H,psnH=psnH,prbH=prbH,plnreH=plnreH,accu=accu,reltol=reltol,verbose=verbose,tol_gtHg=tol_gtHg))
    }
    res$se = sqrt(diag(res$var))
    res$z = res$par/res$se
    res$p = 1 - pchisq(res$z^2, 1)

    # 3.2 Force hetergeneity terms to positive
    if('gamma' %in% names(res$par) && res$par['gamma']<0){
        res$par['gamma'] = -res$par['gamma']
        if('tau' %in% names(res$par)) res$par['tau'] = -res$par['tau']
    }
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
    res$partial = Partial_CRE_SS(res,w,colnames(x),rule)
    wavg = t(colMeans(w)) # convert vector to matrix with a single row (names kept)
    res$partialAvgObs = Partial_CRE_SS(res,wavg,colnames(x),rule)
    writeLines("---Now working on predictions---")
    res$predict = predict_CRE_SS(res$par,y,z,x,w,group,rule)
    # 5. Meta data
    res$iter = tmp.env$iter
    res$input = list(sel_form=sel_form, out_form=out_form, par=par, killed_par=killed_par, par_files=par_files, delta=delta, sigma=sigma, gamma=gamma,max_delta=max_delta,max_sigma=max_sigma,max_gamma=max_gamma, rho=rho, tau=tau, lower=lower, upper=upper, method=method,H=H,psnH=psnH,prbH=prbH,plnreH=plnreH,accu=accu,reltol=reltol,verbose=verbose,tol_gtHg=tol_gtHg)
    writeLines(paste0("\n *** Estimation of CRE_SS model finished, LL=",res$LL," ***"))
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
