paf = function(formula, data, cov)
  {
    attr(formula, ".Environment") = attr(as.formula("aa ~ bb"), ".Environment")
    fit = coxph(formula, data=data, method='breslow')

    #########################################
    ## obtain parameters for calculating paf#
    #########################################
    response= model.response(model.frame(formula, data=data))
    time = response[, attr(response, "dimnames")[[2]][1]]
    status = response[, attr(response, "dimnames")[[2]][2]]
    beta = fit$coefficients
    invI = fit$var
    ## note centered=F
    fit2 = basehaz(fit, centered=F)
    fail = sort(unique(time[status==1]))
    n1 = length(fail)
    ## note hazard>0
    Lambda = unique(fit2$hazard[fit2$hazard>0])
    lambda = Lambda - c(0, Lambda[1:(n1-1)])
    
    data0 = data
    data0[, cov] = 0
    ## Note that model.matrix function has been used for more general model formula such as interaction
    Z = as.matrix(model.matrix(fit$formula, data)[, -1])
    Z0 = as.matrix(model.matrix(fit$formula, data0)[, -1])
    ## Z excludes missing 
    ## Note the following code needed to allow missing covariates
    Z0 = as.matrix(Z0[attr(Z0, "dimnames")[[1]]%in%attr(Z, "dimnames")[[1]], ])
    ## note w, w0 is exponentiated risk score
    w = as.vector(exp(Z%*%beta))
    w0 = as.vector(exp(Z0%*%beta))

    #########################################
    ## calculating paf estimator and other #
    #########################################
    #### obtain paf point estimate
    temp0 = w0%o%Lambda
    temp = w%o%Lambda
    s0= exp(-temp0)
    s= exp(-temp)
    mean.s0 = apply(s0, 2, mean)
    mean.s = apply(s, 2, mean)
    est = 1-(1- mean.s0)/(1- mean.s)

    #### obtain se
    ## obtain invscore1 and invscore2 first
    indicate = (outer(time, fail, ">="))
    index = apply(indicate, 1, sum)    
    w.risk = t(indicate)%*%w # w.risk: n1*1 matrix
    m = dim(Z)[2]
    em = rep(1, m)
    z.wrisk = (t(Z)*(em%*%t(w)))%*%indicate
    E = t(z.wrisk)/(w.risk%*%t(em))#E:n1*m
    Elambda = E *(lambda%*%t(em))# n1*m
    Elambda = apply(Elambda, 2, cumsum)
    wlambda = lambda/as.vector(w.risk)
    wlambda = cumsum(wlambda)

    ## invscore=0 for individual with index=0
    n = length(w)
    invIscore1 = matrix(0, nrow=n, ncol=m)
    invIscore2 = matrix(0, nrow=n, ncol=n1)
    ind = which(index!=0)
    ZZ = Z[ind, ]
    sstatus = status[ind]
    ww = w[ind]
    index = index[ind]
    ## calculate invIscore for individual wiht index!=0
    scoreU = (ZZ-E[index, ])*((sstatus==1)%*%t(em)) - (ww%*%t(em))*(ZZ*(Lambda[index]%*%t(em))-Elambda[index, ])

    invIscore1[ind, ] = scoreU%*%invI
    en1 = rep(1, n1)
    J = c(1:n1)
    invIscore2[ind, ] = ((sstatus==1)%*%t(en1))*(1/w.risk[index, ]%*%t(en1))*outer(index, J, "<=") - invIscore1[ind, ]%*%t(Elambda)
    for(j in 1:n1)
      {
        invIscore2[ind, j] = invIscore2[ind, j] -ww*wlambda[pmin(index, j)]
      }

    ## then obtain se
    h10 = -t(Z0)%*%(s0 * temp0)
    h1 = -t(Z)%*%(s * temp)
    h20 = -as.vector(t(s0)%*%w0)
    h2 = -as.vector(t(s)%*%w)
    n = length(w)
    en = rep(1, n)
    item = s0 - en%*%t(mean.s0) + invIscore2*(en%*%t(h20)) + invIscore1%*%h10 - (en%*%t((1-mean.s0)/(1-mean.s)))*(s - en%*%t(mean.s) + invIscore2*(en%*%t(h2)) + invIscore1%*%h1)
    se = sqrt(apply(item^2, 2, sum))/n/(1-mean.s)
    low = 1-(1-est)*exp(1.96*se/(1-est))
    upp = 1-(1-est)*exp(-1.96*se/(1-est))

    #########################################
    ## output results #
    #########################################
    result = list(time= fail, est = est, se=se, low=low, upp=upp)
    result$fit.cox=fit
    class(result)='paf'
    result
  }


plot.paf = function(x, conf.int=TRUE, lty=1, col=1, ylim=NULL, xlab='Time', ylab='Attributable Fraction Function', ...)
  {
    if(is.null(ylim))
      {
        ylim=c(min(x$low), max(x$upp))
      }
    plot(x$time, x$est, type='s', lty=lty, col=col, ylim=ylim, xlab=xlab, ylab=ylab, ...)
    if(conf.int)
      {
        lines(x$time, x$low, lty=2, col=col, type='s')
        lines(x$time, x$upp, lty=2, col=col, type ='s')
      }
  }

