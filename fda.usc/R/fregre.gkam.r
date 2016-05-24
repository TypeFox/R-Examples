fregre.gkam=function (formula,family = gaussian(),data, weights= rep(1,nobs),
     par.metric = NULL,par.np=NULL,offset=NULL,
     control = list(maxit = 100,epsilon = 0.001, trace = FALSE, inverse="solve"),...)
{
 tf <- terms.formula(formula)
 terms <- attr(tf, "term.labels")
 nt <- length(terms)
 if (attr(tf, "response") > 0) {
        response <- as.character(attr(tf, "variables")[2])
        pf <- rf <- paste(response, "~", sep = "")
    } else pf <- rf <- "~"
 vtab<-rownames(attr(tf,"factors"))
 vnf=intersect(terms,names(data$df))
 vnf2=intersect(vtab[-1],names(data$df)[-1])
 vfunc2=setdiff(terms,vnf)
 vint=setdiff(terms,vtab)
 vfunc=setdiff(vfunc2,vint)
 vnf=c(vnf2,vint)
 off<-attr(tf,"offset")
 name.coef=nam=par.fregre=beta.l=list()
 kterms=1
 if (attr(tf,"intercept")==0) intercept=FALSE
 else intercept=TRUE
 if (length(vnf)>0) {
 XX=data[[1]][,c(response,vnf2)] #data.frame el 1er elemento de la lista
 for ( i in 1:length(vnf)){
   print(paste("Non functional covariate:",vnf[i]))
   print(paste("The procedure considers only functional covariates and therefore the variable",vnf[i]," is not used."))
     if (kterms > 1)   pf <- paste(pf, "+", vnf[i], sep = "")
     else pf <- paste(pf, vnf[i], sep = "")
     kterms <- kterms + 1
     }
if   (!intercept) {
     pf<- paste(pf,-1,sep="")
     }
}
else {
 XX=data.frame(data[[1]][,response])
 names(XX)=response
}
if (is.null(control$maxit))  control$maxit<-100
if (is.null(control$epsilon))  control$epsilon= 0.001
if (is.null(control$trace))  control$trace = FALSE
if (is.null(control$inverse))  control$inverse = "solve"
############################################################
    xlist<-data[-1] #datos funcionales menos el df!
    y0<-y<-data[["df"]][,response]
if (family$family=="binomial") {
   y<-as.numeric(factor(y,labels=c(0,1)))-1
}
####################    eps <- 0.001
#    if (is.null(control$epsilon)) control$epsilon<-0.00
    eps<-control$epsilon
    namesx<- vfunc
    nvars <- length(vfunc)
    ynames <- if (is.matrix(y))  rownames(y)
              else names(y)
    conv <- FALSE
    nobs <- if (is.matrix(y))   nrow(y)
            else length(y)
    if (is.null(offset))        offset <- rep(0,nobs) ##
    EMPTY <- nvars == 0
#    if (is.null(weights))        weights <- rep(1,nobs)
    variance <- family$variance
    linkinv <- family$linkinv
    linkfun <- family$linkfun
    if (!is.function(variance) || !is.function(linkinv))
        stop("'family' argument seems not to be a valid family object",
            call. = FALSE)
    dev.resids <- family$dev.resids
    aic <- family$aic
    mu.eta <- family$mu.eta
    unless.null <- function(x, if.null) {
        if (is.null(x))             if.null
        else x
    }
    valideta <- unless.null(family$valideta, function(eta) TRUE)
    validmu <- unless.null(family$validmu, function(mu) TRUE)
    X = matrix(0, nrow = nobs, ncol = nvars + intercept)
    colnames(X) = c(namesx, "Intercept")
    eqrank = c(rep(0, nvars), 1)
#    result = vector("list", nvars)
    metric<-metric2<-result<-list()
#    metric2 = metric = vector("list", nvars)
    par.np2=par.np
#    names(eqrank)=names(metric2) = names(metric) = namesx
    names(eqrank)<- c(namesx,"Intercept")
    X[, nvars + intercept] = rep(linkfun(mean(y)), nobs)
    if (control$trace)   cat("----Computing the distance matrix ----\n")
    for (i in 1:nvars) {
#        metric2[[namesx[i]]] = metric.lp(xlist[[namesx[i]]], xlist[[namesx[i]]])
        if (is.null(par.metric)) {
            metric[[namesx[i]]] = metric.lp(xlist[[namesx[i]]], xlist[[namesx[i]]])
        }
        else {
           par.metric[[namesx[i]]]$fdata1 <- xlist[[namesx[i]]]
            metric[[namesx[i]]] = do.call(par.metric[[namesx[i]]]$metric,
                par.metric[[namesx[i]]][-1])
        }
       if (is.null(par.np)) {
          par.np2[[namesx[i]]] =list(Ker=AKer.norm,type.S="S.NW",par.S=list(w=weights))
        }
#       if (is.null(par.np[[namesx[i]]]$par.S)) par.np[[namesx[i]]]$par.S=list(w=weights)
       if (is.null(par.np[[namesx[i]]]$Ker)) par.np2[[namesx[i]]]$Ker=AKer.norm
       if (is.null(par.np[[namesx[i]]]$type.S)) par.np2[[namesx[i]]]$type.S="S.NW"
       if (is.null(par.np[[namesx[i]]]$h)) {
              par.np2[[namesx[i]]]$h = h.default(xlist[[namesx[i]]], len = 51,
              prob = c(0.01,0.66),metric =  metric[[namesx[i]]])
              }
    }
#    eta = apply(X, 1, sum)
     eta = rowSums(X)
    mu = linkinv(eta)
    
    conv <- FALSE
    for (iter in 1L:control$maxit) {
        dev <- devold <- sum(dev.resids(y, mu, weights))
        good <- weights > 0
        varmu <- variance(mu)[good]
        if (any(is.na(varmu))) {stop("NAs in V(mu)")}
        if (any(varmu == 0)) {stop("0s in V(mu)")}
        mu.eta.val <- mu.eta(eta)
        if (any(is.na(mu.eta.val[good])))   stop("NAs in d(mu)/d(eta)")
        good <- (weights > 0) & (mu.eta.val != 0)
        if (all(!good)) {
            conv <- FALSE
            warning("no observations informative at iteration ",iter)
            break
        }

        Xold = X
        ngoodobs <- as.integer(nobs - sum(abs(mu.eta.val) <= eps))
        if (ngoodobs < min(0.2 * nobs, 20)) {
            conv = TRUE
            break
        }
   #  z <- (eta - offset)[good] + (y - mu)[good]/mu.eta.val[good]        #GLM
     ytilde <- (eta - offset)[good] + (y - mu)[good]/mu.eta.val[good]
#        ytilde <- eta[good] + (y - mu)[good]/mu.eta.val[good] #ANTES
        w <- sqrt((weights[good] * mu.eta.val[good]^2)/variance(mu)[good])##################################################
        X[, nvars + intercept] = rep(mean(ytilde - apply(X[good,1:nvars, drop = FALSE], 1, sum)), nobs)
        if (control$trace)
            cat("#------------------------------------------------\n")
        if (control$trace)
            cat("Iter:", iter, ngoodobs, "/", nobs, "alpha:",
                X[1, nvars + intercept], "Dev:", dev, "\n")
        for (i in 1:nvars) {
            off = apply(X[good, -i, drop = FALSE], 1, sum)
            z = ytilde - off   #################################################
            if (control$trace)                print(summary(ytilde))
            offdf = sum(eqrank[-i])
            xfunc = xlist[[namesx[i]]][good]
            mgood<- metric[[namesx[i]]]
#            class( mgood)<-"matrix"
#            mgood<- mgood[good,good]
#            attributes(mgood)<-c(attributes(mgood),attributes(metric[[namesx[i]]])[-1])
            h<-par.np2[[namesx[i]]]$h
             if (control$trace)   cat("Range h:", range(h), length(h), "\n")
#   res2 = fregre.np.cv(xfunc, z, h = h, type.CV = dev.S,
#                metric = mgood, par.CV = list(obs = y[good],
#                  family = family, off = off, offdf = offdf,weights = diag(weights)))
           Ker=par.np2[[namesx[i]]]$Ker
           type.S=par.np2[[namesx[i]]]$type.S
           parS=par.np2[[namesx[i]]]$par.S
           parS$w=w
           if (is.function(type.S)) ty<-deparse(substitute(type.S))
           else ty<-type.S
           res = fregre.np.cv(xfunc, z, h = h, type.CV = "dev.S", Ker=Ker,
           type.S=ty,par.S=parS,metric = mgood, par.CV = list(obs = y[good],
           family = family, off = off, offdf = offdf,W = diag(w)))
           if (control$trace)
            cat("Var:",namesx[[i]]," h.opt:", res$h.opt," df:",res$df,"\n")           
           eqrank[namesx[i]] <- res$df
           X[good,namesx[i]] <- res$fitted.values
           result[[namesx[i]]] <- res
        }
        X[, nvars + intercept] = rep(mean(ytilde - apply(X[good,
            1:nvars, drop = FALSE], 1, sum)), nobs)
#        eta <- apply(X, 1, sum)
        eta <- rowSums(X)
        mu <- linkinv(eta <- eta + offset)
#        mu <- linkinv(eta)
        mu.eta.val <- mu.eta(eta)
        good <- (weights > 0) & (mu.eta.val != 0)
        ngoodobs <- as.integer(nobs - sum(mu.eta.val <= eps))
        dev <- sum(dev.resids(y,mu,weights))
        if (control$trace)   {
             par(mfrow = c(1, nvars + 1))
             if (length(table(y) > 10)) {colores = 1}
             else { colores = y + 1 }
             plot(eta, mu, col = colores)
             points(eta, y, col = colores, pch = 2)
             for (i in 1:nvars) {
               plot(X[, i], eta, col = colores, ylab = "Linear Predictor",
                xlab = paste("f(", namesx[i], ")", sep = ""),
                main = paste(namesx[i], "EqPar:", round(eqrank[i],1)))
               if (length(table(y)) == 2) {
                abline(h = 0)
                abline(v = 0)
                }
              }
        }
#        cambio = apply((X - Xold)^2, 2, mean)
        cambio = colMeans((X - Xold)^2)
        if (control$trace) {
            cat("Shift Iter:", iter, "EqRank:", sum(eqrank),
                ngoodobs, "/", nobs, "\n")
             print(cambio)
             }
        if (any(cambio[1:nvars] > control$epsilon)) {conv <- FALSE}
        else{conv <- TRUE;break  }
        if (sum(eqrank) > ngoodobs) {conv <- TRUE;break  }
        if (!(valideta(eta) && validmu(mu))) {
            warning("step size truncated: out of bounds", call. = FALSE)
            ii <- 1
            while (!(valideta(eta) && validmu(mu))) {
                if (ii > control$maxit)
                  stop("inner loop 2; cannot correct step size",
                    call. = FALSE)
                ii <- ii + 1
                X <- (X + Xold)/2
#                eta <- apply(X, 1, sum)
	         eta <- rowSums(X)
                mu <- linkinv(eta)
            }
            boundary <- TRUE
            dev <- sum(dev.resids(y, mu, weights))
            if (control$trace)  cat("Step halved: new deviance =", dev, "\n")
        }
        if (all(!good)) {
            conv <- FALSE
            warning("no observations informative at iteration ",iter)
            break
        }
        if (control$trace)
            cat("Deviance =", dev, "Iterations -", iter, "\n")
        if (abs(dev - devold)/(0.1 + abs(dev)) < control$epsilon |
            dev < control$epsilon) {
            conv <- TRUE
            break
        }
    }
    if (!conv)
        warning("kgam.fit: algorithm did not converge", call. = FALSE)
    if (family$family == "binomial") {
        if (any(mu > 1 - eps) || any(mu < eps))
            warning("kgam.fit: fitted probabilities numerically 0 or 1 occurred",
                call. = FALSE)
    }
    if (family$family == "poisson") {
        if (any(mu < eps))
            warning("kgam.fit: fitted rates numerically 0 occurred",
                call. = FALSE)
    }

    residuals <- (y - mu)#/     (eta) ##### ok??
    nr <- min(sum(good), nvars)
    names(residuals) <- ynames
    names(mu) <- ynames
    names(eta) <- ynames
    wt <- rep.int(0,nobs)
    wt[good] <- w^2
    names(wt) <- ynames
    names(weights) <- ynames
    names(y) <- ynames
    wtdmu <- if (intercept)   sum(weights * y)/sum(weights)
             else linkinv(0)
    nulldev <- sum(dev.resids(y, wtdmu, weights))
    n.ok <- nobs - sum(weights == 0)
    nulldf <- n.ok - as.integer(intercept)
    aic.model <- aic(y, nobs, mu, weights, dev) + 2 * sum(eqrank)
    names(result)<-namesx       
    H<-kgam.H(result,control$inverse)
    
    sr2<-sum(residuals^2)/(nobs - sum(eqrank))
#    if (family$family=="binomial" & !is.factor(y)) y<-factor(y)
    res <- list(result = result, residuals = residuals, fitted.values = mu,
        effects = X, alpha = mean(X[, nvars + intercept]), family = family,
        linear.predictors = eta, deviance = dev, aic = aic.model,
        null.deviance = nulldev, iter = iter, weights = wt, eqrank = eqrank,
        prior.weights = weights, y = y0, converged = conv,H=H,sr2=sr2)
    class(res) <- "fregre.gkam"
    res
}


