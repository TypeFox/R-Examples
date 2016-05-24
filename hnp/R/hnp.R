hnp <-
function(object, sim=99, conf=.95, resid.type, maxit,
                halfnormal=T, scale=F, plot.sim=T, verb.sim=F, warn=F,
                how.many.out=F, print.on=F, paint.out=F, col.paint.out, 
                newclass=F, diagfun, simfun, fitfun, ...) {
  if(!warn) options("warn"=-1)
    
  # use method suitable for object class
  if(newclass) {newhnp(object=object, sim=sim, conf=conf, halfnormal=halfnormal,
                       plot.sim=plot.sim, verb.sim=verb.sim, how.many.out=how.many.out, 
                       print.on=print.on, paint.out=paint.out, col.paint.out=col.paint.out, 
                       diagfun=diagfun, simfun=simfun, fitfun=fitfun, ...) 
  } else {
    UseMethod("hnp")
  }
}

#setGeneric("hnp", def=hnp)

hnp.aodml <-
  function(object, sim=99, conf=.95, resid.type, maxit,
           halfnormal=T, scale=F, plot.sim=T, verb.sim=F, warn=F,
           how.many.out=F, print.on=F, paint.out=F, col.paint.out, ...) {
    
    # preparation and extraction of residuals
    if(missing(maxit)) maxit <- 3000
    if(object$family=="bb") {
      if(missing(resid.type)) resid.type <- "pearson"
      if(resid.type=="deviance") warning("Deviance residuals for beta-binomial models in aods3 do not work properly")
      fam <- "betabin_aods3"
    } else {
      if(missing(resid.type)) resid.type <- "deviance"
      fam <- "negbin_aods3"
    }
    if(halfnormal) {res <- sort(abs(aods3::residuals.aodml(object, type=resid.type)))
    } else {res <- sort(aods3::residuals.aodml(object, type=resid.type))}
    
    # producing the envelope bands
    ## negative binomial
    if(fam=="negbin_aods3") {
      cat("Negative binomial model (using aods3 package)", '\n')
      phi.sc <- object$phi.scale
      if(phi.sc!="inverse") stop("Note that simulation is performed assuming the variance function is parameterized as Var(Y) = mu + mu^2/phi.", "\n", "Please set phi.scale='inverse' in your aodml call.", "\n")
      dat <- object$dat  
      f.nb <- aods3::fitted.aodml(object)
      phi <- aods3::fitted.aodml(object, what="phi")
      rnb <- function(n, mu, phi) {
        k <- if (length(n) > 1L) 
          length(n)
        else n
        rpois(k, (mu * rgamma(k, phi))/phi)
      }
      y. <- lapply(rep(length(f.nb), sim), rnb, f.nb, phi)
      fmla <- as.formula(paste("y.[[i]] ~", object$formula[3]))  
      phi.fmla <- object$phi.formula
      phi.st <- object$call$phi.start
      fixp <- as.list(object$call$fixpar[-1])
      lnk <- object$link
      meth <- object$call$method
      if(halfnormal) {for(i in 1:sim) {
        if(verb.sim) cat("Simulation", i, "out of", sim, "\n")
        res <- cbind(res,sort(abs(aods3::residuals.aodml(
          aods3::aodml(formula=fmla, phi.formula=phi.fmla, phi.scale=phi.sc,
                       phi.start=phi.st, fixpar=fixp, link=lnk, family="nb", 
                       data=dat, method=meth, control=list(maxit=maxit)),type=resid.type))))}
      } else {for(i in 1:sim) {
        if(verb.sim) cat("Simulation", i, "out of", sim, "\n")  
        res <- cbind(res,sort(aods3::residuals.aodml(
          aods3::aodml(formula=fmla, phi.formula=phi.fmla, phi.scale=phi.sc,
                       phi.start=phi.st, fixpar=fixp, link=lnk, family="nb", 
                       data=dat, method=meth, control=list(maxit=maxit)),type=resid.type)))}}
    }
    ## beta-binomial
    if(fam=="betabin_aods3") {
      cat("Beta-binomial model (using aods3 package)", '\n')
      dat <- object$dat  
      m. <- apply(object$resp, 1, sum)
      f.bb <- aods3::fitted.aodml(object)
      phi <- aods3::fitted.aodml(object, what="phi")
      rbb <- function(n, m, p, phi) {
        a <- (1-phi)/phi*p
        b <- (1-phi)/phi*(1-p)
        P <- rbeta(n, a, b)
        Y <- rbinom(n, m, P)
        return(Y)
      }
      y. <- lapply(rep(length(f.bb), sim), rbb, m., f.bb, phi)
      fmla <- as.formula(paste("cbind(y.[[i]], m. - y.[[i]]) ~", object$formula[3]))  
      phi.fmla <- object$phi.formula
      phi.sc <- object$phi.scale
      phi.st <- object$call$phi.start
      fixp <- as.list(object$call$fixpar[-1])
      lnk <- object$link
      meth <- object$call$method
      if(halfnormal) {for(i in 1:sim) {
        if(verb.sim) cat("Simulation", i, "out of", sim, "\n")
        res <- cbind(res,sort(abs(aods3::residuals.aodml(
          aods3::aodml(formula=fmla, phi.formula=phi.fmla, phi.scale=phi.sc,
                       phi.start=phi.st, fixpar=fixp, link=lnk, family="bb", 
                       data=dat, method=meth, control=list(maxit=maxit)),type=resid.type))))}
      } else {for(i in 1:sim) {
        if(verb.sim) cat("Simulation", i, "out of", sim, "\n")  
        res <- cbind(res,sort(aods3::residuals.aodml(
          aods3::aodml(formula=fmla, phi.formula=phi.fmla, phi.scale=phi.sc,
                       phi.start=phi.st, fixpar=fixp, link=lnk, family="bb", 
                       data=dat, method=meth, control=list(maxit=maxit)),type=resid.type)))}}
    }
    
    # now run .makehnp
    .makehnp(obj=res, conf=conf, halfnormal=halfnormal, how.many.out=how.many.out, 
            paint.out=paint.out, col.paint.out=col.paint.out, print.on=print.on, plot.sim=plot.sim, ...)
  }

##setMethod("hnp", "aodml", hnp.aodml)

hnp.aov <-
  function(object, sim=99, conf=.95, resid.type, maxit,
           halfnormal=T, scale=F, plot.sim=T, verb.sim=F, warn=F,
           how.many.out=F, print.on=F, paint.out=F, col.paint.out, ...) {
    
    # preparation and extraction of residuals
    if(missing(resid.type)) resid.type="student"
    if(resid.type=="partial") stop("resid.type should be one of 'deviance', 'pearson', 'response', 'working', 'student', 'standard'")
    get.residuals <- function(obj, type="deviance") {
      if(type=="student") {
        rstudent(obj)
      } else if(type=="standard") {
        rstandard(obj)
      } else resid(obj, type=type)
    }  
    if(halfnormal) {res <- sort(abs(get.residuals(object, type=resid.type)))
    } else {res <- sort(get.residuals(object, type=resid.type))}
    
    # producing the envelope bands
    cat("Gaussian model (aov object)", '\n')
    X <- model.matrix(object)
    dp <- sqrt(anova(object)$"Mean Sq"[length(anova(object)$"Mean Sq")])
    y. <- lapply(rep(length(object$fit), sim), rnorm, object$fit, dp)
    if(is.null(object$offset)) {
      if(halfnormal) {for(i in 1:sim) {
        if(verb.sim) cat("Simulation", i, "out of", sim, "\n")
        res <- cbind(res,sort(abs(get.residuals(aov(y.[[i]] ~ X - 1), type=resid.type))))}
      } else {for(i in 1:sim) {
        if(verb.sim) cat("Simulation", i, "out of", sim, "\n")  
        res <- cbind(res,sort(get.residuals(aov(y.[[i]] ~ X - 1), type=resid.type)))}}
    } else {
      data <- eval(eval(object)$call$data)
      fmla <- as.formula(paste("y.[[i]] ~", paste(as.formula(paste(object$call)[2]))[3]))
      if(halfnormal) {for(i in 1:sim) {
        if(verb.sim) cat("Simulation", i, "out of", sim, "\n")
        res <- cbind(res,sort(abs(get.residuals(aov(formula=fmla, data=data), type=resid.type))))}
      } else {for(i in 1:sim) {
        if(verb.sim) cat("Simulation", i, "out of", sim, "\n")  
        res <- cbind(res,sort(get.residuals(aov(formula=fmla, data=data), type=resid.type)))}}
    }
    
    # now run .makehnp
    .makehnp(obj=res, conf=conf, halfnormal=halfnormal, how.many.out=how.many.out, 
            paint.out=paint.out, col.paint.out=col.paint.out, print.on=print.on, plot.sim=plot.sim, ...)
  }

###setMethod("hnp", "aov", hnp.aov)

hnp.aovlist <-
  function(object, sim=99, conf=.95, resid.type, maxit,
           halfnormal=T, scale=F, plot.sim=T, verb.sim=F, warn=F,
           how.many.out=F, print.on=F, paint.out=F, col.paint.out, 
           newclass=F, diagfun, simfun, fitfun, ...) {
    
    # preparation and extraction of residuals
    if(!missing(resid.type)) warning("resid.type not implemented for 'aovlist' class object")
    fam <- "gaussian.aovlist"
    residuals.aovlist <- function(object, error.term = NULL, ...) {
      aov.proj <- proj(object)
      if (is.null(error.term)) 
        res <- aov.proj[[length(object)]][, "Residuals"]
      else res <- aov.proj[[error.term]][, "Residuals"]
      res
    }
    fitted.aovlist <- function(object, error.term = NULL, ...) {
      aov.proj <- proj(object)
      if(is.null(error.term)) 
        no.strata <- length(object)
      else no.strata <- which(names(aov.proj) == error.term)
      fit <- aov.proj[["(Intercept)"]][, 1]
      for(i in 2:no.strata) {
        nterms <- ncol(aov.proj[[i]])
        if(dimnames(aov.proj[[i]])[[2]][nterms] == "Residuals") 
          nterms <- nterms - 1
        if(nterms > 0) 
          if(nterms == 1) 
            fit <- fit + aov.proj[[i]][, 1]
          else fit <- fit + rowSums(aov.proj[[i]][, 1:nterms])
      }
      fit
    }
    res <- residuals(object)
    if(halfnormal) {res <- sort(abs(res))} else {res <- sort(res)}
    
    # producing the envelope bands
    cat("Gaussian model (aovlist object)", '\n')
    data <- eval(attr(object, "call")$data)
    msq <- summary(object)[[length(summary(object))]][[1]]$"Mean Sq"
    dp <- sqrt(msq[length(msq)])
    obj.fv <- fitted(object)
    y. <- lapply(rep(length(obj.fv), sim), rnorm, obj.fv, dp)
    xnam <- as.character(terms(object))[3]
    fmla <- as.formula(paste("y.[[i]] ~", paste(xnam, collapse="+")))
    if(halfnormal) {for(i in 1:sim) {
      if(verb.sim) cat("Simulation", i, "out of", sim, "\n")
      res <- cbind(res,sort(abs(residuals(aov(formula=fmla, data=data)))))}
    } else {for(i in 1:sim) {
      if(verb.sim) cat("Simulation", i, "out of", sim, "\n")  
      res <- cbind(res,sort(residuals(aov(formula=fmla, data=data))))}}
    
    # now run .makehnp
    .makehnp(obj=res, conf=conf, halfnormal=halfnormal, how.many.out=how.many.out, 
            paint.out=paint.out, col.paint.out=col.paint.out, print.on=print.on, plot.sim=plot.sim, ...)
  }

##setMethod("hnp", "aovlist", hnp.aovlist)

hnp.gamlss <-
  function(object, sim=99, conf=.95, resid.type, maxit,
           halfnormal=T, scale=F, plot.sim=T, verb.sim=F, warn=F,
           how.many.out=F, print.on=F, paint.out=F, col.paint.out, ...) {
    
    # preparation and extraction of residuals
    if(missing(maxit)) maxit <- 25
    if(missing(resid.type)) resid.type <- "simple"
    if(object$family[1]!="ZIBI"&object$family[1]!="ZIBB"&object$family[1]!="BB") stop("This function has been implemented for gamlss objects with family=ZIBI, ZIBB and BB")
    fam <- object$family[1]
    if(halfnormal) {res <- sort(abs(resid(object, type=resid.type)))
    } else {res <- sort(resid(object, type=resid.type))}
    
    # producing the envelope bands
    ## beta-binomial
    if(fam=="BB") {
      cat("Beta-binomial model (using gamlss)", '\n')
      data <- eval(eval(object)$call$data)
      m. <- object$bd
      f.bb <- fitted(object)
      sig <- object$sigma.fv
      y. <- lapply(rep(length(f.bb), sim), gamlss.dist::rBB, f.bb, sig, m.)
      fmla <- as.formula(paste("cbind(y.[[i]],m.-y.[[i]]) ~ ", paste(object$mu.formula)[3]))
      if(halfnormal) {for(i in 1:sim) {
        if(verb.sim) cat("Simulation", i, "out of", sim, "\n")
        res <- cbind(res, sort(abs(resid(gamlss::gamlss(formula=fmla, family=gamlss.dist::BB, control=gamlss::gamlss.control(n.cyc=maxit), data=data), type=resid.type))))}
      } else {for(i in 1:sim) {
        if(verb.sim) cat("Simulation", i, "out of", sim, "\n")  
        res <- cbind(res, sort(resid(gamlss::gamlss(formula=fmla, family=gamlss.dist::BB, control=gamlss::gamlss.control(n.cyc=maxit), data=data), type=resid.type)))}}
    }   
    ## zero-inflated binomial
    if(fam=="ZIBI") {
      cat("Zero-inflated binomial model (using gamlss)", '\n')
      data <- eval(eval(object)$call$data)
      m. <- object$bd
      sig <- object$sigma.fv
      f.zibi <- fitted(object)
      y. <- lapply(rep(length(f.zibi), sim), gamlss.dist::rZIBI, m., f.zibi, sig)
      fmla <- as.formula(paste("cbind(y.[[i]],m.-y.[[i]]) ~ ", paste(object$mu.formula)[3]))
      if(halfnormal) {for(i in 1:sim) {
        if(verb.sim) cat("Simulation", i, "out of", sim, "\n")
        res <- cbind(res, sort(abs(resid(gamlss::gamlss(formula=fmla, nu.formula=object$nu.formula, family=gamlss.dist::ZIBI, control=gamlss::gamlss.control(n.cyc=maxit), data=data), type=resid.type))))}
      } else {for(i in 1:sim) {
        if(verb.sim) cat("Simulation", i, "out of", sim, "\n")  
        res <- cbind(res, sort(resid(gamlss::gamlss(formula=fmla, nu.formula=object$nu.formula, family=gamlss.dist::ZIBI, control=gamlss::gamlss.control(n.cyc=maxit), data=data), type=resid.type)))}}
    } 
    ## zero-inflated beta-binomial
    if(fam=="ZIBB") {
      cat("Zero-inflated beta-binomial model (using gamlss)", '\n')
      data <- eval(eval(object)$call$data)
      m. <- object$bd
      f.zibb <- fitted(object)
      sig <- object$sigma.fv
      nu <- object$nu.fv
      y. <- lapply(rep(length(f.zibb), sim), gamlss.dist::rZIBB, f.zibb, sig, nu, m.)
      fmla <- as.formula(paste("cbind(y.[[i]],m.-y.[[i]]) ~ ", paste(object$mu.formula)[3]))
      if(halfnormal) {for(i in 1:sim) {
        if(verb.sim) cat("Simulation", i, "out of", sim, "\n")
        res <- cbind(res, sort(abs(resid(gamlss::gamlss(formula=fmla, nu.formula=object$nu.formula, family=gamlss.dist::ZIBB, control=gamlss::gamlss.control(n.cyc=maxit), data=data), type=resid.type))))}
      } else {for(i in 1:sim) {
        if(verb.sim) cat("Simulation", i, "out of", sim, "\n")  
        res <- cbind(res, sort(resid(gamlss::gamlss(formula=fmla, nu.formula=object$nu.formula, family=gamlss.dist::ZIBB, control=gamlss::gamlss.control(n.cyc=maxit), data=data), type=resid.type)))}}
    }   
    
    # now run .makehnp
    .makehnp(obj=res, conf=conf, halfnormal=halfnormal, how.many.out=how.many.out, 
            paint.out=paint.out, col.paint.out=col.paint.out, print.on=print.on, plot.sim=plot.sim, ...)
  }

#setMethod("hnp", "gamlss", hnp.gamlss)

hnp.glm <-
  function(object, sim=99, conf=.95, resid.type, maxit,
           halfnormal=T, scale=F, plot.sim=T, verb.sim=F, warn=F,
           how.many.out=F, print.on=F, paint.out=F, col.paint.out, ...) {
    
    # preparation and extraction of residuals
    if(missing(maxit)) maxit <- 25
    if(missing(resid.type)) resid.type <- "deviance"
    if(resid.type=="partial") stop("resid.type should be one of 'deviance', 'pearson', 'response', 'working', 'student', 'standard'")
    wg <- weights(object)
    get.residuals <- function(obj, type="deviance") {
      if(type=="student") {
        rstudent(obj)
      } else if(type=="standard") {
        rstandard(obj)
      } else resid(obj, type=type)
    }  
    fam <- object$family$family
    if(!is.na(grep("Negative", object$family$family)[1])) stop("Please use the glm.nb function available in package MASS or the aodml function available in package aods3 to fit negative binomial models.", "\n", "If you wish to use family=negative.binomial in the glm function, please check the package vignette for detailed instructions on how to supply your codes to the hnp function.")
    link <- object$family$link
    if(halfnormal) {res <- sort(abs(get.residuals(object, type=resid.type)))
    } else {res <- sort(get.residuals(object, type=resid.type))}
    
    # producing the envelope bands
    ## binomial
    if(fam=="binomial") {
      cat("Binomial model", '\n')
      X <- model.matrix(object)
      m. <- object$prior.weights
      m2 <- ifelse(m.==0, 1, m.)
      y. <- lapply(rep(length(object$fit), sim), rbinom, m., object$fit)
      if(is.null(object$offset)) {
        if(halfnormal) {for(i in 1:sim) {
          if(verb.sim) cat("Simulation", i, "out of", sim, "\n")
          res <- cbind(res,sort(abs(get.residuals(glm(cbind(y.[[i]],m.-y.[[i]]) ~ X - 1, family=binomial(link=link), weights=wg/m2, control=glm.control(maxit=maxit)), type=resid.type))))}
        } else {for(i in 1:sim) {
          if(verb.sim) cat("Simulation", i, "out of", sim, "\n")
          res <- cbind(res,sort(get.residuals(glm(cbind(y.[[i]],m.-y.[[i]]) ~ X - 1, family=binomial(link=link), weights=wg/m2, control=glm.control(maxit=maxit)), resid.type)))}}
      } else {
        fmla <- as.formula(paste("cbind(y.[[i]], m. - y.[[i]]) ~", paste(as.formula(paste(object$call)[2]))[3]))
        if(is.data.frame(object$data)) {
          if(halfnormal) {for(i in 1:sim) {
            if(verb.sim) cat("Simulation", i, "out of", sim, "\n")  
            res <- cbind(res,sort(abs(get.residuals(glm(formula=fmla, family=binomial(link=link), weights=wg/m2, control=glm.control(maxit=maxit), data=object$data), type=resid.type))))}
          } else {for(i in 1:sim) {
            if(verb.sim) cat("Simulation", i, "out of", sim, "\n")
            res <- cbind(res,sort(get.residuals(glm(formula=fmla, family=binomial(link=link), weights=wg/m2, control=glm.control(maxit=maxit), data=object$data), type=resid.type)))}}
        } else {
          if(halfnormal) {for(i in 1:sim) {
            if(verb.sim) cat("Simulation", i, "out of", sim, "\n")
            res <- cbind(res,sort(abs(get.residuals(glm(formula=fmla, family=binomial(link=link), weights=wg/m2, control=glm.control(maxit=maxit)), type=resid.type))))}
          } else {for(i in 1:sim) {
            if(verb.sim) cat("Simulation", i, "out of", sim, "\n")  
            res <- cbind(res,sort(get.residuals(glm(formula=fmla, family=binomial(link=link), weights=wg/m2, control=glm.control(maxit=maxit)), type=resid.type)))}}
        }
      }
    }
    ## quasibinomial
    if(fam=="quasibinomial") {
      cat("Quasi-binomial model", '\n')
      X <- model.matrix(object)
      m. <- object$prior.weights
      m2 <- ifelse(m.==0, 1, m.)
      phi <- summary(object)$dispersion
      y. <- lapply(rep(length(object$fit), sim), rbinom, m., object$fit)
      y. <- lapply(y., function(x) x*phi)
      if(halfnormal) {res <- sort(abs((1/sqrt(phi))*get.residuals(object, type=resid.type)))
      } else {res <- sort((1/sqrt(phi))*get.residuals(object, type=resid.type))}  
      if(is.null(object$offset)) {
        if(halfnormal) {for(i in 1:sim) {
          if(verb.sim) cat("Simulation", i, "out of", sim, "\n")  
          res <- cbind(res,sort(abs((1/sqrt(phi))*get.residuals(glm(cbind(y.[[i]],m.*phi-y.[[i]]) ~ X - 1, family=binomial(link=link), weights=wg/m2, control=glm.control(maxit=maxit)), type=resid.type))))}
        } else {for(i in 1:sim) {
          if(verb.sim) cat("Simulation", i, "out of", sim, "\n")  
          res <- cbind(res,sort((1/sqrt(phi))*get.residuals(glm(cbind(y.[[i]],m.*phi-y.[[i]]) ~ X - 1, family=binomial(link=link), weights=wg/m2, control=glm.control(maxit=maxit)), type=resid.type)))}}
      } else {
        fmla <- as.formula(paste("cbind(y.[[i]], m.*phi - y.[[i]]) ~", paste(as.formula(paste(object$call)[2]))[3]))
        if(is.data.frame(object$data)) {
          if(halfnormal) {for(i in 1:sim) {
            if(verb.sim) cat("Simulation", i, "out of", sim, "\n")
            res <- cbind(res,sort(abs((1/sqrt(phi))*get.residuals(glm(formula=fmla, family=quasibinomial(link=link), weights=wg/m2, control=glm.control(maxit=maxit), data=object$data), type=resid.type))))}
          } else {for(i in 1:sim) {
            if(verb.sim) cat("Simulation", i, "out of", sim, "\n")  
            res <- cbind(res,sort((1/sqrt(phi))*get.residuals(glm(formula=fmla, family=quasibinomial(link=link), weights=wg/m2, control=glm.control(maxit=maxit), data=object$data), type=resid.type)))}}
        } else {
          if(halfnormal) {for(i in 1:sim) {
            if(verb.sim) cat("Simulation", i, "out of", sim, "\n")
            res <- cbind(res,sort(abs((1/sqrt(phi))*get.residuals(glm(formula=fmla, family=quasibinomial(link=link), weights=wg/m2, control=glm.control(maxit=maxit)), type=resid.type))))}
          } else {for(i in 1:sim) {
            if(verb.sim) cat("Simulation", i, "out of", sim, "\n")  
            res <- cbind(res,sort((1/sqrt(phi))*get.residuals(glm(formula=fmla, family=quasibinomial(link=link), weights=wg/m2, control=glm.control(maxit=maxit)), type=resid.type)))}}
        }
      }
    }
    ## Poisson
    if(fam=="poisson") {
      cat("Poisson model", '\n')
      X <- model.matrix(object)
      y. <- lapply(rep(length(object$fit), sim), rpois, object$fit)
      if(is.null(object$offset)) {
        if(halfnormal) {for(i in 1:sim) {
          if(verb.sim) cat("Simulation", i, "out of", sim, "\n")
          res <- cbind(res,sort(abs(get.residuals(glm(y.[[i]] ~ X - 1, family=poisson(link=link), weights=wg, control=glm.control(maxit=maxit)), type=resid.type))))}
        } else {for(i in 1:sim) {
          if(verb.sim) cat("Simulation", i, "out of", sim, "\n")  
          res <- cbind(res,sort(get.residuals(glm(y.[[i]] ~ X - 1, family=poisson(link=link), weights=wg, control=glm.control(maxit=maxit)), type=resid.type)))}}
      } else {
        fmla <- as.formula(paste("y.[[i]] ~", paste(as.formula(paste(object$call)[2]))[3]))
        if(is.data.frame(object$data)) {
          if(halfnormal) {for(i in 1:sim) {
            if(verb.sim) cat("Simulation", i, "out of", sim, "\n")
            res <- cbind(res,sort(abs(get.residuals(glm(formula=fmla, family=poisson(link=link), weights=wg, control=glm.control(maxit=maxit), data=object$data), type=resid.type))))}
          } else {for(i in 1:sim) {
            if(verb.sim) cat("Simulation", i, "out of", sim, "\n")  
            res <- cbind(res,sort(get.residuals(glm(formula=fmla, family=poisson(link=link), weights=wg, control=glm.control(maxit=maxit), data=object$data), type=resid.type)))}}
        } else {
          if(halfnormal) {for(i in 1:sim) {
            if(verb.sim) cat("Simulation", i, "out of", sim, "\n")
            res <- cbind(res,sort(abs(get.residuals(glm(formula=fmla, family=poisson(link=link), weights=wg, control=glm.control(maxit=maxit)), type=resid.type))))}
          } else {for(i in 1:sim) {
            if(verb.sim) cat("Simulation", i, "out of", sim, "\n")  
            res <- cbind(res,sort(get.residuals(glm(formula=fmla, family=poisson(link=link), weights=wg, control=glm.control(maxit=maxit)), type=resid.type)))}}    
        }
      }
    }
    ## quasi-Poisson
    if(fam=="quasipoisson") {
      cat("Quasi-Poisson model", '\n')
      X <- model.matrix(object)
      phi <- summary(object)$dispersion
      rqpois <- function(n, mu, theta) rnbinom(n = n, mu = mu, size = mu/(theta-1))
      if(summary(object)$dispersion < 1) {stop("Data is underdispersed. Consider fitting another type of model.")}
      y. <- lapply(rep(length(object$fit), sim), rqpois, object$fit, summary(object)$dispersion)
      if(halfnormal) {res <- sort(abs((1/sqrt(phi))*get.residuals(object, type=resid.type)))
      } else {res <- sort((1/sqrt(phi))*get.residuals(object, type=resid.type))}  
      if(is.null(object$offset)) {
        if(halfnormal) {for(i in 1:sim) {
          if(verb.sim) cat("Simulation", i, "out of", sim, "\n")
          res <- cbind(res,sort(abs((1/sqrt(phi))*get.residuals(glm(y.[[i]] ~ X - 1, family=quasipoisson(link=link), weights=wg, control=glm.control(maxit=maxit)), type=resid.type))))}
        } else {for(i in 1:sim) {
          if(verb.sim) cat("Simulation", i, "out of", sim, "\n")  
          res <- cbind(res,sort((1/sqrt(phi))*get.residuals(glm(y.[[i]] ~ X - 1, family=quasipoisson(link=link), weights=wg, control=glm.control(maxit=maxit)), type=resid.type)))}}
      } else {
        fmla <- as.formula(paste("y.[[i]] ~", paste(as.formula(paste(object$call)[2]))[3]))
        if(is.data.frame(object$data)) {
          if(halfnormal) {for(i in 1:sim) {
            if(verb.sim) cat("Simulation", i, "out of", sim, "\n")
            res <- cbind(res,sort(abs((1/sqrt(phi))*get.residuals(glm(formula=fmla, family=quasipoisson(link=link), weights=wg, control=glm.control(maxit=maxit), data=object$data), type=resid.type))))}
          } else {for(i in 1:sim) {
            if(verb.sim) cat("Simulation", i, "out of", sim, "\n")   
            res <- cbind(res,sort((1/sqrt(phi))*get.residuals(glm(formula=fmla, family=quasipoisson(link=link), weights=wg, control=glm.control(maxit=maxit), data=object$data), type=resid.type)))}}
        } else {
          if(halfnormal) {for(i in 1:sim) {
            if(verb.sim) cat("Simulation", i, "out of", sim, "\n")
            res <- cbind(res,sort(abs((1/sqrt(phi))*get.residuals(glm(formula=fmla, family=quasipoisson(link=link), weights=wg, control=glm.control(maxit=maxit)), type=resid.type))))}
          } else {for(i in 1:sim) {
            if(verb.sim) cat("Simulation", i, "out of", sim, "\n")  
            res <- cbind(res,sort((1/sqrt(phi))*get.residuals(glm(formula=fmla, family=quasipoisson(link=link), weights=wg, control=glm.control(maxit=maxit)), type=resid.type)))}}    
        }
      }
    }  
    ## gamma
    if(fam=="Gamma") {
      cat("Gamma model", '\n')
      X <- model.matrix(object)
      mu <- predict(object, type="response")
      fi <- (nrow(X)-ncol(X))/sum((resid(object,type="response")/mu)^2)
      y. <- lapply(rep(length(object$fit), sim), rgamma, fi, fi/mu)
      if(is.null(object$offset)) {
        if(halfnormal) {for(i in 1:sim) {
          if(verb.sim) cat("Simulation", i, "out of", sim, "\n")
          res <- cbind(res,sort(abs(get.residuals(glm(y.[[i]] ~ X - 1, family=Gamma(link=link), weights=wg, start=coef(object), control=glm.control(maxit=maxit)), type=resid.type))))}
        } else {for(i in 1:sim) {
          if(verb.sim) cat("Simulation", i, "out of", sim, "\n")  
          res <- cbind(res,sort(get.residuals(glm(y.[[i]] ~ X - 1, family=Gamma(link=link), weights=wg, start=coef(object), control=glm.control(maxit=maxit)), type=resid.type)))}}
      } else {
        fmla <- as.formula(paste("y.[[i]] ~", paste(as.formula(paste(object$call)[2]))[3]))
        if(is.data.frame(object$data)) {
          if(halfnormal) {for(i in 1:sim) {
            if(verb.sim) cat("Simulation", i, "out of", sim, "\n")
            res <- cbind(res,sort(abs(get.residuals(glm(formula=fmla, family=Gamma(link=link), weights=wg, control=glm.control(maxit=maxit), data=object$data), type=resid.type))))}
          } else {for(i in 1:sim) {
            if(verb.sim) cat("Simulation", i, "out of", sim, "\n")  
            res <- cbind(res,sort(get.residuals(glm(formula=fmla, family=Gamma(link=link), weights=wg, control=glm.control(maxit=maxit), data=object$data), type=resid.type)))}}
        } else {
          if(halfnormal) {for(i in 1:sim) {
            if(verb.sim) cat("Simulation", i, "out of", sim, "\n")
            res <- cbind(res,sort(abs(get.residuals(glm(formula=fmla, family=Gamma(link=link), weights=wg, control=glm.control(maxit=maxit)), type=resid.type))))}
          } else {for(i in 1:sim) {
            if(verb.sim) cat("Simulation", i, "out of", sim, "\n")  
            res <- cbind(res,sort(get.residuals(glm(formula=fmla, family=Gamma(link=link), weights=wg, control=glm.control(maxit=maxit)), type=resid.type)))}}    
        }
      }
    }
    ## inverse gaussian
    if(fam=="inverse.gaussian") {
      cat("Inverse gaussian model", '\n')
      X <- model.matrix(object)
      rig <- function (n, mean, scale) {
        if(length(n)>1) n <- length(n)
        y <- rnorm(n)^2
        mu2 <- 0 * y + mean^2
        x <- mean + 0.5*scale*(mu2*y - mean*sqrt(4*mean*y/scale + mu2*y^2))
        ind <- runif(n) > mean/(mean + x)
        x[ind] <- mu2[ind]/x[ind]
        x
      }
      n <- nrow(X)
      p <- ncol(X)
      y <- object$y
      m <- predict(object, type="response")
      lam <- (n-p)/sum((resid(object, type="response")/(m*sqrt(y)))^2)
      y. <- lapply(rep(length(object$fit), sim), rig, fitted(object), 1/lam)
      if(is.null(object$offset)) {
        if(halfnormal) {for(i in 1:sim) {
          if(verb.sim) cat("Simulation", i, "out of", sim, "\n")
          res <- cbind(res,sort(abs(get.residuals(glm(y.[[i]] ~ X - 1, family=inverse.gaussian(link=link), weights=wg, control=glm.control(maxit=maxit)), type=resid.type))))}
        } else {for(i in 1:sim) {
          if(verb.sim) cat("Simulation", i, "out of", sim, "\n")  
          res <- cbind(res,sort(get.residuals(glm(y.[[i]] ~ X - 1, family=inverse.gaussian(link=link), weights=wg, control=glm.control(maxit=maxit)), type=resid.type)))}}
      } else {
        fmla <- as.formula(paste("y.[[i]] ~", paste(as.formula(paste(object$call)[2]))[3]))
        if(is.data.frame(object$data)) {
          if(halfnormal) {for(i in 1:sim) {
            if(verb.sim) cat("Simulation", i, "out of", sim, "\n")
            res <- cbind(res,sort(abs(get.residuals(glm(formula=fmla, family=inverse.gaussian(link=link), weights=wg, control=glm.control(maxit=maxit), data=object$data), type=resid.type))))}
          } else {for(i in 1:sim) {
            if(verb.sim) cat("Simulation", i, "out of", sim, "\n")  
            res <- cbind(res,sort(get.residuals(glm(formula=fmla, family=inverse.gaussian(link=link), weights=wg, control=glm.control(maxit=maxit), data=object$data), type=resid.type)))}}
        } else {
          if(halfnormal) {for(i in 1:sim) {
            if(verb.sim) cat("Simulation", i, "out of", sim, "\n")
            res <- cbind(res,sort(abs(get.residuals(glm(formula=fmla, family=inverse.gaussian(link=link), weights=wg, control=glm.control(maxit=maxit)), type=resid.type))))}
          } else {for(i in 1:sim) {
            if(verb.sim) cat("Simulation", i, "out of", sim, "\n")  
            res <- cbind(res,sort(get.residuals(glm(formula=fmla, family=inverse.gaussian(link=link), weights=wg, control=glm.control(maxit=maxit)), type=resid.type)))}}    
        }
      }
    }
    ## gaussian
    if(fam=="gaussian") {
      cat("Gaussian model (glm object)", '\n')
      X <- model.matrix(object)
      dp <- sqrt(summary(object)$deviance/summary(object)$df.residual)
      y. <- lapply(rep(length(object$fit), sim), rnorm, object$fit, dp)
      if(is.null(object$offset)) {
        if(halfnormal) {for(i in 1:sim) {
          if(verb.sim) cat("Simulation", i, "out of", sim, "\n")
          res <- cbind(res,sort(abs(get.residuals(glm(y.[[i]] ~ X - 1, family=gaussian(link=link), control=glm.control(maxit=maxit)), type=resid.type))))}
        } else {for(i in 1:sim) {
          if(verb.sim) cat("Simulation", i, "out of", sim, "\n")  
          res <- cbind(res,sort(get.residuals(glm(y.[[i]] ~ X - 1, family=gaussian(link=link), control=glm.control(maxit=maxit)), type=resid.type)))}}
      } else {
        fmla <- as.formula(paste("y.[[i]] ~", paste(as.formula(paste(object$call)[2]))[3]))
        if(is.data.frame(object$data)) {
          if(halfnormal) {for(i in 1:sim) {
            if(verb.sim) cat("Simulation", i, "out of", sim, "\n")
            res <- cbind(res,sort(abs(get.residuals(glm(formula=fmla, family=gaussian(link=link), control=glm.control(maxit=maxit), data=object$data), type=resid.type))))}
          } else {for(i in 1:sim) {
            if(verb.sim) cat("Simulation", i, "out of", sim, "\n")  
            res <- cbind(res,sort(get.residuals(glm(formula=fmla, family=gaussian(link=link), control=glm.control(maxit=maxit), data=object$data), type=resid.type)))}}
        } else {
          if(halfnormal) {for(i in 1:sim) {
            if(verb.sim) cat("Simulation", i, "out of", sim, "\n")
            res <- cbind(res,sort(abs(get.residuals(glm(formula=fmla, family=gaussian(link=link), control=glm.control(maxit=maxit)), type=resid.type))))}
          } else {for(i in 1:sim) {
            if(verb.sim) cat("Simulation", i, "out of", sim, "\n")  
            res <- cbind(res,sort(get.residuals(glm(formula=fmla, family=gaussian(link=link), control=glm.control(maxit=maxit)), type=resid.type)))}}    
        }
      }
    }
    
    # now run .makehnp
    .makehnp(obj=res, conf=conf, halfnormal=halfnormal, how.many.out=how.many.out, 
            paint.out=paint.out, col.paint.out=col.paint.out, print.on=print.on, plot.sim=plot.sim, ...)
  }

#setMethod("hnp", "glm", hnp.glm)

hnp.glmerMod <-
  function(object, sim=99, conf=.95, resid.type, maxit,
           halfnormal=T, scale=F, plot.sim=T, verb.sim=F, warn=F,
           how.many.out=F, print.on=F, paint.out=F, col.paint.out, ...) {
    
    # preparation and extraction of residuals
    if(missing(maxit)) maxit <- 300
    if(missing(resid.type)) resid.type <- "deviance"
    if(resid.type=="partial") stop("resid.type should be one of 'deviance', 'pearson', 'response', 'working'")
    .fam <- object@resp$family$family 
    if(.fam=="..1") stop("Please do not use update() when fitting another link function")
    if(.fam!="poisson" & .fam!="binomial") stop("Function implemented for family='poisson' or family='binomial'")  
    if(.fam=="binomial") fam <- "binomial-normal"
    if(.fam=="poisson") fam <- "poisson-normal"
    link <- object@resp$family$link
    if(halfnormal) {res <- sort(abs(residuals(object, type=resid.type)))
    } else {res <- sort(residuals(object, type=resid.type))}
    
    # producing the envelope bands
    ## binomial-normal
    if(fam=="binomial-normal") {
      cat("Binomial-normal model", '\n')
      data <- eval(eval(object)@call$data)
      rbn <- function(x) simulate(object)[,1][,1]
      m. <- object@resp$weights
      y. <- lapply(1:sim, rbn)
      fmla <- as.formula(paste("cbind(y.[[i]], m.-y.[[i]]) ~", paste(as.formula(paste(object@call)[2]))[3]))
      if(halfnormal) {for(i in 1:sim) {
        if(verb.sim) cat("Simulation", i, "out of", sim, "\n")
        res <- cbind(res,sort(abs(residuals(lme4::glmer(fmla, family=binomial(link=link), control=lme4::glmerControl(optCtrl=list(maxfun=maxit)), data=data), type=resid.type))))}
      } else {for(i in 1:sim) {
        if(verb.sim) cat("Simulation", i, "out of", sim, "\n")  
        res <- cbind(res,sort(residuals(lme4::glmer(fmla, family=binomial(link=link), control=lme4::glmerControl(optCtrl=list(maxfun=maxit)), data=data), type=resid.type)))}}
    }
    ## Poisson-normal
    if(fam=="poisson-normal") {
      cat("Poisson-normal model", '\n')
      data <- eval(eval(object)@call$data)
      rpn <- function(x) simulate(object)[,1]
      y. <- lapply(1:sim, rpn)
      fmla <- as.formula(paste("y.[[i]] ~", paste(as.formula(paste(object@call)[2]))[3]))
      if(halfnormal) {for(i in 1:sim) {
        if(verb.sim) cat("Simulation", i, "out of", sim, "\n")
        res <- cbind(res,sort(abs(residuals(lme4::glmer(fmla, family=poisson, control=lme4::glmerControl(optCtrl=list(maxfun=maxit)), data=data), type=resid.type))))}
      } else {for(i in 1:sim) {
        if(verb.sim) cat("Simulation", i, "out of", sim, "\n")  
        res <- cbind(res,sort(residuals(lme4::glmer(fmla, family=poisson, control=lme4::glmerControl(optCtrl=list(maxfun=maxit)), data=data), type=resid.type)))}}
    }
    
    # now run .makehnp
    .makehnp(obj=res, conf=conf, halfnormal=halfnormal, how.many.out=how.many.out, 
            paint.out=paint.out, col.paint.out=col.paint.out, print.on=print.on, plot.sim=plot.sim, ...)
  }

#setMethod("hnp", "glmerMod", hnp.glmerMod)

hnp.glmmadmb <-
  function(object, sim = 99, conf = 0.95, resid.type, maxit,
           halfnormal = T, scale = F, plot.sim = T, verb.sim = F,
           warn = F, how.many.out = F, print.on = F, paint.out = F,
           col.paint.out, newclass = F, diagfun, simfun, fitfun,
           ...) {
    
    # preparation and extraction of residuals
    if(missing(maxit)) maxit <- 300
    if(missing(resid.type)) resid.type <- "response"
    if(object$family[1]!="binom"&object$family[1]!="betabinom"|object$zeroInflation==F) stop("This function has been implemented for glmmadmb objects with family='binomial' and 'betabinomial' with zeroInflation=TRUE")
    if(object$family[1]=="binom") fam <- "ZIBI_admb"
    if(object$family[1]=="betabinom") fam <- "ZIBB_admb"
    link <- object$link
    if(halfnormal) {res <- sort(abs(resid(object, type=resid.type)))
    } else {res <- sort(resid(object, type=resid.type))}
    
    # producing the envelope bands
    ## zero-inflated binomial
    if(fam=="ZIBI_admb") {  
      cat("Zero-inflated binomial model (using glmmADMB)", '\n')
      data <- eval(eval(object)$call$data)
      if(is.null(data)) stop("Please use the argument data within your glmmadmb call so that hnp() can use the dataset internally to refit the model with simulated samples")
      rzibinom <- function(n, m, p, p0) {
        Z <- rbinom(n, 1, p0)
        B <- rbinom(n, m, p)
        r <- Z * 0 + (1 - Z) * B
        return(r)
      }
      m. <- apply(object$frame[,1], 1, sum)
      f.zibi <- fitted(object)
      pz <- object$pz
      y. <- lapply(rep(length(f.zibi), sim), rzibinom, m., f.zibi, pz)
      fmla <- as.formula(paste("cbind(y.,m.-y.) ~ ", paste(object$terms)[3]))
      if(halfnormal) {for(i in 1:sim) {
        if(verb.sim) cat("Simulation", i, "out of", sim, "\n")  
        datax <- data.frame(data, y.=y.[[i]], m.)
        res <- cbind(res, sort(abs(resid(glmmADMB::glmmadmb(formula=fmla, family="binomial", link=link, zeroInflation=T, admb.opts=glmmADMB::admbControl(maxfn=maxit), data=datax), type=resid.type))))}
      } else {for(i in 1:sim) {
        if(verb.sim) cat("Simulation", i, "out of", sim, "\n")
        datax <- data.frame(data, y.=y.[[i]], m.)  
        res <- cbind(res, sort(resid(glmmADMB::glmmadmb(formula=fmla, family="binomial", link=link, zeroInflation=T, admb.opts=glmmADMB::admbControl(maxfn=maxit), data=datax), type=resid.type)))}}
    }   
    ## zero-inflated beta-binomial
    if(fam=="ZIBB_admb") {  
      cat("Zero-inflated beta-binomial model (using glmmADMB)", '\n')
      data <- eval(eval(object)$call$data)
      if(is.null(data)) stop("Please use the argument data within your glmmadmb call so that hnp() can use the dataset internally to refit the model with simulated samples")
      rbb <- function(n, m, p, phi) {
        a <- (1-phi)/phi*p
        b <- (1-phi)/phi*(1-p)
        P <- rbeta(n, a, b)
        Y <- rbinom(n, m, P)
        return(Y)
      }
      rzibetabin <- function(n, m, p, p0, phi) {
        Z <- rbinom(n, 1, p0)
        B <- rbb(n, m, p, phi)
        r <- Z * 0 + (1 - Z) * B
        return(r)
      }
      m. <- apply(object$frame[,1], 1, sum)
      pz <- object$pz
      f.zibb <- fitted(object)/(1-pz)
      phi <- 1/object$alpha
      y. <- lapply(rep(length(f.zibb), sim), rzibetabin, m., f.zibb, pz, phi)
      fmla <- as.formula(paste("cbind(y.,m.-y.) ~ ", paste(object$terms)[3]))
      if(halfnormal) {for(i in 1:sim) {
        if(verb.sim) cat("Simulation", i, "out of", sim, "\n")  
        datax <- data.frame(data, y.=y.[[i]], m.)
        res <- cbind(res, sort(abs(resid(glmmADMB::glmmadmb(formula=fmla, family="betabinomial", link=link, zeroInflation=T, admb.opts=glmmADMB::admbControl(maxfn=maxit), data=datax), type=resid.type))))}
      } else {for(i in 1:sim) {
        if(verb.sim) cat("Simulation", i, "out of", sim, "\n")
        datax <- data.frame(data, y.=y.[[i]], m.)
        res <- cbind(res, sort(resid(glmmADMB::glmmadmb(formula=fmla, family="betabinomial", link=link, zeroInflation=T, admb.opts=glmmADMB::admbControl(maxfn=maxit), data=datax), type=resid.type)))}}
    } 
    
    # now run .makehnp
    .makehnp(obj=res, conf=conf, halfnormal=halfnormal, how.many.out=how.many.out, 
            paint.out=paint.out, col.paint.out=col.paint.out, print.on=print.on, plot.sim=plot.sim, ...)
  }

#setMethod("hnp", "glmmadmb", hnp.glmmadmb)

hnp.hurdle <-
  function(object, sim=99, conf=.95, resid.type, maxit,
           halfnormal=T, scale=F, plot.sim=T, verb.sim=F, warn=F,
           how.many.out=F, print.on=F, paint.out=F, col.paint.out, ...) {
    
    # preparation and extraction of residuals
    if(missing(maxit)) maxit <- 10000
    if(missing(resid.type)) resid.type <- "pearson"
    if(resid.type!="pearson"&resid.type!="response") stop("resid.type should be one of 'pearson', 'response'")
    if(object$dist$zero!="binomial") stop("Zero distribution should be specified as 'binomial'")
    if(object$dist$count=="poisson") fam <- "HP"
    if(object$dist$count=="negbin") fam <- "HNB"
    if(object$dist$count=="geometric") stop("This function has not been implemented yet for hurdle geometric models")
    if(halfnormal) {res <- sort(abs(residuals(object, type=resid.type)))
    } else {res <- sort(residuals(object, type=resid.type))}
    
    # producing the envelope bands
    ## hurdle Poisson  
    if(fam=="HP") {
      cat("Hurdle Poisson model", '\n')
      X <- model.matrix(object, type="count")
      Z <- model.matrix(object, type="zero")
      obj.link <- object$link
      rtruncpois <- function(n, mu) {
        if(length(n) > 1) n <- length(n)
        if(any(mu==0)) mu[mu==0] <- min(mu[mu!=0])
        r <- rep(0, n)
        for(i in 1:n) {
          while(r[i]==0) r[i] <- rpois(1, mu[i])
        }
        return(r)
      }  
      rHP <- function(n, mu, p) {
        Z <- rbinom(n, 1, p)
        P <- rtruncpois(n, mu)
        r <- Z * 0 + (1 - Z) * P
        return(r)
      }
      pr.zero <- predict(object, type="zero")
      pr.zero[pr.zero > 1] <- 1
      pr.zero[pr.zero < 0] <- 0
      nsim <- 1
      while(nsim < (sim+1)) {
        y. <- lapply(rep(length(object$fit), sim), rHP, predict(object, type="count"), pr.zero)
        if(is.null(object$offset$count) & is.null(object$offset$zero)) {
          if(halfnormal) {try(for(i in 1:sim) {
            if(verb.sim) cat("Simulation", nsim + i - 1, "out of", sim, "\n")
            res <- cbind(res,sort(abs(residuals(pscl::hurdle(y.[[i]] ~ X - 1 | Z - 1, link=obj.link, control=pscl::hurdle.control(maxit=maxit)), type=resid.type))))}, silent=TRUE)
          } else {try(for(i in 1:sim) {
            if(verb.sim) cat("Simulation", nsim + i - 1, "out of", sim, "\n")  
            res <- cbind(res,sort(residuals(pscl::hurdle(y.[[i]] ~ X - 1 | Z - 1, link=obj.link, control=pscl::hurdle.control(maxit=maxit)), type=resid.type)))}, silent=TRUE)}
        } else {
          data <- eval(eval(object)$call$data)
          fmla <- as.formula(paste("y.[[i]] ~", object$formula[3]))
          if(length(grep("offset", fmla)) > 0) {
            if(halfnormal) {try(for(i in 1:sim) {
              if(verb.sim) cat("Simulation", nsim + i - 1, "out of", sim, "\n")
              res <- cbind(res,sort(abs(residuals(pscl::hurdle(formula=fmla, link=obj.link, control=pscl::hurdle.control(maxit=maxit), data=data), type=resid.type))))}, silent=TRUE)
            } else {try(for(i in 1:sim) {
              if(verb.sim) cat("Simulation", nsim + i - 1, "out of", sim, "\n")  
              res <- cbind(res,sort(residuals(pscl::hurdle(formula=fmla, link=obj.link, control=pscl::hurdle.control(maxit=maxit), data=data), type=resid.type)))}, silent=TRUE)}
          } else {
            obj.offset <- object$offset$count
            if(halfnormal) {try(for(i in 1:sim) {
              if(verb.sim) cat("Simulation", nsim + i - 1, "out of", sim, "\n")
              res <- cbind(res,sort(abs(residuals(pscl::hurdle(formula=fmla, link=obj.link, offset=obj.offset, control=pscl::hurdle.control(maxit=maxit), data=data), type=resid.type))))}, silent=TRUE)
            } else {try(for(i in 1:sim) {
              if(verb.sim) cat("Simulation", nsim + i - 1, "out of", sim, "\n")  
              res <- cbind(res,sort(residuals(pscl::hurdle(formula=fmla, link=obj.link, offset=obj.offset, control=pscl::hurdle.control(maxit=maxit), data=data), type=resid.type)))}, silent=TRUE)}
          }
        }
        nsim <- ncol(cbind(res))
      }
      res <- res[,1:(sim+1)]
    }
    ## hurdle negative binomial
    if(fam=="HNB") {
      cat("Hurdle negative binomial model", '\n')
      X <- model.matrix(object, type="count")
      Z <- model.matrix(object, type="zero")
      obj.link <- object$link
      rtruncnegbin <- function(n, mu, theta) {
        if(length(n) > 1) n <- length(n)
        if(any(mu==0)) mu[mu==0] <- min(mu[mu!=0])
        r <- rep(0, n)
        for(i in 1:n) {
          while(r[i]==0) r[i] <- rnegbin(1, mu[i], theta)
        }
        return(r)
      }  
      rHNB <- function(n, mu, theta, p) {
        Z <- rbinom(n, 1, p)
        NB <- rtruncnegbin(n, mu, theta)
        r <- Z * 0 + (1 - Z) * NB
        return(r)
      }
      pr.zero <- predict(object, type="zero")
      pr.zero[pr.zero > 1] <- 1
      pr.zero[pr.zero < 0] <- 0
      nsim <- 1
      while(nsim < (sim+1)) {
        y. <- lapply(rep(length(object$fit), sim), rHNB, predict(object, type="count"), object$theta, pr.zero)
        if(is.null(object$offset$count) & is.null(object$offset$zero)) {
          if(halfnormal) {try(for(i in 1:sim) {
            if(verb.sim) cat("Simulation", nsim + i - 1, "out of", sim, "\n")
            res <- cbind(res,sort(abs(residuals(pscl::hurdle(y.[[i]] ~ X - 1 | Z - 1, dist="negbin", link=obj.link, control=pscl::hurdle.control(maxit=maxit)), type=resid.type))))}, silent=TRUE)
          } else {try(for(i in 1:sim) {
            if(verb.sim) cat("Simulation", nsim + i - 1, "out of", sim, "\n")  
            res <- cbind(res,sort(residuals(pscl::hurdle(y.[[i]] ~ X - 1 | Z - 1, dist="negbin", link=obj.link, control=pscl::hurdle.control(maxit=maxit)), type=resid.type)))}, silent=TRUE)}
        } else {
          data <- eval(eval(object)$call$data)
          fmla <- as.formula(paste("y.[[i]] ~", object$formula[3]))
          if(length(grep("offset", fmla)) > 0) {
            if(halfnormal) {try(for(i in 1:sim) {
              if(verb.sim) cat("Simulation", nsim + i - 1, "out of", sim, "\n")
              res <- cbind(res,sort(abs(residuals(pscl::hurdle(formula=fmla, dist="negbin", link=obj.link, control=pscl::hurdle.control(maxit=maxit), data=data), type=resid.type))))}, silent=TRUE)
            } else {try(for(i in 1:sim) {
              if(verb.sim) cat("Simulation", nsim + i - 1, "out of", sim, "\n")  
              res <- cbind(res,sort(residuals(pscl::hurdle(formula=fmla, dist="negbin", link=obj.link, control=pscl::hurdle.control(maxit=maxit), data=data), type=resid.type)))}, silent=TRUE)}
          } else {
            obj.offset <- object$offset$count
            if(halfnormal) {try(for(i in 1:sim) {
              if(verb.sim) cat("Simulation", nsim + i - 1, "out of", sim, "\n")
              res <- cbind(res,sort(abs(residuals(pscl::hurdle(formula=fmla, dist="negbin", link=obj.link, offset=obj.offset, control=pscl::hurdle.control(maxit=maxit), data=data), type=resid.type))))}, silent=TRUE)
            } else {try(for(i in 1:sim) {
              if(verb.sim) cat("Simulation", nsim + i - 1, "out of", sim, "\n")  
              res <- cbind(res,sort(residuals(pscl::hurdle(formula=fmla, dist="negbin", link=obj.link, offset=obj.offset, control=pscl::hurdle.control(maxit=maxit), data=data), type=resid.type)))}, silent=TRUE)}
          }
        }
        nsim <- ncol(cbind(res))
      }
      res <- res[,1:(sim+1)]
    }
    
    # now run .makehnp
    .makehnp(obj=res, conf=conf, halfnormal=halfnormal, how.many.out=how.many.out, 
            paint.out=paint.out, col.paint.out=col.paint.out, print.on=print.on, plot.sim=plot.sim, ...)
  }

#setMethod("hnp", "hurdle", hnp.hurdle)

hnp.integer <-
  function(object, sim=99, conf=.95, resid.type, maxit,
           halfnormal=T, scale=F, plot.sim=T, verb.sim=F, warn=F,
           how.many.out=F, print.on=F, paint.out=F, col.paint.out, ...) {
    
    hnp.numeric(object=object, sim=sim, conf=conf, resid.type=resid.type, maxit=maxit,
                halfnormal=halfnormal, scale=scale, plot.sim=plot.sim, verb.sim=verb.sim, warn=warn,
                how.many.out=how.many.out, print.on=print.on, paint.out=paint.out, col.paint.out=col.paint.out, ...)
  }

#setMethod("hnp", "integer", hnp.integer)

hnp.lm <-
  function(object, sim=99, conf=.95, resid.type, maxit,
           halfnormal=T, scale=F, plot.sim=T, verb.sim=F, warn=F,
           how.many.out=F, print.on=F, paint.out=F, col.paint.out, ...) {
    
    # preparation and extraction of residuals
    if(missing(resid.type)) resid.type="student"
    if(resid.type=="partial") stop("resid.type should be one of 'deviance', 'pearson', 'response', 'working', 'student', 'standard'")
    get.residuals <- function(obj, type="deviance") {
      if(type=="student") {
        rstudent(obj)
      } else if(type=="standard") {
        rstandard(obj)
      } else resid(obj, type=type)
    }  
    if(halfnormal) {res <- sort(abs(get.residuals(object, type=resid.type)))
    } else {res <- sort(get.residuals(object, type=resid.type))}
    
    # producing the envelope bands
    cat("Gaussian model (lm object)", '\n')
    X <- model.matrix(object)
    dp <- sqrt(anova(object)$"Mean Sq"[length(anova(object)$"Mean Sq")])
    y. <- lapply(rep(length(object$fit), sim), rnorm, object$fit, dp)
    if(is.null(object$offset)) {
      if(halfnormal) {for(i in 1:sim) {
        if(verb.sim) cat("Simulation", i, "out of", sim, "\n")
        res <- cbind(res,sort(abs(get.residuals(lm(y.[[i]] ~ X - 1), type=resid.type))))}
      } else {for(i in 1:sim) {
        if(verb.sim) cat("Simulation", i, "out of", sim, "\n")  
        res <- cbind(res,sort(get.residuals(lm(y.[[i]] ~ X - 1), type=resid.type)))}}
    } else {
      data <- eval(eval(object)$call$data)
      fmla <- as.formula(paste("y.[[i]] ~", paste(as.formula(paste(object$call)[2]))[3]))
      if(halfnormal) {for(i in 1:sim) {
        if(verb.sim) cat("Simulation", i, "out of", sim, "\n")
        res <- cbind(res,sort(abs(get.residuals(lm(formula=fmla, data=data), type=resid.type))))}
      } else {for(i in 1:sim) {
        if(verb.sim) cat("Simulation", i, "out of", sim, "\n")  
        res <- cbind(res,sort(get.residuals(lm(formula=fmla, data=data), type=resid.type)))}}
    }
    
    # now run .makehnp
    .makehnp(obj=res, conf=conf, halfnormal=halfnormal, how.many.out=how.many.out, 
            paint.out=paint.out, col.paint.out=col.paint.out, print.on=print.on, plot.sim=plot.sim, ...)
  }

#setMethod("hnp", "lm", hnp.lm)

hnp.lmerMod <-
  function(object, sim=99, conf=.95, resid.type, maxit,
           halfnormal=T, scale=F, plot.sim=T, verb.sim=F, warn=F,
           how.many.out=F, print.on=F, paint.out=F, col.paint.out, ...) {
    
    # preparation and extraction of residuals
    if(missing(maxit)) maxit <- 300
    if(missing(resid.type)) resid.type <- "deviance"
    if(resid.type=="partial") stop("resid.type should be one of 'deviance', 'pearson', 'response', 'working'")
    if(halfnormal) {res <- sort(abs(residuals(object, type=resid.type)))
    } else {res <- sort(residuals(object, type=resid.type))}
    
    # producing the envelope bands
    cat("Linear mixed-effects model (using lme4)", '\n')
    data <- eval(eval(object)@call$data)
    rlmm <- function(x) simulate(object)[,1]
    y. <- lapply(1:sim, rlmm)
    fmla <- as.formula(paste("y.[[i]] ~", paste(as.formula(paste(object@call)[2]))[3]))
    if(halfnormal) {for(i in 1:sim) {
      if(verb.sim) cat("Simulation", i, "out of", sim, "\n")
      res <- cbind(res,sort(abs(residuals(lme4::lmer(fmla, control=lme4::lmerControl(optCtrl=list(maxfun=maxit)), data=data), type=resid.type))))}
    } else {for(i in 1:sim) {
      if(verb.sim) cat("Simulation", i, "out of", sim, "\n")
      res <- cbind(res,sort(residuals(lme4::lmer(fmla, control=lme4::lmerControl(optCtrl=list(maxfun=maxit)), data=data), type=resid.type)))}}
    
    # now run .makehnp
    .makehnp(obj=res, conf=conf, halfnormal=halfnormal, how.many.out=how.many.out, 
            paint.out=paint.out, col.paint.out=col.paint.out, print.on=print.on, plot.sim=plot.sim, ...)
  }

#setMethod("hnp", "lmerMod", hnp.lmerMod)

hnp.multinom <-
  function(object, sim=99, conf=.95, resid.type, maxit,
           halfnormal=T, scale=F, plot.sim=T, verb.sim=F, warn=F,
           how.many.out=F, print.on=F, paint.out=F, col.paint.out, ...) {
    
    # preparation and extraction of residuals
    if(!missing(maxit)) warning("multinom function does not let users modify maxit")
    if(!missing(resid.type)) warning("resid.type not implemented for multinomial models")
    fam <- "multinomial"
    if(halfnormal) {res <- sort(abs(as.numeric(resid(object))))
    } else {res <- sort(as.numeric(resid(object)))}
    
    # producing the envelope bands
    cat("Multinomial model", '\n')
    X <- model.matrix(object)
    ncat <- length(object$weights)
    FN <- function(n, obj) {
      a <- list()
      for(i in 1:ncat) a[[i]] <- rmultinom(n, weights(obj)[i], fitted(obj)[i,])
      return(unlist(a))
    }
    y. <- lapply(rep(1, sim), FN, object)
    if(is.null(attr(object$terms, "offset"))) {
      if(halfnormal) {for(i in 1:sim) {
        if(verb.sim) cat("Simulation", i, "out of", sim, "\n")
        res <- cbind(res,sort(abs(as.numeric(resid(nnet::multinom(matrix(y.[[i]], nrow=ncat, byrow=T) ~ X - 1))))))}
      } else {for(i in 1:sim) {
        if(verb.sim) cat("Simulation", i, "out of", sim, "\n")  
        res <- cbind(res,sort(as.numeric(resid(nnet::multinom(matrix(y.[[i]], nrow=ncat, byrow=F) ~ X - 1)))))}}
    } else {
      data <- eval(eval(object)$call$data)
      fmla <- as.formula(paste("matrix(y.[[i]], nrow=ncat, byrow=F) ~", as.formula(paste(object$call[2]))[3]))
      if(halfnormal) {for(i in 1:sim) {
        if(verb.sim) cat("Simulation", i, "out of", sim, "\n")
        res <- cbind(res,sort(abs(as.numeric(resid(nnet::multinom(formula=fmla, data=data))))))}
      } else {for(i in 1:sim) {
        if(verb.sim) cat("Simulation", i, "out of", sim, "\n")  
        res <- cbind(res,sort(as.numeric(resid(nnet::multinom(formula=fmla, data=data)))))}}
    }
    
    # now run .makehnp
    .makehnp(obj=res, conf=conf, halfnormal=halfnormal, how.many.out=how.many.out, 
            paint.out=paint.out, col.paint.out=col.paint.out, print.on=print.on, plot.sim=plot.sim, ...)
  }

#setMethod("hnp", "multinom", hnp.multinom)

hnp.negbin <-
  function(object, sim=99, conf=.95, resid.type, maxit,
           halfnormal=T, scale=F, plot.sim=T, verb.sim=F, warn=F,
           how.many.out=F, print.on=F, paint.out=F, col.paint.out, ...) {
    
    # preparation and extraction of residuals
    if(missing(maxit)) maxit <- 25
    if(missing(resid.type)) resid.type <- "deviance"
    if(resid.type=="partial") stop("resid.type should be one of 'deviance', 'pearson', 'response', 'working', 'student', 'standard'")
    wg <- weights(object)
    get.residuals <- function(obj, type="deviance") {
      if(type=="student") {
        rstudent(obj)
      } else if(type=="standard") {
        rstandard(obj)
      } else resid(obj, type=type)
    }  
    if(halfnormal) {res <- sort(abs(get.residuals(object, type=resid.type)))
    } else {res <- sort(get.residuals(object, type=resid.type))}
    
    # producing the envelope bands
    cat("Negative binomial model (using MASS package)", '\n')
    X <- model.matrix(object)
    nsim <- 1
    while(nsim < (sim+1)) {
      y. <- lapply(rep(length(object$fit), sim), rnegbin, object$fit, object$theta)
      if(is.null(object$offset)) {
        if(halfnormal) {try(for(i in 1:sim) {
          if(verb.sim) cat("Simulation", nsim + i - 1, "out of", sim, "\n")
          res <- cbind(res,sort(abs(get.residuals(glm.nb(y.[[i]] ~ X - 1, weights=wg, control=glm.control(maxit=maxit)), type=resid.type))))}, silent=TRUE)
        } else {try(for(i in 1:sim) {
          if(verb.sim) cat("Simulation", nsim + i - 1, "out of", sim, "\n")  
          res <- cbind(res,sort(get.residuals(glm.nb(y.[[i]] ~ X - 1, weights=wg, control=glm.control(maxit=maxit)), type=resid.type)))}, silent=TRUE)}
      } else {
        data <- eval(eval(object)$call$data)
        fmla <- as.formula(paste("y.[[i]] ~", paste(as.formula(paste(object$call)[2]))[3]))
        if(halfnormal) {try(for(i in 1:sim) {
          if(verb.sim) cat("Simulation", nsim + i - 1, "out of", sim, "\n")
          res <- cbind(res,sort(abs(get.residuals(glm.nb(formula=fmla, weights=wg, control=glm.control(maxit=maxit), data=data), type=resid.type))))}, silent=TRUE)
        } else {try(for(i in 1:sim) {
          if(verb.sim) cat("Simulation", nsim + i - 1, "out of", sim, "\n")  
          res <- cbind(res,sort(get.residuals(glm.nb(formula=fmla, weights=wg, control=glm.control(maxit=maxit), data=data), type=resid.type)))}, silent=TRUE)}
      }
      nsim <- ncol(cbind(res))
    }
    res <- res[,1:(sim+1)]
    
    # now run .makehnp
    .makehnp(obj=res, conf=conf, halfnormal=halfnormal, how.many.out=how.many.out, 
            paint.out=paint.out, col.paint.out=col.paint.out, print.on=print.on, plot.sim=plot.sim, ...)
  }

#setMethod("hnp", "negbin", hnp.negbin)

hnp.numeric <-
  function(object, sim=99, conf=.95, resid.type, maxit,
           halfnormal=T, scale=F, plot.sim=T, verb.sim=F, warn=F,
           how.many.out=F, print.on=F, paint.out=F, col.paint.out, ...) {
    
    # producing envelope bands
    cat("Half-normal plot with simulated envelope generated assuming the residuals are 
        normally distributed under the null hypothesis.", "\n")
    if(scale) {
      ex <- mean(object)
      sdx <- sd(object)
      cat("Estimated mean:", ex, "\n")
      cat("Estimated variance:", sdx^2, "\n")
    }
    FN1 <- function(x) return(sort(abs(rnorm(x, ex, sdx))))
    FN2 <- function(x) return(sort(abs(rnorm(x))))
    FN3 <- function(x) return(sort(rnorm(x, ex, sdx)))
    FN4 <- function(x) return(sort(rnorm(x)))
    if(halfnormal) {
      res1 <- sort(abs(object)) 
      if(scale) {
        res <- cbind(res1, sapply(rep(length(res1), sim), FN1))
      } else {
        res <- cbind(res1, sapply(rep(length(res1), sim), FN2))}
    } else {
      res1 <- sort(object)
      if(scale) {
        res <- cbind(res1, sapply(rep(length(res1), sim), FN3))
      } else {
        res <- cbind(res1, sapply(rep(length(res1), sim), FN4))}
    }
    
    # now run .makehnp
    .makehnp(obj=res, conf=conf, halfnormal=halfnormal, how.many.out=how.many.out, 
            paint.out=paint.out, col.paint.out=col.paint.out, print.on=print.on, plot.sim=plot.sim, ...)
  }

#setMethod("hnp", "numeric", hnp.numeric)

hnp.vglm <-
  function(object, sim=99, conf=.95, resid.type, maxit,
           halfnormal=T, scale=F, plot.sim=T, verb.sim=F, warn=F,
           how.many.out=F, print.on=F, paint.out=F, col.paint.out, ...) {
    
    # preparation and extraction of residuals
    if(missing(maxit)) maxit <- 25
    if(missing(resid.type)) resid.type <- "response"
    if(resid.type!="response") stop("Only resid.type='response' residuals implemented")
    if(object@family@vfamily=="betabinomial") fam <- "beta-binomial_VGAM"
    if(object@family@vfamily=="zibinomial") fam <- "ZIB"
    if(object@family@vfamily[1]!="betabinomial"&object@family@vfamily[1]!="zibinomial") stop("This function has been implemented for VGLM objects with family='betabinomial' and 'zibinomial' only")
    if(sum(grep("zero", object@call, invert=T)==3)==1) {
      if(fam=="beta-binomial_VGAM") zero <- 2 else zero <- NULL
    } else {
      if(length(grep(2, object@call[3]))==1) {
        zero <- 2
      } else {
        if(length(grep(1, object@call[3]))==1) {
          zero <- 1
        } else {
          zero <- NULL
        }
      }}
    if(halfnormal) {res <- sort(abs(VGAM::resid(object, type=resid.type)))
    } else {res <- sort(VGAM::resid(object, type=resid.type))}
    
    # producing the envelope bands
    ## beta-binomial
    if(fam=="beta-binomial_VGAM") {
      cat("Beta-binomial model (using VGAM package)", '\n')
      data <- eval(eval(object)@call$data)
      m. <- object@prior.weights
      f.bb <- object@fitted.values
      rho <- object@misc$rho
      y. <- lapply(rep(length(f.bb), sim), VGAM::rbetabinom, m., f.bb, rho)
      fmla <- as.formula(paste("cbind(y.[[i]],m.-y.[[i]]) ~ ", paste(object@terms$terms)[3]))
      if(halfnormal) {for(i in 1:sim) {
        if(verb.sim) cat("Simulation", i, "out of", sim, "\n")
        res <- cbind(res, sort(abs(resid(VGAM::vglm(formula=fmla, family=VGAM::betabinomial(zero=zero), control=VGAM::vglm.control(maxit=maxit), data=data, subset=eval(eval(object)@call$subset)), type=resid.type))))}
      } else {for(i in 1:sim) {
        if(verb.sim) cat("Simulation", i, "out of", sim, "\n")  
        res <- cbind(res, sort(resid(VGAM::vglm(formula=fmla, family=VGAM::betabinomial(zero=zero), control=VGAM::vglm.control(maxit=maxit), data=data, subset=eval(eval(object)@call$subset)), type=resid.type)))}}
    }  
    ## zero-inflated binomial
    if(fam=="ZIB") {
      cat("Zero-inflated binomial model (using VGAM)", '\n')
      data <- eval(eval(object)@call$data)
      m. <- object@prior.weights
      f.zib <- object@fitted.values
      phi <- object@misc$pstr0
      y. <- lapply(rep(length(f.zib), sim), VGAM::rzibinom, m., f.zib, phi)
      fmla <- as.formula(paste("cbind(y.[[i]],m.-y.[[i]]) ~ ", paste(object@terms$terms)[3]))
      if(halfnormal) {for(i in 1:sim) {
        if(verb.sim) cat("Simulation", i, "out of", sim, "\n") 
        res <- cbind(res, sort(abs(resid(VGAM::vglm(formula=fmla, family=VGAM::zibinomial(zero=zero), control=VGAM::vglm.control(maxit=maxit), data=data, subset=eval(eval(object)@call$subset)), type=resid.type))))}
      } else {for(i in 1:sim) {
        if(verb.sim) cat("Simulation", i, "out of", sim, "\n")  
        res <- cbind(res, sort(resid(VGAM::vglm(formula=fmla, family=VGAM::zibinomial(zero=zero), control=VGAM::vglm.control(maxit=maxit), data=data, subset=eval(eval(object)@call$subset)), type=resid.type)))}}
    }  
    
    # now run .makehnp
    .makehnp(obj=res, conf=conf, halfnormal=halfnormal, how.many.out=how.many.out, 
            paint.out=paint.out, col.paint.out=col.paint.out, print.on=print.on, plot.sim=plot.sim, ...)
  }

#setMethod("hnp", "vglm", hnp.vglm)

hnp.zeroinfl <-
  function(object, sim=99, conf=.95, resid.type, maxit,
           halfnormal=T, scale=F, plot.sim=T, verb.sim=F, warn=F,
           how.many.out=F, print.on=F, paint.out=F, col.paint.out, ...) {
    
    # preparation and extraction of residuals
    ## zero-inflated Poisson
    if(missing(maxit)) maxit <- 10000
    if(missing(resid.type)) resid.type <- "pearson"
    if(resid.type!="pearson"&resid.type!="response") stop("resid.type should be one of 'pearson', 'response'")
    if(object$dist=="poisson") fam <- "ZIP"
    if(object$dist=="negbin") fam <- "ZINB"
    if(object$dist=="geometric") stop("This function has not been implemented yet for zero-inflated geometric models")
    if(halfnormal) {res <- sort(abs(residuals(object, type=resid.type)))
    } else {res <- sort(residuals(object, type=resid.type))}
    
    # producing the envelope bands
    if(fam=="ZIP") {
      cat("Zero-inflated Poisson model", '\n')
      X <- model.matrix(object, type="count")
      Z <- model.matrix(object, type="zero")
      obj.link <- object$link
      rZIP <- function(n, mu, p) {
        Z <- rbinom(n, 1, p)
        P <- rpois(n, mu)
        r <- Z * 0 + (1 - Z) * P
        return(r)
      }
      nsim <- 1
      while(nsim < (sim+1)) {
        y. <- lapply(rep(length(object$fit), sim), rZIP, predict(object, type="count"), predict(object, type="zero"))
        if(is.null(object$offset$count) & is.null(object$offset$zero)) {
          if(halfnormal) {try(for(i in 1:sim) {
            if(verb.sim) cat("Simulation", nsim + i - 1, "out of", sim, "\n")
            res <- cbind(res,sort(abs(residuals(pscl::zeroinfl(y.[[i]] ~ X - 1 | Z - 1, link=obj.link, control=pscl::zeroinfl.control(maxit=maxit)), type=resid.type))))}, silent=TRUE)
          } else {try(for(i in 1:sim) {
            if(verb.sim) cat("Simulation", nsim + i - 1, "out of", sim, "\n")  
            res <- cbind(res,sort(residuals(pscl::zeroinfl(y.[[i]] ~ X - 1 | Z - 1, link=obj.link, control=pscl::zeroinfl.control(maxit=maxit)), type=resid.type)))}, silent=TRUE)}
        } else {
          data <- eval(eval(object)$call$data)
          fmla <- as.formula(paste("y.[[i]] ~", object$formula[3]))
          if(length(grep("offset", fmla)) > 0) {
            if(halfnormal) {try(for(i in 1:sim) {
              if(verb.sim) cat("Simulation", nsim + i - 1, "out of", sim, "\n")
              res <- cbind(res,sort(abs(residuals(pscl::zeroinfl(formula=fmla, link=obj.link, control=pscl::zeroinfl.control(maxit=maxit), data=data), type=resid.type))))}, silent=TRUE)
            } else {try(for(i in 1:sim) {
              if(verb.sim) cat("Simulation", nsim + i - 1, "out of", sim, "\n")  
              res <- cbind(res,sort(residuals(pscl::zeroinfl(formula=fmla, link=obj.link, control=pscl::zeroinfl.control(maxit=maxit), data=data), type=resid.type)))}, silent=TRUE)}
          } else {
            obj.offset <- object$offset$count
            if(halfnormal) {try(for(i in 1:sim) {
              if(verb.sim) cat("Simulation", nsim + i - 1, "out of", sim, "\n")
              res <- cbind(res,sort(abs(residuals(pscl::zeroinfl(formula=fmla, link=obj.link, offset=obj.offset, control=pscl::zeroinfl.control(maxit=maxit), data=data), type=resid.type))))}, silent=TRUE)
            } else {try(for(i in 1:sim) {
              if(verb.sim) cat("Simulation", nsim + i - 1, "out of", sim, "\n")  
              res <- cbind(res,sort(residuals(pscl::zeroinfl(formula=fmla, link=obj.link, offset=obj.offset, control=pscl::zeroinfl.control(maxit=maxit), data=data), type=resid.type)))}, silent=TRUE)}
          }
        }
        nsim <- ncol(cbind(res))
      }
      res <- res[,1:(sim+1)]
    }
    ## zero-inflated negative binomial
    if(fam=="ZINB") {
      cat("Zero-inflated negative binomial model", '\n')
      X <- model.matrix(object, type="count")
      Z <- model.matrix(object, type="zero")
      obj.link <- object$link
      rZINB <- function(n, mu, theta, p) {
        Z <- rbinom(n, 1, p)
        NB <- rnegbin(n, mu, theta)
        r <- Z * 0 + (1 - Z) * NB
        return(r)
      }
      nsim <- 1
      while(nsim < (sim+1)) {
        y. <- lapply(rep(length(object$fit), sim), rZINB, predict(object, type="count"), object$theta, predict(object, type="zero"))
        if(is.null(object$offset$count) & is.null(object$offset$zero)) {
          if(halfnormal) {try(for(i in 1:sim) {
            if(verb.sim) cat("Simulation", nsim + i - 1, "out of", sim, "\n")
            res <- cbind(res,sort(abs(residuals(pscl::zeroinfl(y.[[i]] ~ X - 1 | Z - 1, dist="negbin", link=obj.link, control=pscl::zeroinfl.control(maxit=maxit)), type=resid.type))))}, silent=TRUE)
          } else {try(for(i in 1:sim) {
            if(verb.sim) cat("Simulation", nsim + i - 1, "out of", sim, "\n")  
            res <- cbind(res,sort(residuals(pscl::zeroinfl(y.[[i]] ~ X - 1 | Z - 1, dist="negbin", link=obj.link, control=pscl::zeroinfl.control(maxit=maxit)), type=resid.type)))}, silent=TRUE)}
        } else {
          data <- eval(eval(object)$call$data)
          fmla <- as.formula(paste("y.[[i]] ~", object$formula[3]))
          if(length(grep("offset", fmla)) > 0) {
            if(halfnormal) {try(for(i in 1:sim) {
              if(verb.sim) cat("Simulation", nsim + i - 1, "out of", sim, "\n")
              res <- cbind(res,sort(abs(residuals(pscl::zeroinfl(formula=fmla, dist="negbin", link=obj.link, control=pscl::zeroinfl.control(maxit=maxit), data=data), type=resid.type))))}, silent=TRUE)
            } else {try(for(i in 1:sim) {
              if(verb.sim) cat("Simulation", nsim + i - 1, "out of", sim, "\n")  
              res <- cbind(res,sort(residuals(pscl::zeroinfl(formula=fmla, dist="negbin", link=obj.link, control=pscl::zeroinfl.control(maxit=maxit), data=data), type=resid.type)))}, silent=TRUE)}
          } else {
            obj.offset <- object$offset$count
            if(halfnormal) {try(for(i in 1:sim) {
              if(verb.sim) cat("Simulation", nsim + i - 1, "out of", sim, "\n")
              res <- cbind(res,sort(abs(residuals(pscl::zeroinfl(formula=fmla, offset=obj.offset, dist="negbin", link=obj.link, control=pscl::zeroinfl.control(maxit=maxit), data=data), type=resid.type))))}, silent=TRUE)
            } else {try(for(i in 1:sim) {
              if(verb.sim) cat("Simulation", nsim + i - 1, "out of", sim, "\n")  
              res <- cbind(res,sort(residuals(pscl::zeroinfl(formula=fmla, offset=obj.offset, dist="negbin", link=obj.link, control=pscl::zeroinfl.control(maxit=maxit), data=data), type=resid.type)))}, silent=TRUE)}
          }
        }
        nsim <- ncol(cbind(res))
      }
      res <- res[,1:(sim+1)]
    }
    
    # now run .makehnp
    .makehnp(obj=res, conf=conf, halfnormal=halfnormal, how.many.out=how.many.out, 
            paint.out=paint.out, col.paint.out=col.paint.out, print.on=print.on, plot.sim=plot.sim, ...)
  }

#setMethod("hnp", "zeroinfl", hnp.zeroinfl)

hnp.default <-
  function(object, ...) {
    stop("This function has not been implemented for objects of class '", class(object)[1], 
         "'. If you wish to supply your own fitting, simulation and diagnostic extration codes, 
       see ?hnp for details.")
  }