svygnm<-function(formula, design, ...){
#   .svycheck(design)
  UseMethod("svygnm",design)
}

svygnm.svyrep.design<-function(formula, design, subset=NULL, data.fun=NULL,
                              rho=NULL, return.replicates=FALSE, keep.weights=FALSE,
                              na.action, eliminate, ncpus=getOption("boot.ncpus"), ...){
      if(!inherits(design, "svyrep.design"))
        stop("'design' must be an object of class \"svyrep.design\".")

  subset<-substitute(subset)
  subset<-eval(subset, design$variables, parent.frame())
  if (!is.null(subset))
    design<-design[subset,]

  if(!is.null(data.fun))
    data <- data.fun(design, ...)
  else
    data<-design$variables

  g<-match.call()
  formula<-eval.parent(formula)
  environment(formula)<-environment()
  g$formula<-formula
  g$data<-quote(data)
  g$design<-NULL
  g$var<-g$rho<-g$return.replicates<-NULL
  g[[1]]<-quote(gnm)      
  g$model<-TRUE
  g$x<-TRUE
  g$y<-TRUE

      args <- list(formula=formula, data=quote(data), eliminate=substitute(eliminate),
                   model=TRUE, x=TRUE, y=TRUE, ...)
  
      scale<-design$scale
      rscales<-design$rscales
      if (!is.null(rho)) .NotYetUsed(rho)

      pwts<-design$pweights/sum(design$pweights)
      if (is.data.frame(pwts)) pwts<-pwts[[1]]

      if(is.null(data.fun)) {
        args$weights<-pwts
      }
      
      if (!all(all.vars(formula) %in% names(data))) 
        stop("all variables must be in design= argument")

      full<-do.call(gnm::gnm, args)
  
      if(!is.null(args$method) && args$method != "gnmFit")
          return(full)

      full$naive.cov<-vcov(full)
  
      nas<-attr(full$model, "na.action")

      if(getOption("survey.drop.replicates") && !is.null(design$selfrep) && all(design$selfrep)){

          v<-matrix(0,ncol=length(coef(full)),nrow=length(coef(full)))
          betas<-NULL

      } else {
          betas<-matrix(ncol=length(coef(full)),
                        nrow=ncol(design$repweights))
          
          if (!design$combined.weights)
              pw1<-pwts
          else
              pw1<-rep(1,length(pwts))
          wts<-design$repweights
          if (length(nas)){
              wts<-wts[-nas,]
              pw1<-pw1[-nas]
          }
          beta0<-gnm::parameters(full)

          args$start <- beta0
          args$etastart <- as.numeric(predict(full))
          args$offset <- full$offset
          args$trace <- FALSE
          args$verbose <- FALSE

          if(is.null(ncpus))
            ncpus <- if(requireNamespace("parallel")) min(parallel::detectCores(), 5)
                     else 1

          if (ncpus > 1 && requireNamespace("parallel")){
            cl <- parallel::makePSOCKcluster(rep("localhost", ncpus), outfile="", methods=FALSE)
            on.exit(parallel::stopCluster(cl))

            libpaths <- .libPaths()
            parallel::clusterExport(cl, "libpaths", env=environment())
            # Required for some formulas
            parallel::clusterEvalQ(cl, library(gnm, lib.loc=libpaths,
                                               warn.conflicts=FALSE, quietly=TRUE, verbose=FALSE))
            # Required for to acces wts via `[.repweights_compressed`
            parallel::clusterEvalQ(cl, library(survey, lib.loc=libpaths,
                                               warn.conflicts=FALSE, quietly=TRUE, verbose=FALSE))
  
            betas<-do.call(rbind, parallel::parLapply(cl, 1:ncol(wts), function(i, design, wts, pw1, args, data.fun, ...) {
              cat(".")
              wi<-as.vector(wts[,i])*pw1

              if(is.null(data.fun)) {
                args$weights <- wi/sum(wi)
              }
              else {
                design$pweights <- wi
                args$data <- data.fun(design, ...)
              }

              # FIXME: gnm does not appear to like 0 weights in some models (Dref) when eliminate is set
              # Investigate why...
              args$data <- subset(data, args$weights > 0)
              args$offset <- subset(args$offset, args$weights > 0)
              args$etastart <- subset(args$etastart, args$weights > 0)
              args$weights <- subset(args$weights, args$weights > 0)

              gnm::parameters(do.call(gnm::gnm, args))
            }, subset(design, TRUE, select=all.vars(formula)), wts, pw1, args, data.fun, ...))
          } else {
            args2 <- args
            if(!is.null(data.fun))
              design2 <- design

            for(i in 1:ncol(wts)){
              cat(".")

              wi<-as.vector(wts[,i])*pw1

              if(is.null(data.fun)) {
                args$weights <- wi/sum(wi)
              }
              else {
                design2$pweights <- wi
                args$data <- data.fun(design2, ...)
              }

              # FIXME: gnm does not appear to like 0 weights in some models (Dref) when eliminate is set
              # Investigate why...
              args2$data <- subset(data, args$weights > 0)
              args2$offset <- subset(args$offset, args$weights > 0)
              args2$etastart <- subset(args$etastart, args$weights > 0)
              args2$weights <- subset(args$weights, args$weights > 0)

              betas[i,]<-gnm::parameters(do.call(gnm::gnm, args2))
            }
          }

          cat("\n")

          colnames(betas) <- names(coef(full))

          v<-survey::svrVar(betas,scale, rscales,mse=design$mse,coef=beta0)
  }

  full$x<-NULL
  full$df.residual<-survey::degf(design)+1-full$rank
  
  if (length(nas))
      design<-design[-nas,]

  full$cov.unscaled<-v
  if (return.replicates){
    attr(betas,"scale")<-design$scale
    attr(betas,"rscales")<-design$rscales
    attr(betas,"mse")<-design$mse
    full$replicates<-betas
  }  
  class(full)<-c("svrepgnm", "svygnm", class(full))
  full$call<-match.call()
  if(!("formula" %in% names(full$call))) {
    if (is.null(names(full$call)))
      i<-1
    else
      i<-min(which(names(full$call)[-1]==""))
    names(full$call)[i+1]<-"formula"
  }

  # Needed by residuals.svrepgnm()
  full$survey.design$pweights<-design$pweights
  presid<-residuals.svrepgnm(full,"pearson")

  if(is.null(data.fun))
    full$dispersion<-sum(design$pweights * presid^2,na.rm=TRUE)/sum(design$pweights)
  else
    full$dispersion<-mean(presid^2, na.rm=TRUE)

  full$survey.design$call<-design$call
  full$survey.design$type<-design$type

  if(!keep.weights)
      full$survey.design$pweights<-NULL

  full
}


print.svygnm<-function(x,...){
  print(x$survey.design, varnames=FALSE, design.summaries=FALSE,...)
  NextMethod()

}

vcov.svygnm<-function(object,...) {
  v<-object$cov.unscaled
  dimnames(v)<-list(names(coef(object)),names(coef(object)))
  v
}

summary.svrepgnm<-function (object, correlation = FALSE, df.resid=NULL,...) 
{
    est.disp <- TRUE
    if (is.null(df.resid))
      df.r <- object$df.residual
    else
      df.r<-df.resid

    coef.p <- gnm::parameters(object)
    covmat<-vcov(object)
    dimnames(covmat) <- list(names(coef.p), names(coef.p))
    var.cf <- diag(covmat)
    s.err <- sqrt(var.cf)
    tvalue <- coef.p/s.err
    dn <- c("Estimate", "Std. Error")
    if (!est.disp) {
        pvalue <- 2 * pnorm(-abs(tvalue))
        coef.table <- cbind(coef.p, s.err, tvalue, pvalue)
        dimnames(coef.table) <- list(names(coef.p), c(dn, "z value", 
            "Pr(>|z|)"))
    }
    else if (df.r > 0) {
        pvalue <- 2 * pt(-abs(tvalue), df.r)
        coef.table <- cbind(coef.p, s.err, tvalue, pvalue)
        dimnames(coef.table) <- list(names(coef.p), c(dn, "t value", 
            "Pr(>|t|)"))
    }
    else {
        coef.table <- cbind(coef.p, Inf)
        dimnames(coef.table) <- list(names(coef.p), dn)
    }
    ans <- c(object[c("call", "terms", "family", "deviance", 
        "aic", "contrasts", "df.residual", "null.deviance", "df.null", 
        "iter")], list(deviance.resid = residuals(object, type = "deviance"), 
        aic = object$aic, coefficients = coef.table, dispersion = object$dispersion, 
        df = c(object$rank, df.r), cov.unscaled = covmat, 
        cov.scaled = covmat))
    if (correlation) {
        dd <- sqrt(diag(covmat))
        ans$correlation <- covmat/outer(dd, dd)
    }

    ans$aliased<-is.na(ans$coef)
    ans$survey.design<-list(call=object$survey.design$call,
                            type=object$survey.design$type)
    class(ans) <- c("summary.svygnm","summary.gnm")
    return(ans)
}


print.summary.svygnm<-function (x, digits = max(3, getOption("digits") - 3),
                                symbolic.cor = x$symbolic.cor, 
    signif.stars = getOption("show.signif.stars"), ...) 
{
  cat("\nCall:\n")
    cat(paste(deparse(x$call), sep = "\n", collapse = "\n"), 
        "\n\n", sep = "") 

    cat("Survey design:\n")
    print(x$survey.design$call)
   
        cat("\nCoefficients:\n")

        printCoefmat(x$coefficients, digits = digits, signif.stars = signif.stars, 
            na.print = "NA", ...)
    
    cat("\n(Dispersion parameter for ", x$family$family, " family taken to be ", 
        format(x$dispersion), ")\n\n",  "Number of iterations: ", 
        x$iter, "\n", sep = "")
    correl <- x$correlation
    if (!is.null(correl)) {
        p <- NCOL(correl)
        if (p > 1) {
            cat("\nCorrelation of Coefficients:\n")
            if (is.logical(symbolic.cor) && symbolic.cor) {
                print(symnum(correl, abbr.colnames = NULL))
            }
            else {
                correl <- format(round(correl, 2), nsmall = 2, 
                  digits = digits)
                correl[!lower.tri(correl)] <- ""
                print(correl[-1, -p, drop = FALSE], quote = FALSE)
            }
        }
    }
    cat("\n")
    invisible(x)
}

residuals.svrepgnm<-function(object,type = c("deviance", "pearson", "working", 
    "response", "partial"),...){
        type<-match.arg(type)
        if (type=="pearson"){
           y <- object$y
           mu <- object$fitted.values

           wts <- object$prior.weights
           if(isTRUE(object$data.fun))
             r<-(y - mu) * sqrt(wts)/(sqrt(object$family$variance(mu)))
           else
             r<-(y - mu) * sqrt(wts)/(sqrt(object$family$variance(mu))*sqrt(object$survey.design$pweights/sum(object$survey.design$pweights)))

           if (is.null(object$na.action))
                r
           else 
                naresid(object$na.action, r)
        } else 
                NextMethod()

}


confint.svygnm<-function(object,parm,level=0.95,method=c("Wald","likelihood"),ddf=Inf,...){
  method<-match.arg(method)
  if(method=="Wald"){
    tlevel<-pt(qnorm(level),df=ddf)
    return(confint.default(object,parm=parm,level=tlevel,...))
  }
  pnames <- names(coef(object))
  if (missing(parm)) 
    parm <- seq_along(pnames)
  else if (is.character(parm))
    parm <- match(parm, pnames, nomatch = 0)
  lambda<-diag(object$cov.unscaled[parm,parm,drop=FALSE])/diag(object$naive.cov[parm,parm,drop=FALSE])
  if(is.null(ddf)) ddf<-object$df.residual
  if (ddf==Inf)
    level<- 1-2*pnorm(qnorm((1-level)/2)*sqrt(lambda))
  else {
    level<- 1-2*pnorm(qt((1-level)/2,df=ddf)*sqrt(lambda))
  }
  rval<-vector("list",length(parm))
  for(i in 1:length(parm)){
    rval[[i]]<-NextMethod(object=object,parm=parm[i],level=level[i],...)
  }
  names(rval)<-pnames[parm]
  if (length(rval)==1)
    rval<-rval[[1]]
  else
    rval<-do.call(rbind,rval)
  attr(rval,"levels")<-level
  rval
}

logLik.svrepgnm<-function(object,...){
   stop("svrepglm not fitted by maximum likelihood.")
}

extractAIC.svrepgnm<-function(fit,...){
    stop("svrepglm not fitted by maximum likelihood")
}

model.matrix.svygnm <- function (object, coef = NULL, ...) {
    if (!"x" %in% names(object) || !is.null(coef)) {
        xcall <- object$call
        xcall$method <- "model.matrix"
        xcall$constrain <- object$constrain
        xcall$constrainTo <- object$constrainTo
        xcall[c("weights", "offset")] <- NULL
        xcall$verbose <- FALSE
        if (!is.null(coef)) 
            xcall$start <- coef
        else xcall$start <- coef(object)
        eval(xcall)
    }
    else object[[match("x", names(object))]]
}
