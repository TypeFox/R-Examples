
## Authors 
## Martin Schlather, schlather@math.uni-mannheim.de
##
##
## Copyright (C) 2015 Martin Schlather
##
## This program is free software; you can redistribute it and/or
## modify it under the terms of the GNU General Public License
## as published by the Free Software Foundation; either version 3
## of the License, or (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.  


print.RFratiotest <- function(x, ...) {  
  if (!is.null(x$simu.ratios)) {
    ## MC ratio test
     cat("\nMonte Carlo likelihood ratio test",
         "\n=================================",
         "\nnull model:", rfConvertRMmodel2string(x$model.list$nullmodel), 
         "\nalt. model:", rfConvertRMmodel2string(x$model.list$alternative),
         "\n",x$msg)
 # } else if (is.null(x$model1.df) && is.null(x$loglik)) {
#    ## wann wird diese if-Anweiung verwendet?
 #     print.default(t(as.matrix(x)))
  } else {
    cat("\nApprox. likelihood ratio test\n=============================\n")
    if (is.null(x$model1.df)) {
      cat("null model: df=", x$df[1], "  loglik=", x$loglik[1], "\n", sep="")
      cat("alt. model: df=", x$df[2], "  loglik=", x$loglik[2], "\n", sep="")
      cat("p=", x$p,"\n\n")
    } else {
      len <- length(x$model1.df)
      for (i in 1:len) {
        if (len > 1) cat("Test", i, "\n")
        cat("Null model:", as.character(x$model1.name)[i], "\n")
        cat("Alt. model:", as.character(x$model2.name)[i], "\n")
        cat(as.character(x$txt[i]))
      }
      cat("\n")
    }
  }
}



mess <- function(alpha, p, df, nullmodelname=0, altmodelname=1:length(df)) {
  if (is.na(p)) return("NA") else p <- formatC(p, digits=4)
  nullmodelname <- if (is.numeric(nullmodelname)) paste("model", nullmodelname)
  else paste("'", nullmodelname, "'", sep="")
  altmodelname <- if (is.numeric(altmodelname)) paste("model", altmodelname)
  else paste("'", altmodelname, "'", sep="")
  if (missing(df)) ## simulated
    paste("The p-value equals", p,
          "and the hypothesis that the two models significantly differ at",
          "the level alpha=", alpha, "is",
          ifelse(p <= alpha, "accepted.", "rejected."),
          "\n")
  else if (missing(alpha)) 
    paste("loglikelihood test: ", altmodelname,
          " against ", nullmodelname, ": p=", p,
          " (df=", df, ")\n", sep="")
  else {
    if (length(alpha) != 1) stop("alpha should be a scalar")
    paste("loglikelihood test: ", altmodelname,
          " against ", nullmodelname,
          ": p=", p,
          " (df=", df, ")",
          ifelse(p <= alpha, "Null hypothesis accepted.",
                 "Null hypothesis rejected."),
          "\n",
          sep="")
  }
  ## paste(p, df) ## to do
}


approx_test <- function(modellist, alpha) {
  n <- length(modellist)
  loglik <- df <- numeric()
  
  for (m in 1:n) {
    if (class(modellist[[m]]) == "RF_fit") {
      df[m] <- modellist[[m]]$number.of.parameters
      loglik[m] <- modellist[[m]]$ml$likelihood
    } else if (class(modellist[[m]]) == "RFfit") {
      df[m] <- modellist[[m]]@number.of.parameters
      loglik[m] <- modellist[[m]]["ml"]@likelihood
    } else stop("wrong class ('", class(modellist), "') in approx_test.")
  }
  if (length(df) > 2 && (!missing(alpha) || length(alpha) > 1))
    stop("'alpha' must be a scalar")
  p <- pchisq(diff(loglik),  diff(df), lower.tail = FALSE)   
  txt <- mess(alpha=alpha, p=p, df=diff(df))
  return(invisible(list(df=df, loglik=loglik, p=p, txt=txt)))
}


approx_test_single <- function(model, method, alpha, modelinfo) {
  if (class(model) == "RF_fit") {
    submodels <- model$submodels
    df <-  model$number.of.parameters
    loglik <- model[[method]]$likelihood
    report <- model$report
    p.proj <- model$p.proj
    v.proj <- model$v.proj
    x.proj <- model$x.proj
    true.tsdim <- model$true.tsdim
    true.vdim <- model$true.vdim
    AIC <- model[[method]]$AIC
    BIC <- model[[method]]$BIC
    number.of.data <-  model$number.of.data
    if (missing(modelinfo)) modelinfo <- model$modelinfo
    fixed <- model$fixed
    fitted.model <- model[[method]]$model
  } else { # "RFfit"     
    submodels <- model@submodels
    df <- model@number.of.parameters
    loglik <- model[method]@likelihood
    report <- model@report
    p.proj <- model@p.proj
    v.proj <- model@v.proj
    x.proj <- model@x.proj
    true.tsdim <- model@true.tsdim
    true.vdim <- model@true.vdim
    AIC <- model[method]@AIC
    BIC <- model[method]@BIC
    number.of.data <- model@number.of.data
    if (missing(modelinfo)) modelinfo <- model@modelinfo
    fixed <- NULL
    fitted.model <- PrepareModel2(model[method])
  }
  
  if (!is.logical(x.proj) && length(x.proj)  != true.tsdim) ## todo x.proj!=NULL streichen
    stop("space-time projection can't be evaluated yet. Please contact author.")
  
  nm <- rownames(modelinfo)

  proj.txt <- (if (length(p.proj) == 0)  "user's model" else
               paste("(", paste(nm[p.proj], collapse=", "),  sep=""))
  if (length(fixed$zero) > 0)
    proj.txt <- paste(proj.txt, ", ", sep="",
                      paste(nm[fixed$zero], "=0", collapse=", ", sep=""))
  if (length(fixed$one) > 0)
    proj.txt <- paste(proj.txt, ", ", sep="",
                      paste(nm[fixed$one], "=0", collapse=", "))
 if (length(p.proj) > 0) proj.txt <- paste(proj.txt, ")", sep="")
  
  
  if (length(submodels) > 0) {
    sub.df <- sub.loglik <- sub.report <- sub.p.proj <-
      sub.proj.txt <- sub.fixed <- NULL
    for (i in 1:length(submodels)) {
      ret <- approx_test_single(submodels[[i]], method, alpha, modelinfo)
      sub.df <- c(sub.df, ret$df)
      sub.loglik <- c(sub.loglik, ret$loglik)
      sub.report <- c(sub.report, ret$report)
      sub.fixed <- c(if (i > 1) sub.fixed, list(ret$fixed))                     
      sub.p.proj <- c(if (i > 1) sub.p.proj, list(ret$p.proj))    
      sub.v.proj <- c(if (i > 1) sub.v.proj, list(ret$v.proj))
      sub.proj.txt <-  c(sub.proj.txt, ret$proj.txt)
      sub.fitted.models <-
        c(if (i > 1) sub.fitted.models, ret$fitted.model)
    }
    
    len <- length(sub.df)      
    i <- 1
    
    result <- NULL
    result.model <- list()
    result.n <- 0
    while(i <= len) {
      result.n <- result.n + 1
      j <- i
      if (j > len) stop("Error. Please contact author")
      while(j <= len && sub.report[i] == sub.report[j]) {
        j <- j + 1;
      }
      j <- j - 1
      tot.loglik <- sum(sub.loglik[i:j])
      tot.df <- sum(sub.df[i:j])
      tot.proj.txt <- paste(sub.proj.txt[i:j], collapse=" * ")
      tot.p.proj <- sub.p.proj[i:j]
      tot.v.proj <- sub.v.proj[i:j]
      tot.models <- sub.fitted.models[i:j]
      if (length(sub.fixed) == 0) {
        tot.fixed.zero <- tot.fixed.one <- NULL
      } else {
        tot.fixed.zero <- unlist(lapply(sub.fixed[i:j], function(x) x$zero))
        tot.fixed.one <- unlist(lapply(sub.fixed[i:j], function(x) x$one))
      }
      if (true.vdim != length(tot.v.proj[[1]])) {
        aux.model <- list("+")        
        for (k in 1:length(tot.models)) {
          aux.model[[k+1]] <-
            list("M", M=diag(true.vdim)[, tot.v.proj[[k]], drop=FALSE],
                 tot.models[[k]])
        }
        result.model[[result.n]] <- aux.model
       } else {
        result.model[[result.n]] <- tot.models[[1]]
        if (i!=j) stop("model mismatch. Please contact author")
      }
     
      ## oder einfach nur AIC der submodels addieren
      tot.AIC <-  2 * tot.df - 2 * tot.loglik 
      tot.BIC <- log(number.of.data) * tot.df - 2 * tot.loglik
      
      subpproj <- unlist(tot.p.proj)

      if (length(subpproj) == length(unique(subpproj))) {
        delta.df <- df - tot.df
        if (delta.df <= 0) stop("negative df -- please contact author")
        p <-  pchisq(2 * (loglik - tot.loglik), df=delta.df, lower.tail = FALSE)
      } else {
        delta.df <- p <- NA
      }
      result <-
        rbind(result,
              data.frame(model1.name = tot.proj.txt,
                         model1.loglik = tot.loglik,
                         model1.df = tot.df,
                         model1.AIC = tot.AIC,
                         model1.BIC = tot.BIC,
                         model1.zero = paste(nm[tot.fixed.zero], collapse=","),
                         model1.one = paste(nm[tot.fixed.one], collapse=","),
                         model2.name = proj.txt,
                         model2.loglik=loglik,
                         model2.df = df,
                         model2.AIC = AIC,
                         model2.BIC = BIC,
                         model2.zero =  paste(nm[fixed$zero], collapse=","),
                         model2.one = paste(nm[fixed$one], collapse=","),
                         delta.df = delta.df,                      
                         p = p,
                         txt = mess(alpha=alpha, p=p, df=delta.df,
                           nullmodelname="Null model",
                           altmodelname="Alt. model")
                         )
              )
      i <- j + 1
    }

    return(list(df=df, loglik=loglik, report=report,
                p.proj = p.proj, v.proj=v.proj, AIC=AIC, BIC=BIC,
                fixed = fixed,
                fitted.model=c(result.model, list(fitted.model)),
                proj.txt=proj.txt, result = result))
    
  } else {
    return(list(df=df, loglik=loglik, report=report,
                p.proj = p.proj, v.proj=v.proj, AIC=AIC, BIC=BIC,
                fixed = fixed, fitted.model=list(fitted.model),
                proj.txt = proj.txt,               
                result = data.frame(
                  model1.name = "",
                  model1.loglik = NA,
                  model1.df = -1,
                  model1.AIC = NA,
                  model1.BIC = NA,
                  model1.zero = "",
                  model1.one = "",
                  model2.name = proj.txt,
                  model2.loglik=loglik,
                  model2.df = df,
                  model2.AIC = AIC,
                  model2.BIC = BIC,
                  model2.zero = if (is.null(fixed)) "" else
                         paste(nm[fixed$zero], collapse=","),
                  model2.one = if (is.null(fixed)) "" else
                         paste(nm[fixed$one], collapse=","),
                  delta.df = NA,                      
                  p = NA,
                  txt = "NA"
                  )))
  }
}
 


RFratiotest <-
  function(nullmodel, alternative,
           x, y=NULL, z=NULL, T=NULL,  grid=NULL, data,
           alpha,
           n = 5 / alpha, ## number of simulations to do
           seed = 0,
           lower=NULL, upper=NULL, 
           methods, # "reml", "rml1"),
           sub.methods,
           ## "internal" : name should not be changed; should always be last
           ##              method!
           optim.control=NULL,
           users.guess=NULL,  
           distances=NULL, dim,
           transform=NULL,
           ##type = c("Gauss", "BrownResnick", "Smith", "Schlather",
           ##             "Poisson"),
           ... ) {
       
  classes <- c("RF_fit", "RFfit")
    

  RFoptOld <- internal.rfoptions(#general.modus_operandi="normal",
                                  ..., general.seed=NA)
  on.exit(RFoptions(LIST=RFoptOld[[1]]))
  RFopt <- RFoptOld[[2]]
  printlevel <- RFopt$general$printlevel
  if (RFopt$general$modus_operandi == "neurotic")
    stop("crossvalidation is not a precise method")
  

  if ((!RFopt$fit$ratiotest_approx && (missing(alpha) || n < 1 / alpha)) ||
      (!missing(alpha) && (alpha < 0 && alpha > 1))
      )
    stop("alpha is not given or outside [0,1] or to small")
 
  if (class(nullmodel) %in% classes) {
    if (!missing(alternative) && !(class(alternative) %in% classes))
      stop("alternative model not of the class 'RFfit'")
    if (!RFopt$fit$ratiotest_approx)
      stop("for models of class 'RFfit' the parameter 'ratiotest_approx' must be'TRUE'")
    if (missing(alternative)) {
      ats <- approx_test_single(nullmodel, "ml", alpha)$result
      ats <- ats[!is.na(ats$delta.df) ,
                 c("model1.name", "model1.loglik", "model1.df", "model1.zero",
                   "model1.one",
                   "model2.name", "model2.loglik", "model2.df", "model2.zero",
                   "model2.one",
                   "delta.df", "p", "txt"
                   ), drop=FALSE]
      class(ats) <- "RFratiotest"
      if (RFopt$general$returncall)
        attr(ats, "call") <- as.character(deparse(match.call()))
      attr(ats, "coord_system") <- c(orig=RFopt$coords$coord_system,
                                     model=RFopt$coords$new_coord_system)
    return(ats)
    } else {
      ats <- approx_test(list(nullmodel, alternative), alpha)
      class(ats) <- "RFratiotest"
      if (RFopt$general$returncall)
        attr(ats, "call") <- as.character(deparse(match.call()))
      attr(ats, "coord_system") <- c(orig=RFopt$coords$coord_system,
                                     model=RFopt$coords$new_coord_system)
     return(ats)
    }
  } else if (missing(alternative) || (class(alternative) %in% classes))
    stop("alternative model is not given or not of model type")

  
  if (!is.null(seed) && !is.na(seed)) set.seed(seed)
  else if (!is.na(RFopt$general$seed)) {
    if (printlevel >= PL_IMPORTANT)
      message("NOTE: 'RFratiotest' is performed with fixed random seed ",
              RFopt$general$seed,
              ".\nSet RFoptions(seed=NA) to make the seed arbitrary.")
    set.seed(RFopt$general$seed)
  }

  nullmodel <- PrepareModel2(nullmodel, ...)
  alternative <- PrepareModel2(alternative, ...)

  Z <- StandardizeData(x=x, y=y, z=z, T=T, grid=grid, data=data,
                       distances=distances, dim=dim, RFopt=RFopt)
  values <- try(GetValuesAtNA(NAmodel=nullmodel, valuemodel=alternative,
                              spatialdim=Z$spatialdim, Time=Z$Zeit,
                              shortnamelength=3, skipchecks=FALSE),
                silent=TRUE)
  remove("Z")

  isSubmodel <- is.numeric(values) && all(is.na(values))
  if (!isSubmodel && printlevel >= PL_IMPORTANT)
    message("'nullmodel' cannot be automatically detected as being a nullmodel of 'alternative'")
  
  model.list <- list(nullmodel=nullmodel, alternative=alternative)
  data.fit <- list()
  guess <- users.guess

  for (m in 1:length(model.list)) {
    data.fit[[m]] <-
      RFfit(model.list[[m]], x=x, y=y, z=z, T=T, grid=grid, data=data,
            lower=lower, upper=upper,
             methods=methods,
            sub.methods=sub.methods, optim.control=optim.control,
            users.guess=guess,
            distances=distances, dim=dim,
            transform=transform,
            ..., spConform = FALSE)
    guess <- if (isSubmodel) data.fit[[m]]$ml$model else NULL
  }


  if (RFopt$fit$ratiotest_approx) {
    ats <- approx_test(data.fit)
    class(ats) <- "RFratiotest"
    if (RFopt$general$returncall)
      attr(ats, "call") <- as.character(deparse(match.call()))
    attr(ats, "coord_system") <- c(orig=RFopt$coords$coord_system,
                                   model=RFopt$coords$new_coord_system)
   return(ats)
  }
  
  model <- data.fit[[1]]$ml$model
  data.ratio <- -diff(sapply(data.fit, function(x) x$ml$likelihood))
  stopifnot(!isSubmodel || data.ratio <= 0) # should never appear

  simu.n <- n - 1
  ratio <- numeric(simu.n)
  fit <- numeric(2)
  Z <- StandardizeData(x=x, y=y, z=z, T=T, grid=grid, data=data,
                       distances=distances, dim=dim, RFopt=RFopt)
  dist.given <- Z$dist.given
  if (length(Z$coord) > 1) stop("multisets of data cannot be considered yet")

  newx <- if (!dist.given) Z$coord[[1]]$x # lapply(Z$coord, function(x) x$x)
  newT <- if (!dist.given) Z$coord[[1]]$T # lapply(Z$coord, function(x) x$T),
  
  pch <- if (RFopt$general$pch=="") "" else '@'
  for (i in 1:simu.n) {
    if (printlevel>=PL_SUBIMPORTANT)
      cat("\n ", i, "th simulation out of", simu.n)
    else cat(pch)
    simu <- RFsimulate(model, x=newx, T=newT, grid=grid, 
                       distances=if (dist.given) Z$coord[[1]], dim=dim, spC=FALSE)
    guess <- users.guess   
    for (m in 1:length(model.list)) {
      simufit <-
        RFfit(model.list[[m]], x=newx, T=newT, grid=grid,
              data=simu,
              lower=lower, upper=upper, 
              methods=methods,
              sub.methods=sub.methods, optim.control=optim.control,
              users.guess=guess,
              distances=if (dist.given) Z$coord[[1]], dim=dim,
              transform=transform,
              ..., spConform=FALSE)
      fit[m] <- simufit$ml$ml
      guess <- if (isSubmodel) simufit$ml$model else NULL
    }

    ratio[i] <- -diff(fit)

    stopifnot(!isSubmodel || ratio[i] <= 0)# should never appear
    
    if (printlevel > PL_SUBIMPORTANT)
      Print(c(data.ratio, ratio), fit, rank(c(data.ratio, ratio))[1])#
  }
  
  r <- rank(c(data.ratio, ratio))[1]
  
  p <- r / n

  msg <-
    paste("\nThe likehood ratio test ranks the likelihood of the data on rank",
          r, "among", simu.n, "simulations:", mess(alpha=alpha, p=p))

  res <- list(p=p, n=n, data.ratio=data.ratio, simu.ratios=ratio,
              data.fit=data.fit, msg=msg, model.list=model.list)
  class(res) <- "RFratiotest"
  if (RFopt$general$returncall)
    attr(res, "call") <- as.character(deparse(match.call()))
  attr(res, "coord_system") <- c(orig=RFopt$coords$coord_system,
                                 model=RFopt$coords$new_coord_system)
  return(res)

}
