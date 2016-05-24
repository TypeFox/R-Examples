# created  April 8, 2015 
# author: Mikis Stasinopoulos
# based on the function extractAIC() dropterm() addterm() stepAIC() of the MASS
# and stepGAIC of gamlss
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
extractTGD <- function (fit, newdata, ...) 
{
    if (is.gamlss(fit)) 
       {ls
         out<-getTGD(fit, newdata=newdata)
        edf <-  fit$df.fit #  the df used in the fit 
     # noObs <- out$newGD/out$newPE # the number of new observarions
        tgd <-  out$TGD  #  the test GD
        return( c(edf, tgd))  #  return( c(edf, tgd, noObs))
       }
    else stop(paste("this is not a gamlss object"))
}
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
drop1TGD<-function (object,
                           scope, 
                           newdata,
                           parameter = c("mu", "sigma", "nu", "tau"), 
                           sorted = FALSE, 
                           trace = FALSE, 
                        parallel = c("no", "multicore", "snow"), 
                           ncpus = 1L, 
                              cl = NULL, 
                           ...) 
{
#-------------------------------------------------------------------------------
# like the drop.scope() + parameter
  drop1.scope<-function (terms1, terms2, parameter = c("mu", "sigma", "nu", "tau")) 
  {
    parameter <- match.arg(parameter)
    terms1 <- terms(terms1, parameter)
    f2 <- if (missing(terms2)) 
        numeric(0)
    else attr(terms(terms2, parameter), "factor")
    factor.scope(attr(terms1, "factor"), list(drop = f2))$drop
  }
#-------------------------------------------------------------------------------
# Main function starts here
#-------------------------------------------------------------------------------
    parameter <- match.arg(parameter)
    if (!parameter %in% object$par) 
        stop(paste(parameter, "is not a parameter in the object", "\n"))
    tl <- attr(terms(object, parameter ), "term.labels")
    if (missing(scope)) 
        scope <- drop1.scope(object, parameter = parameter)
    else 
       {
        if (!is.character(scope)) 
            scope <- attr(terms(update.formula(formula(object, parameter=parameter), scope), parameter), 
                "term.labels")
        if (!all(match(scope, tl, FALSE))) 
            stop("scope is not a subset of term labels")
       }
    ns <- length(scope)
    ans <- matrix(nrow = ns + 1, ncol = 2, dimnames = list(c("<none>", 
        scope), c("df", "AIC")))
    ans[1, ] <- extractTGD(object, newdata=newdata,  ...)
#--------------- PARALLEL-------------------------------------------------------
#----------------SET UP PART----------------------------------------------------
if (missing(parallel)) 
  parallel <- "no"
parallel <- match.arg(parallel)
have_mc <- have_snow <- FALSE
if (parallel != "no" && ncpus > 1L) 
{
  if (parallel == "multicore") 
    have_mc <- .Platform$OS.type != "windows"
  else if (parallel == "snow") 
    have_snow <- TRUE
  if (!have_mc && !have_snow) 
    ncpus <- 1L
  loadNamespace("parallel")
}
# -------------- finish parallel------------------------------------------------
#-------------------------------------------------------------------------------
#  function for parallel apply
fn <- function(term)
{
  if (trace) 
    cat("trying -", term, "\n")
  nfit <- update(object,  as.formula(paste("~ . -", term)), parameter = parameter, # MS 
                 evaluate = FALSE, trace=FALSE)
  nfit <- try(eval.parent(nfit), silent=TRUE)
  if (any(class(nfit)%in%"try-error"))
  { 
     cat("Model with term ", term, "has failed \n")       
    c(NA,NA)#
       
  }
   else{
      tgd=try(extractTGD(nfit, newdata=newdata,  ...),silent=TRUE)
      if (any(class(tgd)%in%"try-error"))
      {
          cat("TGD for term ", term, "has failed \n")
          c(NA,NA)#
      }
      else
      {tgd}
  }
}
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# --------  parallel -----------------------------------------------------------
ans[-1,] <- if (ncpus > 1L && (have_mc || have_snow)) 
{
  if (have_mc) 
  {# sapply(scope, fn)
    matrix(unlist(parallel::mclapply(scope, fn, mc.cores = ncpus)), ncol=2, 
           byrow = T)   
  }
  else if (have_snow) 
  {
    list(...)
    if (is.null(cl)) 
    {
      # make the cluster
      # cl <- parallel::makePSOCKcluster(rep("localhost", ncpus))
      cl <- parallel::makeForkCluster(ncpus)
      if (RNGkind()[1L] == "L'Ecuyer-CMRG") 
        parallel::clusterSetRNGStream(cl)
      res <- t(parallel::parSapply(cl, scope, fn))
      parallel::stopCluster(cl)
      res
    }
    else t(parallel::parSapply(cl, scope, fn))
  }
} # end parallel ---------------------------------------------------------------
else t(sapply(scope, fn)) 
     dfs <- ans[1, 1] - ans[, 1]
    dfs[1] <- NA
    aod <- data.frame(Df = dfs, TGD = ans[, 2])
    o <- if (sorted) 
        order(aod$TGD)
    else seq(along = aod$TGD)
    aod <- aod[o, ]
    head <- c("Single term deletions for", parameter, "\nModel:", deparse(as.vector(formula(object, parameter))))
    class(aod) <- c("anova", "data.frame")
    attr(aod, "heading") <- head
    aod
}
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
add1TGD<- function (object, 
                           scope,
                           newdata,
                           parameter = c("mu", "sigma", "nu", "tau"),  
                         sorted = FALSE, 
                          trace = FALSE, 
                       parallel = c("no", "multicore", "snow"), 
                          ncpus = 1L, 
                             cl = NULL, 
                           ...) 
{
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
 add.scope <- function (terms1, terms2, parameter = c("mu", "sigma", "nu", "tau") ) 
    {
      parameter <- match.arg(parameter)
    terms1 <- terms(terms1, parameter)
    terms2 <- terms(terms2, parameter)
    factor.scope(attr(terms1, "factor"), list(add = attr(terms2, 
        "factor")))$add  
   }
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
   parameter <- match.arg(parameter)
    if (missing(scope) || is.null(scope)) 
        stop("no terms in scope")
    if (!is.character(scope)) 
        scope <- add.scope(object, terms(update.formula(formula(object, parameter=parameter), scope)), parameter = parameter)
    if (!length(scope)) 
        stop("no terms in scope for adding to object")
    ns <- length(scope)
    ans <- matrix(nrow = ns + 1, ncol = 2, dimnames = list(c("<none>", 
        scope), c("df", "TGD")))
    ans[1, ] <- extractTGD(object, newdata=newdata, ...)
#--------------- PARALLEL-------------------------------------------------------
#----------------PART-----------------------------------------------------------
if (missing(parallel)) 
  parallel <- "no"
parallel <- match.arg(parallel)
have_mc <- have_snow <- FALSE
if (parallel != "no" && ncpus > 1L) 
{
  if (parallel == "multicore") 
    have_mc <- .Platform$OS.type != "windows"
  else if (parallel == "snow") 
    have_snow <- TRUE
  if (!have_mc && !have_snow) 
    ncpus <- 1L
  loadNamespace("parallel")
}
# -------------- finish parallel------------------------------------------------
#-------------------------------------------------------------------------------
#     for (i in seq(ns)) 
#       {
#         tt <- scope[i]
#         if (trace) 
#             cat("trying +", tt, "\n")
#         nfit <- update(object, as.formula(paste("~ . +", tt)), parameter = parameter, trace=FALSE, 
#             evaluate = FALSE)
#         nfit <- eval.parent(nfit)
#         ans[i + 1, ] <- extractTGD(nfit, newdata=newdata, ...)
#       }
#  function for parallel apply
fn <- function(term)
{
  if (trace) 
    cat("trying -", term, "\n")
  nfit <- update(object,  as.formula(paste("~ . +", term)), parameter = parameter, 
                 evaluate = FALSE, trace=FALSE) 
  nfit <- try(eval.parent(nfit), silent=TRUE)
  if (any(class(nfit)%in%"try-error"))
  { 
     cat("Model with term ", term, "has failed \n")       
    c(NA,NA)        
  }
    else{
      tgd=try(extractTGD(nfit, newdata=newdata,  ...),silent=TRUE)
      if (any(class(tgd)%in%"try-error"))
      {
          cat("TGD for term ", term, "has failed \n")
          c(NA,NA)#
      }
      else {tgd}
  }
}
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# --------  parallel ----------------------------------------------------------- 
ans[-1,] <- if (ncpus > 1L && (have_mc || have_snow)) 
{
  if (have_mc) 
  {# sapply(scope, fn)
    matrix(unlist(parallel::mclapply(scope, fn, mc.cores = ncpus)), ncol=2, 
           byrow = T)   
  }
  else if (have_snow) 
  {
    list(...)
    if (is.null(cl)) 
    {
      # make the cluster
      # cl <- parallel::makePSOCKcluster(rep("localhost", ncpus))
      cl <- parallel::makeForkCluster(ncpus)
      if (RNGkind()[1L] == "L'Ecuyer-CMRG") 
        parallel::clusterSetRNGStream(cl)
      res <- t(parallel::parSapply(cl, scope, fn))
      parallel::stopCluster(cl)
      res
    }
    else t(parallel::parSapply(cl, scope, fn))
  }
} # end parallel ----------------------------------------------------------
else t(sapply(scope, fn)) 
    dfs <- ans[, 1] - ans[1, 1]
    dfs[1] <- NA
    aod <- data.frame(Df = dfs, TGD = ans[, 2])
    o <- if (sorted) 
        order(aod$TGD)
    else seq(along = aod$TGD)
    aod <- aod[o, ]
    head <- c("Single term additions for", parameter,"\nModel:", deparse(as.vector(formula(object,parameter))))
    class(aod) <- c("anova", "data.frame")
    attr(aod, "heading") <- head
    aod
}
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Venable and Ripley step  function
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
stepTGD <-function(object, 
                    scope,
                    newdata, 
                direction = c("both", "backward", "forward"), 
                    trace = TRUE, 
                     keep = NULL, 
                    steps = 1000,
                     parameter = c("mu", "sigma", "nu", "tau"),
                 parallel = c("no", "multicore", "snow"), 
                    ncpus = 1L, 
                       cl = NULL, 
                       ...)                    
                    
{
#-------------------------------------------------------------------------------
    mydeviance <- function(x,newdata, ...) 
      {
        dev <- getTGD(x, newdata=newdata)$TGD
        if (!is.null(dev)) 
            dev
        else extractTGD(x, newdata=newdata)[2]  
      }
#------------------------------------------------------------------------------
    cut.string <- function(string) 
      {
        if (length(string) > 1) 
            string[-1] <- paste("\n", string[-1], sep = "")
        string
      }
#------------------------------------------------------------------------------
    re.arrange <- function(keep) 
      {
        namr <- names(k1 <- keep[[1]])
        namc <- names(keep)
          nc <- length(keep)
          nr <- length(k1)
        array(unlist(keep, recursive = FALSE), c(nr, nc), list(namr, 
            namc))
      }
#------------------------------------------------------------------------------
    step.results<- function(models, fit, object) #
      {
        change <- sapply(models, "[[", "change")
            rd <- sapply(models, "[[", "deviance")
            dd <- c(NA, abs(diff(rd)))
           rdf <- sapply(models, "[[", "df.resid")
           ddf <- c(NA, abs(diff(rdf)))
          # TGD <- sapply(models, "[[", "TGD")
       heading <- c("Stepwise Model Path \nTest Global Deviance Table", 
            "\nInitial", parameter," Model:", deparse(as.vector(formula(object, parameter=parameter))), 
            "\nFinal", parameter, " Model:", deparse(as.vector(formula(fit, parameter=parameter))), 
            "\n")
           aod <-  data.frame(Step = change, Df = ddf, "TGD diff." = dd, 
            "Resid. Df" = rdf, "Test G.Dev" = rd,  
            check.names = FALSE)
        attr(aod, "heading") <- heading
        class(aod) <- c("Anova", "data.frame")
        fit$anova <- aod
        fit
      }
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
drop1TGDP<-function (object, 
                      scope, 
                      newdata,
                      parameter = c("mu", "sigma", "nu", "tau"), 
                      sorted = FALSE, 
                      trace = FALSE, 
                      ...) 
{
  #-------------------------------------------------------------------------------
  # like the drop.scope() + parameter
  drop1.scope<-function (terms1, terms2, parameter = c("mu", "sigma", "nu", "tau")) 
  {
    parameter <- match.arg(parameter)
    terms1 <- terms(terms1, parameter)
    f2 <- if (missing(terms2)) 
      numeric(0)
    else attr(terms(terms2, parameter), "factor")
    factor.scope(attr(terms1, "factor"), list(drop = f2))$drop
  }
  #-------------------------------------------------------------------------------
  # Main function starts here
  #-------------------------------------------------------------------------------
  parameter <- match.arg(parameter)
  if (!parameter %in% object$par) 
    stop(paste(parameter, "is not a parameter in the object", "\n"))
  tl <- attr(terms(object, parameter ), "term.labels")
  if (missing(scope)) 
    scope <- drop1.scope(object, parameter = parameter)
  else 
  {
    if (!is.character(scope)) 
      scope <- attr(terms(update.formula(formula(object, parameter=parameter), scope), parameter), 
                    "term.labels")
    if (!all(match(scope, tl, FALSE))) 
      stop("scope is not a subset of term labels")
  }
  ns <- length(scope)
  ans <- matrix(nrow = ns + 1, ncol = 2, dimnames = list(c("<none>", 
                                                           scope), c("df", "AIC")))
  ans[1, ] <- extractTGD(object, newdata=newdata,  ...)
  #--------------- PARALLEL-------------------------------------------------------
  #-------------------------------------------------------------------------------
  #  function for parallel apply
  fn <- function(term)
  {
    if (trace) 
      cat("trying -", term, "\n")
    nfit <- update(object,  as.formula(paste("~ . -", term)), parameter = parameter, # MS 
                   evaluate = FALSE, trace=FALSE)
    nfit <- try(eval.parent(nfit), silent=TRUE)
    if (any(class(nfit)%in%"try-error"))
    { 
      cat("Model with term ", term, "has failed \n")       
      c(NA,NA)#         
    }
    else{
        tgd=try(extractTGD(nfit, newdata=newdata,  ...),silent=TRUE)
        if (any(class(tgd)%in%"try-error"))
        {
            cat("TGD for term ", term, "has failed \n")
            c(NA,NA)#
        }
        else
        {tgd}
    }
  }
  #-------------------------------------------------------------------------------
  #-------------------------------------------------------------------------------
  #-------------------------------------------------------------------------------
  # --------  parallel -----------------------------------------------------------
  ans[-1,] <- if (ncpus > 1L && (have_mc || have_snow)) 
  {
    if (have_mc) 
    {# sapply(scope, fn)
      matrix(unlist(parallel::mclapply(scope, fn, mc.cores = ncpus)), ncol=2, 
             byrow = T)   
    }
    else if (have_snow) 
    {
      list(...)
      if (is.null(cl)) 
      {
        # make the cluster
        # cl <- parallel::makePSOCKcluster(rep("localhost", ncpus))
        res <- t(parallel::parSapply(cl, scope, fn))
        res
      }
      else t(parallel::parSapply(cl, scope, fn))
    }
  } # end parallel ---------------------------------------------------------------
  else t(sapply(scope, fn)) 
  dfs <- ans[1, 1] - ans[, 1]
  dfs[1] <- NA
  aod <- data.frame(Df = dfs, TGD = ans[, 2])
  o <- if (sorted) 
    order(aod$TGD)
  else seq(along = aod$TGD)
  aod <- aod[o, ]
  head <- c("Single term deletions for", parameter, "\nModel:", deparse(as.vector(formula(object, parameter))))
  class(aod) <- c("anova", "data.frame")
  attr(aod, "heading") <- head
  aod
}
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# end drop1TGDP
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
add1TGDP<- function (object, 
                      scope,
                      newdata,
                      parameter = c("mu", "sigma", "nu", "tau"),  
                      sorted = FALSE, 
                      trace = FALSE, 
                      ...) 
{
  #-------------------------------------------------------------------------------
  #-------------------------------------------------------------------------------
  add.scope <- function (terms1, terms2, parameter = c("mu", "sigma", "nu", "tau") ) 
  {
    parameter <- match.arg(parameter)
    terms1 <- terms(terms1, parameter)
    terms2 <- terms(terms2, parameter)
    factor.scope(attr(terms1, "factor"), list(add = attr(terms2, 
                                                         "factor")))$add  
  }
  #-------------------------------------------------------------------------------
  #-------------------------------------------------------------------------------
  parameter <- match.arg(parameter)
  if (missing(scope) || is.null(scope)) 
    stop("no terms in scope")
  if (!is.character(scope)) 
    scope <- add.scope(object, terms(update.formula(formula(object, parameter=parameter), scope)), parameter = parameter)
  if (!length(scope)) 
    stop("no terms in scope for adding to object")
  ns <- length(scope)
  ans <- matrix(nrow = ns + 1, ncol = 2, dimnames = list(c("<none>", 
                                                           scope), c("df", "TGD")))
  ans[1, ] <- extractTGD(object, newdata=newdata, ...)
  #--------------- PARALLEL-------------------------------------------------------
  fn <- function(term)
  {
    if (trace) 
      cat("trying -", term, "\n")
    nfit <- update(object,  as.formula(paste("~ . +", term)), parameter = parameter, 
                   evaluate = FALSE, trace=FALSE) 
    nfit <- try(eval.parent(nfit), silent=TRUE)
    if (any(class(nfit)%in%"try-error"))
    { 
    cat("Model with term ", term, "has failed \n")       
      c(NA,NA)#        
    }
    else{
        tgd=try(extractTGD(nfit, newdata=newdata,  ...),silent=TRUE)
        if (any(class(tgd)%in%"try-error"))
        {
            cat("TGD for term ", term, "has failed \n")
            c(NA,NA)#
        }
        else
        {tgd}
    }
  }
  #-------------------------------------------------------------------------------
  #-------------------------------------------------------------------------------
  # --------  parallel ----------------------------------------------------------- 
  ans[-1,] <- if (ncpus > 1L && (have_mc || have_snow)) 
  {
    if (have_mc) 
    {# sapply(scope, fn)
      matrix(unlist(parallel::mclapply(scope, fn, mc.cores = ncpus)), ncol=2, 
             byrow = T)   
    }
    else if (have_snow) 
    {
      list(...)
      if (is.null(cl)) 
      {
        res <- t(parallel::parSapply(cl, scope, fn))
        res
      }
      else t(parallel::parSapply(cl, scope, fn))
    }
  } # end parallel ----------------------------------------------------------
  else t(sapply(scope, fn)) 
  dfs <- ans[, 1] - ans[1, 1]
  dfs[1] <- NA
  aod <- data.frame(Df = dfs, TGD = ans[, 2])
  o <- if (sorted) 
    order(aod$TGD)
  else seq(along = aod$TGD)
  aod <- aod[o, ]
  head <- c("Single term additions for", parameter,"\nModel:", deparse(as.vector(formula(object,parameter))))
  class(aod) <- c("anova", "data.frame")
  attr(aod, "heading") <- head
  aod
}
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# main fuction starts here
     parameter <- match.arg(parameter) 
#--------------- PARALLEL-------------------------------------------------------
#----------------PART-----------------------------------------------------------
#--------------- PARALLEL-------------------------------------------------------
#----------------SET UP PART---------------------------------------------------
if (missing(parallel)) 
  parallel <- "no"
parallel <- match.arg(parallel)
have_mc <- have_snow <- FALSE
if (parallel != "no" && ncpus > 1L) 
{
  if (parallel == "multicore") 
    have_mc <- .Platform$OS.type != "windows"
  else if (parallel == "snow") 
    have_snow <- TRUE
  if (!have_mc && !have_snow) 
    ncpus <- 1L
  loadNamespace("parallel")
}
if (have_snow)
{
  cl <- parallel::makeForkCluster(ncpus)
  if (RNGkind()[1L] == "L'Ecuyer-CMRG") 
    parallel::clusterSetRNGStream(cl)
  on.exit(parallel::stopCluster(cl))
}         
# -------------- finish parallel------------------------------------------------
#-------------------------------------------------------------------------------
    Terms <- terms(object, parameter)
    if (parameter=="mu")
     {
           object$formula <- Terms 
      object$call$formula <- Terms
     } 
    else
     {
      object[[paste(parameter,"formula",sep=".")]] <- Terms
      object[[paste(parameter,"formula",sep=".")]][[2]]<-NULL
      if (paste(parameter, "formula", sep=".")%in%names(object$call)) 
        object$call[[paste(parameter,"formula",sep=".")]] <- formula(Terms)[-2]
     else ##this is when the sigma formula is not defined
        {object$call[[paste(parameter,"formula",sep=".")]] <- formula(Terms)[-2]
        names(object$call)[length(names(object$call))]<-paste(parameter,"formula",sep=".")
        }
     }     
           md <- missing(direction)
    direction <- match.arg(direction)
     backward <- direction == "both" | direction == "backward"
      forward <- direction == "both" | direction == "forward"
    if (missing(scope)) 
      {
        fdrop <- numeric(0)
         fadd <- attr(Terms, "factors")
        if (md) 
            forward <- FALSE
      }
    else 
      {
        if (is.list(scope)) 
          {
            fdrop <- if (!is.null(fdrop <- scope$lower)) 
                attr(terms(update.formula(formula(object, parameter=parameter), fdrop), parameter = parameter), "factors")
            else numeric(0)
             fadd <- if (!is.null(fadd <- scope$upper)) 
                attr(terms(update.formula(formula(object, parameter=parameter), fadd), parameter = parameter), "factors")
          }
        else 
          {
             fadd <- if (!is.null(fadd <- scope)) 
                attr(terms(update.formula(formula(object, parameter=parameter), scope), parameter = parameter ), "factors")
            fdrop <- numeric(0)
          }
       }
    models <- vector("list", steps)
    if (!is.null(keep)) 
      keep.list <- vector("list", steps)
    if (is.list(object) && (nmm <- match("nobs", names(object), 
        0)) > 0) 
         n <- object[[nmm]]
    else n <- length(residuals(object))
     fit <- object
    bTGD <- extractTGD(fit, newdata, ...)
     edf <- bTGD[1]
    bTGD <- bTGD[2]
    if (is.na(bTGD)) 
        stop("TGD is not defined for this model, so stepTGD cannot proceed")
      nm <- 1
   Terms <- terms(fit, parameter)
    if (trace)
        cat("Distribution parameter: ", parameter, "\n") 
        cat("Start:  TGD=", format(round(bTGD, 2)), "\n", cut.string(deparse(as.vector(formula(fit, parameter=parameter)))), 
            "\n\n")
    models[[nm]] <- list(deviance = mydeviance(fit, newdata=newdata), df.resid = n - 
        edf, change = "", TGD = bTGD)
    if (!is.null(keep)) 
        keep.list[[nm]] <- keep(fit, bTGD)
   usingCp <- FALSE
    while (steps > 0) 
     {
      steps <- steps - 1
        TGD <- bTGD
       ffac <- attr(Terms, "factors")
        if (!is.null(sp <- attr(Terms, "specials")) && !is.null(st <- sp$strata)) 
            ffac <- ffac[-st, ]
       scope <- factor.scope(ffac, list(add = fadd, drop = fdrop))
         aod <- NULL
      change <- NULL
        if (backward && length(scope$drop)) 
           {
            aod <- drop1TGDP(fit, scope=scope$drop, newdata=newdata, parameter = parameter, 
                           trace = max(0,trace - 1), ...)
             rn <- row.names(aod)
            row.names(aod) <- c(rn[1], paste("-", rn[-1], sep = " "))
            if (any(aod$Df == 0, na.rm = TRUE)) 
               {
                zdf <- aod$Df == 0 & !is.na(aod$Df)
                 nc <-  match(c("Cp", "TGD"), names(aod))
                 nc <- nc[!is.na(nc)][1]
                 ch <- abs(aod[zdf, nc] - aod[1, nc]) > 0.01
                if (any(ch)) 
                  {
                  warning("0 df terms are changing TGD")
                  zdf <- zdf[!ch]
                  }
                if (length(zdf) > 0) 
                  change <- rev(rownames(aod)[zdf])[1]
               }
           }
        if (is.null(change)) 
           {
            if (forward && length(scope$add)) 
              {
                aodf <- add1TGDP(fit, scope=scope$add, newdata=newdata, parameter =parameter,  
                  trace = max(0, trace - 1),  ...)
                  rn <- row.names(aodf)
                row.names(aodf) <- c(rn[1], paste("+", rn[-1], 
                  sep = " "))
                 aod <- if (is.null(aod)) 
                  aodf
                else rbind(aod, aodf[-1, , drop = FALSE])
              }
            attr(aod, "heading") <- NULL
            if (is.null(aod) || ncol(aod) == 0) 
                break
            nzdf <- if (!is.null(aod$Df)) 
                aod$Df != 0 | is.na(aod$Df)
             aod <- aod[nzdf, ]
            if (is.null(aod) || ncol(aod) == 0) 
                break
            nc <- match(c("Cp", "TGD"), names(aod))
            nc <- nc[!is.na(nc)][1]
             o <- order(aod[, nc])
            if (trace) 
                print(aod[o, ])
            if (o[1] == 1) 
                break
            change <- rownames(aod)[o[1]]
           }
        usingCp <- match("Cp", names(aod), 0) > 0
        fit <- update(fit, paste("~ .", change), parameter = parameter,  evaluate = FALSE, trace = FALSE) #MS
        fit <- eval.parent(fit)
        if (is.list(fit) && (nmm <- match("nobs", names(fit), 
            0)) > 0) 
            nnew <- fit[[nmm]]
        else nnew <- length(residuals(fit))
        if (nnew != n) 
            stop("number of rows in use has changed: remove missing values?")
        Terms <- terms(fit, parameter)
         bTGD <- extractTGD(fit, newdata=newdata,  ...)
          edf <- bTGD[1]
         bTGD <- bTGD[2]
        if (trace) 
            cat("\nStep:  TGD=", format(round(bTGD, 2)), "\n", 
                cut.string(deparse(as.vector(formula(fit, parameter)))), 
                "\n\n")
        if (bTGD >= TGD + 1e-07) 
            break
        nm <- nm + 1
        models[[nm]] <- list(deviance = mydeviance(fit, newdata=newdata), df.resid = n - 
            edf, change = change, TGD = bTGD)
        if (!is.null(keep)) 
            keep.list[[nm]] <- keep(fit, bTGD)
    }
    if (!is.null(keep)) 
        fit$keep <- re.arrange(keep.list[seq(nm)])
    step.results(models = models[seq(nm)], fit, object)
}
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
