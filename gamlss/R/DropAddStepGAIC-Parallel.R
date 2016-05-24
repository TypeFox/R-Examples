# This is an attempt make parallel the computation in dropterm(),
# addterm() and stepGAIC()
# Mikis Stasinopoulos  wrote the original with some help from Fernanda de Bastiani
# Aendment by Daniil Kosie to avoid the failure of the funcion 
# This file contains
# i)   dropterm.gamlss()
# ii)  addterm.gamlss()
# iii) stepGAIC()
# 
#-------------------------------------------------------------------------------
# this refers to how the function histrorically was developled and can be ingnore now 
# The approach is as follows
#  i) STEP 1: trying to rewrite the dropterm.gamlss() function having 
#      instead of the loop 
#      for (i in seq(ns)) 
#      a function fn() which is call with lapply() to produce the updated fitted
#      gamlss objects. 
#      the function dropterm1()  is doing just that
#  ii) STEP 2 Take now the new function and add the parallel stuff from bootT
#      this is function dropterm2() 
# iii) Fernanda will try to do the same with addterm - it is done
#      it is omplemented in addterm2()  
#  iv) After that the idea is to imlement this in stepGAIC 
#       stepGAIC1(): uses the functions dropterm2() and addterm2() 
#                    without any modification (not faster than stepGAIC)
#       stepGAIC2(): uses the functions local droptermP() and addtermP() but with 
#                    parallel called every time the fuctions are called 
#                   (not faster than stepGAIC)
#       stepGAIC.P(): uses modified local droptermP() and addtermP() but the creation #   #         of parallel is done in the begining it looks that that is the best  in that 
#                   with snow reduce the time    
#     v) we have to implement stepGAICAll A and B now 
#            stepGAICAll.Bp is ready 
# 
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
dropterm.gamlss <- function (object, 
                     scope, 
                     what = c("mu", "sigma", "nu", "tau"),
                     parameter = NULL,
                    scale = 0, 
                     test = c("none", "Chisq"), 
                        k = 2, 
                   sorted = FALSE, 
                    trace = FALSE, 
                 parallel = c("no", "multicore", "snow"), #The type of parallel operation to be used (if any). If missing, the default is taken from the option "boot.parallel" (and if that is not set, "no")
                    ncpus = 1L, #integer: number of processes to be used in parallel operation: typically one would chose this to the number of available CPUs
                       cl = NULL, # An optional parallel or snow cluster for use if parallel = "snow". If not supplied, a cluster on the local machine is created for the duration of the boot call.
                     ...) 
{
  #-----------------------------------------------------------------------------
  #-----------------------------------------------------------------------------
  drop1.scope<-function (terms1, terms2, what = c("mu", "sigma", "nu", "tau"), parameter= NULL) 
  {
     what <- if (!is.null(parameter))  {
    match.arg(parameter, choices=c("mu", "sigma", "nu", "tau"))} else  match.arg(what)
    terms1 <- terms(terms1, what)
    f2 <- if (missing(terms2)) 
      numeric(0)
    else attr(terms(terms2, what), "factor")
    factor.scope(attr(terms1, "factor"), list(drop = f2))$drop
  }
  #-----------------------------------------------------------------------------
  #-----------------------------------------------------------------------------
  safe_pchisq <- function (q, df, ...) 
  {
    df[df <= 0] <- NA
    pchisq(q = q, df = df, ...)
  }
#-------------------------------------------------------------------------------
#   main functions starts here 
#-------------------------------------------------------------------------------
  what <- if (!is.null(parameter))  {
    match.arg(parameter, choices=c("mu", "sigma", "nu", "tau"))} else  match.arg(what)
  if (!what %in% object$par) 
    stop(paste(what, "is not a parameter in the object", "\n"))
  tl <- attr(terms(object, what ), "term.labels")
  if (missing(scope)) 
    scope <- drop1.scope(object, what = what)
  else 
  {
    if (!is.character(scope)) 
      scope <- attr(terms(update.formula(formula(object, what=what), scope), what), 
                    "term.labels")
    if (!all(match(scope, tl, FALSE))) 
      stop("scope is not a subset of term labels")
  }
  ns <- length(scope)
  ans <- matrix(nrow = ns + 1, ncol = 2, dimnames = list(c("<none>", scope),
            c("df", "AIC")))
  ans[1, ] <- extractAIC(object, scale, k = k,   ...)
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
    nfit <- update(object,  as.formula(paste("~ . -", term)), what = what, evaluate = FALSE, trace=FALSE) 
    nfit <- try(eval.parent(nfit), silent=TRUE)
    if (any(class(nfit)%in%"try-error"))
    { 
      cat("Model with term ", term, "has failed \n")       
      c(NA,NA) # Daniil: prevents execution to stop when fitting of term fails, returns c(df=NA, GAIC = NA)    
    }
       else
    {
        gaic=try(extractAIC(nfit, scale, k = k,   ...),silent=TRUE)
        if (any(class(gaic)%in%"try-error"))
        {
            cat("GAIC for term ", term, "has failed \n")
            c(NA,NA)#
        }
        else
        {gaic}
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
} # end parallel ---------------------------------------------------------------
else t(sapply(scope, fn)) 
  dfs <- ans[1, 1] - ans[, 1]
  dfs[1] <- NA
  aod <- data.frame(Df = dfs, AIC = ans[, 2])
  o <- if (sorted) 
    order(aod$AIC)
  else seq(along = aod$AIC)
  test <- match.arg(test)
  if (test == "Chisq") 
  {
    dev <- ans[, 2] - k * ans[, 1]
    dev <- dev - dev[1]
    dev[1] <- NA
    nas <- !is.na(dev)
    P <- dev
    P[nas] <- safe_pchisq(dev[nas], dfs[nas], lower.tail = FALSE)
    aod[, c("LRT", "Pr(Chi)")] <- list(dev, P)
  }
  aod <- aod[o, ]
  head <- c("Single term deletions for", what, "\nModel:", deparse(as.vector(formula(object, what))))
  if (scale > 0) 
    head <- c(head, paste("\nscale: ", format(scale), "\n"))
  class(aod) <- c("anova", "data.frame")
  attr(aod, "heading") <- head
  aod
}
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#  addterm
#-------------------------------------------------------------------------------
addterm.gamlss <- function (object, 
                     scope,
                     what = c("mu", "sigma", "nu", "tau"), 
                     parameter= NULL,
                    scale = 0, 
                     test = c("none", "Chisq"), 
                        k = 2, 
                   sorted = FALSE, 
                    trace = FALSE, 
                 parallel = c("no", "multicore", "snow"), 
                    ncpus = 1L, 
                       cl = NULL, 
                       ...) 
{
  #-----------------------------------------------------------------------------
  add.scope <- function (terms1, terms2, what = c("mu", "sigma", "nu", "tau"), parameter= NULL ) 
  {
    what <- if (!is.null(parameter))  {
    match.arg(parameter, choices=c("mu", "sigma", "nu", "tau"))} else  match.arg(what)
    terms1 <- terms(terms1, what)
    terms2 <- terms(terms2, what)
    factor.scope(attr(terms1, "factor"), list(add = attr(terms2, 
                                                         "factor")))$add  
  }
  #-----------------------------------------------------------------------------
  safe_pchisq <- function (q, df, ...) 
  {
    df[df <= 0] <- NA
    pchisq(q = q, df = df, ...)
  }
  #-----------------------------------------------------------------------------
  what <- if (!is.null(parameter))  {
    match.arg(parameter, choices=c("mu", "sigma", "nu", "tau"))} else  match.arg(what)
  if (!what %in% object$par) 
    stop(paste(what, "is not a parameter in the object", "\n"))
  if (missing(scope) || is.null(scope)) 
    stop("no terms in scope")
  if (!is.character(scope)) 
    scope <- add.scope(object, terms(update.formula(formula(object, what=what), scope)), what = what)
  if (!length(scope)) 
    stop("no terms in scope for adding to object")
  ns <- length(scope)
  ans <- matrix(nrow = ns + 1, ncol = 2, dimnames = list(c("<none>", scope), 
                c("df", "AIC")))
  ans[1, ] <- extractAIC(object, scale, k = k,   ...)
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
#  function for parallel apply
fn <- function(term)
{
  if (trace) 
    cat("trying -", term, "\n")
  nfit <- update(object,  as.formula(paste("~ . +", term)), what = what, evaluate = FALSE, trace=FALSE) 
  nfit <- try(eval.parent(nfit), silent=TRUE)
  if (any(class(nfit)%in%"try-error"))
  { 
    cat("Model with term ", term, "has failed \n")       
    c(NA,NA) # Daniil: prevents execution to stop when fitting of term fails, returns c(df=NA, GAIC = NA)       
  }
    else
  {
      gaic=try(extractAIC(nfit, scale, k = k,   ...),silent=TRUE)
      if (any(class(gaic)%in%"try-error"))
      {
          cat("GAIC for term ", term, "has failed \n")
          c(NA,NA)#
      }
      else
      {gaic}
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
} # end parallel ---------------------------------------------------------------
else t(sapply(scope, fn)) 

#   for (i in seq(ns)) 
#   {
#     tt <- scope[i]
#     if (trace) 
#       cat("trying +", tt, "\n")
#     nfit <- update(object, as.formula(paste("~ . +", tt)), what = what, trace=FALSE, 
#                    evaluate = FALSE)
#     nfit <- try(eval.parent(nfit), silent=TRUE)
#     if (any(class(nfit)%in%"try-error"))
#     { 
#       cat("Model with term ", tt, "has failed \n")       
#       ans[i + 1, ] <- NA# extractAIC(object, scale, k = k, ...)          
#     }
#     else ans[i + 1, ] <- extractAIC(nfit, scale, k = k,  ...)
#   }
#-------------------------------------------------------------------------------
  dfs <- ans[, 1] - ans[1, 1]
  dfs[1] <- NA
  aod <- data.frame(Df = dfs, AIC = ans[, 2])
  o <- if (sorted) 
    order(aod$AIC)
  else seq(along = aod$AIC)
  test <- match.arg(test)
  if (test == "Chisq") 
  {
    dev <- ans[, 2] - k * ans[, 1]
    dev <- dev[1] - dev
    dev[1] <- NA
    nas <- !is.na(dev)
    P <- dev
    P[nas] <- safe_pchisq(dev[nas], dfs[nas], lower.tail = FALSE)
    aod[, c("LRT", "Pr(Chi)")] <- list(dev, P)
  }
  aod <- aod[o, ]
  head <- c("Single term additions for", what,"\nModel:", deparse(as.vector(formula(object,what))))
  if (scale > 0) 
    head <- c(head, paste("\nscale: ", format(scale), "\n"))
  class(aod) <- c("anova", "data.frame")
  attr(aod, "heading") <- head
  aod
}
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# modification of Venable and Ripley stepAIC() function
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
stepGAIC <-function(object, 
                     scope, 
                     direction = c("both", "backward", "forward"), 
                     trace = T, 
                     keep = NULL, 
                     steps = 1000,
                     scale = 0,
                     what = c("mu", "sigma", "nu", "tau"),
                     parameter= NULL,
                     k = 2,
                     parallel = c("no", "multicore", "snow"), 
                     ncpus = 1L, 
                     cl = NULL, 
                     ...)                    

{
#-----------------------------------------------------------------------------
# local functions 
#  1)    mydevianc()
#  2)    cut.string()
#  3)    re.arrange()
#  4)    step.results() 
#  5)    droptermP()
#  6)    addtermP()
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#  1 
mydeviance <- function(x, ...) 
  {
    dev <- deviance(x)
    if (!is.null(dev)) 
      dev
    else extractAIC(x, k = 0)[2]
  }
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
# 2 
  cut.string <- function(string) 
  {
    if (length(string) > 1) 
      string[-1] <- paste("\n", string[-1], sep = "")
    string
  }
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
# 3 
 re.arrange <- function(keep) 
  {
    namr <- names(k1 <- keep[[1]])
    namc <- names(keep)
    nc <- length(keep)
    nr <- length(k1)
    array(unlist(keep, recursive = FALSE), c(nr, nc), list(namr, 
                                                           namc))
  }
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
# 4 
  step.results <- function(models, fit, object, usingCp = FALSE) #
  {
    change <- sapply(models, "[[", "change")
    rd <- sapply(models, "[[", "deviance")
    dd <- c(NA, abs(diff(rd)))
    rdf <- sapply(models, "[[", "df.resid")
    ddf <- c(NA, abs(diff(rdf)))
    AIC <- sapply(models, "[[", "AIC")
    heading <- c("Stepwise Model Path \nAnalysis of Deviance Table", 
                 "\nInitial", what," Model:", deparse(as.vector(formula(object, what=what))), 
                 "\nFinal", what, " Model:", deparse(as.vector(formula(fit, what=what))), 
                 "\n")
    aod <- if (usingCp) 
      data.frame(Step = change, Df = ddf, Deviance = dd, 
                 "Resid. Df" = rdf, "Resid. Dev" = rd, Cp = AIC, 
                 check.names = FALSE)
    else data.frame(Step = change, Df = ddf, Deviance = dd, 
                    "Resid. Df" = rdf, "Resid. Dev" = rd, AIC = AIC, 
                    check.names = FALSE)
    attr(aod, "heading") <- heading
    class(aod) <- c("Anova", "data.frame")
    fit$anova <- aod
    fit
  }
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
# 5
droptermP <- function (object, 
                       scope, 
                       what = c("mu", "sigma", "nu", "tau"),
                       parameter= NULL,
                       scale = 0, 
                       test = c("none", "Chisq"), 
                       k = 2, 
                       sorted = FALSE, 
                       trace = FALSE, 
                       parallel = c("no", "multicore", "snow"), #The type of parallel operation to be used (if any). If missing, the default is taken from the option "boot.parallel" (and if that is not set, "no")
                       ncpus = 1L, #nteger: number of processes to be used in parallel operation: typically one would chose this to the number of available CPUs
                       cl = NULL, # An optional parallel or snow cluster for use if parallel = "snow". If not supplied, a cluster on the local machine is created for the duration of the boot call.
                       ...) 
{
#-----------------------------------------------------------------------------
    drop1.scope<-function (terms1, terms2, what = c("mu", "sigma", "nu", "tau"), parameter= NULL) 
    {
      what <- if (!is.null(parameter))  {
    match.arg(parameter, choices=c("mu", "sigma", "nu", "tau"))} else  match.arg(what)
      terms1 <- terms(terms1, what)
      f2 <- if (missing(terms2)) 
        numeric(0)
      else attr(terms(terms2, what), "factor")
      factor.scope(attr(terms1, "factor"), list(drop = f2))$drop
    }
 #-----------------------------------------------------------------------------
    safe_pchisq <- function (q, df, ...) 
    {
      df[df <= 0] <- NA
      pchisq(q = q, df = df, ...)
    }
#-----------------------------------------------------------------------------
#  droptermP  starts here 
#----------------------------------------------------------------------------
    what <- if (!is.null(parameter))  {
    match.arg(parameter, choices=c("mu", "sigma", "nu", "tau"))} else  match.arg(what)
    if (!what %in% object$par) 
      stop(paste(what, "is not a parameter in the object", "\n"))
    tl <- attr(terms(object, what ), "term.labels")
    if (missing(scope)) 
      scope <- drop1.scope(object, what = what)
    else 
    {
      if (!is.character(scope)) 
        scope <- attr(terms(update.formula(formula(object, what=what), scope), what), 
                      "term.labels")
      if (!all(match(scope, tl, FALSE))) 
        stop("scope is not a subset of term labels")
    }
    ns <- length(scope)
    ans <- matrix(nrow = ns + 1, ncol = 2, dimnames = list(c("<none>", scope),
                                                           c("df", "AIC")))
    ans[1, ] <- extractAIC(object, scale, k = k,   ...)
    #  function for parallel apply
    fn <- function(term)
    {
      if (trace) 
        cat("trying -", term, "\n")
      nfit <- update(object,  as.formula(paste("~ . -", term)), what = what, evaluate = FALSE, trace=FALSE) 
      nfit <- try(eval.parent(nfit), silent=TRUE)
      if (any(class(nfit)%in%"try-error"))
      { 
        cat("Model with term ", term, "has failed \n")       
        c(NA, NA)# Daniil        
      }
     else
      {
          gaic=try(extractAIC(nfit, scale, k = k,   ...),silent=TRUE)
          if (any(class(gaic)%in%"try-error"))
          {
              cat("GAIC for term ", term, "has failed \n")
              c(NA,NA)#
          }
          else
          {gaic}
      }
    }
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
#           cl <- parallel::makeForkCluster(ncpus)
#           if (RNGkind()[1L] == "L'Ecuyer-CMRG") 
#             parallel::clusterSetRNGStream(cl)
          res <- t(parallel::parSapply(cl, scope, fn))
#          parallel::stopCluster(cl)
          res
        }
        else t(parallel::parSapply(cl, scope, fn))
      }
    } # end parallel ----------------------------------------------------------
    else t(sapply(scope, fn)) 
    dfs <- ans[1, 1] - ans[, 1]
    dfs[1] <- NA
    aod <- data.frame(Df = dfs, AIC = ans[, 2])
    o <- if (sorted) 
      order(aod$AIC)
    else seq(along = aod$AIC)
    test <- match.arg(test)
    if (test == "Chisq") 
    {
      dev <- ans[, 2] - k * ans[, 1]
      dev <- dev - dev[1]
      dev[1] <- NA
      nas <- !is.na(dev)
      P <- dev
      P[nas] <- safe_pchisq(dev[nas], dfs[nas], lower.tail = FALSE)
      aod[, c("LRT", "Pr(Chi)")] <- list(dev, P)
    }
    aod <- aod[o, ]
    head <- c("Single term deletions for", what, "\nModel:", deparse(as.vector(formula(object, what))))
    if (scale > 0) 
      head <- c(head, paste("\nscale: ", format(scale), "\n"))
    class(aod) <- c("anova", "data.frame")
    attr(aod, "heading") <- head
    aod
  }
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#  6 
#  addterm
addtermP<- function (object, 
                     scope,
                     what = c("mu", "sigma", "nu", "tau"), 
                     parameter= NULL,
                     scale = 0, 
                     test = c("none", "Chisq"), 
                     k = 2, 
                     sorted = FALSE, 
                     trace = FALSE, 
                     parallel = c("no", "multicore", "snow"), 
                     ncpus = 1L, 
                     cl = NULL, 
                     ...) 
{
  #-----------------------------------------------------------------------------
  add.scope <- function (terms1, terms2, what = c("mu", "sigma", "nu", "tau"), parameter= NULL ) 
  {
    what <- if (!is.null(parameter))  {
    match.arg(parameter, choices=c("mu", "sigma", "nu", "tau"))} else  match.arg(what)
    terms1 <- terms(terms1, what)
    terms2 <- terms(terms2, what)
    factor.scope(attr(terms1, "factor"), list(add = attr(terms2, 
                                                         "factor")))$add  
  }
  #-----------------------------------------------------------------------------
  safe_pchisq <- function (q, df, ...) 
  {
    df[df <= 0] <- NA
    pchisq(q = q, df = df, ...)
  }
  #-----------------------------------------------------------------------------
  #-----------------------------------------------------------------------------
  what <- if (!is.null(parameter))  {
  match.arg(parameter, choices=c("mu", "sigma", "nu", "tau"))} else  match.arg(what)  
    if (!what %in% object$par) 
    stop(paste(what, "is not a parameter in the object", "\n"))
  if (missing(scope) || is.null(scope)) 
    stop("no terms in scope")
  if (!is.character(scope)) 
    scope <- add.scope(object, terms(update.formula(formula(object, what=what), scope)), what = what)
  if (!length(scope)) 
    stop("no terms in scope for adding to object")
  ns <- length(scope)
  ans <- matrix(nrow = ns + 1, ncol = 2, dimnames = list(c("<none>", scope), 
                                                         c("df", "AIC")))
  ans[1, ] <- extractAIC(object, scale, k = k,   ...)
  #  function for parallel apply------------------------------------------------
  fn <- function(term)
  {
    if (trace) 
      cat("trying -", term, "\n")
    nfit <- update(object,  as.formula(paste("~ . +", term)), what = what, evaluate = FALSE, trace=FALSE) 
    nfit <- try(eval.parent(nfit), silent=TRUE)
    if (any(class(nfit)%in%"try-error"))
    { 
      cat("Model with term ", term, "has failed \n")       
      c(NA, NA)# extractAIC(object, scale, k = k, ...)          
    }
       else
    {
        gaic=try(extractAIC(nfit, scale, k = k,   ...),silent=TRUE)
        if (any(class(gaic)%in%"try-error"))
        {
            cat("GAIC for term ", term, "has failed \n")
            c(NA,NA)#
        }
        else
        {gaic}
    }
  }
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
#        cl <- parallel::makeForkCluster(ncpus)
#        if (RNGkind()[1L] == "L'Ecuyer-CMRG") 
#          parallel::clusterSetRNGStream(cl)
        res <- t(parallel::parSapply(cl, scope, fn))
#        parallel::stopCluster(cl)
        res
      }
      else t(parallel::parSapply(cl, scope, fn))
    }
  } # end parallel -------------------------------------------------------------
  else t(sapply(scope, fn)) 
#-----------------------------------------------------------------------------
  dfs <- ans[, 1] - ans[1, 1]
  dfs[1] <- NA
  aod <- data.frame(Df = dfs, AIC = ans[, 2])
  o <- if (sorted) 
    order(aod$AIC)
  else seq(along = aod$AIC)
  test <- match.arg(test)
  if (test == "Chisq") 
  {
    dev <- ans[, 2] - k * ans[, 1]
    dev <- dev[1] - dev
    dev[1] <- NA
    nas <- !is.na(dev)
    P <- dev
    P[nas] <- safe_pchisq(dev[nas], dfs[nas], lower.tail = FALSE)
    aod[, c("LRT", "Pr(Chi)")] <- list(dev, P)
  }
  aod <- aod[o, ]
  head <- c("Single term additions for", what,"\nModel:", deparse(as.vector(formula(object,what))))
  if (scale > 0) 
    head <- c(head, paste("\nscale: ", format(scale), "\n"))
  class(aod) <- c("anova", "data.frame")
  attr(aod, "heading") <- head
  aod
}
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# main stepGAIC function starts here
#-------------------------------------------------------------------------------
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
      what <- if (!is.null(parameter))  {
    match.arg(parameter, choices=c("mu", "sigma", "nu", "tau"))} else  match.arg(what)
     Terms <- terms(object, what)
if (what=="mu")
{
       object$formula <- Terms 
  object$call$formula <- Terms
} 
else
{
  object[[paste(what,"formula",sep=".")]] <- Terms
  object[[paste(what,"formula",sep=".")]][[2]]<-NULL
  if (paste(what, "formula", sep=".")%in%names(object$call)) 
    object$call[[paste(what,"formula",sep=".")]] <- formula(Terms)[-2]
  else ##this is when the sigma formula is not defined
  {object$call[[paste(what,"formula",sep=".")]] <- formula(Terms)[-2]
   names(object$call)[length(names(object$call))]<-paste(what,"formula",sep=".")
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
      attr(terms(update.formula(formula(object, what=what), fdrop), what = what), "factors")
    else numeric(0)
    fadd <- if (!is.null(fadd <- scope$upper)) 
      attr(terms(update.formula(formula(object, what=what), fadd), what = what), "factors")
  }
  else 
  {
    fadd <- if (!is.null(fadd <- scope)) 
      attr(terms(update.formula(formula(object, what=what), scope), what = what ), "factors")
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
bAIC <- extractAIC(fit, scale, k = k,  ...)
edf <- bAIC[1]
bAIC <- bAIC[2]
if (is.na(bAIC)) 
  stop("AIC is not defined for this model, so stepAIC cannot proceed")
nm <- 1
Terms <- terms(fit, what)
if (trace)
  cat("Distribution parameter: ", what, "\n") 
cat("Start:  AIC=", format(round(bAIC, 2)), "\n", cut.string(deparse(as.vector(formula(fit, what=what)))), 
    "\n\n")
models[[nm]] <- list(deviance = mydeviance(fit), df.resid = n - 
                       edf, change = "", AIC = bAIC)
if (!is.null(keep)) 
  keep.list[[nm]] <- keep(fit, bAIC)
usingCp <- FALSE
while (steps > 0) 
{
  steps <- steps - 1
  AIC <- bAIC
  ffac <- attr(Terms, "factors")
  if (!is.null(sp <- attr(Terms, "specials")) && !is.null(st <- sp$strata)) 
    ffac <- ffac[-st, ]
  scope <- factor.scope(ffac, list(add = fadd, drop = fdrop))
  aod <- NULL
  change <- NULL
  if (backward && length(scope$drop)) 
  {
    aod <- droptermP(fit, scope$drop, what = what, scale = scale, 
                     trace = max(0,trace - 1), k = k, parallel = parallel, 
                     ncpus = ncpus, cl = cl)
    rn <- row.names(aod)
    row.names(aod) <- c(rn[1], paste("-", rn[-1], sep = " "))
    if (any(aod$Df == 0, na.rm = TRUE)) 
    {
      zdf <- aod$Df == 0 & !is.na(aod$Df)
      nc <- match(c("Cp", "AIC"), names(aod))
      nc <- nc[!is.na(nc)][1]
      ch <- abs(aod[zdf, nc] - aod[1, nc]) > 0.01
      if (any(ch)) 
      {
        warning("0 df terms are changing AIC")
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
      aodf <- addtermP(fit, scope$add, what =what, scale = scale, 
                       trace = max(0, trace - 1), k = k,  parallel = parallel, 
                       ncpus = ncpus, cl = cl)
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
    nc <- match(c("Cp", "AIC"), names(aod))
    nc <- nc[!is.na(nc)][1]
    o <- order(aod[, nc])
    if (trace) 
      print(aod[o, ])
    if (o[1] == 1) 
      break
    change <- rownames(aod)[o[1]]
  }
  usingCp <- match("Cp", names(aod), 0) > 0
  fit <- update(fit, paste("~ .", change), what = what,  evaluate = FALSE, trace = FALSE) #MS
  # the final at this stage fit
  fit <- eval.parent(fit)
  if (is.list(fit) && (nmm <- match("nobs", names(fit), 0)) > 0) 
    nnew <- fit[[nmm]]
  else nnew <- length(residuals(fit))
  if (nnew != n) 
    stop("number of rows in use has changed: remove missing values?")
  Terms <- terms(fit, what)
  bAIC <- extractAIC(fit, scale, k = k, ...)
  edf <- bAIC[1]
  bAIC <- bAIC[2]
  if (trace) 
    cat("\nStep:  AIC=", format(round(bAIC, 2)), "\n", 
        cut.string(deparse(as.vector(formula(fit, what)))), 
        "\n\n")
  if (bAIC >= AIC + 1e-07) 
    break
  nm <- nm + 1
  models[[nm]] <- list(deviance = mydeviance(fit), df.resid = n - 
                         edf, change = change, AIC = bAIC)
  if (!is.null(keep)) 
    keep.list[[nm]] <- keep(fit, bAIC)
}
if (!is.null(keep)) 
  fit$keep <- re.arrange(keep.list[seq(nm)])
step.results(models = models[seq(nm)], fit, object, usingCp)
}
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
