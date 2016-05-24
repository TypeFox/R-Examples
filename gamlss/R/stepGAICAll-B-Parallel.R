# require(parallel)
# created Thusday 7-4-2015
# author: Mikis Stasinopoulos
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
drop1All <- function (object, 
                           scope,  
                            test = c("Chisq", "none"), 
                               k = 2, 
                          sorted = FALSE, 
                           trace = FALSE, 
                        parallel = c("no", "multicore", "snow"), 
                           ncpus = 1L, 
                              cl = NULL, 
                           ...) 
{
#-----------------------
  drop1.scope<-function (terms1, terms2) 
  {
    terms1 <- terms(terms1, "mu")
    f2 <- if (missing(terms2)) 
        numeric(0)
    else attr(terms(terms2, "mu"), "factor")
    factor.scope(attr(terms1, "factor"), list(drop = f2))$drop
  }
 #-----------------------------------------------------------------------------
  safe_pchisq <- function (q, df, ...) 
  {
    df[df <= 0] <- NA
    pchisq(q = q, df = df, ...)
  }
#-------------------------------------------------------------------------------
#   main functions starts here 
#-------------------------------------------------------------------------------
    #if (!what %in% object$par) 
    #    stop(paste(what, "is not a parameter in the object", "\n"))
    tl <- attr(terms(object, "mu"), "term.labels")
    if (missing(scope)) 
        {scope <- drop1.scope(object)}
    else 
       {
        if (!is.character(scope)) 
            scope <- attr(terms(update.formula(formula(object, "mu"), scope), "mu"), 
                "term.labels")
        if (!all(match(scope, tl, FALSE))) 
            stop("scope is not a subset of term labels")
       }
    ns <- length(scope)
    ans <- matrix(nrow = ns + 1, ncol = 2, dimnames = list(c("<none>", 
        scope), c("df", "AIC")))
    ans[1, ] <- extractAIC(object, scale, k = k, ...)
#-------------------------------------------------------------------------------
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
#------------------------------------------------------------------------------- 
#  function for parallel apply
  fn <- function(term)
  {
    if (trace) 
      cat("trying -", term, "\n")
          nfit <- update(object,  as.formula(paste("~ . -", term)), what="All", evaluate = FALSE, trace=FALSE)
    nfit <- try(eval.parent(nfit), silent=TRUE)
    if (any(class(nfit)%in%"try-error"))
    { 
      cat("Model with term ", term, "has failed \n")       
      NA# extractAIC(object, scale, k = k, ...)          
    }
    else  extractAIC(nfit, scale, k = k,   ...)
  }
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# --------  parallel -----------------------------------------------------------

 #   for (i in seq(ns)) 
  #    {
#        tt <- scope[i]
 #       if (trace) 
#            cat("trying -", tt, "\n")
#        nfit <- update(object,  as.formula(paste("~ . -", tt)), what="All", evaluate = FALSE, #trace=FALSE)
        
       # nfit <- eval.parent(nfit)
#        nfit <- try(eval.parent(nfit), silent=TRUE)
 #             if (any(class(nfit)%in%"try-error"))
  #               { 
 #                 cat("Model with term ", tt, "has failed \n")       
 #                ans[i + 1, ] <- NA# extractAIC(object, scale, k = k, ...)          
 #                }
 #           else ans[i + 1, ] <- extractAIC(nfit, scale, k = k, ...)
 #     }
   #-------------------------------------------------------------------------------------   
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
    head <- c("Single term deletions", "\nModel:", deparse(as.vector(formula(object))))
    class(aod) <- c("anova", "data.frame")
    attr(aod, "heading") <- head
    aod
}
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#  addterm
add1All   <- function (object, 
                           scope,
                           test = c("Chisq", "none" ), 
                              k = 2, 
                         sorted = FALSE, 
                          trace = FALSE,
                       parallel = c("no", "multicore", "snow"), 
                          ncpus = 1L, 
                             cl = NULL,                       
                           ...) 
{
#-------------------------------------------------------------------------------
 add.scope <- function (terms1, terms2 ) 
    {
    terms1 <- terms(terms1)
    terms2 <- terms(terms2)
    factor.scope(attr(terms1, "factor"), list(add = attr(terms2, "factor")))$add  
   }
#------------------------------------------------------------------------------
#-----------------------------------------------------------------------------
safe_pchisq <- function (q, df, ...) 
{
  df[df <= 0] <- NA
  pchisq(q = q, df = df, ...)
}
#-----------------------------------------------------------------------------
    if (missing(scope) || is.null(scope)) 
        stop("no terms in scope")
    if (!is.character(scope))
        scope <- add.scope(object, terms(update.formula(formula(object, "mu"), scope)))

    if (!length(scope)) 
        stop("no terms in scope for adding to object")
    ns <- length(scope)
    ans <- matrix(nrow = ns + 1, ncol = 2, dimnames = list(c("<none>", 
        scope), c("df", "AIC")))
    ans[1, ] <- extractAIC(object, scale, k = k, ...)
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
  nfit <- update(object, as.formula(paste("~ . +", term)), what="All", trace=FALSE, evaluate = FALSE)
  nfit <- try(eval.parent(nfit), silent=TRUE)
  if (any(class(nfit)%in%"try-error"))
  { 
    cat("Model with term ", term, "has failed \n")       
    NA# extractAIC(object, scale, k = k, ...)          
  }
  else  extractAIC(nfit, scale, k = k,   ...)
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

#     for (i in seq(ns)) 
#       {
#         tt <- scope[i]
#         if (trace) 
#             cat("trying +", tt, "\n")
#         nfit <- update(object, as.formula(paste("~ . +", tt)), what="All", trace=FALSE, evaluate = FALSE)
#         nfit <- try(eval.parent(nfit), silent=TRUE)
#               if (any(class(nfit)%in%"try-error"))
#                  { 
#                   cat("Model with term ", tt, "has failed \n")       
#                   ans[i + 1, ] <- NA# extractAIC(object, scale, k = k, ...)          
#                  }
#             else ans[i + 1, ] <- extractAIC(nfit, scale, k = k, ...)
#       }
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
    head <- c("Single term additions for", "\nModel:", deparse(as.vector(formula(object))))
    class(aod) <- c("anova", "data.frame")
    attr(aod, "heading") <- head
    aod
}
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
# Venable and Ripley modification of the stepAIC function
#----------------------------------------------------------------------------------------
stepGAICAll.B <-function(object, 
                    scope, 
                direction = c("both", "backward", "forward"), 
                    trace = T, 
                     keep = NULL, 
                    steps = 1000,
                    scale = 0,
                        k = 2,
                 parallel = c("no", "multicore", "snow"), 
                    ncpus = 1L, 
                       cl = NULL, 
                       ...)                    
                    
{
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
    mydeviance <- function(x, ...) 
      {
        dev <- deviance(x)
        if (!is.null(dev)) 
            dev
        else extractAIC(x, k = 0)[2]
      }
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
    cut.string <- function(string) 
      {
        if (length(string) > 1) 
            string[-1] <- paste("\n", string[-1], sep = "")
        string
      }
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
    re.arrange <- function(keep) 
      {
        namr <- names(k1 <- keep[[1]])
        namc <- names(keep)
          nc <- length(keep)
          nr <- length(k1)
        array(unlist(keep, recursive = FALSE), c(nr, nc), list(namr, 
            namc))
      }
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
    step.results <- function(models, fit, object, usingCp = FALSE) #
      {
        change <- sapply(models, "[[", "change")
            rd <- sapply(models, "[[", "deviance")
            dd <- c(NA, abs(diff(rd)))
           rdf <- sapply(models, "[[", "df.resid")
           ddf <- c(NA, abs(diff(rdf)))
           AIC <- sapply(models, "[[", "AIC")
       heading <- c("Stepwise Model Path \nAnalysis of Deviance Table", 
            "\nInitial  Model:", deparse(as.vector(formula(object))), 
            "\nFinal   Model:", deparse(as.vector(formula(fit))), 
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
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
droptermAllP <- function (object, 
                       scope,  
                       test = c("Chisq", "none"), 
                          k = 2, 
                     sorted = FALSE, 
                      trace = FALSE, 
                   parallel = c("no", "multicore", "snow"), 
                      ncpus = 1L, 
                         cl = NULL, 
                       ...) 
{
  #-----------------------
  drop1.scope<-function (terms1, terms2) 
  {
    terms1 <- terms(terms1, "mu")
    f2 <- if (missing(terms2)) 
      numeric(0)
    else attr(terms(terms2, "mu"), "factor")
    factor.scope(attr(terms1, "factor"), list(drop = f2))$drop
  }
  #-----------------------------------------------------------------------------
  safe_pchisq <- function (q, df, ...) 
  {
    df[df <= 0] <- NA
    pchisq(q = q, df = df, ...)
  }
  #-------------------------------------------------------------------------------
  #   main function droptermAllP  starts here 
  #-------------------------------------------------------------------------------
  tl <- attr(terms(object, "mu"), "term.labels")
  if (missing(scope)) 
  {scope <- drop1.scope(object)}
  else 
  {
    if (!is.character(scope)) 
      scope <- attr(terms(update.formula(formula(object, "mu"), scope), "mu"), 
                    "term.labels")
    if (!all(match(scope, tl, FALSE))) 
      stop("scope is not a subset of term labels")
  }
  ns <- length(scope)
  ans <- matrix(nrow = ns + 1, ncol = 2, dimnames = list(c("<none>", 
                                                           scope), c("df", "AIC")))
  ans[1, ] <- extractAIC(object, scale, k = k, ...)
  #-------------------------------------------------------------------------------
  #--------------- PARALLEL-------------------------------------------------------
  #------------------------------------------------------------------------------- 
  #  function for parallel apply
  fn <- function(term)
  {
    if (trace) 
      cat("trying -", term, "\n")
    nfit <- update(object,  as.formula(paste("~ . -", term)), what="All", evaluate = FALSE, trace=FALSE)
    nfit <- try(eval.parent(nfit), silent=TRUE)
    if (any(class(nfit)%in%"try-error"))
    { 
      cat("Model with term ", term, "has failed \n")       
      NA# extractAIC(object, scale, k = k, ...)          
    }
    else  extractAIC(nfit, scale, k = k,   ...)
  }
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# --------  parallel -----------------------------------------------------------
#-------------------------------------------------------------------------------------   
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
  head <- c("Single term deletions", "\nModel:", deparse(as.vector(formula(object))))
  class(aod) <- c("anova", "data.frame")
  attr(aod, "heading") <- head
  aod
}
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
addtermAllP   <- function (object, 
                        scope,
                        test = c("Chisq", "none" ), 
                        k = 2, 
                        sorted = FALSE, 
                        trace = FALSE,
                        parallel = c("no", "multicore", "snow"), 
                        ncpus = 1L, 
                        cl = NULL,                       
                        ...) 
{
  #-------------------------------------------------------------------------------
  add.scope <- function (terms1, terms2 ) 
  {
    terms1 <- terms(terms1)
    terms2 <- terms(terms2)
    factor.scope(attr(terms1, "factor"), list(add = attr(terms2, "factor")))$add  
  }
  #------------------------------------------------------------------------------
  #-----------------------------------------------------------------------------
  #-----------------------------------------------------------------------------
  safe_pchisq <- function (q, df, ...) 
  {
    df[df <= 0] <- NA
    pchisq(q = q, df = df, ...)
  }
  #-----------------------------------------------------------------------------
  #-----------------------------------------------------------------------------
  if (missing(scope) || is.null(scope)) 
    stop("no terms in scope")
  if (!is.character(scope))
    scope <- add.scope(object, terms(update.formula(formula(object, "mu"), scope)))
  
  if (!length(scope)) 
    stop("no terms in scope for adding to object")
  ns <- length(scope)
  ans <- matrix(nrow = ns + 1, ncol = 2, dimnames = list(c("<none>", 
                                                           scope), c("df", "AIC")))
  ans[1, ] <- extractAIC(object, scale, k = k, ...)

  #-------------------------------------------------------------------------------
  #  function for parallel apply
  fn <- function(term)
  {
    if (trace) 
      cat("trying -", term, "\n")
    nfit <- update(object, as.formula(paste("~ . +", term)), what="All", trace=FALSE, evaluate = FALSE)
    nfit <- try(eval.parent(nfit), silent=TRUE)
    if (any(class(nfit)%in%"try-error"))
    { 
      cat("Model with term ", term, "has failed \n")       
      NA# extractAIC(object, scale, k = k, ...)          
    }
    else  extractAIC(nfit, scale, k = k,   ...)
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
#        cl <- parallel::makeForkCluster(ncpus)
#        if (RNGkind()[1L] == "L'Ecuyer-CMRG") 
#          parallel::clusterSetRNGStream(cl)
        res <- t(parallel::parSapply(cl, scope, fn))
#        parallel::stopCluster(cl)
        res
      }
      else t(parallel::parSapply(cl, scope, fn))
    }
  } # end parallel ----------------------------------------------------------
  else t(sapply(scope, fn)) 
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
  head <- c("Single term additions for", "\nModel:", deparse(as.vector(formula(object))))
  class(aod) <- c("anova", "data.frame")
  attr(aod, "heading") <- head
  aod
}

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# main function starts here
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
    Terms <- terms(object)
           object$formula <- Terms 
      object$call$formula <- Terms 
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
                attr(terms(update.formula(formula(object, what="mu"), fdrop), what = "mu"), "factors")
            else numeric(0)
             fadd <- if (!is.null(fadd <- scope$upper)) 
                attr(terms(update.formula(formula(object, what="mu"), fadd), what = "mu"), "factors")
          }
        else 
          {
             fadd <- if (!is.null(fadd <- scope)) 
                attr(terms(update.formula(formula(object, what="mu"), scope), what = "mu" ), "factors")
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
    bAIC <- extractAIC(fit, scale, k = k, ...)
     edf <- bAIC[1]
    bAIC <- bAIC[2]
    if (is.na(bAIC)) 
        stop("AIC is not defined for this model, so stepAIC cannot proceed")
      nm <- 1
   Terms <- terms(fit, "mu")
    if (trace)
        #cat("Distribution parameter: ", what, "\n") 
        cat("Start:  AIC=", format(round(bAIC, 2)), "\n", cut.string(deparse(as.vector(formula(fit, what="mu")))), 
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
            aod <- droptermAllP(fit, scope$drop, trace = max(0,trace - 1), k = k,
                    test="none",  parallel = parallel, ncpus = ncpus, cl = cl)
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
                aodf <- addtermAllP(fit, scope$add, trace = max(0, trace - 1), k = k,
                        test="none",parallel = parallel, ncpus = ncpus, cl = cl)
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
        fit <- update(fit, paste("~ .", change), evaluate = FALSE, what="All", trace = FALSE) #MS
        fit <- eval.parent(fit)
        if (is.list(fit) && (nmm <- match("nobs", names(fit), 
            0)) > 0) 
            nnew <- fit[[nmm]]
        else nnew <- length(residuals(fit))
        if (nnew != n) 
            stop("number of rows in use has changed: remove missing values?")
        Terms <- terms(fit, "mu")
         bAIC <- extractAIC(fit, scale, k = k, ...)
          edf <- bAIC[1]
         bAIC <- bAIC[2]
        if (trace) 
            cat("\nStep:  AIC=", format(round(bAIC, 2)), "\n", 
                cut.string(deparse(as.vector(formula(fit, "mu")))), 
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
