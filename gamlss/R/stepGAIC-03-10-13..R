# created Wednesday, September 15, 2004 
# author: Mikis Stasinopoulos
# based on functions in MASS 
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
# extractAIC  needs MASS
extractAIC.gamlss<-function (fit, scale, k = 2, c = FALSE,  ...) 
{
    if (is.gamlss(fit)) 
       {
        edf <- fit$df.fit
          N <- fit$N
        Cor <- if ((k == 2)&&(c==TRUE)) (2*edf*(edf+1))/(N-edf-1) else 0 
        aic <- fit$G.dev + fit$df.fit * k + Cor
        return( c(edf, aic))
       }
    else stop(paste("this is not a gamlss object"))
}
#----------------------------------------------------------------------------------------
#library(MASS)
# these are not needed with NAMESPACE
#dropterm <- function (object, ...) 
#UseMethod("dropterm")
#addterm <- function (object, ...) 
#UseMethod("addterm")
# dropterm extractAIC  needs MASS
# dropterm.gamlss<-function (object, 
#                            scope, 
#                            what = c("mu", "sigma", "nu", "tau"),
#                            scale = 0, 
#                            test = c("none", "Chisq"), 
#                            k = 2, 
#                            sorted = FALSE, 
#                            trace = FALSE, 
#                            ...) 
# {
# #-------------------------------------------------------------------------------
#   drop1.scope<-function (terms1, terms2, what = c("mu", "sigma", "nu", "tau")) 
#   {
#     what <- match.arg(what)
#     terms1 <- terms(terms1, what)
#     f2 <- if (missing(terms2)) 
#         numeric(0)
#     else attr(terms(terms2, what), "factor")
#     factor.scope(attr(terms1, "factor"), list(drop = f2))$drop
#   }
# #-------------------------------------------------------------------------------
# safe_pchisq <- function (q, df, ...) 
# {
#   df[df <= 0] <- NA
#   pchisq(q = q, df = df, ...)
# }
# #-------------------------------------------------------------------------------
#     what <- match.arg(what)
#     if (!what %in% object$par) 
#         stop(paste(what, "is not a parameter in the object", "\n"))
#     tl <- attr(terms(object, what ), "term.labels")
#     if (missing(scope)) 
#         scope <- drop1.scope(object, what = what)
#     else 
#        {
#         if (!is.character(scope)) 
#             scope <- attr(terms(update.formula(formula(object, what=what), scope), what), 
#                 "term.labels")
#         if (!all(match(scope, tl, FALSE))) 
#             stop("scope is not a subset of term labels")
#        }
#     ns <- length(scope)
#     ans <- matrix(nrow = ns + 1, ncol = 2, dimnames = list(c("<none>", 
#         scope), c("df", "AIC")))
#     ans[1, ] <- extractAIC(object, scale, k = k,   ...)
#     for (i in seq(ns)) 
#       {
#         tt <- scope[i]
#         if (trace) 
#             cat("trying -", tt, "\n")
#             
#         nfit <- update(object,  as.formula(paste("~ . -", tt)), what = what, evaluate = FALSE, trace=FALSE) 
#         nfit <- try(eval.parent(nfit), silent=TRUE)
#               if (any(class(nfit)%in%"try-error"))
#                  { 
#                   cat("Model with term ", tt, "has failed \n")       
#                   ans[i + 1, ] <- NA# extractAIC(object, scale, k = k, ...)          
#                  }
#             else ans[i + 1, ] <- extractAIC(nfit, scale, k = k,   ...)
#       }
#     dfs <- ans[1, 1] - ans[, 1]
#     dfs[1] <- NA
#     aod <- data.frame(Df = dfs, AIC = ans[, 2])
#     o <- if (sorted) 
#         order(aod$AIC)
#     else seq(along = aod$AIC)
#     test <- match.arg(test)
#     if (test == "Chisq") 
#       {
#         dev <- ans[, 2] - k * ans[, 1]
#         dev <- dev - dev[1]
#         dev[1] <- NA
#         nas <- !is.na(dev)
#         P <- dev
#         P[nas] <- safe_pchisq(dev[nas], dfs[nas], lower.tail = FALSE)
#         aod[, c("LRT", "Pr(Chi)")] <- list(dev, P)
#       }
#     aod <- aod[o, ]
#     head <- c("Single term deletions for", what, "\nModel:", deparse(as.vector(formula(object, what))))
#     if (scale > 0) 
#         head <- c(head, paste("\nscale: ", format(scale), "\n"))
#     class(aod) <- c("anova", "data.frame")
#     attr(aod, "heading") <- head
#     aod
# }
#----------------------------------------------------------------------------------------
#  addterm
# addterm.gamlss <- function (object, 
#                            scope,
#                            what = c("mu", "sigma", "nu", "tau"), 
#                           scale = 0, 
#                            test = c("none", "Chisq"), 
#                               k = 2, 
#                          sorted = FALSE, 
#                           trace = FALSE, 
#                            ...) 
# {
# #-------------------------------------------------------------------------------
#  add.scope <- function (terms1, terms2, what = c("mu", "sigma", "nu", "tau") ) 
#     {
#       what <- match.arg(what)
#     terms1 <- terms(terms1, what)
#     terms2 <- terms(terms2, what)
#     factor.scope(attr(terms1, "factor"), list(add = attr(terms2, 
#         "factor")))$add  
#    }
# #-------------------------------------------------------------------------------
# safe_pchisq <- function (q, df, ...) 
# {
#   df[df <= 0] <- NA
#   pchisq(q = q, df = df, ...)
# }
# #-------------------------------------------------------------------------------
#    what <- match.arg(what)
#     if (missing(scope) || is.null(scope)) 
#         stop("no terms in scope")
#     if (!is.character(scope)) 
#         scope <- add.scope(object, terms(update.formula(formula(object, what=what), scope)), what = what)
#     if (!length(scope)) 
#         stop("no terms in scope for adding to object")
#     ns <- length(scope)
#     ans <- matrix(nrow = ns + 1, ncol = 2, dimnames = list(c("<none>", 
#         scope), c("df", "AIC")))
#     ans[1, ] <- extractAIC(object, scale, k = k,   ...)
#     for (i in seq(ns)) 
#       {
#         tt <- scope[i]
#         if (trace) 
#             cat("trying +", tt, "\n")
#         nfit <- update(object, as.formula(paste("~ . +", tt)), what = what, trace=FALSE, 
#             evaluate = FALSE)
#         nfit <- try(eval.parent(nfit), silent=TRUE)
#              if (any(class(nfit)%in%"try-error"))
#                  { 
#                   cat("Model with term ", tt, "has failed \n")       
#                   ans[i + 1, ] <- NA# extractAIC(object, scale, k = k, ...)          
#                  }
#             else ans[i + 1, ] <- extractAIC(nfit, scale, k = k,  ...)
#       }
#     dfs <- ans[, 1] - ans[1, 1]
#     dfs[1] <- NA
#     aod <- data.frame(Df = dfs, AIC = ans[, 2])
#     o <- if (sorted) 
#         order(aod$AIC)
#     else seq(along = aod$AIC)
#     test <- match.arg(test)
#     if (test == "Chisq") 
#       {
#         dev <- ans[, 2] - k * ans[, 1]
#         dev <- dev[1] - dev
#         dev[1] <- NA
#         nas <- !is.na(dev)
#         P <- dev
#         P[nas] <- safe_pchisq(dev[nas], dfs[nas], lower.tail = FALSE)
#         aod[, c("LRT", "Pr(Chi)")] <- list(dev, P)
#       }
#     aod <- aod[o, ]
#     head <- c("Single term additions for", what,"\nModel:", deparse(as.vector(formula(object,what))))
#     if (scale > 0) 
#         head <- c(head, paste("\nscale: ", format(scale), "\n"))
#     class(aod) <- c("anova", "data.frame")
#     attr(aod, "heading") <- head
#     aod
# }
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
# Venable and Ripley step AIC function
#----------------------------------------------------------------------------------------
stepGAIC.VR <-function(object, 
                    scope, 
                direction = c("both", "backward", "forward"), 
                    trace = T, 
                     keep = NULL, 
                    steps = 1000,
                    scale = 0,
                     what = c("mu", "sigma", "nu", "tau"),
                     parameter= NULL, 
                        k = 2,
                       ...)                    
                    
{
#--------------------------------------------------------------------
    mydeviance <- function(x, ...) 
      {
        dev <- deviance(x)
        if (!is.null(dev)) 
            dev
        else extractAIC(x, k = 0)[2]
      }
#--------------------------------------------------------------------
    cut.string <- function(string) 
      {
        if (length(string) > 1) 
            string[-1] <- paste("\n", string[-1], sep = "")
        string
      }
#--------------------------------------------------------------------
    re.arrange <- function(keep) 
      {
        namr <- names(k1 <- keep[[1]])
        namc <- names(keep)
          nc <- length(keep)
          nr <- length(k1)
        array(unlist(keep, recursive = FALSE), c(nr, nc), list(namr, 
            namc))
      }
#--------------------------------------------------------------------
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
#----------------------------------------------------------------------
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
            aod <- dropterm(fit, scope$drop, what = what, scale = scale, 
                           trace = max(0,trace - 1), k = k, ...)
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
                aodf <- addterm(fit, scope$add, what =what, scale = scale, 
                  trace = max(0, trace - 1), k = k, ...)
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
#----------------------------------------------------------------------------------------
# created Wednesday, September 15, 2004 
# author: Mikis Stasinopoulos
# based on the step.gam() function of Trevor Hastie
#----------------------------------------------------------------------------------------
stepGAIC.CH <-function(object, 
                    scope = gamlss.scope(model.frame(object)), 
                direction = c("both", "backward", "forward"), 
                    trace = T, 
                     keep = NULL, 
                    steps = 1000,
                     what = c("mu", "sigma", "nu", "tau"),
                     parameter= NULL,
                        k = 2,
                      ...)
{
# -------------------------------------------------------------------
scope.char <- function(form)
     {
     form <- as.character(form)
     form <- unlist(strsplit(form[2], "\\+"))
       f1 <- lapply(form, sub, pattern=" ",replacement= "")
     form <- lapply(f1, sub, pattern=" ",replacement= "")
     form
     } 
# -------------------------------------------------------------------
re.arrange <- function(keep)
       {
        namr <- names(k1 <- keep[[1]])
        namc <- names(keep)
          nc <- length(keep)
          nr <- length(k1)
        array(unlist(keep, recursive = F), c(nr, nc), list(namr, namc))
       }
# -------------------------------------------------------------------
untangle.scope <- function(terms, regimens)
       {
               a <- attributes(terms)
        response <- deparse(a$variables[[2]])
     term.labels <- a$term.labels
              nt <- length(regimens)
          select <- integer(nt)
        for(i in seq(nt)) 
          {
            j <- match(regimens[[i]], term.labels, 0)
            if(any(j)) 
            {
                if(sum(j > 0) > 1)
                stop(paste("The elements of a regimen", i, "appear more than once in the initial model", sep = " "))
                select[i] <- seq(j)[j > 0]
                term.labels <- term.labels[ - sum(j)]
            }
            else 
            {
                if(!(j <- match("1", regimens[[i]], 0)))
                stop(paste("regimen", i, "does not appear in the initial model", sep = " "))
                select[i] <- j
            }
          }
        if(length(term.labels))
            term.labels <- paste(term.labels, "+")
        return(list( response = paste(response, term.labels, sep = " ~ "), select=select))
        }
# -------------------------------------------------------------------
make.step <- function(models, fit, scale, object)
    {
      as.anova <- function(df, heading)
         {
           if(!inherits(df, "data.frame"))
           stop("df must be a data frame")
           attr(df, "heading") <- heading
           #class(df) <- c("anova", class(df))
           df
         }  
        chfrom <- sapply(models, "[[", "from")
        chfrom[chfrom == "1"] <- ""
        chto <- sapply(models, "[[", "to")
        chto[chto == "1"] <- ""
        dev <- sapply(models, "[[", "deviance")
        df <- sapply(models, "[[", "df.resid")
        ddev <- c(NA, diff(dev))
        ddf <- c(NA, diff(df))
        AIC <- sapply(models, "[[", "AIC")
        heading <- c("Stepwise Model Path \nAnalysis of Deviance Table", "\nInitial Model:", deparse(as.vector(formula(object))),
            "\nFinal Model:", deparse(as.vector(formula(fit))), paste("\nScale: ", format(scale), "\n", sep = ""))
        aod <- data.frame(From = chfrom, To = chto, Df = ddf, Deviance = ddev, "Resid. Df" = df, "Resid. Dev" = dev, AIC = AIC, 
            check.names = F)
        fit$anova <- as.anova(aod, heading)
        fit
    }
# -------------------------------------------------------------------
# the function starts here
what <- if (!is.null(parameter))  {
    match.arg(parameter, choices=c("mu", "sigma", "nu", "tau"))} else  match.arg(what)    
direction <- match.arg(direction)
    #if(missing(scope))
    #    stop("you must supply a scope argument to step.gam(); the gamlss.scope() function might be useful")
    if(!is.character(scope[[1]]))
        scope <- lapply(scope, scope.char)
    response <- untangle.scope(terms(object,what), scope)
    form.y <- response$response
   # if (what=="mu") 
   # { form.y <- response$response }# MS Monday, October 25, 2004 
   # else  
   # {
    #    form.y <- if (length(response$response)==2) response$response
    #            else response$response[-2]
    #}
    backward <- direction == "both" | direction == "backward"
    forward <- direction == "both" | direction == "forward"
    items <- response$select
    Call <- object$call
    term.lengths <- sapply(scope, length)
    n.items <- length(items)
    visited <- array(F, term.lengths)
    visited[array(items, c(1, n.items))] <- T
    if(!is.null(keep)) 
    {
        keep.list <- vector("list", length(visited))
        nv <- 1
    }
    models <- vector("list", length(visited))
    nm <- 2
    form.vector <- character(n.items)
    for(i in seq(n.items))
        form.vector[i] <- scope[[i]][items[i]]
    form <- deparse(formula(object, what))
    if(trace)
        cat("Distribution parameter: ", what, "\n")
        cat("Start: ", form)
    fit <- object
    n <- length(fit$N)
    bAIC <- AIC(object, k=k)
    if(trace)
        cat("; AIC=", format(round(bAIC, 4)), "\n")
    models[[1]] <- list(deviance = deviance(fit), df.resid = fit$df.resid, AIC = bAIC, from = "", to = "")
    if(!is.null(keep)) 
    { 
        keep.list[[nv]] <- keep(fit, bAIC)
        nv <- nv + 1
    }
    AIC <- bAIC + 1
    while(bAIC < AIC & steps > 0) 
    {
        steps <- steps - 1
        AIC <- bAIC
        bitems <- items
        bfit <- fit
        for(i in seq(n.items)) 
        {
# try go down a level
            if(backward) 
            
            {
                trial <- items
                trial[i] <- trial[i] - 1
                if(trial[i] > 0 && !visited[array(trial, c(1, n.items))]) 
                {
                  visited[array(trial, c(1, n.items))] <- T
                  tform.vector <- form.vector
                  tform.vector[i] <- scope[[i]][trial[i]]
                  form <- paste(form.y, paste(tform.vector, collapse = " + "))
                  if(trace)
                    cat("Trial: ", form)
                  tfit <- update(object, eval(parse(text = form)), what = what, trace = F, ...)
                  tAIC <- IC(tfit, k=k)
                  if(!is.null(keep)) 
                  {
                    keep.list[[nv]] <- keep(tfit, tAIC)
                    nv <- nv + 1
                  }
                  if(tAIC < bAIC) 
                  {
                    bAIC <- tAIC
                    bitems <- trial
                    bfit <- tfit
                    bform.vector <- tform.vector
                    bfrom <- form.vector[[i]]
                    bto <- tform.vector[[i]]
                  }
                  if(trace)
                    cat("; AIC=", format(round(tAIC, 4)), "\n")
                }
            }
            if(forward) 
            {
                trial <- items
                trial[i] <- trial[i] + 1
                if(trial[i] <= term.lengths[i] && !visited[array(trial, c(1, n.items))]) 
                {
                  visited[array(trial, c(1, n.items))] <- T
                  tform.vector <- form.vector
                  tform.vector[i] <- scope[[i]][trial[i]]
                  form <- paste(form.y, paste(tform.vector, collapse = " + "))
                  if(trace)
                    cat("Trial: ", form)
                  tfit <- update(object, eval(parse(text = form)), what = what, trace = F, ...)
                  tAIC <- IC(tfit,k=k) 
                  if(!is.null(keep)) 
                  {
                    keep.list[[nv]] <- keep(tfit, tAIC)
                    nv <- nv + 1
                  }
                  if(tAIC < bAIC) 
                  {
                    bAIC <- tAIC
                    bitems <- trial
                    bfit <- tfit
                    bform.vector <- tform.vector
                    bfrom <- form.vector[[i]]
                    bto <- tform.vector[[i]]
                  }
                  if(trace)
                    cat("; AIC=", format(round(tAIC, 4)), "\n")
                }
            }
        }
        if(bAIC >= AIC | steps == 0) 
        {
            if(!is.null(keep))
                fit$keep <- re.arrange(keep.list[seq(nv - 1)])
            return(make.step(models[seq(nm - 1)], fit, scale, object))  #time to quit
        }
        else 
        {
            if(trace)
                cat("Step : ", deparse(formula(bfit,what)), "; AIC=", format(round(bAIC, 4)), "\n\n")
            items <- bitems
            models[[nm]] <- list(deviance = deviance(bfit), df.resid = bfit$df.resid, AIC = bAIC, from = bfrom, to = bto)
            nm <- nm + 1
            fit <- bfit
            form.vector <- bform.vector
        }
    }
}
#----------------------------------------------------------------------------------------
# copy of the gam.scope S-plus function
# its needs or a model.frame or a data frame 
gamlss.scope <- function(frame, 
                      response = 1, 
                      smoother = "cs", 
                      arg = NULL, 
                      form = TRUE)
{
    vnames <- names(frame)
    vnames <- vnames[ - response]
    step.list <- as.list(vnames)
    names(step.list) <- vnames
    for(vname in vnames) {
        junk <- c("1", vname)
        if(is.vector(frame[[vname]]))
            junk <- c(junk, paste(smoother, "(", vname, if(is.null(arg)) ")" else paste(",", arg, ")", sep = ""), sep = ""))
        if(form)
            junk <- eval(parse(text = paste("~", paste(junk, collapse = "+"))))
        step.list[[vname]] <- junk
    }
    step.list
}
#----------------------------------------------------------------------------------------
# stepGAIC <-function(object, 
#                     scope = gamlss.scope(model.frame(object)), 
#                 direction = c("both", "backward", "forward"), 
#                     trace = T, 
#                      keep = NULL, 
#                     steps = 1000,
#                      what = c("mu", "sigma", "nu", "tau"),
#                         k = 2,
#                  additive = FALSE,
#                     ...)
# {
# stepGAICcall <- match.call()  
#     mnames <- c("", "object", "scope", "direction", "trace", "keep", "steps", "what", "k", "..." ) #  
#     cnames <- names(stepGAICcall)  # 
#     cnames <- cnames[match(mnames,cnames,0)] 
#    newcall <- stepGAICcall[cnames] #   
# newcall[[1]] <- if (additive == TRUE)  as.name("stepGAIC.CH") else as.name("stepGAIC.VR")
# finalmodel <-eval(newcall)
# finalmodel 
# }
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
drop1.gamlss<-function(object, ...) dropterm(object, test="Chisq", ...)
 add1.gamlss<-function(object, ...)  addterm(object, test="Chisq", ...)
#----------------------------------------------------------------------------------------