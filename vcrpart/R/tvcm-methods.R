##' -------------------------------------------------------- #
##' Author:          Reto Buergin
##' E-Mail:          rbuergin@gmx.ch
##' Date:            2016-02-16
##'
##' Description:
##' S3 methods for tvcm objects
##'
##' References:
##' party:           http://CRAN.R-project.org/package=party
##' partykit:        http://CRAN.R-project.org/package=partykit
##'
##' Methods:
##' coef, coefficients:
##' depth:               depth of trees
##' extract:
##' fitted:
##' formula:
##' getCall:             extract original call
##' logLik:              extract log Likelihood
##' model.frame:         extract the total model frame including model
##'                      and partitioning variables
##' nobs:                extract the number of observations
##' predict:             predict responses (see prediction of 'olmm' class)
##' print:               print tvcm objects
##' prune:               prune the tree
##' ranef:               extract the predicted random effects
##' resid, residuals:    extract residuals
##' splitpath:           show information of the splitting procedure
##' weights:             extract the weights
##' width:               width of trees
##' 
##' Modifications:
##' 2016-02-16: changed titles for varying coefficients in 'tvcm_print'.
##' 2016-01-10: changed 'neglokLik2.tvcm' to allow for subsets
##' 2015-08-21: implemented changes to 'tvcm_formula' in 'prune.tvcm'.
##' 2014-10-18: added 'depth' and 'width' methods
##' 2014-10-14: adapt print.splitpath to new dev-grid structure
##' 2014-10-03: add option 'cv' to 'extract.tvcm'
##' 2014-09-17: prune.tvcm:
##'             - 'keepdev' argument in 'prune.tvcm' dropped (to complicated
##'               to explain)
##'             - add 'control' argument 
##' 2014-09-15: 'tab', 'ntab' and 'otab' are not matrices
##' 2014-09-02: added 'prunepath.tvcm' and 'print.prunepath.tvcm' functions
##' 2014-09-02: various modifications in 'prune.tvcm' for accelerations:
##'             - 'do.call' was replaced by 'eval'
##'             - new option 'keepdev' (reuses information of the previous
##'               pruning step)
##'             - new option 'papply' (even though the accrelation is not
##'               as efficient as expected)
##' 2014-08-13: modified option 'type = "coef"' of predict.tvcm
##'             to deal with multiple components
##' -------------------------------------------------------- #

coef.tvcm <- function(object, ...) tvcm_get_estimates(object, ...)


coefficients.tvcm <- coef.tvcm 


depth.tvcm <- function(x, root = FALSE, ...) {
  rval <- sapply(x$info$node, depth, root = root)
  names(rval) <- LETTERS[seq_along(rval)]
  return(rval)
}


extract.tvcm <- function(object, what = c("control", "model", 
                                   "nodes", "sctest", "p.value",
                                   "devgrid", "cv",  "selected", 
                                   "coef", "sd", "var"),
                         steps = NULL, ...) {
  
  what <- match.arg(what)
  splitpath <- object$info$splitpath
  if (length(splitpath) > 0 && is.null(steps))
    steps <- seq(1L, object$info$nstep)
  rval <- NULL
  sctest <- object$info$control$sctest
  
  if (what == "control") {

    return(object$info$control)
    
  } else if (what == "sctest" && sctest) {
    
    rval <- lapply(splitpath[steps], function(x) x$sctest)
    return(rval)
    
  } else if (what == "devgrid" && !is.null(splitpath)) {
    
    rval <- lapply(splitpath[steps], function(x) x$grid)
    return(rval)

  } else if (what == "cv") {
    if (is.null(object$info$cv)) {
      warning("no information on cross-validation")
    } else {
      return(object$info$cv)
    }
    
  } else if (what == "model") {
    
    return(object$info$model)
    
  } else if (what == "selected") {

    rval <-
      lapply(object$info$node,
             function(node) {
               if (depth(node) > 0L) {
                 ids <- setdiff(nodeids(node), nodeids(node, terminal = TRUE))
                 rval <- unlist(nodeapply(node, ids, function(node) node$split$varid))
                 if (length(rval) > 0L) rval <- unique(colnames(object$data)[rval])
               } else {
                 rval <- NULL
               }
               return(rval)
             })
    return(rval)
    
  } else if (sctest && what == "p.value" && !is.null(splitpath)) {
    
    rval <- unlist(sapply(splitpath[steps],
                          function(x) {
                            if (is.null(x$sctest)) return(NA)
                            rval <- na.omit(unlist(x$sctest))
                            if (length(rval) == 0) return(NA)
                            return(min(rval, na.rm = TRUE))
                          }))
    return(rval)
    
  } else if (what %in% c("coef", "sd", "var")) {

    return(tvcm_get_estimates(object, what = what, ...))   

  } else if (what == "nodes") {

    return(object$info$node)
    
  }
  return(rval)
}


extractAIC.tvcm <- function(fit, scale, k = 2, ...) {
  extractAIC(extract(fit, "model"), scale, k, ...)
}


fitted.tvcm <- function(object, ...) vcrpart_fitted(object, ...)


formula.tvcm <- function(x, ...) return(x$info$formula$original)


getCall.tvcm <- function(x, ...) return(x$info$call)


logLik.tvcm <- function(object, ...)
  return(logLik(extract(object, "model")))


model.frame.tvcm <- function(formula, ...) {
  rval <- cbind(model.frame(formula$info$model), formula$data)
  attr(rval, "terms") <- attr(formula$data, "terms")
  attr(rval, "na.action") <- attr(formula$data, "na.action")
  return(rval)
}


neglogLik2.tvcm <- function(object, ...)
    return(neglogLik2(extract(object, "model"), ...))


nobs.tvcm <- function(object, ...) nobs(extract(object, "model"), ...)


predict.tvcm <- function(object, newdata = NULL,
                         type = c("link", "response", "prob", "class",
                           "node", "coef", "ranef"),
                         ranef = FALSE, na.action = na.pass, ...) {

  ## match type
  type <- match.arg(type)
  if (type == "prob") type = "response"
  
  ## resolve conflicts with the 'ranef' argument
  if (!is.null(newdata) && is.logical(ranef) && ranef)
    stop("'ranef' should be 'FALSE' or a 'matrix' if 'newdata' is not 'NULL'.")
  if (type == "ranef" & (!is.logical(ranef) | is.logical(ranef) && ranef))
    stop("for 'type = 'ranef'' the argument 'ranef' must be 'FALSE'.")
  if (type == "ranef" & !is.null(newdata))
    stop("prediction for random effect for 'newdata' is not implemented.")
    
  if (type == "ranef") return(ranef(object$info$model, ...))
  
  ## the terminal node identifiers
  ids <- lapply(object$info$node, nodeids, terminal = TRUE)
  
  ## substitute 'newdata' by learning sample
  if (is.null(newdata)) newdata <- model.frame(object)

  ## extract node ids
  fitted <- as.data.frame(tvcm_get_node(object, newdata, formList = object$info$formula))
  newdata[, names(fitted)] <- fitted

  ## return fitted node ids
  if (type == "node") {

    ## return fitted nodes
    return(fitted)
    
  } else if (type == "coef") {
      
    ## extract individual effects
    what <- list(...)$what
    if (is.null(what)) what <- "coef"
    coef <- extract(object, what)
    fe <- vc <- re <- NULL
    
    if (length(coef$fe) > 0L)
      fe <- matrix(coef$fe, nrow = nrow(fitted),
                   ncol = length(coef$fe), byrow = TRUE,
                   dimnames = list(rownames(fitted), names(coef$fe)))

    if (!is.null(coef$vc)) {
      for (pid in 1:length(coef$vc)) {
        if (ncol(coef$vc[[pid]]) > 0L) {
          vcPid <- lapply(fitted[, pid], function(x) coef$vc[[pid]][as.character(x), ])
          vcPid <- matrix(unlist(vcPid), nrow = nrow(fitted), byrow = TRUE)
          rownames(vcPid) <- rownames(fitted)
          colnames(vcPid) <- colnames(coef$vc[[pid]])
          vc <- cbind(vc, vcPid)
        }
      }
    }
    
    if (length(coef$re) > 0L)
      re <- matrix(coef$re, nrow = nrow(fitted),
                   ncol = length(coef$re), byrow = TRUE,
                   dimnames = list(rownames(fitted), names(coef$re)))
    
    rval <- cbind(fe, vc)
    terms <- unique(colnames(rval))
    rval <-
      sapply(terms, function(x)
          rowSums(rval[, colnames(rval) == x, drop = FALSE], na.rm = TRUE))

    vcparm <- unique(unlist(tvcm_get_vcparm(object)))
    parm <- union(colnames(rval), vcparm)
    if (length(misC <- setdiff(parm, colnames(rval))) > 0L) {
        cn <- c(colnames(rval), misC)
        rval <- cbind(rval, matrix(0, nrow(rval), length(misC)))
        colnames(rval) <- cn
    }
    
    rval <- cbind(rval, re)
    
    rownames(rval) <- rownames(fitted)
    
    return(rval)
  } else {

    ## call predict of fitting function
    return(predict(object$info$model, newdata = newdata, type = type,
                   ranef = ranef, na.action = na.action, ...))
  }
  return(fitted)
}


tvcm_print <- function(x, type = c("print", "summary"),
                       etalab = c("int", "char", "eta"), ...) {

  type <- match.arg(type)
  etalab <- match.arg(etalab)
  coef <- extract(x, "coef")
  sd <- extract(x, "sd")
  yLevs <- if (x$info$fit == "olmm") levels(x$info$model$y) else NULL
  
  cat(x$info$title, "\n")
  cat("  Family:", x$info$family$family, x$info$family$link, "\n")
  cat(" Formula:", deparse(x$info$formula$original, 500L)[1L], "\n")  
  if (length(str <- deparseCall(x$info$call$data)) > 0L)
    cat("    Data:", str, "\n")
  if (length(str <- deparseCall(x$info$call$subset)) > 0L)
    cat("  Subset:", str, "\n")
  if (x$info$control$sctest)
    cat(paste("   Tests: alpha = ", format(x$info$control$alpha, ...),
              if (x$info$control$bonferroni) ", nodewise Bonferroni corrected\n",
              sep = ""))
    
  if (length(coef$re) > 0L) {
    cat("\nRandom Effects:\n")
    VarCorr <- VarCorr(extract(x, "model"))
    if (x$info$fit == "olmm") {
      VarCorr <- olmm_rename(VarCorr, yLevs, x$info$family, etalab)
      print.VarCorr.olmm(VarCorr, ...)
    } else {
      print(VarCorr)
    }
  }
    
  if (length(coef$fe) > 0L) {
    cat("\nFixed Effects:\n")
    if (type == "print") {
      coefMat <- matrix(coef$fe, 1)
      colnames(coefMat) <- names(coef$fe) 
       
    } else {
      coefMat <- cbind("Estimate" = coef$fe,
                       "Std. Error" = sd$fe,
                       "z value" = coef$fe / sd$fe)
    }
    if (x$info$fit == "olmm")
      coefMat <- olmm_rename(coefMat, yLevs, x$info$family, etalab)
    print(coefMat, ...)
  }

  terminal_panel <- function(node, partid) {
    rval <- function(node) {
        nid <- as.character(id_node(node))
        if (nrow(coef$vc[[partid]]) > 0L) {
            if (type == "print") {
                coefMat <- matrix(coef$vc[[partid]][nid, ], 1)
                colnames(coefMat) <- colnames(coef$vc[[partid]])
            } else {
                coefMat <- cbind("Estimate" = coef$vc[[partid]][nid, ],
                                 "Std. Error" = sd$vc[[partid]][nid, ],
                                 "z value" = coef$vc[[partid]][nid, ] / sd$vc[[partid]][nid, ])
                rownames(coefMat) <- colnames(coef$vc[[partid]])
            }
            if (x$info$fit == "olmm")
                coefMat <- olmm_rename(coefMat, yLevs, x$info$family, etalab)
            return(c("", unlist(strsplit(formatMatrix(coefMat, ...), "\n"))))
        } else {
            return(NULL) 
        }
    }
    return(rval)
}
    
  class(terminal_panel) <- "grapcon_generator"

  vcLabs <- tvcm_print_vclabs(x$info$formula, TRUE)

  for (pid in seq_along(coef$vc)) {
    cat(paste("\nVarying Coefficient ", LETTERS[pid],
              ": ", vcLabs[pid], "\n", sep = ""))
    x$node <- x$info$node[[pid]]
    print.party(x, terminal_panel = terminal_panel,
                tp_args = list(partid = pid))
    ## 2015-10-10: no more necessary???
    ## if (depth(x$info$node[[pid]]) == 0L && length(coef$vc[[pid]]) > 0L) {
    ##   if (type == "print") {
    ##     coefMat <- coef$vc[[pid]]
    ##   } else {
    ##     coefMat <- cbind("Estimate" = coef$vc[[pid]],
    ##                      "Std. Error" = as.double(sd$vc[[pid]]),
    ##                      "z value" = as.double(coef$vc[[pid]] / sd$vc[[pid]]))
    ##   }
    ##   print(coefMat, ...)
    ## }
  }
  
  ## print footer
  if (nzchar(mess <- naprint(attr(x$data, "na.action")))) 
    cat(paste("\n(", mess, ")\n", sep = ""))
}


print.tvcm <- function(x, ...)
  tvcm_print(x, type = "print", ...)


prune.tvcm <- function(tree, cp = NULL, alpha = NULL, maxstep = NULL,
                       terminal = NULL, original = FALSE, ...) {

  mc <- match.call()

  ## checking arguments  
  tunepar <- c(!is.null(cp),
               !is.null(alpha),
               !is.null(maxstep),
               !is.null(terminal))

  if (sum(tunepar) < 1L) stop("no tuning parameter specified.")
  if (sum(tunepar) > 1L) stop("only one tuning parameter is allowed.")
  
  tunepar <- c("cp", "alpha", "maxstep", "terminal")[tunepar]
  
  if (tunepar == "cp") {
    stopifnot(is.numeric(cp) && length(cp) == 1L)

  } else if (tunepar == "alpha") {
    if (!tree$info$control$sctest) {
      warning("'alpha' is not a tuning parameter for 'tree'.")
      alpha <- NULL
    }
    stopifnot(is.numeric(alpha) && length(alpha) == 1L &&
              alpha >= 0.0 && alpha <= 1.0)

  } else if (tunepar == "maxstep") {   
    stopifnot(is.numeric(maxstep) && length(maxstep) == 1L &&
              maxstep >= 0L)
    
  } else if (tunepar == "terminal") {
    errMess <- paste("'terminal' must be a list with ",
                     length(tree$info$node), "elements.")
    if (!is.list(terminal)) {
      if (is.numeric(terminal) && length(tree$info$node) == 1L) {
        terminal <- list(terminal)
      } else {
        stop(errMess)
      }
    }
    if (is.list(terminal)) {
      if (length(terminal) != length(tree$info$node))
        stop(errMess)
      if (!all(sapply(terminal, function(x) is.null(x) | is.numeric(x))))
        stop("the elements of 'terminal' must be 'NULL' or 'numeric'")
      terminal <- lapply(terminal, as.integer)
    }
  }

  refit <- function(tree) {

    ## extract new formula
    env <- environment()
    formList <- tree$info$formula
    vcRoot <- sapply(tree$info$node, width) == 1L
    ff <- tvcm_formula(formList, vcRoot, tree$info$family, env)
    
    ## overwrite node predictors
    data <- model.frame(tree)
    data[, paste("Node", LETTERS[seq_along(tree$info$node)], sep = "")] <-
      tvcm_get_node(tree, tree$data, TRUE, tree$fitted[,"(weights)"], formList)
    
    ## refit model
    call <- list(name = as.name(tree$info$fit),
                 formula = quote(ff$full),
                 data = quote(data),
                 family = quote(tree$info$family),
                 weights = tree$fitted[,"(weights)"])
    call[names(tree$info$dotargs)] <- tree$info$dotargs
    mode(call) <- "call"
    tree$info$model <- suppressWarnings(try(eval(call), TRUE))
    
    ## check if the refitting failed
    if (inherits(tree$info$model, "try-error"))
      stop("tree model fitting failed.")
    
    ## heavy terms
    if (inherits(tree$info$model, "glm")) {
      attr(attr(tree$info$model$model, "terms"), ".Environment") <- NULL
      environment(tree$info$model$formula) <- NULL
      attr(tree$info$model$terms, ".Environment") <- NULL
    }
    
    ## update control
    tree$info$control <-
      tvcm_grow_setcontrol(tree$info$control, tree$info$model, formList, vcRoot)
    
    return(tree)
  }
  
  ## get the original tree
  if (original && !identical(tree$info$node, tree$info$grownode)) {
    tree$node <- tree$info$grownode[[1L]]
    tree$info$node <- tree$info$grownode
    tree <- refit(tree)         
  }
  
  if (tunepar %in% c("alpha", "maxstep", "terminal")) {
    
    ## prune the tree structure
    node <- tvcm_prune_node(tree, alpha, maxstep, terminal)
    
    ## refit the model if something has changed
    if (!identical(node, tree$info$node)) {
      
      ## attach nodes
      tree$node <- node[[1L]]
      tree$info$node <- node
      
      ## refit the model
      tree <- refit(tree)     
    }
    
  } else if (tunepar == "cp") {
    
    ## prune as long as there remain 'weak' links   
    run <- 1L; step <- 0L;
    cols <- c("part", "node", "loss", "npar", "nsplit", "dev")
    evalcols <- c("loss", "npar", "nsplit")
    prunepath <- list()
    
    while (run > 0) {
      
      run <- 0 ; step <- step + 1L;
      ncollapse <- NULL
      
      ids <- lapply(tree$info$node, function(node) {
        setdiff(nodeids(node), nodeids(node, terminal = TRUE))
      })
      
      ids00 <- lapply(seq_along(tree$info$node),
                      function(pid) {
                        unlist(nodeapply(tree$info$node[[pid]], ids[[pid]],
                                         function(node) node$info$id$original))
                      })
      
      ## 'call' to evaluate the collapses        
      prStatCall <- list(name = as.name(tree$info$control$papply),
                         X = quote(subs),
                         FUN = quote(prStat))
      prStatCall[names(tree$info$control$papply.args)] <-
        tree$info$control$papply.args
      mode(prStatCall) <- "call"
      
      ntab <- matrix(, length(unlist(ids)), length(cols))
      colnames(ntab) <- cols
      ntab[, "part"] <- rep(seq_along(tree$info$node), sapply(ids, length))
      ntab[, "node"] <- unlist(ids)
      ntab[, "loss"] <- rep.int(Inf, length(unlist(ids)))
      
      if (nrow(ntab) > 0L) {
        
        if (step > 1L) {
          
          ## search for the parents of the previously selecte node
          spart <- otab[ocollapse, "part"]
          snode <- otab[ocollapse, "node"]
          root <- 1L
          npath <- c(root); opath <- c(root);
          while (root != snode) {
            nkids <-
              unlist(nodeapply(tree$info$node[[spart]], root,
                               function(node)
                               sapply(node$kids, function(kids) kids$id)))
            okids <- unlist(nodeapply(tree$info$node[[spart]], nkids,
                                      function(node) node$info$id$last))
            kidkids <- nodeapply(tree$info$node[[spart]], nkids, nodeids)
            kid <- sapply(kidkids, function(x) snode %in% x)
            root <- nkids[kid]
            if (root != snode) {
              npath <- c(npath, nkids[kid])
              opath <- c(opath, okids[kid])
            }
          }
          
          ## insert results of parent models
          ntab[ntab[, "node"] %in% npath & ntab[, "part"] == spart, evalcols] <-
            otab[otab[, "node"] %in% opath & otab[, "part"] == spart, evalcols]
        }
      }
      
      loss0 <- tree$info$control$lossfun(tree)
      npar0 <- extractAIC(tree)[1L]
      nsplit0 <- sum(sapply(tree$info$node, width) - 1L)
      
      subs <- which(is.infinite(ntab[, "loss"]))
      if (length(subs) > 0L) {
        
        ## re-estimate all models which collapse each on inner node
        prStat <- function(i) {
          term <- lapply(seq_along(tree$info$node),
                         function(pid) {
                           if (pid == ntab[i, "part"]) return(ntab[i, "node"])
                           return(NULL)
                         })
          prTree <- try(prune(tree, terminal = term, papply = lapply), TRUE)
          if (!inherits(prTree, "try-error"))
            return(c(prTree$info$control$lossfun(prTree),
                     extractAIC(prTree$info$model)[1L],
                     sum(sapply(prTree$info$node, width) - 1L)))
          return(c(Inf, NA, NA))
        }
        stat <- eval(prStatCall)
        ntab[subs, evalcols] <- t(sapply(stat, function(x) x))
      }
      
      if (nrow(ntab) > 0L) {
        if (any(ntab[, "loss"] < Inf)) {
          
          ## minimum dfsplit such that the smaller model improves the fit
          ntab[, "dev"] <-
            (ntab[, "loss"] - loss0) /
              (tvcm_complexity(npar0, tree$info$control$dfpar,
                               nsplit0, tree$info$control$dfsplit) -
               tvcm_complexity(ntab[, "npar"], tree$info$control$dfpar,
                               ntab[, "nsplit"], tree$info$control$dfsplit))
          
          ## prune selected inner node
          if (any(ntab[, "dev"] <= cp)) {
            ncollapse <- which(!is.na(ntab[, "dev"]) &
                               !is.nan(ntab[, "dev"]) &
                               ntab[, "dev"] == min(ntab[, "dev"]))
            if (length(ncollapse) > 1L) ncollapse <- sample(ncollapse, 1L)
            if (length(ncollapse) > 0L) {
              term <- lapply(seq_along(tree$info$node),
                             function(pid) {
                               if (pid == ntab[ncollapse, "part"])
                                 return(ntab[ncollapse, "node"])
                               return(NULL)
                             })
              tree <- prune(tree, terminal = term)
              run <- 1
            }
          }
          otab <- ntab
          ocollapse <- ncollapse
          
        } else {
          stop("fitting of nested models failed.")
        }
      }
      
      ## create a 'data.frame' for the 'prunepath' method
      
      tab <- matrix(c(NA, NA, loss0, npar0, nsplit0, NA),
                    ncol = length(cols), dimnames = list("<none>", cols))
      if (nrow(ntab) > 0L) {
        ntab[, "node"] <- unlist(ids00)
        rownames(ntab) <- seq(1L, nrow(ntab), length.out = nrow(ntab))
        tab <- rbind(tab, ntab)
      }
      prunepath[[step]] <- list(step = step, tab = tab)
    }
    tree$info$prunepath <- prunepath
  }
  
  ## return pruned model  
  return(tree)
}


prunepath.tvcm <- function(tree, steps = 1L, ...) {
  steps <- intersect(steps, seq_along(tree$info$prunepath))
  rval <- tree$info$prunepath[steps]
  rval <- lapply(rval, function(x) {
    x$tab <- as.data.frame(x$tab)
    x$tab$part <- LETTERS[x$tab$part]
    return(x)
  })
  class(rval) <- "prunepath.tvcm"
  return(rval)
}


print.prunepath.tvcm <- function(x, na.print = "", ...) {
  for (i in seq_along(x)) {
    if (!is.null(x[[i]])) {
      if (i != 1L) cat("\n")
      cat("Step:", x[[i]]$step, "\n")
      print(as.matrix(x[[i]]$tab), na.print = na.print, quote = FALSE, ...)
    }
  }
}


summary.tvcm <- function(object, ...)
  tvcm_print(object, type = "summary", ...)


ranef.tvcm <- function(object, ...)
  return(ranef(object$info$model, ...))


resid.tvcm <- function(object, ...)
  return(resid(object = object$info$model, ...))


residuals.tvcm <- resid.tvcm


splitpath.tvcm <- function(tree, steps = 1L,
                           details = FALSE, ...) {

    steps <- intersect(steps, seq_along(tree$info$splitpath))
    rval <- tree$info$splitpath[steps]
    if (!details) {
        for (i in seq_along(steps)) {
            rval[[i]]$sctest <- NULL
            rval[[i]]$grid <- NULL
        }
    }
    class(rval) <- "splitpath.tvcm"
    return(rval)
}


print.splitpath.tvcm <- function(x, ...) {  
  for (i in seq_along(x)) {
    if (i != 1L) cat("\n")
    cat("Step:", x[[i]]$step)
    if (is.null(unlist(x[[i]]$varid))) {
      cat(" (no splitting processed)\n")
    } else {
      cat("\n\nSelected Split:")
      cat("\nPartition:", LETTERS[x[[i]]$partid])
      cat("\nNode:", x[[i]]$node)
      cat("\nVariable:", x[[i]]$var)
      cat("\nCutpoint: ")
      cat(paste("{",paste(x[[i]]$cutpoint, collapse = "}, {"), "}", sep = ""))
  }
    
    if (!is.null(x[[i]]$sctest)) {
      cat("\nCoefficient Constancy Tests (p-value):\n")
      for (pid in seq_along(x[[i]]$sctest)) {
        cat(paste("\nPartition ", LETTERS[pid], ":\n", sep = ""))
        print(as.data.frame(x[[i]]$sctest[[pid]]), ...)
      }
    } else cat("\n")
    
    if (!is.null(x[[i]]$grid)) {
      cat("\nLoss Reduction Statistics:")
      for (pid in seq_along(x[[i]]$grid)) {
        for (nid in seq_along(x[[i]]$grid[[pid]])) {
          for (vid in seq_along(x[[i]]$grid[[pid]][[nid]])) {
            if (is.null(x[[i]]$sctest) | !is.null(x[[i]]$sctest) &&
                any(x[[i]]$grid[[pid]][[nid]][[vid]][[2L]] > -Inf)) {
              cat("\nPartition:", LETTERS[pid],
                  "Node:", sub("Node", "",names(x[[i]]$grid[[pid]])[nid]),
                  "Variable:", names(x[[i]]$grid[[pid]][[nid]])[vid], "\n")
              print(as.data.frame(
                      cbind(x[[i]]$grid[[pid]][[nid]][[vid]][[1L]],
                            "dev" = x[[i]]$grid[[pid]][[nid]][[vid]][[2L]],
                            "npar" = x[[i]]$grid[[pid]][[nid]][[vid]][[3L]])), ...)
            }
          }
        }
      }
    }
  }
}


weights.tvcm <- function(object, ...) {
  weights(extract(object, "model"))
}


width.tvcm <- function(x, ...) {
  rval <- sapply(x$info$node, width)
  names(rval) <- LETTERS[seq_along(rval)]
  return(rval)
}
