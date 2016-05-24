##' -------------------------------------------------------- #
##' Author:          Reto Buergin
##' E-Mail:          rbuergin@gmx.ch
##' Date:            2016-02-08
##'
##' Description:
##' Workhorse functions for the 'tvcm' function.
##'
##' Overview:
##'
##' Workhorse functions for partitioning:
##' tvcm_complexity:          computes the complexity of the model
##' tvcm_grow:                main function for growing the trees
##' tvcm_grow_fit:            fits the current model
##' tvcm_grow_update:          refits the model (with utility function
##'                           'glm.doNotFit')
##' tvcm_getNumSplits         get cutpoints of numeric variables
##' tvcm_getOrdSplits         get cutpoints of ordinal variables
##' tvcm_getNomSplits         get cutpoints of nominal variables
##' tvcm_grow_setsplits:      get current splits
##' tvcm_setsplits_sctest:    update splits after tests
##' tvcm_setsplits_splitnode: 
##' tvcm_setsplits_rselect:   randomly select partitions, variables and nodes
##' tvcm_grow_sctest:         run coefficient constancy tests
##' tvcm_grow_exsearch:       compute the 'dev' statistics
##' tvcm_grow_splitnode:      split in variable x.
##' tvcm_formula:             extract separate formulas for
##'                           model and partitioning from
##'                           input formula.
##' tvcm_grow_setcontrol:     update the control argument before
##'                           fitting the tree.
##' tvcm_grow_setparm:        update the 'parm' slot
##'
##' Utility functions used by various functions:
##' tvcm_get_node:            extract node vectors and assign the contrasts
##' tvcm_get_terms:           creates a list which assigns coefficients
##'                           to the corresponding type, partition etc.
##' tvcm_get_vcparm:          extracts the names of the predictors on 'vc' terms
##' tvcm_get_estimates:       extracts the estimates from a fitted
##'                           'tvcm' object and creates a list with
##'                           an entry for each different type of
##'                           estimate ('fe', 'vc' or 're')
##' tvcm_print_vclabs:        creates short labels for 'vc' terms
##'
##' Functions for pruning:
##' tvcm_prune_node:          main function for pruning 'partynode'
##'                           objects
##' tvcm_prune_maxstep:       recursive function for pruning
##' tvcm_prune_terminal:      prunes branches
##' tvcm_grow_splitpath:      creates a 'splitpath.tvcm' object
##'
##' Last modifications:
##' 2016-01-08: change in 'tvcm_grow_fit' to allow fitting the approximate
##'             search modell locally. For now only for 'glm' fits!
##' 2015-11-31: enable the setting 'mtry <- Inf'
##' 2015-10-15: add function 'tvcm_grow_setparm'
##' 2015-08-25: replace 'fit' argument in 'tvcm_formula' by 'family'.
##' 2015-08-21: - small changes in 'tvcm_grow_fit'.
##'             - replace 'family' argument in 'tvcm_formula' by 'fit'.
##' 2015-02-25: add check for fixed effects model matrix in 'tvcm_grow_update'.
##' 2015-02-24: - improved 'tvcm_getNumSplits' (bugs for the upper limits)
##' 2014-12-10: - added 'drop = FALSE' commands in 'tvcm_exsearch_nomToOrd'
##'               which produced errors
##'             - 'tvcm_getNumSplits' yielded sometimes more than 'maxnumsplit'
##'               values. Now a random selection is applied for these cases
##' 2014-12-09: implemented accurate search model. Involves changes in
##'             'tvcm_formula', 'tvcm_grow_exsearch', 'tvcm_exsearch_dev'
##'             and 'tvcm_control'.
##' 2014-11-11: modified transformation of nominal into ordinal variables
##'             to accelerate exhaustive search. There is now a function
##'             'tvcm_exsearch_nomToOrd'.
##' 2014-11-05: parallelized 'estfun.olmm' call in 'tvcm_grow_sctest'
##' 2014-10-14: - modify dev-grid structure obtained from 'tvcm_grow_exsearch'
##'               each combination of part/node/var has now a list of three
##'               elements where the first contains the cuts, the second the
##'               the loss reduction and the third the difference in the number
##'               of parameters. Modifications concerned:
##'               - tvcm_grow_setsplits
##'               - tvcm_setsplits_sctest
##'               - tvcm_setsplits_rselect
##'               - print.splitpath
##'             - tvcm_setsplits_splitnode allocates for the splitted
##'               node a list structure (before it was a empty list
##'               on the node level)
##'             - small modifications in getSplits function
##'             - deleted 'tvcm_setsplits_validcats'
##'             - added 'tvcm_getNumSplits', 'tvcm_getOrdSplits' and
##'               'tvcm_getNomSplits'
##' 2014-09-22: deleted unnecessary 'subs' object in 'tvcm_grow_exsearch'
##'             which I didn't remove when removing the 'keepdev'
##'             option
##' 2014-09-17: - delete 'keepdev' argument (also for prune.tvcm)
##'             - add function 'tvcm_complexity'
##' 2014-09-15: changed 'dev' labels to 'dev' etc.
##' 2014-09-10: - add 'control' argument for 'tvcm_grow_update'
##'               to allow the control of variable centering
##'             - add variable centering in 'tvcm_grow_update'
##'               (which was not implemented for some curious reasons)
##' 2014-09-08: substitute 'rep' function by 'rep.int' or 'rep_len'
##' 2014-09-07: - added 'tvcm_get_vcparm' function
##'             - set default values in 'glm.doNotFit'
##' 2014-09-06: modified function names for 'tvcm_fit_model' and
##'             'tvcm_refit_model' for consistency reasons. The
##'             new names are 'tvcm_grow_fit' and 'tvcm_grow_update'
##' 2014-09-06: added new function 'tvcm_grow', which was formerly
##'             in 'tvcm'
##' 2014-09-04: added new function 'tvcm_print_vclabs'
##' 2014-09-02: modifications on 'tvcm_get_node' to accelerate
##'             the code
##' 2014-08-10: modifications to speed-up the code
##'             - update formulas of 'tvcm_formula' are now
##'               always identical
##'             - if 'doFit = FALSE', the call of 'glm.fit'
##'               is avoided
##' 2014-08-08: correct bug in 'tvcm_get_terms' for cases where
##'             multiple vc() terms with equal 'by' arguments
##'             are present
##' 2014-08-08: correct bug in 'tvcm_grow_setsplits' regarding
##'             'keepdev'
##' 2014-08-08: add suppressWarnings in tvcm_grow_fit
##' 2014-07-22: the list of splits is now of structure
##'             partitions-nodes-variables
##' 2014-07-22: AIC and BIC are no longer criteria and therefore
##'             multiple functions were adjusted
##' 2014-07-22: modified some function names
##' 2014-07-06: implement method to deal with many nominal categories
##' 2014-06-30: implement random selection if split is not unique
##' 2014-06-23: correct bug for 'start' argument in 'tvcm_grow_exsearch' 
##' 2014-06-17: modify documentation style
##' 2014-06-16: deleted several 'tvcm_prune_XXX' functions
##' 2014-06-03: modify 'tvcm_formula' to allow partition-wise
##'             trees
##' 2014-04-27: complete revision and improved documentation
##' 2014-04-01: rename 'fluctest' to 'sctest'
##' 2013-12-02: remove 'tvcm_grow_setupnode'
##' 2013-11-01: modify 'restricted' and 'terms' correctly in
##'             'tvcm_modify_modargs'
##'
##' Bottleneck functions:
##' - tvcm_grow_sctest
##' - tvcm_grow_setsplits
##' - tvcm_grow_exsearch
##'
##' To do:
##' - fitting local models when 'fast = TRUE' for 'olmm'
##'   objects
##' -------------------------------------------------------- #

##' -------------------------------------------------------- #
##' Compute the complexity of the model.
##'
##' @param npar    the number of coefficients
##' @param dfpar   the degree of freedom per parameter
##' @param nsplit  the number of splits
##' @param dfsplit the degree of freedom per split
##' 
##' @return a scalar giving the complexity of the model
##' -------------------------------------------------------- #

tvcm_complexity <- function(npar, dfpar, nsplit, dfsplit)
    return(dfsplit * nsplit + dfpar * npar)


##' -------------------------------------------------------- #
##' Growing a 'tvcm' tree.
##'
##' @param object  a 'tvcm' object
##' @param subset  a vector indicating the subset on which
##'    the model is to be fitted
##' @param weights a vector of weights corresponding to the
##'    subset entries
##' 
##' @return A 'tvcm' object.
##'
##' @details Used in 'cvloss.tvcm' 
##' -------------------------------------------------------- #

tvcm_grow <- function(object, subset = NULL, weights = NULL) {

  mcall <- object$info$mcall
  environment(mcall) <- environment()
  formList <- object$info$formula
  model <- object$info$model
  mf <- model.frame(model)
  partData <- object$data
  control <- object$info$control
  family <- model$family
  
  if (!is.null(subset)) {
    mf <- mf[subset,,drop = FALSE]
    partData <- partData[subset,, drop = FALSE]
  }

  if (!is.null(weights)) {
    mcall$weights <- weights
  } else {
    weights <- weights(model)
  }
      
  ## get number of partitions
  nPart <- length(formList$vc)

  ## get partitioning variables
  partVars <- lapply(formList$vc, function(x) attr(terms(x$cond), "term.labels"))
  varid <- lapply(partVars, function(x) {
    as.integer(sapply(x, function(x) which(colnames(partData) == x))) })
     
  ## set the root node
  nodes <-
    replicate(nPart, partynode(id = 1L, info = list(dims = nobs(model), depth = 0L)))
  names(nodes) <- names(formList$vc)
  where <- vector("list", length = nPart)
  
  partid <- seq(1, nPart, length.out = nPart)
  spart <- 0 # pseudo value

  ## allocate 'splits'
  splits <- lapply(seq_along(partid), function(pid) {
    lapply(1L, function(nid, pid) {
      lapply(seq_along(varid[[pid]]), function(x) {
        vector("list", 3)
      })
    }, pid = pid)
  })

  ## allocate 'splitpath'
  splitpath <- list()
  
  run <- 1L
  step <- 0L
  
  while (run > 0L) {

    step <- step + 1L; nstep <- step;
    test <- NULL; dev <- NULL;

    ## get current partitions and add them to the model data
    for (pid in seq_along(nodes)) {
      where[[pid]] <- factor(fitted_node(nodes[[pid]], partData))
      if (nlevels(where[[pid]]) > 1L)
      contrasts(where[[pid]]) <- contr.wsum(where[[pid]], weights)      
      mf[, paste0("Node", LETTERS[pid])] <- where[[pid]]
    }
    
    nodeid <- lapply(nodes, function(x) 1:width(x))
    
    if (control$verbose) cat("\n* starting step", step, "...")
    
    ## --------------------------------------------------- #
    ## Step 1: fit the current model
    ## --------------------------------------------------- #

    vcRoot <- sapply(nodeid, length) == 1L
    ff <- tvcm_formula(formList, vcRoot, family,
                       environment(formList$original))
    model <- try(tvcm_grow_fit(mcall))
    
    if (inherits(model, "try-error")) stop(model)

    control <- tvcm_grow_setcontrol(control, model, formList, vcRoot, TRUE)

    if (control$verbose) {
      cat("\n\nVarying-coefficient(s) of current model:\n")
      if (length(unlist(control$parm)) > 0L) {
        print(data.frame(Estimate = coef(model)[unique(unlist(control$parm))]),
              digits = 2)
      } else {
        cat("<no varying-coefficients>\n")
      }
    }
    
    ## compute / update splits
    splits <- tvcm_grow_setsplits(splits, partid, nodeid, varid, model,
                                  nodes, where, partData, control)
    
    ## check if there is at least one admissible split
    if (length(unlist(splits)) == 0L | step > control$maxstep |
        length(control$parm) == 0L) {
      run <- 0L
      if (step > control$maxstep) {
        stopinfo <- "maximal number of steps reached"
      } else if (length(control$parm) == 0L) {
        stopinfo <- "no varying coefficients"
      } else {
        stopinfo <- "no admissible splits (exceeded tree size parameters)"
      }
      nstep <- step - 1L
    }

    ## random selection (used by 'fvcm')
    if (control$mtry < .Machine$integer.max)
        splits <- tvcm_setsplits_rselect(splits, control)
    
    if (run > 0L && control$sctest) {
      
      ## --------------------------------------------------- #
      ## Step 2: variable selection via coefficient constancy tests
      ## --------------------------------------------------- #
      
      ## get raw p-values
      test <- try(tvcm_grow_sctest(model, nodes, where, partid, nodeid, varid, 
                                   splits, partData, control), TRUE)
      
      ## return error if test failed
      if (inherits(test, "try-error")) {
        run <- 0L
        stopinfo <- test
        
      } else {      
        testAdj <-
          tvcm_sctest_bonf(test,ifelse(control$bonferroni,"nodewise", "none"))
        run <- 1L * (min(c(1.0 + .Machine$double.eps, unlist(testAdj)),
                         na.rm = TRUE) <= control$alpha)
      }
      
      if (run > 0L) {
        
        ## extract the selected partition
        testAdjPart <-
          tvcm_sctest_bonf(test,ifelse(control$bonferroni,"partitionwise","none"))
        minpval <- min(unlist(testAdjPart), na.rm = TRUE)
        spart <-
          which(sapply(testAdjPart, function(x)any(sapply(x,identical,minpval))))
        if (length(spart) > 1L) spart <- sample(spart, 1L)
        
        ## select variable and node
        minsubs <- which(sapply(test[[spart]], identical,
                                min(test[[spart]], na.rm = TRUE)))
        if (length(minsubs) > 1L) minsubs <- sample(minsubs, 1L)
        svar <- ceiling(minsubs / nrow(test[[spart]]))
        snode <- minsubs - (svar - 1L) * nrow(test[[spart]])
        
        ## print results
        if (control$verbose) {
          
          ## tests
          cat("\nCoefficient constancy tests (p-value):\n")   
          for (pid in seq_along(nodes)) {
            cat(paste0("\nPartition ", LETTERS[pid], ":\n"))
            print(data.frame(format(testAdj[[pid]], digits = 2L)))              
          }
          
          ## selections
          cat("\nPartition:", LETTERS[spart])
          cat("\nNode:", levels(where[[spart]])[snode])
          cat("\nVariable:", names(partData)[varid[[spart]][svar]], "\n")
          
        }

        ## set dev statistic of not to selected nodes to 'Inf' to avoid
        ## model evaluations
        splits <- tvcm_setsplits_sctest(splits, partid, spart,
                                        nodeid, snode, varid, svar)
        
      } else {
        stopinfo <- "p-values of coefficient constancy tests exceed alpha"
      }  
    }
      
    if (run > 0L) {
      
        ## ------------------------------------------------- #
        ## Step 3: search a cutpoint
        ## ------------------------------------------------- #
      
        ## compute the dev of all candidate splits and extract the best split
        dev <- try(tvcm_grow_exsearch(splits, partid, nodeid, varid, 
                                      model, nodes, where, partData,
                                      control, mcall, formList, step),
                   silent = TRUE)
        
        ## handling stops
        if (inherits(dev, "try-error")) {
            run <- 0L
            stopinfo <- dev
            nstep <- step - 1L
            
        } else {
            splits <- dev$grid
            spart <- dev$partid
        
            if (is.null(dev$cut)) {
                run <- 0L
                stopinfo <- "found no admissible split"
                nstep <- step - 1L
            }
        
            if (run > 0L) {
                if (dev$pendev < control$mindev) {
                    run <- 0
                    stopinfo <- paste("no split with",
                                      if (control$cp > 0) "penalized",
                                      "loss reduction > mindev")
                    nstep <- nstep - 1L
                }
            }
        }   
    }
    
    ## incorporate the split into 'nodes'
    if (run > 0L)
        nodes <- tvcm_grow_splitnode(nodes, where, dev, partData,
                                     step, weights)

    if (run > 0L)
      splits <- tvcm_setsplits_splitnode(splits, dev$partid, dev$nodeid, nodeid)
      
    ## update 'splitpath' to make the splitting process traceable
    if (run >= 0L)
      splitpath[[step]] <-
        list(step = step,
             dev = control$lossfun(model),
             npar = extractAIC(model)[1L],
             nsplit = step - 1L)

    if (!inherits(test, "try-error"))
      splitpath[[step]]$sctest <- test
    
    if (!inherits(dev, "try-error") && run > 0L) 
      splitpath[[step]] <- append(splitpath[[step]], dev)
    
    ## print the split
    if (control$verbose) {
      if (run > 0L) {
        if (!control$sctest) {
            cat("\n\nPartition:", LETTERS[dev$partid])
            cat("\nNode:", levels(where[[dev$partid]])[dev$nodeid])
            cat("\nVariable:", names(partData)[dev$varid])
        } else {
            cat("\n")
        }

        cat("\nCutpoint:\n")
        print(as.data.frame(matrix(dev$cut, 1L,
                                   dimnames = list(dev$cutid,
                                     names(dev$cut)))))
        
        cat("Model comparison:\n")
        print(data.frame(
                cbind("loss" = c(
                        round(control$lossfun(model), 2),
                        round(control$lossfun(model) - dev$dev, 2L)),
                      ## if 'cp == 0'
                      "dev" = if (control$cp == 0)
                      c("", round(dev$dev, 2L)),
                      ## if 'cp > 0'
                      "penalized dev" = if (control$cp > 0)
                      c("", round(dev$pendev, 2L)),
                      deparse.level = 2),
                      row.names = paste("step", step + c(-1, 0)),
                      check.names = FALSE))
        
      } else {
        cat("\n\nStopping partitioning.\nMessage:", as.character(stopinfo), "\n")
        if (inherits("try-error", stopinfo)) warning(as.character(stopinfo))
        
      }
    }
  }
  
  ## inscribe original node names for later pruning
  for (pid in seq_along(nodes)) {
    nodes[[pid]] <- as.list(nodes[[pid]])
    for (nid in 1:length(nodes[[pid]])) {
      nodes[[pid]][[nid]]$info$id$original <- nodes[[pid]][[nid]]$id
      nodes[[pid]][[nid]]$info$id$last <- nodes[[pid]][[nid]]$id
    }
    nodes[[pid]] <- as.partynode(nodes[[pid]])
  }
  
  ## prepare the title
  title <- c("Tree-Based Varying Coefficients Model")
  
  ## modify splitpath    
  splitpath <- tvcm_grow_splitpath(splitpath, varid, nodes, partData, control)
  
  ## the output object
  if (nPart == 0L) {
    tree <- model
    
  } else {

    ## delete environments
    environment(mcall) <- NULL
    environment(object$info$call) <- NULL
    formList <- vcrpart_formula_delEnv(formList)
    attr(attr(partData, "terms"), ".Environment") <- NULL
    
    ## build 'tvcm' object
    tree <- party(nodes[[1L]],
                  data = partData,
                  fitted = data.frame(
                    "(response)" = model.response(model.frame(model)),
                    "(weights)" = weights(model),
                    check.names = FALSE),
                  terms = terms(formList$original, keep.order = TRUE),
                  info = list(
                    title = title,
                    call = object$info$call,
                    mcall = mcall,
                    formula = formList,
                    direct = object$info$direct,
                    fit = object$info$fit,
                    family = family,
                    control = control,
                    info = stopinfo,
                    model = model,
                    node = nodes,
                    grownode = nodes,
                    nstep = nstep,
                    splitpath = splitpath,
                    dotargs = object$info$dotargs))
    class(tree) <- c("tvcm", "party")
  }

  return(tree)
}


##' -------------------------------------------------------- #
##' Avoid calling glm.fit if 'fit == FALSE'
##'
##' @param call    an object of class call
##' @param doFit   a logical indicating whether the parameters
##'    have to be optimized.
##'
##' @return A list of formulas ('root', 'tree' and 'original').
##'
##' @details Used in 'tvcm' and 'tvcm_grow_sctest'. 'glm.doNotFit'
##' is just an utility function to skip \code{\link{glm.fit}}
##' if \code{doFit = FALSE}
##' -------------------------------------------------------- #

glm.doNotFit <- function(x, y, weights = NULL, start = NULL, etastart = NULL,
                         mustart = NULL, offset = NULL, family = gaussian(),
                         control = list(), intercept = TRUE) {
  coefficients <- rep.int(0, NCOL(x))
  names(coefficients) <- colnames(x)
  if (is.null(weights)) weights <- rep.int(1.0, NROW(x))
  if (is.null(offset)) offset <- rep.int(0.0, NROW(x))
  if (!is.null(start)) {
    if (!is.null(names(start))) {
      start <- start[intersect(names(coefficients), names(start))]
      coefficients[names(start)] <- start
    } else {
      if (length(start) > length(coefficients))
        start <- start[seq_along(coefficients)]
      coefficients[seq_along(start)] <- start
    }
  } 
  return(list(coefficients = coefficients, residuals = NULL,
              effect = NULL, R = NULL, rank = NULL,
              qr = NULL, family = family,
              linear.predictor = etastart,
              deviance = NULL, aic = NULL, null.deviance = NULL,
              iter = 0, weights = NULL, prior.weights = weights,
              df.residual = NULL, df.null = NULL, y = y,
              converged = TRUE, boundary = TRUE))
}

tvcm_grow_fit <- function(mcall, doFit = TRUE) {
    
  ## extract information from 'mcall'
  env <- environment(mcall)
  
  ## set mcall if coefficients are not to optimized
  if (!doFit) {
      fit <- deparse(mcall$name)
      if (fit == "olmm") {
          mcall$doFit <- FALSE
      } else if (fit == "glm") {        
          mcall$method <- glm.doNotFit # skips glm.fit
      } 
  }
  
  ## fit model
  object <- suppressWarnings(eval(mcall, env))

  ## return error if fitting failed
  if (inherits(object, "try-error")) stop("model fitting failed.")
  if (doFit && !is.null(object$converged) && !object$converged)
      stop("no convergence")

  ## delete heavy objects
  if (inherits(object, "glm")) {
    attr(attr(object$model, "terms"), ".Environment") <- NULL
    environment(object$formula) <- NULL
    attr(object$terms, ".Environment") <- NULL
  }
  
  ## return model
  return(object)
}


##' -------------------------------------------------------- #
##' Updates the model matrix and re-fits the current node
##' model. Used for the grid-search in 'tvcm_grow_exsearch'.
##'
##' @param object a prototype model
##' @param mcall   the mcall for the prototype model
##'
##' @return A list of formulas ('root', 'tree' and 'original').
##'
##' @details Used in 'tvcm' and 'tvcm_grow_exsearch'. Note that the
##' function will modify the slots of the original object as well!
##'
##' To do:
##' Improve performance for non 'olmm' objects
##' -------------------------------------------------------- #

tvcm_grow_update <- function(object, control, subs = NULL) {
  
  if (inherits(object, "olmm")) {
  
    ## set new partition
    data <- model.frame(object)
    nodeVars <- grep("Node[A-Z]", colnames(data), value = TRUE)
    object$frame[, nodeVars] <- data[, nodeVars]

    ## get terms
    termsFeCe <- terms(object, "fe-ce")
    termsFeGe <- terms(object, "fe-ge")

    ## get constrasts
    contrasts <- object$contrasts
    conCe <- contrasts[intersect(names(contrasts), all.vars(termsFeCe))]
    if (length(conCe) == 0L) conCe <- NULL
    conGe <- contrasts[intersect(names(contrasts), all.vars(termsFeGe))]    
    if (length(conGe) == 0L) conGe <- NULL
    
    ## update model matrix  
    object$X <-
      olmm_merge_mm(x = model.matrix(termsFeCe, object$frame, conCe),
                    y = model.matrix(termsFeGe, object$frame, conGe), TRUE)
    object$X <- olmm_check_mm(object$X)
    
    if (control$center) {
      
      ## extract interaction predictors to be centered
      ## (the ones with 'Left' and 'Right')
      cColsCe <- which(rownames(attr(termsFeCe, "factors")) %in% c("Left", "Right"))
      if (length(cColsCe) > 0L) {
        cColsCe <-
          which(colSums(attr(termsFeCe, "factors")[cColsCe,,drop = FALSE]) > 0 &
                !colnames(attr(termsFeCe, "factors")) %in% c("Left", "Right"))
        cColsCe <-
          which(attr(object$X, "assign") %in% cColsCe & attr(object$X, "merge") == 1)
      }
      cColsGe <- which(rownames(attr(termsFeGe, "factors")) %in% c("Left", "Right"))
      if (length(cColsGe) > 0L) {
        cColsGe <-
          which(colSums(attr(termsFeGe, "factors")[cColsGe,,drop = FALSE]) > 0 &
                !colnames(attr(termsFeGe, "factors")) %in% c("Left", "Right"))
        cColsGe <-
          which(attr(object$X, "assign") %in% cColsGe & attr(object$X, "merge") == 2)
      }
      
      ## center the predictors
      object$X[,  c(cColsCe, cColsGe)] <-
        scale(object$X[,  c(cColsCe, cColsGe)], center = TRUE, scale = FALSE)
    }
    
    ## prepare optimization
    optim <- object$optim
    optim[[1L]] <- object$coefficients
    optim[[4L]] <- object$restricted
    environment(optim[[2L]]) <- environment()
    if (!object$dims["numGrad"]) environment(optim[[3L]]) <- environment()
    FUN <- optim$fit
    subs <- which(names(optim) == "fit")
    optim <- optim[-subs]

    ## run optimization
    output <- try(suppressWarnings(do.call(FUN, optim)), TRUE)
    
    ## check optimized model 
    if (!inherits(output, "try-error")) {
      object$output <- output
      object$conv <- switch(object$optim$fit,
                            optim = object$output$convergence == 0,
                            nlminb = object$output$convergence == 0,
                            ucminf = object$output$convergence %in% c(1, 2, 4))
      if (!object$conv) object <- try(stop("not converged"), TRUE)
    } else {
      object <- output
    }
    
  } else {

    ## extract interaction predictors to be centered
    ## (the ones with 'Left' and 'Right')
    X <- model.matrix(object$formula, model.frame(object))

    if (control$center) {

        ## get the columns of 'X' to be centered
        terms <- terms(object$formula)       
        cCols <- which(rownames(attr(terms, "factors")) %in% c("Left", "Right"))
        if (length(cCols) > 0L) {
            cCols <- which(colSums(attr(terms, "factors")[cCols,,drop = FALSE]) > 0 &
                          !colnames(attr(terms, "factors")) %in% c("Left", "Right"))
            cCols <- which(attr(X, "assign") %in% cCols)
        }

        ## centering
        X[, cCols] <- scale(X[, cCols], center = TRUE, scale = TRUE)
    }

    ## if 'fast = TRUE' model is fitted locally (thereby nuisance parameter is 
    ## left as a free parameter)
    if (!control$fast | is.null(subs)) subs <- rep(TRUE, nrow(X))
    start <- if (control$fast) object$coefficients else NULL
    
    ##subs <- rep(TRUE, nrow(X))
    
    object <- try(suppressWarnings(
                    glm.fit(x = X[subs,,drop=FALSE], 
                            y = object$y[subs], 
                            weights = object$prior.weights[subs],
                            start = start, 
                            offset = object$offset[subs],
                            family = object$family, 
                            control = object$control,
                            intercept = TRUE)), TRUE)
    
    if (!inherits(object, "try-error")) {
      class(object) <- c("glm", "lm")
      if (!object$conv) object <- try(stop("not converged"), TRUE)
    }
  }
  
  ## return model
  return(object)
}


tvcm_grow_gefp <- gefp.olmm # see 'olmm-methods.R'


##'-------------------------------------------------------- #
##' Computes cutpoints of 'numeric' partitioning variables
##'
##' @param z a numeric vector
##' @param w a numeric vector of weights
##' @param minsize numerical scalar. The minimum node size.
##' @param maxnumsplit integer scalar. The maximum number
##'    of cutpoints.
##' 
##' @return A matrix with one column and one row for each
##'    cutpoint
##'
##' @details Used in 'tvcm_grow_setsplits'.
##'-------------------------------------------------------- #

tvcm_getNumSplits <- function(z, w, minsize, maxnumsplit) {
  
  ## order the partitioning variable
  ord <- order(z)
  z <- z[ord]; w <- w[ord];
  cw <- cumsum(w)

  ## result if there is no split
  rval0 <- matrix(numeric(), ncol = 1L)
  colnames(rval0) <- "cut"
  attr(rval0, "type") <- "dev"
  
  ## get the first index
  subsL <- which(cw >= minsize)
  if (length(subsL) < 1L) return(rval0)
  subsL <- min(subsL)

  ## get the last index
  subsR <- cw < (cw[length(cw)] - minsize + 1)
  if (!any(subsR)) return(rval0)
  if (any(!subsR)) subsR <- subsR & z < min(z[!subsR])
  if (!any(subsR)) return(rval0)
  subsR <- max(which(subsR))

  if (subsL <= subsR && z[subsL] <= z[subsR]) {

    ## valid splits available
    z <- z[subsL:subsR]
    cw <- cw[subsL:subsR] - cw[subsL] + w[subsL]
    
    ## apply a cutpoint reduction if necessary
    if (length(unique(z)) > maxnumsplit) {
      nq <- maxnumsplit - 1
      rval <- c()
      cw <- cw / cw[length(cw)]
      while (length(unique(rval)) < maxnumsplit) {
        nq <- nq + 1L
        q <- (1:nq) / (nq + 1L)
        rval <- unique(sapply(q, function(p) z[max(which(cw <= p))]))
      }
    } else {
      rval <- unique(z)
    }

    ## sometimes the while loop yields too many values ...
    if (length(rval) > maxnumsplit)
      rval <- sort(sample(rval, 9))
    
  } else {

    ## no valid splits
    rval <- numeric()
  }

  ## prepare return value
  rval <- matrix(rval, ncol = 1L)
  colnames(rval) <- "cut"
  attr(rval, "type") <- "dev"

  ## return matrix with cutpoints
  return(rval)
}


##'-------------------------------------------------------- #
##' Computes cutpoints of 'ordinal' partitioning variables
##'
##' @param z a vector of class 'ordered'
##' @param w a numeric vector of weights
##' @param minsize numerical scalar. The minimum node size.
##' @param maxordsplit integer scalar. The maximum number
##'    of cutpoints.
##' 
##' @return A matrix with one column for each category and
##'    one row for each cutpoint
##'
##' @details Used in 'tvcm_grow_setsplits'.
##'-------------------------------------------------------- #

tvcm_getOrdSplits <- function(z, w, minsize, maxordsplit) {

  ## get integer cutpoints using 'tvcm_getNumSplits'
  nl <- nlevels(z) # all levels
  cuts <- tvcm_getNumSplits(as.integer(z), w, minsize, maxordsplit)
  zdlev <- 1:nl %in% as.integer(cuts)

  ## create a matrix to apply categorical splits
  rval <- diag(nl)
  rval[lower.tri(rval)] <- 1L
  rval <- rval[zdlev,, drop = FALSE]

  ## prepare return value
  colnames(rval) <- levels(z)
  attr(rval, "type") <- "dev"

  ## return matrix with cutpoints
  return(rval)
}


##'-------------------------------------------------------- #
##' Computes cutpoints of 'nominal' partitioning variables
##'
##' @param z a vector of class 'factor'
##' @param w a numeric vector of weights
##' @param minsize numerical scalar. The minimum node size.
##' @param maxnomsplit integer scalar. The maximum number
##'    of cutpoints.
##' 
##' @return A matrix with one column for each category and
##'    one row for each cutpoint
##'
##' @details Used in 'tvcm_grow_setsplits'.
##'-------------------------------------------------------- #

tvcm_getNomSplits <- function(z, w, minsize, maxnomsplit) {

  zdlev <- 1 * (levels(z) %in% levels(droplevels(z)))
  if (sum(zdlev) < maxnomsplit) {
  
    ## exhaustive search
    rval <- .Call("tvcm_nomsplits",
                  as.integer(zdlev),
                  PACKAGE = "vcrpart")
    type <- "dev"
    
  } else {

    ## Heuristic reduction of splits: in tvcm_grow_exsearch,
    ## the 'isolated' coefficients of each category are
    ## computed. The coefficients are used for ordering
    ## the categories and finally the variable is treated
    ## as ordinal. See tvcm_grow_exsearch

    rval <- diag(nlevels(z))
    type <- "coef"
  }

  ## remove splits according to 'minsize'
  sumWTot <- sum(w)
  sumWCat <- tapply(w, z, sum)
  sumWCat[is.na(sumWCat)] <- 0
  valid <- apply(rval, 1, function(x) {
    all(c(sum(sumWCat[x > 0]), sumWTot - sum(sumWCat[x > 0])) > minsize)
  })
  rval <- rval[valid,, drop = FALSE]
  
  ## return matrix with cutpoints
  colnames(rval) <- levels(z)
  attr(rval, "type") <- type
  return(rval)
}


##'-------------------------------------------------------- #
##' Computes candidate splits for the current step
##'
##' @param splits   a list. The former list of splits
##' @param partid   a vector of candidate partitions for splitting.
##' @param spart    integer scalar. The partition in which the
##'                 last split was employed
##' @param varid    a \code{list} with a vector for each partition that
##'    that specifies candidate variables for splitting.
##' @param nodeid   a \code{list} with a vector for each partition that
##'    that specifies candidate nodes for splitting.
##' @param model    a fitted model of class \code{\link{olmm}} or
##'    \code{\link{glm}}.
##' @param nodes    a \code{list} with a \code{\link{partynode}}
##'    object for each partition.
##' @param where    a \code{list} with a factor vector for each
##'    partition that the assigns observations to nodes.
##' @param partData a data frame with variables for
##'    splitting.
##' @param control    a \code{list} of control parameters as produced
##'    by 'tvcm_control.'
##' 
##' @return A list of splits. Entries for splits that
##'    exceed the tuning parameters are a vector of length
##'    zero.
##'
##' @details Used in 'tvcm'.
##'-------------------------------------------------------- #

tvcm_grow_setsplits <- function(splits, partid, nodeid, varid,
                                model, nodes, where,
                                partData, control) {
  
  ## get tree size criteria of current tree(s)
  width <- sapply(nodes, width)  
  depth <- lapply(nodes, function(node) {
    unlist(nodeapply(node, nodeids(node, terminal = TRUE), function(node) {
      info_node(node)$depth }))
  })
  w <- weights(model)
  
  getSplits <- function(pid, nid, vid) {

    ## prepare data
    subs <- where[[pid]] == levels(where[[pid]])[nid]
    z <- partData[subs, vid]
    w <- w[subs]
    
    if (width[pid] >= control$maxwidth[pid] |
        depth[[pid]][nid] >= control$maxdepth[pid] |
        sum(subs) < 1L |
        sum(w) < 2 * control$minsize[pid]) {

      ## return an empty matrix if control parameters exceeded
      rval <- matrix(, 0, ifelse(is.numeric(z), 1L, nlevels(z)))
      colnames(rval) <- if (is.numeric(z)) "cut" else levels(z)
      attr(rval, "type") <- "dev"

    } else {

      ## compute matrix according to their scale
      rval <- switch(class(z)[1L],
                     numeric = tvcm_getNumSplits(z, w,
                       control$minsize[pid],
                       control$maxnumsplit),
                     integer = tvcm_getNumSplits(z, w,
                       control$minsize[pid],
                       control$maxnumsplit),
                     ordered = tvcm_getOrdSplits(z, w,
                       control$minsize[pid],
                       control$maxordsplit),
                     factor = tvcm_getNomSplits(z, w,
                       control$minsize[pid],
                       control$maxnomsplit),
                     stop("class of variable '",
                          colnames(partData)[vid],
                          "' not recognized"))
    }
    return(rval)
  }

  ## compute splits in all partitions, partitioning variables and nodes  
  for (pid in seq_along(partid)) {
    for (nid in seq_along(nodeid[[partid[pid]]])) {
      for (vid in seq_along(varid[[partid[pid]]])) {
        split <- splits[[pid]][[nid]][[vid]]
        if (is.null(split[[1L]]) |
            (!is.null(split[[1L]]) && attr(split[[1L]], "type") == "coef") |
            width[pid] >= control$maxwidth[pid]) {
          
          split[[1L]] <- getSplits(partid[pid],
                                   nodeid[[partid[pid]]][nid],
                                   varid[[partid[pid]]][vid])
          split[[2L]] <- split[[3L]] <- rep(NA, nrow(split[[1L]]))
        } else {
          if (nrow(split[[1]]) > 0L)
            split[[2L]][] <- split[[3L]][] <- NA
        }
        splits[[pid]][[nid]][[vid]] <- split
      }
    }
  }

  ## return list of updated 'splits'
  return(splits)
}


##'------------------------------------------------------ #
##' Updates the list of splits after coefficient
##' constancy tests.
##'
##' @param splits  a 'list' as produced by 'tvcm_grow_splits'
##' @param partid     a vector of candidate partitions for splitting.
##' @param spart      the selected partition
##' @param nodeid     a list with a vector for each partition that
##'    that specifies candidate nodes for splitting.
##' @param snode      the selected node
##' @param varid      a list with a vector for each partition that
##'    that specifies candidate variables for splitting.
##' @param svar       the selected variable

##' 
##' @return An modified list of splits.
##'
##' @details Used in 'tvcm'.
##'------------------------------------------------------ #

tvcm_setsplits_sctest <- function(splits, partid, spart,
                                  nodeid, snode, varid, svar) {
  
  ## set loss of not selected parts to -Inf
  for (pid in seq_along(partid))
    for (nid in seq_along(nodeid[[pid]]))
      for (vid in seq_along(varid[[pid]])) {
        if (pid == spart & nid == snode & vid == svar) {
          splits[[pid]][[nid]][[vid]][[2L]][] <- NA
        } else {
          splits[[pid]][[nid]][[vid]][[2L]][] <- -Inf
        }
      }
  ## return updated 'splits'
  return(splits)
}


##'------------------------------------------------------ #
##' Updates the list of splits after grid search
##'
##' @param splits  a 'list' as produced by 'tvcm_grow_splits'
##' @param partid     a vector of candidate partitions for splitting.
##' @param spart      the selected partition
##' @param nodeid     a list with a vector for each partition that
##'    that specifies candidate nodes for splitting.
##' @param snode      the selected node
##' @param varid      a list with a vector for each partition that
##'    that specifies candidate variables for splitting.
##' @param svar       the selected variable
##' 
##' @return An modified list of splits.
##'
##' @details Used in 'tvcm'.
##'------------------------------------------------------ #

tvcm_setsplits_splitnode <- function(splits, spart, snode, nodeid) {
  
  ## expand the splits list
  lnodes <- nodeid[[spart]][nodeid[[spart]] < snode]
  unodes <- nodeid[[spart]][nodeid[[spart]] > snode]
  split0 <- splits[[spart]]
  split <- vector("list", length(split0) + 1L) 
  if (length(lnodes) > 0L) split[lnodes] <- split0[lnodes]
  if (length(unodes) > 0L) split[unodes + 1L] <- split0[unodes]
  nvar <- length(split0[[snode]])
  split[[snode]] <- split[[snode + 1]] <-
    replicate(nvar, vector("list", 3L), simplify = FALSE)
  splits[[spart]] <- split
  
  ## return updated 'splits'
  return(splits)
}


##'------------------------------------------------------ #
##' Randomly select partitions, variables and nodes
##'
##' @param splits  a 'list' as produced by 'tvcm_grow_splits'
##' @param partid     a vector of candidate partitions for splitting.
##' @param spart      the selected partition
##' @param varid      a list with a vector for each partition that
##'    that specifies candidate variables for splitting.
##' @param svar       the selected variable
##' @param nodeid     a list with a vector for each partition that
##'    that specifies candidate nodes for splitting.
##' @param snode      the selected node
##' 
##' @return An modified list of splits.
##'
##' @details Used in 'tvcm'.
##'------------------------------------------------------ #

tvcm_setsplits_rselect <- function(splits, control) {

  ## extract valid partition/node/variable combinations
  idmat <- matrix(, 0, 3)
  for (pid in seq_along(splits))
    for (nid in seq_along(splits[[pid]]))
      for (vid in seq_along(splits[[pid]][[nid]]))
        if (nrow(splits[[pid]][[nid]][[vid]][[1L]]) > 0L)
          idmat <- rbind(idmat, c(pid, nid, vid))
  
  ## randomly select partition/node/variable combinations
  if (nrow(idmat) > control$mtry)
    idmat <- idmat[sort(sample(1:nrow(idmat), control$mtry)),,drop = FALSE]
  
  ## delete unselected nodes
  for (pid in seq_along(splits))
    for (nid in seq_along(splits[[pid]]))
      for (vid in seq_along(splits[[pid]][[nid]]))
        if (!any(apply(idmat, 1, function(x) all(x == c(pid, nid, vid)))))        
          splits[[pid]][[nid]][[vid]][[2L]][] <- -Inf
  
  ## return updated 'splits'
  return(splits)
}


##'------------------------------------------------------ #
##' Processing of nodewise coefficient constancy tests.
##'
##' @param model    the current model.
##' @param nodes    an object of class 'partynode'.
##' @param where    a list of vectors that locate the observations
##'    to their corresponding node(s)
##' @param partid   an integer vector that indicates which
##'    partitions should be tested.
##' @param varid    a list with a vector for each partition that
##'    that specifies the variables to be tested
##' @param partData a 'data.frame' with the partitioning variables.
##' @param control  an object of class 'tvcm_control'.
##' 
##' @return A list with partitions 'statistic' and 'p.value'.
##'
##' @details Used in 'tvcm'.
##'------------------------------------------------------ #

tvcm_grow_sctest <- function(model, nodes, where, partid, nodeid, varid, 
                             splits, partData, control) {
    
  ## get variable types
  functional <- sapply(partData, function(x) {
    switch(class(x)[1],
           factor = control$functional.factor,
           ordered = control$functional.ordered,
           integer = control$functional.numeric,
           numeric = control$functional.numeric)
  })
  
  ## prepare list with arguments for 'sctest'
  rval <- vector("list", length(nodes))
  for (pid in seq_along(nodes)) {
    dim <- c(nlevels(where[[pid]]), length(varid[[pid]]), control$nimpute)
    dn <- list(paste0("Node", LETTERS[pid], levels(where[[pid]])),
               colnames(partData)[varid[[pid]]], 1:control$nimpute)
    rval[[pid]] <- array(, dim = dim, dimnames = dn)
  }

  ## call 'estfun' (eventually parallelized)
  eCall <-
    list(name = as.name(ifelse(inherits(model, "olmm"),"estfun.olmm", "estfun")))
  eCall$x <- quote(model)
  eCall[names(control$estfun.args)] <- control$estfun.args
  mode(eCall) <- "call"
  pCall <- list(name = as.name(control$papply),
               X = quote(seq(1L, control$nimpute, 1L)),
               FUN = quote(function(i) eval(eCall)))  
  pCall[names(control$papply.args)] <- control$papply.args
  mode(pCall) <- "call"
  scores <- simplify2array(eval(pCall))
  
  ## set the 'gefp' call (which is called in each iteration below)
  gCall <- call(name = "tvcm_grow_gefp", object = quote(model),
                scores = quote(sc),
                order.by = quote(z), subs = quote(rows),
                parm = quote(cols), center = TRUE, silent = TRUE)

  terms <- # useful information to identify the coefficients to test
    tvcm_get_terms(
        names = dimnames(scores)[[2L]],
        ids = lapply(nodes, function(node) nodeids(node, terminal = TRUE)),
        parm = control$parm,
        fixed = NULL)
  ## note: fixed is set to 'NULL' for the case where a variable that is at
  ## the same time a fixed i.e. main effect and a 'by' variable in a 'vc'
  ## term without split is not dropped (which is meaningful for the 'coef'
  ## method)
  
  ## apply test for each variable and partition separately
  for (pid in seq_along(partid)) { # loop over partitions
      
    for (vid in seq_along(varid[[pid]])) { # loop over variables
      
      ## get variable to test
      z <- partData[, varid[[pid]][vid]]       
      
      for (nid in seq_along(nodeid[[pid]])) { # loop over nodes
        
        ## check if there is a permitted split
        if (nrow(splits[[pid]][[nid]][[vid]][[1L]]) > 0L &&
            all(is.na(splits[[pid]][[nid]][[vid]][[2L]]))) {
          
          ## observations of the current partition
          rows <- where[[pid]] == levels(where[[pid]])[nid] 
          
          ## columns corresponding to the tested subset
          if (nlevels(where[[pid]]) == 1L) {
            cols <- unlist(control$parm[[pid]])
          } else {
            cols <- terms$partition == LETTERS[pid] & terms$node == 
              levels(where[[pid]])[nid]
            cols <- dimnames(scores)[[2L]][cols]
          }
          
          for (k in 1:control$nimpute) {
            sc <- matrix(scores[,,k,drop=FALSE], dim(scores)[1], dim(scores)[2],
                         dimnames = dimnames(scores)[1L:2L])            
            gefp <- try(eval(gCall), TRUE)
            
            ## extract test statistic
            if (!inherits(gefp, "try-error")) {            
              if (is.character(functional)) {
                functional <- tolower(functional)
                fi <- switch(functional[varid[[pid]][vid]],
                             "suplm" = supLM(from = control$trim),
                             "lmuo" = catL2BB(gefp),
                             stop("Unknown efp functional."))
              } else {
                fi <- functional[varid[[pid]][vid]]
              }              
              test <- try(sctest(x = gefp, functional = fi), TRUE)
              if (!inherits(test, "try-error")) {
                rval[[pid]][nid, vid, k] <- test$p.value
              }
            }
          }
        }
      }
    }
  }
  rval <- lapply(rval, function(x) apply(x, c(1L, 2L), mean, 
                                         na.rm = TRUE))
  rval <- lapply(rval, function(x) {
      x[is.nan(x)] <- NA
      return(x)
  })
  names(rval) <- LETTERS[partid]
  return(rval)
}


tvcm_sctest_bonf <- function(test, type) {
  for (pid in seq_along(test)) {
    for (nid in 1:nrow(test[[pid]])) {
      k <- switch(type,
                  "none" = 1.0,
                  "all" = sum(!is.na(unlist(test))),
                  "nodewise" = sum(!is.na(test[[pid]][nid, ])),
                  "partitionwise" = sum(!is.na(test[[pid]])))
      pval1 <- pmin(1.0, k * test[[pid]][nid, ])
      pval2 <- c(1.0 - (1.0 - test[[pid]][nid, ]) ^ k)
      test[[pid]][nid, ] <-
        ifelse(!is.na(test[[pid]][nid, ]) & (test[[pid]][nid, ] > 0.01),
               pval2, pval1)
    }
  }
  return(test)
}


##'-------------------------------------------------------- #
##' Computes the loss for each possible split.
##'
##' @param varid    an integer vector indicating the partitioning
##'    variables to be evaluated.
##' @param nodeid   an integer vector indicating the nodes to be
##'    evaluated.
##' @param model    the current model
##' @param nodes    an object of class 'partynode'.
##' @param partData a 'data.frame' with the partitioning variables.
##' @param control  an object of class 'tvcm_control'.
##' @param mcall
##' @param step     integer. The current step number.
##'
##' @return A nested list with loss matrices. Partitions of nodes
##'    are nested in partitions for variables. 
##'
##' @details Used in 'tvcm'. 'tvcm_grow_ordnom' and
##'    'tvcm_grow_dev' are help functions for 'tvcm_grow_exsearch'
##'-------------------------------------------------------- #

tvcm_exsearch_nomToOrd <- function(cp, pid, nid, vid, 
                                   mcall, weights,
                                   where, partData,
                                   control) {
  
  ## get partitioning variable
  subs <- where[[pid]] == levels(where[[pid]])[nid]
  z <- partData[, vid]
  
  ## get categories for which the coefficients should be estimated
  levs <- which(colSums(cp) > 0)
    
  ## adjust formula
  ff <- mcall$formula
  aTerms <- attr(terms(ff),"term.labels")
  lTerms <- grep("Left",aTerms,value=TRUE)
  ff <- update(ff,paste(".~.-", paste(aTerms,collapse = "-")))
  lTerms <- unlist(lapply(levs, function(i) {
    lapply(lTerms, function(t) sub("Left", paste0("Left", i), t))
  })) 
  ff <- update(ff,paste(".~+",paste(lTerms,collapse="+"),"+."))
  mcall$formula <- ff
  
  ## create a dummy for each category
  newdata <- as.data.frame(lapply(levs, function(i) 1 * (z == levels(z)[i])))
  colnames(newdata) <- paste0("Left", levs)
  mcall$data <- cbind(mcall$data, newdata)
  
  ## fit the model
  model <- tvcm_grow_fit(mcall, doFit = TRUE)
  
  if (!inherits(model, "try-error")) {
    
    ## extract coefficients
    st <- coef(model)
    stLabs <- strsplit(names(st), ":")
    stLabs <- sapply(stLabs, function(x) x[grep("Left[1-9]+", x)])
    ind <- sapply(paste0("Left", levs), function(x) which(stLabs == x))
    if (!is.matrix(ind)) ind <- matrix(ind, nrow = 1)
    ind <- t(ind)
    st <- matrix(st[c(ind)], nrow(ind), ncol(ind))
    
    ## get ordering
    score <- rep.int(0, nlevels(z))
    score[levs] <- prcomp(st)$x[,1]
    
    ## transform 'z'              
    z <- ordered(z, levels = levels(z)[order(score)])
    cp <- tvcm_getOrdSplits(z[subs], weights[subs],
                            control$minsize[pid],
                            control$maxordsplit)
    cp <- cp[, levels(partData[, vid]),drop=FALSE]
  } else {

    ## avoid splitting (better solution?)
    cp <- cp[-(1:nrow(cp)),,drop=FALSE]
  }
  attr(cp, "type") <- "coef"

  ## return list
  return(list(cp, rep(NA, nrow(cp)), rep(NA, nrow(cp))))
}

tvcm_exsearch_dev <- function(cutpoint, 
                              pid, nid, vid, 
                              model, start,
                              modelNuis, startNuis, 
                              nuisance,
                              where, partData, 
                              control, loss0,
                              mfName) {
      
  ## set node indicator
  subs <- where[[pid]] == levels(where[[pid]])[nid]
  z <- partData[, vid]    
  if (is.numeric(z)) {
    zs <- z <= cutpoint
  } else {
    zs <- z %in% levels(z)[cutpoint > 0L]            
  }

  if (control$fast) {
    model[[mfName]]$Left <- 1 * (subs & zs)
    model[[mfName]]$Right <- 1 * (subs & !zs)
  } else {
    model[[mfName]][subs & zs, paste0("Node", LETTERS[pid])] <- "Left"
    model[[mfName]][subs & !zs, paste0("Node", LETTERS[pid])] <- "Right"
    model[[mfName]][!subs, paste0("Node", LETTERS[pid])] <-
      as.integer(droplevels(where[[pid]][!subs]))
  }
  model$coefficients <- vcrpart_copy(start)
  eta0 <- model$offset
  model <- tvcm_grow_update(model, control, subs)
  rval <- rep(NA, 2L)
  
  if (!inherits(model, "try-error")) { # new from the 2016-01-10

      dev.fast <- function(object, eta0) {
          return(
              sum(object$family$dev.resids(
                  y = object$y,
                  mu = object$family$linkinv(eta0),
                  wt = object$prior.weights)) -
                      sum(object$family$dev.resids(
                          y = object$y,
                          mu = object$family$linkinv(object$linear.predictors),
                          wt = object$prior.weights)))
      }
      
      rval[1] <- ifelse(inherits(model, "glm") & control$fast, 
                        dev.fast(model, eta0[subs]),
                        loss0 - control$lossfun(model))
      rval[2L] <- length(coef(model)[grep("Right", names(coef(model)))])
    if (is.null(modelNuis)) {
      return(rval)
    } else {
      modelNuis[[mfName]]$Left <- 1 * (subs & zs)
      modelNuis[[mfName]]$Right <- 1 * (subs & !zs)
      modelNuis$coefficients <- vcrpart_copy(startNuis)
      modelNuis <- tvcm_grow_update(modelNuis, control, subs)
      rval[1] <- rval[1L] -
          ifelse(inherits(model, "glm") & control$fast,
                 dev.fast(modelNuis, eta0[subs]),
                 loss0 - control$lossfun(modelNuis))
      rval[2L] <- rval[2L] -
        length(coef(modelNuis)[grep("Right", names(coef(modelNuis)))])
      return(rval)
    } 
  } else {
    return(rval)        
  }
}


tvcm_grow_exsearch <- function(splits, partid, nodeid, varid, 
                               model, nodes, where, partData, 
                               control, mcall, formList, step) {

  verbose <- control$verbose; control$verbose <- FALSE;
  loss0 <- control$lossfun(model)
  
  mfName <- switch(deparse(mcall[[1]]), glm = "model", olmm = "frame")
  
  if (verbose) cat("\n* computing splits ")

  mcall$data <- eval(mcall$data, environment(mcall))
  nodeData <- mcall$data[paste0("Node", names(nodes))]
  
  w <- weights(model)

  if (control$fast)
      mcall$offset <- predict(model, type = "link")
  
  if (inherits(model, "glm")) {
    mcall$x <- TRUE
    mcall$y <- TRUE
    mcall$model <- TRUE
  } else if (inherits(model, "olmm")) {
    if (control$fast) {
      mcall$restricted <- grep("ranefCholFac", names(coef(model)), value = TRUE)
      mcall$start <- coef(model)[mcall$restricted]
    }
  }

  if (control$fast) {
    root <- rep.int(FALSE, length(partid))
  } else {
    root <- sapply(nodes, width) == 1
  }
  
  ff <- tvcm_formula(formList, root,
                     eval(mcall$family, environment(mcall)),
                     environment(mcall), full = FALSE,
                     update = TRUE, fast = control$fast)

  Left <- sample(c(0, 1), nobs(model), replace = TRUE)
  Right <- -(Left - 1)
  Node <- lapply(where, function(x) {
    levs <- c("Left", "Right", seq(1, nlevels(x) - 1, length.out = nlevels(x) - 1))
    return(factor(rep(levs, length.out = length(x)), levels = levs))
  })
  
  for (pid in seq_along(partid)) { 
    if (length(unlist(splits[[pid]])) > 0L) {
     
      mcall$formula <- ff$update[[pid]][[1L]]

      ## for 'control$fast = TRUE'
      mcall$data$Left <- Left
      mcall$data$Right <- Right
      ## for 'control$fast = FALSE'
      mcall$data[, paste0("Node", LETTERS[pid])] <- Node[[pid]]

      sModel <- tvcm_grow_fit(mcall, doFit = FALSE)
      sStart <- vcrpart_copy(sModel$coefficients)
      sStart[intersect(names(sStart), names(coef(model)))] <-
        coef(model)[intersect(names(sStart), names(coef(model)))]
      
      sModel$control <- model$control
      
      if (length(control$nuisance[[pid]]) == 0L) {
        sModelN <- NULL
        sNStart <- NULL
      } else {
        mcallN <- mcall
        mcallN$offset <- predict(model, type = "link")
        mcallN$formula <- ff$update[[pid]][[2L]]
        sModelN <- tvcm_grow_fit(mcallN, doFit = FALSE)
        sNStart <- sModelN$coefficients
        sNStart[intersect(names(sNStart), names(coef(model)))] <-
          coef(model)[intersect(names(sNStart), names(coef(model)))]
        sModelN$control <- model$control
      } 
      
      ## run computation
      for (nid in seq_along(splits[[pid]])) {
        if (length(unlist(splits[[pid]][[nid]])) > 0L) {
          for (vid in seq_along(splits[[pid]][[nid]])) {

            ## get reduced cupoints for nominal variables with many categories
            if (attr(splits[[pid]][[nid]][[vid]][[1L]], "type") == "coef" &&
                any(is.na(splits[[pid]][[nid]][[vid]][[2L]])))
              splits[[pid]][[nid]][[vid]] <-
                tvcm_exsearch_nomToOrd(cp = splits[[pid]][[nid]][[vid]][[1L]],
                                       pid = partid[pid],
                                       nid = nodeid[[partid[pid]]][nid],
                                       vid  = varid[[partid[pid]]][vid],
                                       mcall = mcall, weights = w,
                                       where = where, partData = partData,
                                       control = control)
            
            ## extract cutpoints
            cp <- splits[[pid]][[nid]][[vid]][[1L]]
            subs <- is.na(splits[[pid]][[nid]][[vid]][[2L]])
            cp <- cp[subs,, drop = FALSE]
            
            if (nrow(cp) > 0L) {

              ## employ exhaustive search
              st <- apply(cp, 1, tvcm_exsearch_dev, 
                          pid = partid[pid],
                          nid = nodeid[[partid[pid]]][nid],
                          vid = varid[[partid[pid]]][vid],
                          model = sModel, start = sStart,
                          modelNuis = sModelN, startNuis = sNStart,
                          nuisance = control$nuisance[[pid]],
                          where = where, partData = partData,
                          control = control, loss0 = loss0,
                          mfName = mfName)

              if (is.matrix(st)) st <- t(st) else st <- matrix(st, ncol = 1L)
              splits[[pid]][[nid]][[vid]][[2L]][subs] <- st[, 1L]
              splits[[pid]][[nid]][[vid]][[3L]][subs] <- st[, 2L]
            }
          }
        }
      }
  }
    mcall$data[, paste0("Node", LETTERS[pid])] <-
        nodeData[, paste0("Node", names(nodes)[pid])]  
}
  
  ## extracts the penalized loss reduction
  pendev <- lapply(splits, function(part) {
    ## partition level
    lapply(part, function(node) {
      ## node level
      lapply(node, function(var) {
        ## variable level
        if (length(var[[2L]]) > 0L) {
          return(var[[2L]] - control$cp *
                 tvcm_complexity(var[[3]], control$dfpar, 1, control$dfsplit))
        } else {
          return(numeric())
        }
      })
    })
  })

  
  ## function to extract the maximum loss reduction for each part/node/var
  getMaxPenDev <- function(x) {
      x <- unlist(x)
      if (length(x) == 0L) return(-Inf)
      x <- na.omit(x)
      if (length(x) > 0L) return(max(x)) else return(-Inf)
  }

  ## the maximum loss reduction
  maxpendev <- max(c(-Inf, na.omit(unlist(pendev))))
  
  if (maxpendev > -Inf) {
    
    ## select the partition, node and variable
    spart <- which(sapply(sapply(pendev, getMaxPenDev), identical, maxpendev))
    if (length(spart) > 1L) spart <- sample(spart, 1L)
    snode <-
      which(sapply(sapply(pendev[[spart]], getMaxPenDev), identical, maxpendev))
    if (length(snode) > 1L) snode <- sample(snode, 1L)
    svar <- which(sapply(sapply(pendev[[spart]][[snode]], getMaxPenDev),
                         identical, maxpendev))
    if (length(svar) > 1L) svar <- sample(svar, 1L)
    
    ## select the cut
    stat <- splits[[spart]][[snode]][[svar]]
    cutid <- which(stat[[2L]] == max(stat[[2L]], na.rm = TRUE))
    if (length(cutid) > 1L) cutid <- sample(cutid, 1L)
    
    if (verbose) cat("OK")
    
    return(list(partid = partid[spart],
                nodeid = nodeid[[partid[spart]]][snode],
                varid = varid[[partid[spart]]][svar],
                cutid = cutid,
                cut = stat[[1L]][cutid, ],
                dev = as.numeric(stat[[2L]][cutid]),
                pendev = maxpendev,
                npar = as.numeric(stat[[3L]][cutid]),
                grid = splits))
  } else {
    
    if (verbose) cat("failed")
    
    return(list(partid = NULL, nodeid = NULL, varid = NULL, 
                cutid = NULL, cut = NULL, dev = NULL,
                pendev = NULL, grid = splits))
    
  }
}


##'-------------------------------------------------------- #
##' Incorporates a new binary split into an existing
##' tree structure.
##'
##' @param nodes    an object of class 'partynode'.
##' @param loss     a list produced by 'tvcm_grow_exsearch'.
##' @param partData a 'data.frame' with the partitioning variables.
##' @param step     integer. The current algorithm step.
##'
##' @return A list of formulas ('root', 'tree' and 'original').
##'
##' Used in 'tvcm'.
##'-------------------------------------------------------- #

tvcm_grow_splitnode <- function(nodes, where, dev, partData, step, weights) {

  pid <- dev$partid
  nid <- dev$nodeid
  vid <- dev$varid
  nidLab <- nodeids(nodes[[pid]], terminal = TRUE)[nid]
  cut <- dev$cut
  x <- partData[, vid]
  
  ## collect information for the split
  subs <- where[[pid]] == levels(where[[pid]])[nid]
  if (is.numeric(x)) { # numerical variables
    breaks <- as.double(cut)
    index <- NULL
    ordered <- TRUE
    subsL <- subs & x <= breaks
    subsR <- subs & x > breaks
  } else {
    subsL <- subs & x %in% levels(x)[cut == 1L]
    subsR <- subs & x %in% levels(x)[cut == 0L]
    if (is.ordered(x)) { 
      breaks <- as.double(max(which(cut == 1)))
      index <- NULL
      ordered <- TRUE
    } else {
      breaks <- NULL
      index <- as.integer(-cut + 2)
      index[table(x[subs]) == 0] <- NA
      ordered <- FALSE
    }
  }
    
  ## get current nodes
  oldnodes <- as.list(nodes[[pid]])
  
  ## setup 'newnodes' object
  subsN <- which(sapply(oldnodes, function(node) node$id) == nidLab)
  newnodes <- vector("list", length(oldnodes) + 2)
  newnodes[1:subsN] <- oldnodes[1:subsN]
  if (length(oldnodes) > nidLab)
    newnodes[(subsN + 3L):length(newnodes)] <-
      oldnodes[(subsN + 1L):length(oldnodes)]
  
  ## adjust ids of children
  ids <- sapply(newnodes, function(x) if (!is.null(x$id)) x$id else -1)
  for (i in 1L:length(newnodes))
    if (!is.null(newnodes[[i]]$kids))
      newnodes[[i]]$kids <- which(ids %in% newnodes[[i]]$kids)
  
  ## setup new split
  newnodes[[subsN]]$split <-
    partysplit(varid = vid, breaks = breaks, index = index,
               info = list(ordered = ordered, step = step))
  newnodes[[subsN]]$kids <- nidLab + 1L:2L
  newnodes[[subsN]]$info$dims <- c(n = sum(weights[subs]))
  
  ## add new children
  newnodes[[subsN + 1L]] <-
    list(id = nidLab + 1L,
         info = list(
           dims = c(n = sum(weights[subsL])),
           depth = newnodes[[subsN]]$info$depth + 1L))
  
  newnodes[[subsN + 2L]] <-
    list(id = nidLab + 2L,
         info =
         list(dims = c(n = sum(weights[subsR])),
              depth = newnodes[[subsN]]$info$depth + 1L))
  
  ## adjust ids
  for (i in 1L:length(newnodes))
    newnodes[[i]]$id <- i
  
  ## return new nodes
  nodes[[pid]] <- as.partynode(newnodes)
  return(nodes)
}


##'-------------------------------------------------------- #
##' Extracts the formula for the root node and the tree
##' from the output of \code{\link{vcrpart_formula}}.
##'
##' @param formList a list of formulas from 'vcrpart_formula'.
##' @param root     logical vector of the same length as the
##'    'vc' slot of 'formList'.
##' @param family   an object of class 'family' or 'family.olmm'.
##' @param env      the environment where the output formula
##'    is to be evaluated.
##' @param full     whether the full formula should be derived
##' @param update   whether the formula for the update model
##'    should be derived.
##' @return A list of formulas ('root', 'tree' and 'original').
##'
##' @details Used in \code{\link{predict.fvcm}} and
##'    \code{\link{tvcm}}.
##'-------------------------------------------------------- #

tvcm_formula <- function(formList, root, family,
                         env = parent.frame(),
                         full = TRUE, update = FALSE,
                         fast = TRUE) {

    
  yName <- rownames(attr(terms(formList$original), "factors"))[1L]

  ## puts the predictors for fixed effects and varying effects
  ## into one formula
  getTerms <- function(x, effect, root, family, nuisance = NULL) {
    
    ## get 'vc' terms
    vcTerms <- x$vc
    if (!is.null(vcTerms)) {
      vcTerms <-
        lapply(x$vc, function(x) attr(terms(x$eta[[effect]],
                                            keep.order = TRUE), "term.labels"))
      for (i in 1:length(root)) {
        if (root[i]) {
          vcTerms[[i]] <-
            vcTerms[[i]][vcTerms[[i]] != paste0("Node", LETTERS[i])]
          vcTerms[[i]] <- gsub("Node[A-Z]:", "", vcTerms[[i]])    
        }
      }
      vcTerms <- unlist(vcTerms)
    }

    ## get 'fe' terms
    feTerms <- x$fe$eta[[effect]]
    if (!is.null(feTerms))
      feTerms <- attr(terms(feTerms, keep.order = TRUE), "term.labels")
    
    vcTerms <- setdiff(vcTerms, feTerms)
    
    rval <- ""
    if (length(vcTerms) > 0L)
      rval <- paste0(rval, paste(vcTerms, collapse = "+"))
    if (length(vcTerms) > 0L & length(feTerms) > 0L)
      rval <- paste0(rval, "+")
    if (length(feTerms) > 0L)
      rval <- paste0(rval, paste(feTerms, collapse = "+"))
    if (rval != "" && inherits(family, "family.olmm"))
      rval <- paste0(effect, "(", rval, ")")
   
    return(c(vcTerms, feTerms))
  }
  
  ## incorporate fixed effects terms
  feCeTerms <- getTerms(formList, "ce", root, family)
  feGeTerms <- getTerms(formList, "ge", root, family) 
  
  ## intercept
  feInt <- formList$fe$intercept
  vcInt <- unlist(lapply(formList$vc, function(x) x$intercept))
  direct <- sapply(formList$vc, function(x) x$direct)
  if (!is.null(vcInt) && any(direct) && root[direct])
    feInt <- "ce"
  
  ## random effects
  if (!is.null(formList$re)) {
    subjectName <- attr(terms(formList$re$cond), "term.labels")      
    getReTerms <- function(effect) {
      rval <- attr(terms(formList$re$eta[[effect]], keep.order = TRUE), "term.labels")
      if (formList$re$intercept == effect) rval <- c("1", rval)
      if (length(rval) == 0L) return(NULL)
      rval <- paste(rval, collapse = "+")
      if (inherits(family, "family.olmm"))
        rval <- paste0(effect, "(", rval, ")")
      return(rval)
    }
    reForm <- unlist(lapply(c("ce", "ge"), getReTerms))
    reForm <- paste(reForm, collapse = "+")
    if (inherits(family, "family.olmm")) {        
      reForm <- paste0("re(", reForm, "|", subjectName,
                      ",intercept='", formList$re$intercept, "')")
    } else {
        reForm <- NULL
    }
  } else {
    reForm <- NULL
  }

  getForm <- function(yName, feCeTerms, feGeTerms, feInt, reForm, family, env) {

    fTree <- "" # the return value
    
    feCeForm <- if (length(feCeTerms) == 0L) "" else paste(feCeTerms, collapse = "+")
    if (feCeForm != "" & inherits(family, "family.olmm"))
      feCeForm <- paste0("ce(", feCeForm, ")")

    feGeForm <- if (length(feGeTerms) == 0L) "" else paste(feGeTerms, collapse = "+")
    if (feGeForm != "" & inherits(family, "family.olmm"))
      feGeForm <- paste0("ge(", feGeForm, ")")

    if (feCeForm != "") fTree <- paste0(fTree, feCeForm)
    if (feCeForm != "" & feGeForm != "") fTree <- paste0(fTree, " + ")
    if (feGeForm != "") fTree <- paste0(fTree, feGeForm)
    if (fTree == "") fTree <- "1"

    if (inherits(family, "family.olmm")) {
      fTree <- paste0(fTree, ", intercept='", feInt, "'")
    } else {
      if (feInt == "none")
        fTree <- paste("-1", fTree, sep = "+")
    }
    
    if (inherits(family, "family.olmm"))
      fTree <- paste0("fe(", fTree, ")")

    if (!is.null(reForm)) fTree <- paste(fTree, "+", reForm)

    fTree <- paste(yName, "~", fTree)
    
    return(as.formula(fTree, env = env))
  }

  ## full formula
  
  fFull <- NULL
  if (full)
    fFull <- getForm(yName, feCeTerms, feGeTerms, feInt, reForm, family, env)
  
  ## update formulas
  
  fUpdate <- NULL

  feCeTerms <- getTerms(formList, "ce", rep(FALSE, length(root)), family)
  feGeTerms <- getTerms(formList, "ge", rep(FALSE, length(root)), family)
  
  if (update) {

    ## get nuisance terms
    nuisance <- lapply(seq_along(formList$vc), function(pid) {
      return(formList$vc[[pid]]$nuisance)
    })
  
    ## update formulas for tvcm_grow_exsearch
    fUpdate <- vector("list", length(formList$vc))
    for (pid in seq_along(fUpdate)) {
      fUpdate[[pid]] <- vector("list", 2L)
      nLab <- paste0("Node", LETTERS[pid])

      if (fast) {
        feIntTmp <- "none"
        feCeTmp <- feCeTerms[grep(nLab, feCeTerms)]
        feGeTmp <- feGeTerms[grep(nLab, feGeTerms)]
        feCeTmp <- c(gsub(nLab, "Left", feCeTmp), gsub(nLab, "Right", feCeTmp))
        feGeTmp <- c(gsub(nLab, "Left", feGeTmp), gsub(nLab, "Right", feGeTmp))
      } else {
        rootTmp <- root
        rootTmp[pid] <- FALSE
        feIntTmp <- formList$fe$intercept
        if (!is.null(vcInt) && any(direct) && rootTmp[direct])
          feIntTmp <- "ce"
        feCeTmp <- getTerms(formList, "ce", rootTmp, family)
        feGeTmp <- getTerms(formList, "ge", rootTmp, family)         
      }

      ## full formula
      fUpdate[[pid]][[1L]] <-
          getForm(yName,feCeTmp,feGeTmp,feIntTmp,reForm, family, env)
      
      ## null formula (always use approximative model, even if fast = FALSE)
      feCeTmp <- feCeTerms[grep(nLab, feCeTerms)]
      feCeTmp <- intersect(feCeTmp, nuisance[[pid]])
      
      feGeTmp <- feGeTerms[grep(nLab, feGeTerms)]
      feGeTmp <- intersect(feGeTmp, nuisance[[pid]])

      feCeTmp <- c(gsub(nLab, "Left", feCeTmp), gsub(nLab, "Right", feCeTmp))
      feGeTmp <- c(gsub(nLab, "Left", feGeTmp), gsub(nLab, "Right", feGeTmp))
           
      fUpdate[[pid]][[2L]] <-
          getForm(yName,feCeTmp,feGeTmp,"none",reForm,family, env)
    }
  }
  return(list(full = fFull, update = fUpdate))
}


##'-------------------------------------------------------- #
##' Update the 'control_tvcm' object for internal purposes.
##'
##' @param control  an object of class 'tvcm_control'.
##' @param model    a root node regression model, e.g., an 'olmm'
##'    or a 'glm' object
##' @param formList a list of formulas from 'vcrpart_formula'.
##' @param root parm.only
##'
##' @return An updated 'tvcm_control' object.
##'
##' @details Used in 'tvcm' and 'tvcm_grow'.
##'-------------------------------------------------------- #

tvcm_grow_setcontrol <- function(control, model, formList, root, parm.only = TRUE) {

  family <- model$family
  if (is.null(formList$vc)) return(control)

  ## specify the tree size parameters separately for each partition
  
  if (!parm.only) {

    npart <- length(formList$vc)
    
    if (!length(control$maxwidth) %in% c(1L, npart))
      stop("'maxwidth' must be either of length ", 1L, " or ", npart, ".")
    control$maxwidth <- rep_len(control$maxwidth, npart)
    
    if (!is.null(control$minsize) && !length(control$minsize) %in% c(1L, npart))
      stop("'minsize' must be either of length ", 1L, " or ", npart, ".")
    if (!is.null(control$minsize))
      control$minsize <- rep_len(control$minsize, npart)
    
    if (!length(control$maxdepth) %in% c(1L, npart))
      stop("'maxdepth' must be either of length ", 1L, " or ", npart, ".")
    control$maxdepth <- rep_len(control$maxdepth, npart)
  }
    
  ## update the 'parm' slot
  control$parm <- tvcm_grow_setparm(model, formList, root, control$direct)
  
  ## set 'nuisance' slots
  control$nuisance <- lapply(formList$vc, function(x) x$nuisance)
  control$estfun.args$nuisance <-
    unique(c(control$estfun.args$nuisance,
             setdiff(names(coef(model)), names(fixef(model)))))

  return(control)
}


##'-------------------------------------------------------- #
##' Update the 'parm' slot
##'
##' @param model a root node regression model, e.g., an 'olmm'
##'    or a 'glm' object
##' @param formList a list of formulas from 'vcrpart_formula'.
##' @param root a logical vector. Indicates for each varying
##'    coefficient whether there is at least one in the current
##'    model.
##' @param direct logical scalar. Whether there is a direct
##'    effect defined.
##'
##' @return An updated 'tvcm_control' object.
##'
##' @details Used in 'tvcm_grow_setcontrol' and
##'    'tvcm_get_estimates'.
##'-------------------------------------------------------- #

tvcm_grow_setparm <- function(model, formList, root, direct) {
    
    rval <- lapply(formList$vc, function(x) {
        lapply(x$eta, function(x) attr(terms(x), "term.labels"))
    })
    
    for (pid in seq_along(rval)) {
        for (j in names(rval[[pid]])) {
            terms <- rval[[pid]][[j]]
            if (root[pid]) {
                terms <- terms[terms != paste0("Node", LETTERS[pid])]
                terms <- gsub("Node[A-Z]:", "", terms)
            }  
            if (length(terms) > 0L) {
                type <- paste("fe", j, sep = "-")
                X <- model.matrix(model, which = type)
                assign <- attr(X, "assign")
                subs <- which(attr(terms(model, type), "term.labels") %in% terms)
                if (length(subs) == 0L) {
                    terms <- sapply(terms, function(x) {
                        x <- strsplit(x, ":")[[1L]]
                        len <- length(x)
                        x <- c(if (len > 2) x[1:(len-2)], x[len], x[len-1])
                        return(paste(x, collapse = ":"))
                    })
                    subs <- which(attr(terms(model, type), "term.labels") %in% terms)
                }
                if (length(subs) > 0L) terms <- colnames(X)[assign %in% subs]
            }
            rval[[pid]][[j]] <- terms
        }
    }  
    if (inherits(model$family, "family.olmm")) {
        for (pid in seq_along(rval)) {
            if ((len <- length(rval[[pid]][[1L]])) > 0L)
                rval[[pid]][[1L]] <-
                    paste0("Eta", rep(1L:model$dims["nEta"], each = len), ":",
                          rep(rval[[pid]][[1L]], model$dims["nEta"]))
        }
    }

    ## set the 'intercept' slot (which is always in the first 'vc' term)
    if (direct && root[1L]) {
        if (inherits(model$family, "family.olmm")) {
            rval[[1L]]$ce <- c(grep("Eta[1-9]+:\\(Intercept\\)",
                                    names(coef(model)), value = TRUE),
                               rval[[1L]]$ce)
        } else {
            rval[[1L]]$ce <- c("(Intercept)", rval[[1L]]$ce)
        }
    }

    ## return the list
    return(rval)
}


##'-------------------------------------------------------- #
##' Extract the node vector from 'newdata' and assign
##' the contrasts.
##'
##' @param object  a fitted 'tvcm' object.
##' @param newdata a data.frame from which the nodes are to be
##'    extracted.
##' @param weights a numeric vector of weights with the same
##'    size than data.
##' @param formList a list of formulas as produced by
##'    \code{vcrpart_formula}.
##'
##' @return A list with values of the variables in 'data'
##'
##' @details Used in tvcm.predict and tvcm.prune.
##'    \code{tvcm_get_fitted} is merely a help function
##'-------------------------------------------------------- #

tvcm_get_fitted <- function(pid, object, newdata, weights, setContrasts) {
  fitted <- fitted_node(object$info$node[[pid]],
                        newdata[,colnames(object$data), drop = FALSE])
  fitted <- factor(fitted, nodeids(object$info$node[[pid]], terminal = TRUE))
  names(fitted) <- rownames(newdata) 
  if (nlevels(fitted) > 1L) {
    if (setContrasts) {
      contrasts(fitted) <- contr.wsum(fitted, weights)
    } else {
      contrasts(fitted) <-
        object$contrasts[, paste0("Node", LETTERS[pid])]
    }
  }
  return(fitted)
}

tvcm_get_node <- function(object, newdata, setContrasts = FALSE, weights,
                          formList = NULL) {
  if (is.null(formList))
    formList <- vcrpart_formula(object$info$formula, object$info$family)
  fitted <- lapply(seq_along(object$info$node), tvcm_get_fitted, object = object,
                   newdata = newdata, weights = weights, setContrasts = setContrasts)
  names(fitted) <- paste0("Node", LETTERS[seq_along(object$info$node)])
  return(fitted)
}


##'-------------------------------------------------------- #
##' Creates a list which can be used to extract the varying
##' coefficients.
##'
##' @param names character vector. Names of coefficients of the
##'    current model.
##' @param ids   character vector. Names of the current nodes.
##' @param parm  the 'control$parm' slot
##'
##' @return A list with slots
##'    names:     the original coefficient names.
##'    terms:     the names of the terms to which coefficients
##'               belong according to the original formula.
##'    type:      which type of 'fe', 'vc' or 're' the term
##'               belongs to.
##'    node:      the node to which a coefficient belongs to.
##'    partition: the partition to which a coefficients
##'               belongs to.
##'-------------------------------------------------------- #

tvcm_get_terms <- function(names, ids, parm, fixed = NULL) {

    ## modify 'parm':    
    ## modification 1: remove 'by' variables of 'vc' terms
    ## without splits in cases a corresponding fixed i.e. main effect
    ## is specified for the same variable.
    if (!is.null(fixed))
        parm <- lapply(parm, lapply, function(x) setdiff(x, unlist(fixed)))
    ## modification 2: remove duplicated variables that
    ## appear as 'by' variables in multiple 'vc' terms without splits
    ## By definition, such 'by' variables are assigned to the
    ## corresponding first 'vc' term.
    if (any(subs <- duplicated(unlist(parm)))) {
        doubles <- unique(unlist(parm)[subs])
        for (i in seq_along(doubles)) {
            first <- which(sapply(parm, function(x) doubles[i] %in% unlist(x)))[1L]
            parm <- lapply(parm[-first], lapply, function(x) setdiff(x, doubles[i]))
        }
    }
    
    ## change 'names' if no split was applied in some node
    if (any(unlist(ids) == 1L)) {
        getNames <- function(x) {           
            if (!x %in% unlist(parm)) return(x)
            pid <- which(sapply(parm, function(p) x %in% unlist(p)))            
            if (grepl("(Intercept)", x))
                return(gsub("(Intercept)", paste0("Node", LETTERS[pid], 1),
                            x, fixed = TRUE))
            if (!grepl("Node", x)) {
                x <- strsplit(x, ":")
                x <- rep(x, length(pid))
                subs <- ifelse(grepl("Eta[1-9]+", x[[1L]][1L]), 2L, 1L)
                for (i in 1:length(x)) {
                    x[[i]][subs] <-
                        paste0("Node", LETTERS[pid[i]], 1, ":", x[[i]][subs])
                    x[[i]] <- paste(x[[i]], collapse = ":")
                }
            }
            return(x)
        }
        names <- unlist(lapply(names, getNames))
    }
    
    nodes <- unlist(lapply(seq_along(ids), function(i)
        paste0("Node", names(ids)[i], ids[[i]])))
    split <- strsplit(names, ":")
    type <-
        sapply(split, function(x) {
            rval <- "fe"
            if (any(x %in% nodes)) rval <- "vc"
            if (any(substr(x, 1, 12) %in% "ranefCholFac")) rval <- "re"
            return(rval)
        })
    terms <-
        sapply(split, function(x) {
            paste(x[!x %in% nodes], collapse = ":")
        })
    node <- 
        sapply(split, function(x) {
            rval <- ""
            if (any(subs <- x %in% nodes))
                rval <- substr(x[subs], 6, 500)
            return(rval)
        })
    partition <-
        sapply(split, function(x) {
            rval <- ""
            if (any(subs <- x %in% nodes))
                rval <- substr(x[subs], 5, 5)
            return(rval)
        })
    return(list(names = names, 
                terms = terms,
                type = type,
                node = node,
                partition = partition))
}


##'-------------------------------------------------------- #
##' Extracts the names of the predictors on 'vc' terms
##' 
##' @param object a 'tvcm' object.
##'
##' @return A character vector.
##'-------------------------------------------------------- #

tvcm_get_vcparm <- function(object) {
  vcTerms <- terms(object$info$formula$original, specials = "vc")
  vcTerms <- rownames(attr(vcTerms, "factors"))[attr(vcTerms, "specials")$vc]
  parm <- lapply(vcTerms, function(x) {
    eta <- eval(parse(text = x))$eta
    if (inherits(object$info$model, "olmm")) {
      etaList <- vcrpart_formula(eta)
      parmCe <- all.vars(etaList$fe$eta$ce)
      if (length(parmCe) > 0L) {
        parmCe <-
          paste0("Eta", 1:object$info$model$dims["nEta"], ":", parmCe)
      } else {
        parmCe <- NULL
      }
      parmGe <- all.vars(etaList$fe$eta$ge)
      parm <- c(parmCe, parmGe)
    } else {
      parm <- all.vars(eta)
    }
    return(parm)
  })
  for (pid in seq_along(object$info$formula$vc)) {
    int <- object$info$formula$vc[[pid]]$intercept
    if (inherits(object$info$model, "olmm")) {
      if (int == "ce") {
        parm[[pid]] <- c(paste0("Eta", 1:object$info$model$dims["nEta"],
                               ":(Intercept)"), parm[[pid]])
      } else if (int == "ge") {
        parm[[pid]] <- c("(Intercept)", parm)
      }
    } else {
      if (int != "none") parm[[pid]] <- c("(Intercept)", parm[[pid]])
    }
  }
  return(parm)
}


##'-------------------------------------------------------- #
##' Extracts the estimates a fitted 'tvcm' object and
##' creates a list used in various methods.
##' variances of a fitted 'tvcm' object.
##'
##' @param object an object of class 'partynode'.
##' @param what   the type of estimate to be extracted.
##'
##' @return A list of 'fe' (fixed effects), 'vc' (varying
##' coefficients) and 're' (random effects) coefficients. The
##' 'vc' slot consists of separate matrices for each varying
##' coefficient partition.
##'
##' @details Used in 'coef', 'extract'.
##'-------------------------------------------------------- #

tvcm_get_estimates <- function(object, what = c("coef", "sd", "var"), ...) {

    ## extract slots for the readability of the code
    what <- match.arg(what)
    model <- object$info$model
    control <- object$info$control

    ## prepare return value
    rval <- list(fe = numeric(),
                 vc = replicate(length(object$info$node), matrix(,0,0)),
               re = numeric())
    names(rval$vc) <-  LETTERS[seq_along(object$info$node)]
    
    ## extract coefficients
    estimates <- switch(what,
                        coef = coef(model),
                        sd = diag(vcov(model)),
                        var = diag(vcov(model)))

    ## ids of terminal nodes of tree structures
    ids <- lapply(object$info$node, nodeids, terminal = TRUE)

    
    ## the 'vc' terms that should exist in theory
    vcVars <- lapply(object$info$formula$vc, function(x) {
        rval <- lapply(x$eta, function(x) attr(terms(x), "term.labels"))
        if (length(rval$ce) > 0 && inherits(model$family, "family.olmm")) {
            rval$ce <-
                paste0("Eta", rep(1:model$dims["nEta"], each = length(rval$ce)),
                       ":", rep(rval$ce, model$dims["nEta"]))
        }
        return(lapply(rval, function(x) {
            x <- strsplit(x, ":")
            return(lapply(x, function(x) {
                x <- x[-grep("Node[A-Z]", x)]
              return(paste(x, collapse = ":"))
            }))
        }))
    })

    ## get fixed i.e. main effects
    feTerms <- lapply(object$info$formula$fe$eta,
                      function(x) attr(terms(x), "term.labels"))
    
    ## the terms for which estimates for 'type' (really) exist
    termsE <- tvcm_get_terms(names(estimates), ids, control$parm, feTerms)
    
    ## restricted coefficients
    if (any(termsE$type == "fe"))
        rval$fe <- estimates[termsE$type == "fe"]
    
    ## random effects
    if (any(termsE$type == "re"))
        rval$re <- estimates[termsE$type == "re"]
    
    ## varying coefficients
    if (any(termsE$type == "vc")) {
        
        for (pid in seq_along(object$info$node)) {
            
            nnodes <- length(ids[[pid]]) # number of nodes
            rval$vc[[pid]] <-
                matrix(, nnodes, length(unlist(vcVars[[pid]])),
                       dimnames = list(ids[[pid]], unlist(vcVars[[pid]])))
            
            ## fill the matrix
            vcTermsE <- unique(termsE$terms[termsE$partition == LETTERS[pid]])
            for (i in seq_along(vcTermsE)) {
                subs <- termsE$terms == vcTermsE[i] & termsE$partition == LETTERS[pid]
                subsRows <- termsE$node[subs]
                subsCols <- which(colnames(rval$vc[[pid]]) == vcTermsE[i])
                rval$vc[[pid]][subsRows, subsCols] <- estimates[subs]
            }
            
            ## add column names for varying intercepts
            subs <- which(colnames(rval$vc[[pid]]) %in% "")
            if (length(subs) > 0L) colnames(rval$vc[[pid]])[subs] <- "(Intercept)"
            if (ncol(rval$vc[[pid]]) > 0 &
                inherits(object$info$family, "family.olmm")) {
                tmp <- strsplit(colnames(rval$vc[[pid]]), ":")
                subs <- sapply(tmp, length) == 1L &
                    sapply(tmp, function(x) substr(x[1L], 1L, 3L) == "Eta")
                colnames(rval$vc[[pid]])[subs] <-
                    paste(colnames(rval$vc[[pid]])[subs], "(Intercept)", sep = ":")
            }
            
            ## fill the last row if necessary
            if (nrow(rval$vc[[pid]]) > 1L) {
                subs <- is.na(rval$vc[[pid]][nnodes, ])
                if (any(subs)) {
                
                    ## compute the coefficients of the omitted node
                    con <- model$contrasts[[paste0("Node", LETTERS[pid])]][nnodes, ]
                    for (i in which(subs)) {
                        rval$vc[[pid]][nnodes, i] <-
                            switch(what,
                                   coef = sum(con * rval$vc[[pid]][-nnodes, i],
                                       na.rm = TRUE),
                                   sd = sum(con^2 * rval$vc[[pid]][-nnodes, i],
                                       na.rm = TRUE),
                                   var = sum(con^2 * rval$vc[[pid]][-nnodes, i],
                                       na.rm = TRUE))
                    }
                }
            }
        }
    }
    
    if (what == "sd") {
        getSqrt <- function(x) {
            if (is.list(x)) {
                return(lapply(x, getSqrt))
            } else if (length(x) > 0) {
                return(sqrt(x))
            } else {
                return(x)
            }
        }
        rval <- lapply(rval, getSqrt)
    }
    return(rval)
}


##'-------------------------------------------------------- #
##' Create short labes for 'vc' terms
##'
##' @param object a \code{\link{tvcm}} object
##' @param intercept logical scalar. Whether '(Intercept)'
##'    should be added to the label.
##' @return A character vector with a label for each 'vc'
##'    term.
##'
##' @details Used in 'tvcm_print' and 'plot.tvcm'.
##'-------------------------------------------------------- #

tvcm_print_vclabs <- function(formList, intercept = FALSE) {
  
  if (length(formList$vc) == 0) return(NULL)
  
  ## conditioning variables
  cond <- lapply(formList$vc, function(x) {
    rval <- all.vars(x$cond)
    if (length(rval) > 2L) rval <- c(rval[1L], "...")
    return(rval)
  })
  
  ## 'by' terms
  vcLabs <- terms(formList$original, specials = "vc")
  if (length(attr(vcLabs, "specials")$vc) == 0L) return(NULL)
  vcLabs <- rownames(attr(vcLabs, "factors"))[attr(vcLabs, "specials")$vc]
  vcLabs <- paste("getBy", vcLabs, sep = "_")
  getBy_vc <- function(..., by, intercept, nuisance) {
    mc <- match.call()
    rval <- deparse(mc$by, 500L)
    if (rval == "NULL") return("") else return(rval)
  }
  by <- sapply(vcLabs, function(x) eval(parse(text = x)))
  if (intercept) {
      hasInt <- sapply(formList$vc, function(x) x$intercept) != "none"
      for (i in which(hasInt)) 
          by[i] <- paste0("(Intercept)", ifelse(by[i] ==  "", "", " + "), by[i])
  }
  
  ## collapse the short labels
  rval <- rep.int("vc(", length(formList$vc))
  for (pid in seq_along(rval)) {
    if (length(cond) > 0L) 
      rval[pid] <- paste0(rval[pid], paste(cond[[pid]], collapse = ", "))
    if (by[pid] != "")
      rval[pid] <- paste0(rval[pid], ", by = ", by[pid])
  }
  rval <- paste0(rval, ")")
  return(rval)
}


##'-------------------------------------------------------- #
##' The main function to prune an object of class 'partynode'
##' according to some criteria.
##'
##' @param nodes   an object of class 'partynode'.
##' @param alpha   criteria according to which 'nodes' is to be
##'    pruned.
##' @param maxstep the maximal number of step (splits) to be
##'    accomplished
##' 
##' @return A list of class 'partynode'.
##'
##' @details Used in 'tvcm.prune'.
##'-------------------------------------------------------- #

tvcm_prune_node <- function(object, alpha = NULL, maxstep = NULL, terminal = NULL) {

  stopifnot(class(object)[1] %in% c("tvcm", "party", "partynode"))
  
  if ("partynode" %in% class(object)) {
    rval <- list(object)
  } else {
    rval <- object$info$node    
  }

  if (all(c(is.null(alpha), is.null(maxstep), is.null(terminal))))
    return(object$info$node)

  if (!is.null(alpha) && depth(rval[[pid]]) > 0L) {
      splitpath <- object$info$splitpath
      p.value <- extract(object, "p.value")
      ms <- c(0, which(!is.na(p.value) & p.value <= alpha))
      if (length(ms) > 0L)
          ms <- max(ms[c(1, diff(ms)) == 1])
      maxstep <- min(maxstep, ms)
  }
  
  for (pid in seq_along(rval)) {

    ## prune the tree

    if (!is.null(maxstep))
      rval[[pid]] <- tvcm_prune_maxstep(rval[[pid]], maxstep)

    if (!is.null(terminal[[pid]]))
      rval[[pid]] <- tvcm_prune_terminal(rval[[pid]], terminal[[pid]])
    
    ## delete empty nodes and adjust id labeling
    tmp <- as.list(rval[[pid]])
    tmp <- tmp[sapply(tmp, function(x) !is.null(x))] # delete nodes
    ids_old <- unlist(lapply(tmp, function(x) x$id))
    ids_new <- 1L:length(tmp)
    for (i in ids_new) {
        tmp[[i]]$info$id$last <- tmp[[i]]$id
        tmp[[i]]$id <- ids_new[ids_old == tmp[[i]]$id]
        for (j in 1:length(tmp[[i]]$kids))
            tmp[[i]]$kids[j] <- ids_new[ids_old == tmp[[i]]$kids[j]]
    }
    rval[[pid]] <- as.partynode(tmp)
  }
  
  return(rval)
}


##'-------------------------------------------------------- #
##' Recursive function to delete nodes according to
##' the step in the algorithm the nodes were created.
##'
##' @param node    an object of class 'partynode'.
##' @param maxstep an integer scalar. The maximum step allowed.
##' 
##' @return A list of class 'partynode'.
##'
##' @details Used in 'tvcm_prune_node'.
##'-------------------------------------------------------- #

tvcm_prune_maxstep <- function(node, maxstep) {
  
  if (!is.terminal(node)) {
    if (node$split$info$step <= maxstep) {
      kids <- sapply(node$kids, function(kids) kids$id)
      for (i in 1:length(kids)) {
        node$kids[[i]] <- tvcm_prune_maxstep(node$kids[[i]], maxstep)
      }
    } else {
      node$kids <- NULL
      node$split <- NULL
    }
  }
  return(node)
}


##'-------------------------------------------------------- #
##' Recursive function to delete nodes according node labels
##'
##' @param node     an object of class 'partynode'.
##' @param terminal an integer vector. The nodes considered
##'    as terminal nodes
##' 
##' @return A list of class 'partynode'.
##'
##' @details Used in 'tvcm_prune_node'.
##'-------------------------------------------------------- #

tvcm_prune_terminal <- function(node, terminal) {
  
  if (!is.terminal(node)) {
    if (!node$id %in% terminal) {
      kids <- sapply(node$kids, function(kids) kids$id)
      for (i in 1:length(kids)) {
        node$kids[[i]] <- tvcm_prune_terminal(node$kids[[i]], terminal)
      }
    } else {
      node$kids <- NULL
      node$split <- NULL
    }  }
  return(node)
}


##'-------------------------------------------------------- #
##' Creates a 'splitpath.tvcm' object for tracing the
##' fitting process.
##'
##' @param splitpath the splitpath obtained by the fitting process.
##' @param nodes     an object of class 'partynode'.
##' @param partData  a 'data.frame' with the partitioning variables.
##' @param control   an object of class 'tvcm_control'.
##'
##' @return A list of class 'splitpath.tvcm'.
##'
##' @details Used in 'tvcm' and 'prune.tvcm'.
##'-------------------------------------------------------- #

tvcm_grow_splitpath <- function(splitpath, varid, nodes, partData, control) {
    
  if (all(lapply(nodes, width) < 2L)) return(splitpath)  

  ## get the terminal node labels for each varying coefficient partition
  ids <- lapply(nodes, nodeids)
  
  ## assign each terminal node the step in which it was created
  steps <- lapply(nodes, function(node) {
    ids <- nodeids(node)
    rval <- nodeapply(node, ids, function(node) node$split$info$step)
    rval <- sapply(rval, function(x) if (is.null(x)) Inf else x)
    return(rval)
  })
  
  ## upate each step
  for (step in seq_along(splitpath)) {
    
    ## the selected partition
    partid <- splitpath[[step]]$partid
    
    ## get the nodes at disposition at the current step
    kidids <- vector("list", length(nodes))
    for (pid in seq_along(nodes)) {
      parentids <- ids[[pid]][steps[[pid]] < step]
      if (length(parentids) == 0L) {
        kidids[[pid]] <- 1L
      } else {
        kids <- nodeapply(nodes[[pid]], parentids, function(node) node$kids)
        for (i in seq_along(kids))
          kidids[[pid]] <-
            c(kidids[[pid]], unlist(lapply(kids[[i]], function(x) x$id)))
        kidids[[pid]] <- setdiff(sort(kidids[[pid]]), parentids)
      }
    }          
    
    ## re-label the descriptions of splits
    if (!is.null(splitpath[[step]]$varid)) {
      splitpath[[step]]$node <- kidids[[partid]][splitpath[[step]]$nodeid]
      splitpath[[step]]$var <- colnames(partData)[splitpath[[step]]$varid]
      splitpath[[step]]$cutpoint <- character_split(nodeapply(nodes[[partid]], splitpath[[step]]$node, function(node) node$split)[[1]], data = partData)$levels
    }
    
    if (control$sctest) {

      ## change row names of p-values tables
      for (pid in seq_along(nodes))
        if (!is.null(splitpath[[step]]$sctest[[pid]]))
          rownames(splitpath[[step]]$sctest[[pid]]) <-
            paste0("Node", LETTERS[pid], kidids[[pid]])
      
    }

    
    ## change the names of the dev grid elements !!!
    if (!is.null(splitpath[[step]]$grid)) {
      names(splitpath[[step]]$grid) <- LETTERS[seq_along(nodes)]
      for (pid in seq_along(splitpath[[step]]$grid)) {
        names(splitpath[[step]]$grid[[pid]]) <-
          paste0("Node", kidids[[pid]])
        for (nid in seq_along(splitpath[[step]]$grid[[pid]])) { 
          names(splitpath[[step]]$grid[[pid]][[nid]]) <-
            colnames(partData)[varid[[pid]]]
        }
      }
    }
  }
    
  class(splitpath) <- "splitpath.tvcm"
  return(splitpath)
}
