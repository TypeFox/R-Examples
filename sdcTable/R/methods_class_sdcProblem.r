#' @aliases get.sdcProblem,sdcProblem,character-method
#' @rdname get.sdcProblem-method
setMethod(f="get.sdcProblem", signature=c("sdcProblem", "character"),
  definition=function(object, type) {
    if ( !type %in% c("dataObj", "problemInstance", "partition", "elapsedTime", "dimInfo", "indicesDealtWith",
        "startI", "startJ", "innerAndMarginalCellInfo") ) {
      stop("get.sdcProblem:: argument 'type' is not valid!\n")
    }
    if ( type == "dataObj" ) {
      return(g_dataObj(object))
    }
    if ( type == "problemInstance" ) {
      return(g_problemInstance(object))
    }
    if ( type == "partition" ) {
      return(g_partition(object))
    }
    if ( type == "elapsedTime" ) {
      return(g_elapsedTime(object))
    }
    if ( type == "dimInfo" ) {
      return(g_dimInfo(object))
    }
    if ( type == "indicesDealtWith" ) {
      return(g_indicesDealtWith(object))
    }
    if ( type == "startI" ) {
      return(g_startI(object))
    }
    if ( type == "startJ" ) {
      return(g_startJ(object))
    }
    if ( type == "innerAndMarginalCellInfo" ) {
      return(g_innerAndMarginalCellInfo(object))
    }
  }
)

#' @aliases set.sdcProblem,sdcProblem,character,list-method
#' @rdname set.sdcProblem-method
setMethod(f="set.sdcProblem", signature=c("sdcProblem", "character", "list"),
  definition=function(object, type, input) {
    if ( !type %in% c("problemInstance", "partition", "rule.freq", "rule.nk",
        "rule.p", "rule.pk", "startI", "startJ", "indicesDealtWith", "elapsedTime") ) {
      stop("set.sdcProblem:: check argument 'type'!\n")
    }
    if ( type == "problemInstance" ) {
      s_problemInstance(object) <- input[[1]]
    }
    if ( type == "partition" ) {
      s_partition(object) <- input[[1]]
    }
    if ( type == "startI" ) {
      s_startI(object) <- input[[1]]
    }
    if ( type == "startJ" ) {
      s_startJ(object) <- input[[1]]
    }
    if ( type == "indicesDealtWith" ) {
      s_indicesDealtWith(object) <- input[[1]]
    }
    if ( type == "elapsedTime" ) {
      s_elapsedTime(object) <- input[[1]]
    }
    validObject(object)
    return(object)
  }
)

#' @aliases calc.sdcProblem,sdcProblem,character,list-method
#' @rdname calc.sdcProblem-method
setMethod(f="calc.sdcProblem", signature=c("sdcProblem", "character", "list"),
  definition=function(object, type, input) {
    if ( !type %in% c("rule.freq", "rule.nk", "rule.p", "rule.pq", "heuristicSolution",
      "cutAndBranch", "anonWorker", "ghmiter", "preprocess", "cellID",
      "finalize", "ghmiter.diagObj", "ghmiter.calcInformation",
      "ghmiter.suppressQuader", "ghmiter.selectQuader",
      "ghmiter.suppressAdditionalQuader", "contributingIndices",
      "reduceProblem", "genStructuralCuts") ) {
      stop("calc.sdcProblem:: check argument 'type'!\n")
    }
    # frequency-rule
    if ( type == "rule.freq" ) {
      return(c_rule_freq(object, input))
    }
    # nk-dominance rule
    if ( type == "rule.nk" ) {
      return(c_rule_nk(object, input))
    }
    # p%-rule
    if ( type == "rule.p" ) {
      return(c_rule_p(object, input))
    }
    # pq-rule
    if ( type == "rule.pq" ) {
      return(c_rule_pq(object, input))
    }
    if ( type == "heuristicSolution" ) {
      return(c_heuristic_solution(object, input))
    }
    if ( type == "cutAndBranch" ) {
      return(c_cut_and_branch(object, input))
    }
    if ( type == "anonWorker" ) {
      return(c_anon_worker(object, input))
    }
    if ( type == "ghmiter" ) {
      return(c_ghmiter(object, input))
    }
    if ( type == "preprocess" ) {
      return(c_preprocess(object, input))
    }
    if ( type == "cellID" ) {
      return(c_cellID(object, input))
    }
    if ( type == "finalize" ) {
      return(c_finalize(object, input))
    }
    if ( type == "ghmiter.diagObj" ) {
      return(c_ghmiter_diag_obj(object, input))
    }
    if ( type == "ghmiter.calcInformation" ) {
      return(c_ghmiter_calc_info(object, input))
    }
    if ( type == "ghmiter.suppressQuader" ) {
      return(c_ghmiter_suppress_quader(object, input))
    }
    if ( type == "ghmiter.selectQuader" ) {
      return(c_ghmiter_select_quader(object, input))
    }
    if ( type == "ghmiter.suppressAdditionalQuader" ) {
      return(c_ghmiter_supp_additional(object, input))
    }
    if ( type == "contributingIndices" ) {
      return(c_contributing_indices(object, input))
    }
    if ( type == "reduceProblem" ) {
      return(c_reduce_problem(object, input))
    }
    if ( type == "genStructuralCuts" ) {
      return(c_gen_structcuts(object, input))
    }
  }
)

#' summarize object of class \code{\link{sdcProblem-class}} or \code{\link{safeObj-class}}.
#'
#' extract and show relevant information stored in object ofs class \code{\link{sdcProblem-class}} or \code{\link{safeObj-class}}.
#'
#' @aliases summary,sdcProblem-method
#' @rdname summary.sdcProblem-method
#' @param object Objects of either class \code{\link{sdcProblem-class}} or \code{\link{safeObj-class}}.
#' @param ... currently not used.
#' @export
#' @docType methods
setMethod(f="summary", signature="sdcProblem",
  definition=function(object, ...) {
    pI <- g_problemInstance(object)
    dO <- g_dataObj(object)
    dI <- g_dimInfo(object)
    if ( g_is_microdata(dO) ) {
      cat("The raw data contains micro data!")
      if ( length(pI@numVars) > 0 ) {
        cat("--> the use of dominance rules for primary suppressions is possible!")
      }
      cat("\n")
    } else {
      cat("The raw data contain pre-aggregated (tabular) data!\n")
    }

    nrcells <- g_nrVars(pI)
    dim_names <- g_varname(dI)

    cat("\nThe complete table to protect consists of",nrcells,"cells and has",length(dim_names),"spanning variables.")

    cat("\nThe distribution of\n")
    cat("- primary unsafe (u)\n")
    cat("- secondary suppressed (x)\n")
    cat("- forced to publish (z) and\n")
    cat("- selectable for secondary suppression (s) cells is shown below:\n")
    print(table(g_sdcStatus(pI)))

    nr_tables <- g_partition(object)$nrTables
    cat("\nIf this table is protected with heuristic methods, a total of",nr_tables,"has (sub)tables must be considered!\n")
  }
)


#' print objects of class \code{\link{sdcProblem-class}}.
#'
#' print some useful information instead of just displaying the entire object (which may be large)
#'
#' @aliases print,sdcProblem-method
#' @rdname print.sdcProblem-method
#' @param x an objects of class \code{\link{sdcProblem-class}}
#' @param ... currently not used.
#' @export
#' @docType methods
setMethod("print", signature="sdcProblem",
  definition=function(x, ...) {
    dims <- x@dimInfo@dimInfo
    nrDims <- length(dims)
    nrCells <- length(x@problemInstance@strID)
    cat(paste("The object is an 'sdcProblem' with",nrCells,"cells in",nrDims, "dimension(s)!\n"))
    cat("\nThe dimensions are:\n")
    for ( i in 1:nrDims ) {
      nrCodes <- length(dims[[i]]@codesOriginal)
      nrAggregates <- sum(dims[[i]]@codesMinimal==FALSE)
      maxHier <- length(dims[[i]]@structure)
      cat(paste0("\t- ",names(dims)[i]," (",maxHier," levels; ",nrCodes," codes; of these being ",nrAggregates," aggregates)\n"))
    }
    cat("\nCurrent suppression pattern:\n")
    cat("\t- Primary suppressions:",sum(x@problemInstance@sdcStatus=="u"),"\n")
    cat("\t- Secondary suppressions:",sum(x@problemInstance@sdcStatus=="x"),"\n")
    cat("\t- Publishable cells:",sum(x@problemInstance@sdcStatus%in% c("s","z")),"\n")
  }
)

#' show objects of class \code{\link{sdcProblem-class}}.
#'
#' just calls the corresponding print-method
#'
#' @aliases show,sdcProblem-method
#' @rdname show.sdcProblem-method
#' @param object an objects of class \code{\link{sdcProblem-class}}
#' @export
#' @docType methods
setMethod("show", signature="sdcProblem",
  definition=function(object) {
    print(object)
  }
)

setMethod("g_problemInstance", signature="sdcProblem", definition=function(object) {
  object@problemInstance
})

setMethod("g_dimInfo", signature="sdcProblem", definition=function(object) {
  object@dimInfo
})

setMethod("g_partition", signature="sdcProblem", definition=function(object) {
  object@partition
})

setMethod("g_elapsedTime", signature="sdcProblem", definition=function(object) {
  object@elapsedTime
})

setMethod("g_dataObj", signature="sdcProblem", definition=function(object) {
  object@dataObj
})

setMethod("g_startI", signature="sdcProblem", definition=function(object) {
  object@startI
})

setMethod("g_startJ", signature="sdcProblem", definition=function(object) {
  object@startJ
})

setMethod("g_indicesDealtWith", signature="sdcProblem", definition=function(object) {
  object@indicesDealtWith
})

setMethod("g_innerAndMarginalCellInfo", signature="sdcProblem", definition=function(object) {
  pI <- g_problemInstance(object)
  strIDs <- g_strID(pI)
  strInfo <- g_str_info(g_dimInfo(object))

  out <- lapply(1:length(strInfo), function(x) {
    sort(unique(mySplit(strIDs, strInfo[[x]][1]:strInfo[[x]][2])))[-1]
  })

  # deal with 'tot' levels
  ind <- which(sapply(out, length) ==0)
  out[ind] <- lapply(ind, function(x) { "0" } )

  innerCells <- apply(matrix(unlist(expand(out)), ncol=length(out), byrow=FALSE),1,paste, collapse="")
  totCells <- setdiff(strIDs, innerCells)
  indexTotCells <- match(totCells, strIDs)
  indexInnerCells <- match(innerCells, strIDs)
  return(list(innerCells=innerCells, totCells=totCells, indexInnerCells=indexInnerCells, indexTotCells=indexTotCells))
})

setMethod("g_df", signature="sdcProblem", definition=function(object, addDups=FALSE, addNumVars=FALSE) {
  xx <- strID <- NULL
  pI <- g_problemInstance(object)
  dt <- data.table(
    strID=g_strID(pI),
    freq=g_freq(pI),
    sdcStatus=g_sdcStatus(pI)
  )
  if (addNumVars & !is.null(pI@numVars) ) {
    dt <- cbind(dt, as.data.table(pI@numVars))
  }
  dI <- g_dimInfo(object)
  strInfo <- g_str_info(dI)
  dimObj <- g_dim_info(dI)
  vNames <- g_varname(dI)
  res <- as.data.table(cpp_splitByIndices(g_strID(pI), strInfo))
  setnames(res, vNames)
  dt <- cbind(dt, res)
  for ( i in 1:length(strInfo) ) {
    v <- paste0(vNames[i],"_o",sep="")
    dt[[v]] <- c_match_orig_codes(object=dimObj[[i]], input=dt[[vNames[i]]])
  }

  if ( addDups ) {
    dims <- g_dim_info(dI)
    for ( i in seq_along(dims) ) {
      if ( g_has_dups(dims[[i]]) ) {
        dU <- dims[[i]]@dupsUp
        dL <- dims[[i]]@dups
        vName <- paste0(dims[[i]]@vName,"_o")
        for ( j in 1:length(dL) ) {
          cmd <- paste0("xx <- dt[",vName,"=='",dU[j],"']")
          eval(parse(text=cmd))
          if ( !is.numeric(dt[[vName]]) ) {
            cmd <- paste0("xx[,",vName,":='",dL[j],"']")
          } else {
            cmd <- paste0("xx[,",vName,":=",dL[j],"]")
          }
          eval(parse(text=cmd))
          dt <- rbind(dt, xx); rm(xx)
        }
      }
    }
  }
  setkey(dt, strID)
  return(dt)
})

setReplaceMethod("s_problemInstance", signature=c("sdcProblem", "problemInstance"), definition=function(object, value) {
  object@problemInstance <- value
  validObject(object)
  object
})

setReplaceMethod("s_partition", signature=c("sdcProblem"), definition=function(object, value) {
  object@partition <- value
  validObject(object)
  object
})

setReplaceMethod("s_startI", signature=c("sdcProblem", "numeric"), definition=function(object, value) {
  object@startI <- value
  validObject(object)
  object
})

setReplaceMethod("s_startJ", signature=c("sdcProblem", "numeric"), definition=function(object, value) {
  object@startJ <- value
  validObject(object)
  object
})

setReplaceMethod("s_indicesDealtWith", signature=c("sdcProblem"), definition=function(object, value) {
  object@indicesDealtWith <- value
  validObject(object)
  object
})

setReplaceMethod("s_elapsedTime", signature=c("sdcProblem"), definition=function(object, value) {
  object@elapsedTime <- value
  validObject(object)
  object
})

setMethod("c_rule_freq", signature=c("sdcProblem", "list"), definition=function(object, input) {
  pI <- g_problemInstance(object)
  if ( input$allowZeros == TRUE ) {
    suppInd <- which(g_freq(pI) <= input$maxN)
  } else {
    f <- g_freq(pI)
    suppInd <- which(f > 0 & f <= input$maxN)
    zeroInd <- which(g_freq(pI) == 0 )
    if ( length(zeroInd) > 0 ) {
      s_sdcStatus(pI) <- list(index=zeroInd, vals=rep("z", length(zeroInd)))
    }
  }
  if ( length(suppInd) > 0 ) {
    s_sdcStatus(pI) <- list(index=suppInd, vals=rep("u", length(suppInd)))
  }
  s_problemInstance(object) <- pI
  validObject(object)
  return(object)
})

setMethod("c_rule_nk", signature=c("sdcProblem", "list"), definition=function(object, input) {
  nkRule <- function(celltot, sumNcont, k) {
    # if TRUE, cell needs to be suppressed
    (sumNcont) > (k/100*celltot)
  }
  if ( !g_is_microdata(g_dataObj(object)) ) {
    stop("nk-dominance rule can only be applied if micro-data are available!\n")
  }
  pI <- g_problemInstance(object)
  dataObj <- g_dataObj(object)
  strIDs <- g_strID(pI)
  numVarInds <- g_numvar_ind(dataObj)

  numVal <- g_raw_data(dataObj)[[numVarInds[input$numVarInd]]]
  if ( any(numVal < 0 ) ) {
    stop("dominance rules can only be applied to numeric variables with only positive values!\n")
  }

  # calculate contributing indices
  indices <- lapply(1:g_nrVars(pI), function(x) {
    c_contributing_indices(object, input=list(strIDs[x]))
  })

  minContributingUnits <- min(setdiff(unique(sapply(indices, length)), 0))
  if ( input$n < 1 | input$n > minContributingUnits ) {
    stop("set.sdcProblem:: parameter 'n' must be >= 1 and <",minContributingUnits,"!\n")
  }

  # values of contributing units
  valueList <- lapply(1:g_nrVars(pI), function(x) {
    sum(rev(tail(sort(numVal[indices[[x]]]), input$n)))
  })

  cellTotals <- g_numVars(pI)[[input$numVarInd]]
  # suppStatus: TRUE:unsafe, FALSE: safe
  nkState <- sapply(1:g_nrVars(pI), function(x) {
    nkRule(cellTotals[x], valueList[[x]], input$k)
  })

  addSupps <- which(sapply(indices, length) %in% 1:input$n)
  if ( length(addSupps) > 0 ) {
    nkState[addSupps] <- TRUE
  }

  suppIndex <- which(nkState==TRUE)
  if ( length(suppIndex) > 0 ) {
    s_sdcStatus(pI) <- list(index=suppIndex, vals=rep("u", length(suppIndex)))
  }

  if ( input$allowZeros == FALSE ) {
    indZero <- which(g_freq(pI)==0)
    if ( length(indZero) > 0 ) {
      s_sdcStatus(pI) <- list(index=indZero, vals=rep("z", length(indZero)))
    }
  }
  s_problemInstance(object) <- pI
  validObject(object)
  return(object)
})

setMethod("c_rule_p", signature=c("sdcProblem", "list"), definition=function(object, input) {
  pPercRule <- function(celltot, cont1, cont2, p) {
    # if TRUE, cell needs to be suppressed
    (celltot - cont1 - cont2) < (p/100*cont1)
  }

  if ( !g_is_microdata(g_dataObj(object)) ) {
    stop("p-percent rule can only be applied if micro-data are available!\n")
  }

  pI <- g_problemInstance(object)
  dataObj <- g_dataObj(object)
  numVarInds <- g_numvar_ind(dataObj)
  strIDs <- g_strID(pI)

  numVal <- g_raw_data(dataObj)[[numVarInds[input$numVarInd]]]
  if ( any(numVal < 0 ) ) {
    stop("dominance rules can only be applied to numeric variables with only positive values!\n")
  }

  # calculate contributing indices
  indices <- lapply(1:g_nrVars(pI), function(x) {
    c_contributing_indices(object, input=list(strIDs[x]))
  })

  # values of contributing units
  valueList <- lapply(1:g_nrVars(pI), function(x) {
    rev(tail(sort(numVal[indices[[x]]]),2))
  })
  cellTotals <- g_numVars(pI)[[input$numVarInd]]

  # suppStatus: TRUE:unsafe, FALSE: safe
  pState <- sapply(1:g_nrVars(pI), function(x) {
    pPercRule(cellTotals[x], valueList[[x]][1], valueList[[x]][2], input$p)
  })

  suppIndex <- which(pState==TRUE)
  if ( length(suppIndex) > 0 ) {
    s_sdcStatus(pI) <- list(index=suppIndex, vals=rep("u", length(suppIndex)))
  }

  if ( input$allowZeros == FALSE ) {
    indZero <- which(g_freq(pI)==0)
    if ( length(indZero) > 0 ) {
      s_sdcStatus(pI) <- list(index=indZero, vals=rep("u", length(indZero)))
    }
  }
  s_problemInstance(object) <- pI
  validObject(object)
  return(object)
})

setMethod("c_rule_pq", signature=c("sdcProblem", "list"), definition=function(object, input) {
  pqRule <- function(celltot, cont1, cont2, p, q) {
    # if TRUE, cell needs to be suppressed
    (celltot - cont1 - cont2) < (p/q)*cont1
  }
  if ( !g_is_microdata(g_dataObj(object)) ) {
    stop("p-percent rule can only be applied if micro-data are available!\n")
  }

  pI <- g_problemInstance(object)
  dataObj <- g_dataObj(object)
  numVarInds <- g_numvar_ind(dataObj)
  strIDs <- g_strID(pI)

  numVal <- g_raw_data(dataObj)[[numVarInds[input$numVarInd]]]
  if ( any(numVal < 0 ) ) {
    stop("dominance rules can only be applied to numeric variables with only positive values!\n")
  }

  # calculate contributing indices
  indices <- lapply(1:g_nrVars(pI), function(x) {
    c_contributing_indices(object, input=list(strIDs[x]))
  })

  # values of contributing units
  valueList <- lapply(1:g_nrVars(pI), function(x) {
    rev(tail(sort(numVal[indices[[x]]]),2))
  })
  cellTotals <- g_numVars(pI)[[input$numVarInd]]

  # suppStatus: TRUE:unsafe, FALSE: safe
  pState <- sapply(1:g_nrVars(pI), function(x) {
    pqRule(cellTotals[x], valueList[[x]][1], valueList[[x]][2], input$pq[1], input$pq[2])
  })

  suppIndex <- which(pState==TRUE)
  if ( length(suppIndex) > 0 ) {
    s_sdcStatus(pI) <- list(index=suppIndex, vals=rep("u", length(suppIndex)))
  }

  if ( input$allowZeros == FALSE ) {
    indZero <- which(g_freq(pI)==0)
    if ( length(indZero) > 0 ) {
      s_sdcStatus(pI) <- list(index=indZero, vals=rep("u", length(indZero)))
    }
  }
  s_problemInstance(object) <- pI
  validObject(object)
  return(object)
})

setMethod("c_heuristic_solution", signature=c("sdcProblem", "list"), definition=function(object, input) {
  aProb <- input[[1]]
  validCuts <- input[[2]]
  solver <- input[[3]]$solver
  verbose <- input[[3]]$verbose

  ### create incremental attacker problem
  pI <- g_problemInstance(object)
  dimInfoObj <- g_dimInfo(object)
  primSupps <- g_primSupps(pI)
  secondSupps <- g_secondSupps(pI)
  forcedCells <- g_forcedCells(pI)
  nrVars <- g_nrVars(pI)
  weights <- ci <- g_weight(pI)
  ci[primSupps] <- 0
  ci[secondSupps] <- 0

  lb <- g_lb(pI)
  ub <- g_ub(pI)

  # required later in the cleanup-procedure
  LB <- LBdefault <- weights - lb
  UB <- UBdefault <- ub - weights

  m1 <- c_gen_mat_m(input=list(objectA=pI, objectB=dimInfoObj))
  m2 <- m1
  m2@v <- -1* m2@v
  AInc <- c_bind(object=m1, input=list(m2, bindRow=FALSE))

  nrConstraints <- nrow(AInc)
  objective <- rep(ci, 2)

  direction <- rep("==", g_nr_rows(AInc))
  rhs <- rep(0, g_nr_rows(AInc))

  types <- rep("C", g_nr_cols(AInc))
  boundsLower <- list(ind=1:g_nr_cols(AInc), val=rep(0, g_nr_cols(AInc)))
  boundsUpper <- list(ind=1:g_nr_cols(AInc), val=c(UB, LB))

  aProbInc <- new("linProb",
    objective=objective,
    constraints=AInc,
    direction=direction,
    rhs=rhs,
    boundsLower=boundsLower,
    boundsUpper=boundsUpper,
    types=types)

  # make sure that cells that must be published
  # are not part of the heuristic solution
  if ( length(forcedCells) > 0 ) {
    for ( u in 1:length(forcedCells) ) {
      con <- rep(0, g_nr_cols(AInc))
      con[c(forcedCells[u], nrVars+forcedCells[u])] <- c(1,-1)
      aCon <- init.cutList(type='singleCut', input=list(vals=con, dir="==", rhs=0))
      s_add_complete_constraint(aProbInc) <- list(aCon)
    }
  }

  x <- rep(0, g_nr_cols(AInc))
  UPL <- g_UPL(pI)
  LPL <- g_LPL(pI)
  SPL <- g_SPL(pI)

  SUP <- primSupps

  for ( i in seq_along(primSupps) ) {
    cellInd <- primSupps[i]
    if ( verbose ) {
      cat("finding additional cells to protect primSupp",i,"|",length(primSupps),"...\n")
    }
    con1 <- con2 <- x
    con1[cellInd] <- 1
    con2[nrVars+cellInd] <- 1
    con3 <- con1 - con2 # page 1018: fichetti and salazar!! (- is correct!)
    if ( UPL[cellInd] > 0 ) {
      # update and solve: y_ik_minus <- 0 and y_ik_plus <- UPL_ik
      aCon <- init.cutList(type='multipleCuts', input=list(mat=init.simpleTriplet(type='simpleTriplet', input=list(mat=rbind(con1, con2))), dir=rep("==", 2), rhs=c(UPL[cellInd], 0)))
      prob <- aProbInc
      s_add_complete_constraint(prob) <- list(aCon)
      sol <- c_solve_problem(prob, input=list(solver))$solution
      v <- sol[1:nrVars]+sol[(nrVars+1):length(sol)]
      v[which(is.zero(v))] <- 0
      addIndex <- which ( v > 0 )
      if ( length(addIndex) > 0 ) {
        SUP <- unique(c(SUP, addIndex))#
        ci[SUP] <- 0
        s_objective(aProbInc) <- list(rep(ci, 2))
        LB <- ci - lb
        UB <- ub - ci
      }
    }
    if ( LPL[cellInd] > 0 ) {
      # update and solve: y_ik_minus <- LPL_ik and y_ik_plus <- 0
      aCon <- init.cutList(type='multipleCuts', input=list(mat=init.simpleTriplet(type='simpleTriplet', input=list(mat=rbind(con1, con2))), dir=rep("==", 2), rhs=c(0, LPL[cellInd])))
      prob <- aProbInc
      s_add_complete_constraint(prob) <- list(aCon)
      sol <- c_solve_problem(prob, input=list(solver))$solution
      v <- sol[1:nrVars]+sol[(nrVars+1):length(sol)]
      v[which(is.zero(v))] <- 0
      addIndex <- which ( v > 0 )
      if ( length(addIndex) > 0 ) {
        SUP <- unique(c(SUP, addIndex))
        ci[SUP] <- 0
        s_objective(aProbInc) <- list(rep(ci, 2))
        LB <- ci - lb
        UB <- ub - ci
      }
    }
    if ( SPL[cellInd] > 0 ) {
      # update and solve: y_ik_plus + y_ik_minus <- SPL_ik
      aCon <- init.cutList(type='singleCut', input=list(vals=con3, dir="==", rhs=SPL[cellInd]))
      prob <- aProbInc
      s_add_complete_constraint(prob) <- list(aCon)
      sol <- c_solve_problem(prob, input=list(solver))$solution
      v <- sol[1:nrVars]+sol[(nrVars+1):length(sol)]
      v[which(is.zero(v))] <- 0
      addIndex <- which ( v > 0 )
      if ( length(addIndex) > 0 ) {
        SUP <- unique(c(SUP, addIndex))
        ci[SUP] <- 0
        s_objective(aProbInc) <- list(rep(ci, 2))
        LB <- ci - lb
        UB <- ub - ci
      }
    }
  }
  if ( verbose ) {
    cat(length(SUP) - length(primSupps),"additional cells have been suppressed in the heuristic solution!\n")
  }

  ### cleanup: remove redundant suppressions....
  # FIXME: use constraint pool to search for violations in the constraint pool
  # aProb has already been calculated and is an input parameter of this method!
  # aProb <- c_make_att_prob(input=list(objectA=pI, objectB=dimInfoObj))
  nrConstraints <- length(g_objective(aProb)) - 2*length(weights)
  addedSupps <- setdiff(SUP, primSupps)
  orderAddedSupps <- order(weights[addedSupps], decreasing=TRUE)
  xi <- rep(0, length(UPL))
  xi[SUP] <- 1 # we need to check xi

  counter <- 0
  for ( i in orderAddedSupps ) {
    counter <- counter + 1
    if ( verbose ) {
      cat("checking if removing cell",counter,"|",length(addedSupps),"still yields a valid suppression pattern...\n")
    }

    cellInd <- addedSupps[i]
    limits <- c(LPL[cellInd], UPL[cellInd], SPL[cellInd])
    xiWorking <- xi
    xiWorking[cellInd] <- 0 # we need to check if xi without xi[i] is a valid pattern
    UBWorking <- UB
    LBWorking <- LB
    UBWorking[cellInd] <- UBdefault[cellInd]
    LBWorking[cellInd] <- LBdefault[cellInd]

    ######################
    # check if any validCuts are not valid with xiWorking!
    # validCuts needs to be supplemented in function call
    ### get a constraint from validCuts
    if ( g_nr_constraints(validCuts) > 0 ) {
      conMat <- g_constraints(validCuts)
      result <- rep(NA, g_nr_rows(conMat))
      for ( z in 1:g_nr_rows(conMat) ) {
        expr <- paste(sum(xiWorking[g_col_ind(g_row(conMat, input=list(z)))]), g_direction(validCuts)[z], g_rhs(validCuts)[z])
        result[z] <- eval(parse(text=expr))
      }
    } else {
      result <- TRUE
    }

    if ( any(result==FALSE) ) {
      #cat("additionally suppressed cell cannot be removed (violated constraint in the pool found)!\n")
    } else {
      # no constraint was violated, we have to solve the incremental attacker problems
      # we need to solve the attackers problem for each sensitive cell twice
      if ( limits[3] != 0 ) {
        # solveAttackerProblem (upper bound)
        rhs <- rep(0, length(g_rhs(aProb)))
        rhs[cellInd] <- 1
        s_rhs(aProb) <- list(rhs)
        s_objective(aProb) <- list(c(weights + UBWorking*xiWorking, -(weights-xiWorking*LBWorking), rep(0, nrConstraints)))
        up <- c_solve_problem(aProb, input=list(solver))

        # solveAttackerProblem (lower bound)
        s_rhs(aProb) <- list(-1*rhs)
        down <- c_solve_problem(aProb, input=list(solver))

        calcDown <- -down$optimum
        calcUp <- up$optimum
      } else {
        # solve attackers problem (minimize)
        if ( limits[1] != 0 ) {
          rhs <- rep(0, length(g_rhs(aProb)))
          rhs[cellInd] <- -1
          s_rhs(aProb) <- list(rhs)
          s_objective(aProb) <- list(c(weights + UBWorking*xiWorking, -(weights-xiWorking*LBWorking), rep(0, nrConstraints)))
          down <- c_solve_problem(aProb, input=list(solver))
          calcDown <- -down$optimum
        }
        # solve attackers problem (maximize)
        if ( limits[2] != 0 ) {
          rhs <- rep(0, length(g_rhs(aProb)))
          rhs[cellInd] <- 1
          s_rhs(aProb) <- list(rhs)
          s_objective(aProb) <- list(c(weights + UBWorking*xiWorking, -(weights-xiWorking*LBWorking), rep(0, nrConstraints)))
          up <- c_solve_problem(aProb, input=list(solver))
          calcUp <- up$optimum
        }
      }

      # check for feasibility
      valid <- TRUE
      if ( limits[3] > 0 & calcUp - calcDown < SPL[i] ) {
        valid <- FALSE
      } else {
        if ( limits[1] > 0 & calcDown > weights[i] - LPL[i] ) {
          valid <- FALSE
        }
        if ( limits[2] > 0 & calcUp < weights[i] + UPL[i] ) {
          valid <- FALSE
        }
      }
      if ( valid ) {
        xi[cellInd] <- 0
        SUP <- setdiff(SUP, cellInd)
        if ( verbose ) {
          cat("redundant suppression found! --> removing cell!\n")
        }
      }
      # else: additionally suppressed cell cannot be removed!
    }
  }
  return(xi)
})

setMethod("c_anon_worker", signature=c("sdcProblem", "list"), definition=function(object, input) {
  timeLimit <- input$timeLimit
  verbose <- input$verbose
  save <- input$save

  if ( save == TRUE ) {
    files <- NULL
  }

  start.time <- proc.time()
  pI <- g_problemInstance(object)
  sdcStatusBegin <- g_sdcStatus(pI)
  primSupps <- primSuppsOrig <- g_primSupps(pI)

  indexPool <- numeric()
  allStrIDs <- g_strID(pI)

  if ( input$method == 'OPT' ) {
    s_elapsedTime(object) <- g_elapsedTime(object) + (proc.time()-start.time)[3]

    if ( input$useC == TRUE ) {
      result <- csp_cpp(sdcProblem=object, attackonly=FALSE, verbose=input$verbose)
    } else {
      result <- c_cut_and_branch(object, input)
    }

    if ( save==TRUE ) {
      fn <- paste(input$method,"_Object-Final.RData", sep="")
      save(object, file=fn)
    }
    return(result)
  }

  # HITAS or HYPERCUBE
  # check where we should start (saved file)
  partition <- g_partition(object)
  startI <- g_startI(object)
  startJ <- g_startJ(object)

  if ( startI !=1 | startJ != 1 ) {
    maxI <- partition$nrGroups
    if ( startJ < length(partition$indices[[startI]]) ) {
      startJ <- startJ+1
    } else {
      startJ <- 1
      startI <- startI+1
    }
  }

  if ( input$method == 'HITAS' ) {
    for ( i in startI:(partition$nrGroups) ) {
      s_startJ(object) <- 1 # reset j before updating i
      s_startI(object) <- i

      #indexPool <- g_indicesDealtWith(object)
      if ( i == 1 ) {
        indexPool <- c()
      } else {
        indexPool <- sort(unique(unlist(partition$indices[1:(i-1)])))
      }
      ind <- partition$indices[[i]]

      beginJ <- ifelse(i==startI, startJ, 1)
      for ( j in beginJ:(length(ind)) ) {
        s_startJ(object) <- j
        currentIndices <- ind[[j]] # within complete table

        ### cells with status 'u' or 'x' exist
        pI <- g_problemInstance(object)
        if ( any(g_sdcStatus(pI)[currentIndices] %in% c("u","x")) & length(currentIndices) > 1 ) {
          if ( verbose ) {
            cat("Starting to solve problem",j,"/",length(ind),"in group",i,"/",partition$nrGroups,"!\n")
          }
          ### if we have cells with "u" or "x" we need to protect
          ### the corresponding subtable
          ### reduce problemInstance
          probNew <- c_reduce_problem(object, input=list(currentIndices))
          pI.new <- g_problemInstance(probNew)

          ### is it necessary to protect the table??
          currentPrimSupps <- primSupps[!is.na(match(primSupps, currentIndices ))]

          ### indices that have already been inner cells in
          ### tables dealt earlier
          # FIXME: save indexpool somehow
          indicesDealtWith <- which(currentIndices %in% indexPool) #in current current subproblem

          ### fix marginal-cells
          ### --> its-suppression state must not change!
          currentPattern <- g_sdcStatus(g_problemInstance(probNew))

          introducedSupps <- indicesDealtWith[which(currentPattern[indicesDealtWith] == "x")]
          if ( length(introducedSupps) > 0 ) {
            ### secondary suppressions from upper tables
            s_sdcStatus(pI.new) <- list(index=introducedSupps,vals=rep("u", length(introducedSupps)))

            ### temporarily change LPL, UPL, SPL for these cells
            LPL.current <- g_LPL(pI.new)[introducedSupps]
            UPL.current <- g_UPL(pI.new)[introducedSupps]
            SPL.current <- g_SPL(pI.new)[introducedSupps]

            s_LPL(pI.new) <- list(index=introducedSupps, vals=rep(0, length(introducedSupps)))
            s_UPL(pI.new) <- list(index=introducedSupps, vals=rep(0, length(introducedSupps)))
            s_SPL(pI.new) <- list(index=introducedSupps, vals=rep(0.1, length(introducedSupps)))
          }

          ### force non-suppression of cells that have already been dealt with
          indForced <- indicesDealtWith[which(currentPattern[indicesDealtWith] == "s")]
          if ( length(indForced) > 0 ) {
            s_sdcStatus(pI.new) <- list(index=indForced,vals=rep("z", length(indForced)))
          }
          s_problemInstance(probNew) <- pI.new

          ### solving the problem
          if ( input$useC == TRUE ) {
            probNew <- csp_cpp(sdcProblem=probNew, attackonly=FALSE, verbose=input$verbose)
          } else {
            probNew <- c_cut_and_branch(object=probNew, input=input)
          }

          ### update sdcStatus
          status <- g_sdcStatus(g_problemInstance(probNew))

          pI <- g_problemInstance(object)
          if ( length(indForced) > 0 ) {
            status[indForced] <- "s"
          }
          if ( length(introducedSupps) > 0 ) {
            status[introducedSupps] <- "x"
            s_LPL(pI) <- list(index=currentIndices[introducedSupps], vals=LPL.current)
            s_UPL(pI) <- list(index=currentIndices[introducedSupps], vals=UPL.current)
            s_SPL(pI) <- list(index=currentIndices[introducedSupps], vals=SPL.current)
          }
          s_sdcStatus(pI) <- list(index=currentIndices, vals=status)
          s_problemInstance(object) <- pI
        }
        if ( save == TRUE ) {
          if ( verbose ) {
            cat("saving object after i=",i,"and j=",j,"\n")
          }
          fn <- paste(input$method,"_Object_",i,"-",j,".RData", sep="")
          files <- c(files, fn)
          save(object, file=fn)

          # removing old files
          if ( length(files) > 1 ) {
            sapply(rev(files)[-1], file.remove)
            files <- files[length(files)]
          }
        }
      }
      ### update indices that we have already dealt with
      s_indicesDealtWith(object) <- unique(c(indexPool, currentIndices))
    }
  }

  if ( input$method == 'HYPERCUBE' ) {
    runInd <- TRUE
    nrRuns <- 1
    while ( runInd == TRUE ) {
      cat("The algorithm is now starting run",nrRuns,"\n")

      tmpSupps <- c(g_primSupps(g_problemInstance(object)), g_secondSupps(g_problemInstance(object)))
      forcedCells <- g_forcedCells(g_problemInstance(object))

      for ( i in startI:(partition$nrGroups) ) {
        s_startJ(object) <- 1 # reset j before updating i
        s_startI(object) <- i

        ind <- partition$indices[[i]]

        beginJ <- ifelse(i==startI, startJ, 1)
        for ( j in beginJ:(length(ind)) ) {
          s_startJ(object) <- j

          currentIndices <- ind[[j]] # within complete table

          ### cells with status 'u' or 'x' exist
          pI <- g_problemInstance(object)
          # when using HYPERCUBE: we only check primary suppressions because we
          # temporarily set secondary suppressions to "u"
          if ( any(g_sdcStatus(pI)[currentIndices] == "u") & length(currentIndices) > 1 ) {
            if ( verbose ) {
              cat("Starting to solve problem",j,"/",length(ind),"in group",j,"/",partition$nrGroups,"!\n")
            }

            ### if we have cells with "u",  we need to protect
            ### the corresponding subtable
            ### reduce problemInstance
            probNew <- c_reduce_problem(object, input=list(currentIndices))
            pI.new <- g_problemInstance(probNew)

            ### is it necessary to protect the table??
            currentPrimSupps <- primSupps[!is.na(match(primSupps, currentIndices ))]

            s_problemInstance(probNew) <- pI.new

            ### solving the problem
            probNew <- c_ghmiter(object=probNew, input=input)

            ### update sdcStatus
            status <- g_sdcStatus(g_problemInstance(probNew))
            pI <- g_problemInstance(object)
            s_sdcStatus(pI) <- list(index=currentIndices, vals=status)
            s_problemInstance(object) <- pI
          }
          if ( save == TRUE ) {
            if ( verbose ) {
              cat("saving object after i=",i,"and j=",j,"\n")
            }
            fn <- paste(input$method,"_Object_",i,"-",j,".RData", sep="")
            files <- c(files, fn)
            save(object, file=fn)

            # removing old files
            if ( length(files) > 1 ) {
              sapply(rev(files)[-1], file.remove)
              files <- files[length(files)]
            }
          }
        }
      }

      ### protect secondary suppressions ###
      pI <- g_problemInstance(object)
      allSupps <- c(g_primSupps(pI), g_secondSupps(pI))
      newSupps <- setdiff(allSupps, tmpSupps)

      pI <- g_problemInstance(object)

      nrVars <- length(g_freq(pI))
      if ( length(newSupps) == 0 ) {
        runInd <- FALSE
        newSdcStatus <- rep('s', length=nrVars)
        newSdcStatus[forcedCells] <- 'z'
        newSdcStatus[tmpSupps] <- 'x'
        newSdcStatus[primSuppsOrig] <- 'u'
        s_sdcStatus(pI) <- list(index=1:nrVars, vals=newSdcStatus)
        s_problemInstance(object) <- pI
      } else {
        newSdcStatus <- rep('s', length=nrVars)
        newSdcStatus[forcedCells] <- 'z'
        newSdcStatus[tmpSupps] <- 'x'
        newSdcStatus[newSupps] <- 'u'
        s_sdcStatus(pI) <- list(index=1:nrVars, vals=newSdcStatus)
        s_problemInstance(object) <- pI
        nrRuns <- nrRuns + 1
      }
    } # while loop
  }
  return(object)
})

setMethod("c_opt_cpp", signature=c("sdcProblem", "list"), definition=function(object, input) {
  timeLimit <- input$timeLimit
  verbose <- input$verbose
  save <- input$save

  start.time <- proc.time()
  pI <- g_problemInstance(object)
  sdcStatusBegin <- g_sdcStatus(pI)
  primSupps <- primSuppsOrig <- g_primSupps(pI)

  indexPool <- numeric()
  allStrIDs <- g_strID(pI)

  s_elapsedTime(object) <- g_elapsedTime(object) + (proc.time()-start.time)[3]
  invisible(csp_cpp(sdcProblem=object, attackonly=FALSE, verbose=input$verbose))
})

setMethod("c_hitas_cpp", signature=c("sdcProblem", "list"), definition=function(object, input) {
  timeLimit <- input$timeLimit
  verbose <- input$verbose
  save <- input$save

  start.time <- proc.time()
  pI <- g_problemInstance(object)
  sdcStatusBegin <- g_sdcStatus(pI)
  primSupps <- primSuppsOrig <- g_primSupps(pI)

  indexPool <- numeric()
  allStrIDs <- g_strID(pI)

  partition <- g_partition(object)

  # save protection levels
  LPL.start <- g_LPL(pI)
  UPL.start <- g_UPL(pI)
  SPL.start <- g_SPL(pI)

  run_ind <- TRUE
  while ( run_ind ) {
    for ( i in 1:(partition$nrGroups) ) {
      s_startJ(object) <- 1 # reset j before updating i
      s_startI(object) <- i

      indexPool <- NULL
      if ( i > 1 ) {
        indexPool <- sort(unique(unlist(partition$indices[1:(i-1)])))
      }
      ind <- partition$indices[[i]]

      for ( j in 1:(length(ind)) ) {
        is_ok <- TRUE
        s_startJ(object) <- j
        currentIndices <- ind[[j]] # within complete table

        # cells with status "u" or "x" exist
        pI <- g_problemInstance(object)
        if ( any(g_sdcStatus(pI)[currentIndices] %in% c("u","x")) & length(currentIndices) > 1 ) {
          if ( verbose ) {
            cat("Starting to solve problem",j,"/",length(ind),"in group",i,"/",partition$nrGroups,"!\n")
          }
          # if we have cells with "u" or "x" we need to protect
          # the corresponding subtable --> reduce problemInstance
          probNew <- c_reduce_problem(object, input=list(currentIndices))
          pI.new <- g_problemInstance(probNew)

          # is it necessary to protect the table??
          currentPrimSupps <- primSupps[!is.na(match(primSupps, currentIndices))]

          # indices that have already been inner cells in tables dealt earlier
          indicesDealtWith <- which(currentIndices %in% indexPool) # in current current subproblem

          # we have to fix marginal-cells: --> their suppression status must not change!
          currentPattern <- g_sdcStatus(g_problemInstance(probNew))

          introducedSupps <- indicesDealtWith[which(currentPattern[indicesDealtWith] == "x")]
          if ( length(introducedSupps) > 0 ) {
            # secondary suppressions from upper tables
            s_sdcStatus(pI.new) <- list(index=introducedSupps, vals=rep("u", length(introducedSupps)))

            # temporarily change LPL, UPL, SPL for these cells
            LPL.current <- g_LPL(pI.new)[introducedSupps]
            UPL.current <- g_UPL(pI.new)[introducedSupps]
            SPL.current <- g_SPL(pI.new)[introducedSupps]

            s_LPL(pI.new) <- list(index=introducedSupps, vals=rep(0, length(introducedSupps)))
            s_UPL(pI.new) <- list(index=introducedSupps, vals=rep(0, length(introducedSupps)))
            s_SPL(pI.new) <- list(index=introducedSupps, vals=rep(0.1, length(introducedSupps)))
          }

          # force non-suppression of cells that have already been dealt with
          indForced <- indicesDealtWith[which(currentPattern[indicesDealtWith] == "s")]
          if ( length(indForced) > 0 ) {
            s_sdcStatus(pI.new) <- list(index=indForced, vals=rep("z", length(indForced)))
          }
          s_problemInstance(probNew) <- pI.new

          # solve the problem using c++ implementation
          res <- csp_cpp(sdcProblem=probNew, attackonly=FALSE, verbose=input$verbose)
          if ( is.null(res) ) {
            cat("\nWe got a problem and need to relax some conditions!\n\n")
            old.status <- probNew@problemInstance@sdcStatus
            ii <- which(old.status %in% c("z") & probNew@problemInstance@Freq > 0)
            if ( length(ii) == 0 ) {
              stop("This is a really nasty problem. No solution can be computed. Please contact the package maintainer.\n")
            }
            probNew@problemInstance@sdcStatus[ii] <- "s"
            res <- csp_cpp(sdcProblem=probNew, attackonly=FALSE, verbose=input$verbose)
            if ( is.null(res) ) {
              stop("This is a really nasty problem. No solution can be computed. Please contact the package maintainer.\n")
            }
            probNew <- res
            new.status <- probNew@problemInstance@sdcStatus
            xx <- data.frame(currentIndices=currentIndices,old=old.status, new=new.status, freq=probNew@problemInstance@Freq)

            ii <- which(new.status=="x" & old.status=="z")

            updated_status <- rep("s", length(old.status))
            updated_status[which(old.status=="u")] <- "u"
            updated_status[probNew@problemInstance@Freq==0] <- "z"
            updated_status[ii] <- "u"# previously "z", now "u"

            xx <- sdcStatusBegin
            xx[currentIndices] <- updated_status
            pI <- g_problemInstance(object)
            s_LPL(pI) <- list(index=currentIndices[ii], vals=rep(0, length(ii)))
            s_UPL(pI) <- list(index=currentIndices[ii], vals=rep(0, length(ii)))
            s_SPL(pI) <- list(index=currentIndices[ii], vals=rep(0.1, length(ii)))

            s_sdcStatus(pI) <- list(index=1:length(xx), vals=xx)
            s_problemInstance(object) <- pI
            is_ok <- FALSE
          } else {
            probNew <- res
          }
          # break j-loop
          if ( !is_ok ) {
            break
          }

          # update sdcStatus
          status <- g_sdcStatus(g_problemInstance(probNew))

          pI <- g_problemInstance(object)
          if ( length(indForced) > 0 ) {
            status[indForced] <- "s"
          }
          if ( length(introducedSupps) > 0 ) {
            status[introducedSupps] <- "x"
            s_LPL(pI) <- list(index=currentIndices[introducedSupps], vals=LPL.current)
            s_UPL(pI) <- list(index=currentIndices[introducedSupps], vals=UPL.current)
            s_SPL(pI) <- list(index=currentIndices[introducedSupps], vals=SPL.current)
          }
          s_sdcStatus(pI) <- list(index=currentIndices, vals=status)
          s_problemInstance(object) <- pI
        }
      }
      # break i-loop
      if ( !is_ok ) {
        break
      }
      # update indices that we have already dealt with
      s_indicesDealtWith(object) <- unique(c(indexPool, currentIndices))
    }
    if ( is_ok ) {
      run_ind <- FALSE
    }
  }

  pI <- g_problemInstance(object)
  s_LPL(pI) <- list(index=1:g_nrVars(pI), vals=LPL.start)
  s_UPL(pI) <- list(index=1:g_nrVars(pI), vals=UPL.start)
  s_SPL(pI) <- list(index=1:g_nrVars(pI), vals=SPL.start)
  sdcStatus <- g_sdcStatus(pI)

  ii <- which(sdcStatus %in% c("u", "x"))
  sdcStatus[ii] <- "x"
  sdcStatus[primSuppsOrig] <- "u"
  s_sdcStatus(pI) <- list(index=1:g_nrVars(pI), vals=sdcStatus)
  s_problemInstance(object) <- pI
  invisible(object)
})

setMethod("c_quick_suppression", signature=c("sdcProblem", "list"), definition=function(object, input) {
  freq <- id <- sdcStatus <- weights <- NULL
  verbose <- input$verbose
  pI <- g_problemInstance(object)
  indices <- g_partition(object)$indices
  dimInfo <- g_dimInfo(object)
  strInfo <- g_str_info(dimInfo)
  vNames <- g_varname(dimInfo)

  if ( verbose ) {
    cat("calculating subIndices (this may take a while) ...")
  }

  dat <- as.data.table(cpp_splitByIndices(g_strID(pI), strInfo))
  setnames(dat, vNames)
  dat[,id:=1:nrow(dat)]
  dat[,freq:=g_freq(pI)]
  dat[,weights:=g_weight(pI)]
  dat[,sdcStatus:=g_sdcStatus(pI)]
  dimVars <- match(vNames, names(dat))
  nDims <- length(dimVars)
  freqInd <- match("freq", colnames(dat))
  if ( length(vNames)==1 ) {
    combs <- combn(vNames, 1)
  } else {
    combs <- combn(vNames, length(vNames)-1)
  }

  tmpIndices <- rep(NA, length(vNames))

  nrGroups <- length(indices)
  subIndices <- list(); length(subIndices) <- nrGroups

  for ( group in 1:nrGroups ) {
    nrTabs <- length(indices[[group]])
    subIndices[[group]] <- list()
    length(subIndices[[group]]) <- nrTabs
    for (tab in 1:nrTabs) {
      subDat <- dat[indices[[group]][[tab]],]
      # only one dimension!
      if ( ncol(combs) == 1 ) {
        subDat$ind_1_tmp <- 1
        tmpIndices[1] <- ncol(subDat)
      } else {
        for ( i in 1:ncol(combs)) {
          setkeyv(subDat, combs[,i])
          cn <- paste0("ind_",i,"_tmp")
          expr <- parse(text = paste0(cn, ":=.GRP"))
          subDat[,eval(expr), by=key(subDat)]
          tmpIndices[i] <- ncol(subDat)
        }
      }
      setkeyv(subDat, vNames)
      subIndices[[group]][[tab]] <- as.list(subDat[,tmpIndices, with=F])
    }
  }
  if ( verbose ) {
    cat("[done]\n");
  }
  res <- greedyMultDimSuppression(dat,indices,subIndices,dimVars,verbose=verbose)
  if ( verbose ) {
    cat("finishing output...")
  }
  s_sdcStatus(pI) <- list(index=res$id, vals=res$sdcStatus)
  s_problemInstance(object) <- pI
  if ( verbose ) {
    cat("[done]\n")
  }
  invisible(list(object=object, zstatus=res$status_z))
})

setMethod("c_cut_and_branch", signature=c("sdcProblem", "list"), definition=function(object, input) {
  time.start <- proc.time()

  timeLimit <- input$timeLimit
  fixVariables <- input$fixVariables
  maxVars <- input$maxVars
  fastSolution <- input$fastSolution
  approxPerc <- input$approxPerc
  verbose <- input$verbose
  solver <- input$solver

  problemInstance <- g_problemInstance(object)
  dimInfo <- g_dimInfo(object)
  nrVars <- g_nrVars(problemInstance)
  freqs <- g_freq(problemInstance)
  primSupps <- g_primSupps(problemInstance)
  publishVars <- which(g_sdcStatus(problemInstance)=="z")
  noBranchVars <- unique(c(primSupps, publishVars))

  # Nothing to protect here
  if ( !g_hasPrimSupps(problemInstance) ) {
    return(object)
  }

  # returning heuristic solution
  # only if problem size is too large
  if ( is.null(maxVars) ) {
    maxVars <- nrVars + 1
  }
  if ( fastSolution ) {
    maxVars <- 0
  }

  approx <- ifelse(is.null(approxPerc), FALSE, TRUE)

  if ( nrVars >= maxVars ) {
    res <- c_make_att_prob(input=list(objectA=problemInstance, objectB=dimInfo))
    aProb <- res$aProb
    validCuts <- res$newCutsMaster

    heuristicSolution <- c_heuristic_solution(object, input=list(aProb, validCuts, input))

    secondSupps <- setdiff(which(heuristicSolution==1), primSupps)
    if (verbose) {
      cat("Result: we are returning a possibly non-optimal solution with",length(secondSupps),"secondary suppressions because of parameter 'fastSolution' or 'maxVars'!\n")
    }
    if ( length(secondSupps) > 0 ) {
      s_sdcStatus(problemInstance) <- list(index=secondSupps,vals=rep("x", length(secondSupps)))
    }
    out <- new("sdcProblem",
               dataObj=g_dataObj(object),
               dimInfo=dimInfo,
               problemInstance=problemInstance,
               partition=g_partition(object),
               startI=g_startI(object),
               startJ=g_startJ(object),
               indicesDealtWith=g_indicesDealtWith(object),
               elapsedTime=g_elapsedTime(object)+(proc.time()-time.start)[3]
    )
    return(out)
  }

  if ( verbose ) {
    cat("running pre-process procedure...\n")
  }

  resultPreProcess <- c_preprocess(object, input=input)

  object <- resultPreProcess$sdcProblem
  validCuts <- resultPreProcess$validCuts
  aProb <- resultPreProcess$aProb

  # no valid cuts have been generated in preprocessing!
  if ( g_nr_rows(g_constraints(validCuts)) == 0 ) {
    return(object)
  }

  if ( verbose ) {
    cat("calculating a heuristic solution...\n")
  }

  heuristicSolution <- c_heuristic_solution(object, input=list(aProb, validCuts, input))
  ### all solutions found and current best solution
  solutions <- list(); bestSolution <- heuristicSolution
  solutions[[1]] <- heuristicSolution

  startTime <- Sys.time()
  timeStop <- FALSE

  ### cuts due to hierarchical structure
  if ( verbose ) {
    cat("calculating structural cuts...\n")
  }

  structureCuts <- c_gen_structcuts(object, input=list())
  # does heuristicSolution violates any cuts from structure-cuts??
  #c_check_violation(structureCuts, input=list(heuristicSolution, g_weight(problemInstance)))
  validCuts <- c_bind_together(validCuts, input=list(structureCuts))
  #######

  ### create master problem and add constraints derived in pre-processing
  mProb <- c_make_masterproblem(problemInstance, input=list())
  s_add_complete_constraint(mProb) <- list(validCuts)
  if ( verbose ) {
    cat("solving the original master problem (no additional constraints)...\n")
  }
  masterSolution <- c_solve_problem(mProb, input=list(solver))
  masterObj <- masterSolution$optimum
  xi <- masterSolution$solution
  xi[is.zero(xi)] <- 0
  xi[is.one(xi)] <- 1

  ### initialize bounds
  currentBestBoundDown <- masterObj
  currentBestBoundUp <- sum(g_objective(mProb) * heuristicSolution)
  branchedVars <- NULL

  ### check if we have already the optimum solution (without rounding errors)
  runInd <- TRUE
  if ( abs(masterObj-currentBestBoundUp) < 0.1 ) {
    runInd <- FALSE
  } else {
    ### fixing variables
    if ( fixVariables == TRUE & currentBestBoundUp >= currentBestBoundDown ) {
      if ( verbose ) {
        cat("fixing variables...\n")
      }
      fixedVars <- c_fix_variables(mProb, input=list(currentBestBoundDown, currentBestBoundUp, primSupps))
      if (length(fixedVars) > 0 ) {
        if ( verbose ) {
          cat("--> setting",length(fixedVars),"variables to 0!\n")
        }
        bounds <- g_bounds(mProb)
        bounds$upper$val[fixedVars] <- 0
        s_bounds(mProb) <- bounds
      }
    }

    ### constraint pool initialization
    problemPool <- list()
    problemPool[[1]] <- init.cutList(type='empty', input=list(nrCols=nrVars))

    ### solving
    selectFirst <- TRUE
    LPL <- g_LPL(problemInstance)
    UPL <- g_UPL(problemInstance)
    SPL <- g_SPL(problemInstance)

    weights <- g_weight(problemInstance)
    LB <- weights - g_lb(problemInstance)
    UB <- g_ub(problemInstance) - weights
    nrConstraints <- length(g_objective(aProb)) - 2*length(weights)

    ### initialize constants (probably function-parameters later)
    selectFirst <- FALSE

    ### TODO: stop procedure after given time or nr or solutions..
    iter <- 0
  }

  while( runInd ) {
    iter <- iter + 1
    selectInd <- ifelse(selectFirst==TRUE, 1, length(problemPool))
    newCuts <- init.cutList(type='empty', input=list(nrCols=nrVars))
    AttProbDown <- AttProbUp <- rep(NA, length(primSupps))
    status <- NULL
    for ( i in 1:length(primSupps) ) {
      cellInd <- primSupps[i]
      limits <- c(LPL[cellInd], UPL[cellInd], SPL[cellInd])

      # we need to solve the attackers problem for each sensitive cell twice
      if ( limits[3] != 0 ) {
        # solveAttackerProblem (upper bound)
        rhs <- rep(0, length(g_rhs(aProb)))
        rhs[cellInd] <- 1
        s_rhs(aProb) <- list(rhs)
        s_objective(aProb) <- list(c(weights + UB*xi, -(weights-xi*LB), rep(0, nrConstraints)))
        up <- c_solve_problem(aProb, input=list(solver))

        ### solveAttackerProblem (lower bound)
        s_rhs(aProb) <- list(-1*rhs)
        down <- c_solve_problem(aProb, input=list(solver))

        AttProbDown[i] <- -down$optimum
        AttProbUp[i] <- up$optimum
        status <- c(status, down$status, up$status)
        #cat('limits ( origValue=',weights[cellInd],') : [',AttProbDown[i],':',AttProbUp[i],']\n')

        alpha.down <- down$solution[1:nrVars]
        alpha.up <- up$solution[1:nrVars]

        beta.down <- down$solution[(nrVars+1):(2*nrVars)]
        beta.up <- up$solution[(nrVars+1):(2*nrVars)]
      } else {
        # solve attackers problem (minimize)
        if ( limits[1] != 0 ) {
          rhs <- rep(0, length(g_rhs(aProb)))
          rhs[cellInd] <- -1
          s_rhs(aProb) <- list(rhs)
          down <- c_solve_problem(aProb, input=list(solver))
          AttProbDown[i] <- -down$optimum
        }
        # solve attackers problem (maximize)
        if ( limits[2] != 0 ) {
          rhs <- rep(0, length(g_rhs(aProb)))
          rhs[cellInd] <- 1
          s_rhs(aProb) <- list(rhs)
          s_objective(aProb) <- list(c(weights + UB*xi, -(weights-xi*LB), rep(0, nrConstraints)))
          up <- c_solve_problem(aProb, input=list(solver))
          AttProbUp[i] <- up$optimum
        }
      }
      # SPL
      if ( limits[3] != 0 & AttProbUp[i] - AttProbDown[i] < limits[3] ) {
        status <- c(status, down$status, up$status)
        alpha.down <- down$solution[1:nrVars]
        alpha.up <- up$solution[1:nrVars]
        beta.down <- down$solution[(nrVars+1):(2*nrVars)]
        beta.up <- up$solution[(nrVars+1):(2*nrVars)]

        v <- (alpha.down+alpha.up)*UB + (beta.down+beta.up)*LB
        v[which(is.zero(v))] <- 0
        if ( any(v != 0) )
          s_add_complete_constraint(newCuts) <- list(init.cutList(type='singleCut', input=list(vals=v, dir=">=", rhs=limits[3])))
      } else  {
        if ( limits[1] != 0 & freqs[primSupps[i]] - AttProbDown[i] < limits[1] ) { # LPL
          status <- c(status, down$status)
          alpha.down <- down$solution[1:nrVars]
          beta.down <- down$solution[(nrVars+1):(2*nrVars)]

          v <- alpha.down*UB + beta.down*LB
          v[which(is.zero(v))] <- 0
          if ( any(v != 0) )
            s_add_complete_constraint(newCuts) <- list(init.cutList(type='singleCut', input=list(vals=v, dir=">=", rhs=limits[1])))
        }
        if ( limits[2] != 0 & AttProbUp[i] - freqs[primSupps[i]] < limits[2] ) { # UPL
          status <- c(status, up$status)
          alpha.up <- up$solution[1:nrVars]
          beta.up <- up$solution[(nrVars+1):(2*nrVars)]

          v <- alpha.up*UB + beta.up*LB
          v[which(is.zero(v))] <- 0
          if ( any(v != 0) )
            s_add_complete_constraint(newCuts) <- list(init.cutList(type='singleCut', input=list(vals=v, dir=">=", rhs=limits[2])))
        }
      }
      #cat('limits ( origValue=',weights[cellInd],') : [',AttProbDown[i],':',AttProbUp[i],']\n')
    }

    if ( g_nr_constraints(newCuts) > 0 ) {
      # strengthen cuts
      if ( verbose ) {
        cat("strengthening the cuts and adding",g_nr_constraints(newCuts),"new derived cuts to master problem...\n")
      }
      newCuts <- c_strengthen(newCuts)
      s_add_complete_constraint(mProb) <- list(newCuts)
    }
    ### check for duplicated constraints
    indRem <- which(duplicated(cbind(as.matrix(g_constraints(mProb)), g_direction(mProb), g_rhs(mProb))))
    if ( length(indRem) > 0 ) {
      #if ( verbose ) {
      # cat("removing",length(indRem),"duplicated constraints...\n")
      #}
      s_remove_complete_constraint(mProb) <- list(indRem)
    }
    ### bridgeless inequalities only at root-node ####
    bridgelessCandidates <- setdiff(which(xi == 1), primSupps)
    if ( iter == 1 & length(bridgelessCandidates) > 0 ) {
      brCuts <- init.cutList(type='empty', input=list(nrCols=nrVars))
      if ( verbose ) {
        cat("adding",length(bridgelessCandidates),"bridgeless ineqalities!\n")
      }
      for ( z in seq_along(bridgelessCandidates) ) {
        bridgelessInd <- bridgelessCandidates[z]
        ### solveAttackerProblem (upper bound)
        rhs <- rep(0, length(g_rhs(aProb)))
        rhs[bridgelessInd] <- 1
        s_rhs(aProb) <- list(rhs)
        s_objective(aProb) <- list(c(weights + UB*xi, -(weights-xi*LB), rep(0, nrConstraints)))
        up <- c_solve_problem(aProb, input=list(solver))

        ### solveAttackerProblem (lower bound)
        s_rhs(aProb) <- list(-1*rhs)
        down <- c_solve_problem(aProb, input=list(solver))

        alpha.down <- down$solution[1:nrVars]
        alpha.up <- up$solution[1:nrVars]

        beta.down <- down$solution[(nrVars+1):(2*nrVars)]
        beta.up <- up$solution[(nrVars+1):(2*nrVars)]

        brIneq <- (alpha.down+alpha.up)*UB + (beta.down+beta.up)*LB
        brIneq[is.zero(brIneq)] <- 0
        brIneq[brIneq > 0] <- 1
        brIneq[bridgelessInd] <- -1
        s_add_complete_constraint(brCuts) <- list(init.cutList(type='singleCut', input=list(vals=brIneq, dir=">=", rhs=0)))

      }
      if ( g_nr_constraints(brCuts) > 0 ) {
        s_add_complete_constraint(mProb) <- list(brCuts)
      }
    }

    mProbWorking <- mProb
    # eventually update the lower bound...
    tmpSolution <- c_solve_problem(mProbWorking, input=list(solver))
    tmpObj <- tmpSolution$optimum
    if ( tmpObj > currentBestBoundDown & tmpObj <= currentBestBoundUp ) {
      currentBestBoundDown <- tmpObj
    }
    if ( abs(currentBestBoundUp - currentBestBoundDown) < 1 ) {
      # optimal solution found!
      break
    }

    if ( g_nr_constraints(problemPool[[selectInd]]) > 0 ) {
      if ( verbose ) {
        cat("adding",g_nr_constraints(problemPool[[selectInd]]),"constraints to the master problem...\n")
      }
      mProbWorking <- mProb
      s_add_complete_constraint(mProbWorking) <- list(problemPool[[selectInd]])
    }

    if ( verbose ) {
      cat("solving the master problem with",length(g_rhs(mProbWorking)),"constraints...\n")
    }
    masterSolution <- c_solve_problem(mProbWorking, input=list(solver))
    masterObj <- masterSolution$optimum
    xi <- masterSolution$solution
    xi[is.zero(xi)] <- 0
    xi[is.one(xi)] <- 1

    if ( verbose ) {
      cat("best-bounds: [",currentBestBoundDown,":",currentBestBoundUp,"] and objVal =",masterObj,"with sum(xi)=",sum(xi),"\n")
    }
    #cat("current boundUp =",currentBestBoundUp,"and objVal =",masterObj,"with sum(xi)=",sum(xi),"\n")

    ### again fixing variables
    #if ( fixVariables == TRUE ) {
    # newFixedVars <- c_fix_variables(mProb, input=list(currentBestBoundDown, currentBestBoundUp, primSupps))
    # if ( !all(newFixedVars) %in% fixedVars ) {
    #   cat("setting",length(newFixedVars),"variables to 0!\n")
    #   bounds <- g_bounds(mProb)
    #   bounds$upper$val[newFixedVars] <- 0
    #   s_bounds(mProb) <- bounds
    # }
    #}

    ### checking if we can prune the current node
    prune <- FALSE
    pruneReason <- NULL
    # a) valid (protected) integer solution
    if ( all(is.wholenumber(xi)) && c_is_protected_solution(problemInstance, input=list(input1=AttProbDown, input2=AttProbUp)) ) {
      prune <- TRUE
      pruneReason <- c(pruneReason, "V") # valid
    }
    # b) infeasibility
    if ( sum(status) != 0 ) {
      prune <- TRUE
      pruneReason <- c(pruneReason, "I") # infeasible
    }
    # c) bounds
    if ( approx == TRUE ) {
      if ( currentBestBoundUp - masterObj <= currentBestBoundUp * (approxPerc/100) ) {
        prune <- TRUE
        pruneReason <- c(pruneReason, "B") # bounds
      }
      if ( masterObj - currentBestBoundDown < 0 ) {
        prune <- TRUE
        pruneReason <- c(pruneReason, "B") # bounds
      }

    } else {
      if ( abs(masterObj - currentBestBoundUp) <= 0.01 ) {
        prune <- TRUE
        pruneReason <- c(pruneReason, "B") # bounds
      }
      if ( masterObj - currentBestBoundDown < 0 ) {
        prune <- TRUE
        pruneReason <- c(pruneReason, "B") # bounds
      }
    }

    if ( prune == TRUE ) {
      pruneReason <- unique(pruneReason) # remove eventually 2-'Bs'
      if ( length(pruneReason) == 2 ) {
        if ( pruneReason[1] == "V" & pruneReason[2]=="B" ) {
          #cat("found worse-than optimal integer-solution -> pruning by bounds!\n")
          if ( masterObj < currentBestBoundUp ) {
            pruneReason <- "V"
          } else {
            pruneReason <- "B"
          }
        }
      }
      if ( length(pruneReason) > 1 ) {
        stop("Error: only one pruning reason possible!\n")
      }
      if ( pruneReason == "V") {
        solutions[[length(solutions)+1]] <- as.integer(xi)
        if ( masterObj < currentBestBoundUp ) {
          if ( verbose ) {
            cat("new best integer solution (objval=",masterObj,") found!:\n")
          }
          currentBestBoundUp <- masterObj
          bestSolution <- as.integer(xi)
        }
      }
      #if ( pruneReason == "I") {
      # cat("pruning because of infeasibility!\n")
      #}
      #if ( pruneReason == "B") {
      # cat("pruning because of known bounds!\n")
      #}
      problemPool[[selectInd]] <- NULL
      if ( verbose ) {
        cat("pruning the current node: reason=",pruneReason,"!. Still",length(problemPool),"nodes in the pool!\n")
      }
    } else {
      ## 2) Branching: we extend the problemPool and remove the current node afterwards
      branchedVars <- g_col_ind(g_constraints(problemPool[[selectInd]]))
      branchVar <- getBranchingVariable(xi, branchedVars, noBranchVars)

      if ( length(branchVar) == 1 ) {
        cl <- problemPool[[selectInd]]
        v <- rep(0, nrVars)
        v[branchVar] <- 1
        c1 <- c2 <- cl
        s_add_complete_constraint(c1) <- list(init.cutList(type='singleCut', input=list(vals=v, dir="==", rhs=0)))
        s_add_complete_constraint(c2) <- list(init.cutList(type='singleCut', input=list(vals=v, dir="==", rhs=1)))

        problemPool[[length(problemPool)+1]] <- c1
        problemPool[[length(problemPool)+1]] <- c2; rm(cl)

        # now we can prune the current node
        problemPool[[selectInd]] <- NULL
        rm(c1,c2)
        if ( verbose ) {
          cat("branching was required. Problem pool has now",length(problemPool),"nodes!\n")
        }
      } else {
        if ( verbose ) {
          cat("no further branching possible! all branching variables tried!\n")
        }
        problemPool[[selectInd]] <- NULL
      }
    }

    timeSpent <- as.numeric(floor(difftime(Sys.time(), startTime, units="mins")))
    #cat("timeSpent:"); print(timeSpent)

    if ( length(problemPool)==0 ) {
      runInd <- FALSE
    } else {
      if ( !is.null(timeLimit) && timeSpent > timeLimit && length(solutions) > 0 ) {
        runInd <- FALSE
        timeStop <- TRUE
      }
      if ( !is.null(timeLimit) && timeSpent > timeLimit && length(solutions) == 0 ) {
        if ( verbose ) {
          cat("Result: the time-limit was reached and no (heuristic) solution could be generated!\n")
        }
        return(object)
      }
    }
  }

  secondSupps <- setdiff(which(bestSolution==1), primSupps)
  objVarHeuristic <- sum(g_objective(mProb) * heuristicSolution)
  if ( timeStop==FALSE ) {
    if ( currentBestBoundUp == objVarHeuristic ) {
      if ( verbose ) {
        cat('Result: the heuristic solution was already optimal and has',length(secondSupps),'secondary suppressions!\n')
      }
    } else {
      improvement <- 100 - (100 / objVarHeuristic ) * currentBestBoundUp
      if ( verbose ) {
        cat('Result: the heuristic solution was improved by',format(improvement, digits=2, nsmall=2),'% and has',length(secondSupps),'secondary suppressions!!\n ')
      }
    }
  } else {
    if ( verbose ) {
      cat("Result: we are returning a possibly non-optimal solution with",length(secondSupps),"secondary suppressions because of argument 'timeLimit'!\n")
    }
  }
  if ( length(secondSupps) > 0 ) {
    s_sdcStatus(problemInstance) <- list(index=secondSupps,vals=rep("x", length(secondSupps)))
  }
  out <- new("sdcProblem",
    dataObj=g_dataObj(object),
    dimInfo=dimInfo,
    problemInstance=problemInstance,
    partition=g_partition(object),
    startI=g_startI(object),
    startJ=g_startJ(object),
    indicesDealtWith=g_indicesDealtWith(object),
    elapsedTime=g_elapsedTime(object)+(proc.time()-time.start)[3]
  )
  return(out)
})

setMethod("c_ghmiter", signature=c("sdcProblem", "list"), definition=function(object, input) {
  time.start <- proc.time()
  protectionLevel <- input$protectionLevel
  suppMethod <- input$suppMethod
  suppAdditionalQuader <- input$suppAdditionalQuader
  verbose <- input$verbose

  pI <- g_problemInstance(object)
  strIDs <- g_strID(pI)
  strInfo <- g_str_info(g_dimInfo(object))

  sdcStatus <- g_sdcStatus(pI)
  cellsToProtect <- g_primSupps(pI)

  freqs <- g_freq(pI)

  # calc infomation on inner|marginal cells
  cellInfo <- g_innerAndMarginalCellInfo(object)

  # replaces f.recodeIndexVars
  indexList <- lapply(1:length(strInfo), function(x) { mySplit(strIDs, strInfo[[x]][1]:strInfo[[x]][2]) } )

  for ( i in 1:length(cellsToProtect) ) {
    if ( verbose ) {
      cat("--> Cell",i,"|",length(cellsToProtect)," (ID:",strIDs[cellsToProtect[i]],")...")
    }

    diagObj <- c_ghmiter_diag_obj(object, input=list(cellsToProtect[i], indexList, FALSE))
    # calculate required information using diagObj
    infoObj <- c_ghmiter_calc_info(object, input=list(diagObj, indexList, protectionLevel, FALSE))

    if ( !is.null(infoObj) & length(infoObj) == 0 ) {
      diagObj <- c_ghmiter_diag_obj(object, input=list(cellsToProtect[i], indexList, TRUE))

      infoObj <- c_ghmiter_calc_info(object, list(diagObj, indexList, protectionLevel, TRUE))

      if ( !is.null(infoObj) & length(infoObj) == 0 ) {
        stop("ghmiter::: error - could not find sensible cube!\n")
      }
      warning("Cell with Freq=0 has been selected as partner in suppression pattern!\n")
    }

    # cellToProtect used from diagObj$cellToProtect
    suppObj <- c_ghmiter_select_quader(object, input=list(infoObj, input))

    if ( !is.null(suppObj) ) {
      object <- c_ghmiter_suppress_quader(object, input=suppObj)

      # additional quader needs to be found
      # only if it is not a single value in the margins
      # and the cube includes cells with frequency=1
      if ( suppAdditionalQuader==TRUE & suppObj$indikatorSingleItems==TRUE & !(cellsToProtect[i] %in% cellInfo$indexTotCells) ) {
        # find additional cube that does not contain the single cells
        object <- c_ghmiter_supp_additional(object, input=list(diagObj, infoObj, suppObj, input))
      }
    }
    if ( verbose ) {
      cat("[done]\n")
    }
  }
  s_elapsedTime(object) <- g_elapsedTime(object)+(proc.time()-time.start)[3]
  return(object)
})

setMethod("c_preprocess", signature=c("sdcProblem", "list"), definition=function(object, input) {
  solver <- input$solver
  problemInstance <- g_problemInstance(object)
  if ( !g_hasPrimSupps(problemInstance) ) {
    return(object)
  }
  dimInfo <- g_dimInfo(object)
  nrVars <- g_nrVars(problemInstance)
  freqs <- g_freq(problemInstance)
  primSupps <- g_primSupps(problemInstance)

  LPL <- g_LPL(problemInstance)[primSupps]
  UPL <- g_UPL(problemInstance)[primSupps]
  SPL <- g_SPL(problemInstance)[primSupps]

  weights <- g_weight(problemInstance)
  HIGH <- LOW <- weights[primSupps]

  # order of calculations
  myOrder <- order(sapply(1:length(primSupps), function(x) { max(SPL[x], (LPL+UPL)[x]) }), decreasing=TRUE)
  LB <- weights - g_lb(problemInstance)
  UB <- g_ub(problemInstance) - weights

  xi <- g_suppPattern(problemInstance)
  res <- c_make_att_prob(input=list(objectA=problemInstance, objectB=dimInfo))
  aProb <- res$aProb
  validCuts <- res$newCutsMaster

  nrConstraints <- length(g_objective(aProb)) - 2*length(weights)

  #validCuts <- init.cutList(type='empty', input=list(nrCols=nrVars))
  for ( i in myOrder ) {
    if ( i %% 10 == 1 && input$verbose ) {
      cat("preprocessing variable",i,"|",length(myOrder),"...\n")
    }
    cellInd <- primSupps[i]
    limits <- c(LPL[i], UPL[i], SPL[i])

    # solveAttackerProblem (upper bound)
    rhs <- rep(0, length(g_rhs(aProb)))
    rhs[cellInd] <- 1
    s_rhs(aProb) <- list(rhs)
    s_objective(aProb) <- list(c(weights + UB*xi, -(weights-xi*LB), rep(0, nrConstraints)))
    up <- c_solve_problem(aProb, input=list(solver))
    if ( up$status != 0 ) {
      stop("unsolvable problem (up)!\n")
    }
    calcUp <- up$optimum
    HIGH[i] <- max(HIGH[i], calcUp)
    LOW[i] <- min(LOW[i], calcUp)

    # solveAttackerProblem (lower bound)
    s_rhs(aProb) <- list(-1*rhs)
    down <- c_solve_problem(aProb, input=list(solver))
    if ( down$status != 0 ) {
      stop("unsolvable problem (down)!\n")
    }
    calcDown <- -down$optimum
    HIGH[i] <- max(HIGH[i], calcDown)
    LOW[i] <- min(LOW[i], calcDown)

    alpha.down <- down$solution[1:nrVars]
    alpha.up <- up$solution[1:nrVars]

    beta.down <- down$solution[(nrVars+1):(2*nrVars)]
    beta.up <- up$solution[(nrVars+1):(2*nrVars)]

    if ( limits[1] != 0 ) { # LPL
      v <- alpha.down*UB + beta.down*LB
      v[which(is.zero(v))] <- 0
      if ( any(v > 0) )
        s_add_complete_constraint(validCuts) <- list(init.cutList(type='singleCut', input=list(vals=v, dir=">=", rhs=limits[1])))
    }
    if ( limits[2] != 0 ) { # UPL
      v <- alpha.up*UB + beta.up*LB
      v[which(is.zero(v))] <- 0
      if ( any(v > 0) )
        s_add_complete_constraint(validCuts) <- list(init.cutList(type='singleCut', input=list(vals=v, dir=">=", rhs=limits[2])))
    }
    if ( limits[3] != 0 & calcUp - calcDown < limits[3] ) { # SPL
      v <- (alpha.down+alpha.up)*UB + (beta.down+beta.up)*LB
      v[which(is.zero(v))] <- 0
      if ( any(v > 0) )
      s_add_complete_constraint(validCuts) <- list(init.cutList(type='singleCut', input=list(vals=v, dir=">=", rhs=limits[3])))
    }
  }

  if ( g_nr_constraints(validCuts) > 0 ) {
    validCuts <- c_strengthen(validCuts)

    setZeroUPL <- which(freqs[primSupps]+UPL <= HIGH) # -> set UPL = 0
    setZeroLPL <- which(freqs[primSupps]-LPL >= LOW) # -> set LPL = 0
    setZeroSPL <- which(HIGH-LOW >= SPL ) # -> set SPL = 0

    if ( length(setZeroUPL) > 0 ) {
      UPL[setZeroUPL] <- 0
      s_UPL(problemInstance) <- list(index=primSupps, vals=UPL)
    }
    if ( length(setZeroLPL) > 0 ) {
      LPL[setZeroLPL] <- 0
      s_LPL(problemInstance) <- list(index=primSupps, vals=LPL)
    }
    if ( length(setZeroSPL) > 0 ) {
      SPL[setZeroSPL] <- 0
      s_SPL(problemInstance) <- list(index=primSupps, vals=SPL)
    }
  }
  s_problemInstance(object) <- problemInstance
  return(list(sdcProblem=object, aProb=aProb, validCuts=validCuts))
})

setMethod("c_cellID", signature=c("sdcProblem", "list"), definition=function(object, input) {
  para.names <- input[[1]]
  para.codes <- input[[2]]
  para.verbose <- input[[3]]

  pI <- g_problemInstance(object)
  dimInfoObj <- g_dimInfo(object)

  vNames <- g_varname(dimInfoObj)
  vIndex <- g_pos_index(dimInfoObj)

  indexVar <- match(para.names, vNames)
  strInfo <- g_str_info(dimInfoObj)
  dimInfo <- g_dim_info(dimInfoObj)

  # calculate original codes
  codesDefault <- lapply(1:length(strInfo), function(x) {
    mySplit(g_strID(pI), strInfo[[x]][1]:strInfo[[x]][2])
  })
  codesOrig <- list()
  for ( i in 1:length(codesDefault) ) {
    codesOrig[[i]] <- c_match_orig_codes(object=dimInfo[[i]], input=codesDefault[[i]])
  }

  if ( length(input) != 3 ) {
    stop("c_cellID:: length of argument 'input' must equal 3!\n")
  }
  if ( length(para.names) != length(para.codes) ) {
    stop("c_cellID:: check argument 'input'!\n")
  }
  if ( !all(para.names %in% vNames) ) {
    stop("c_cellID:: check variable names in 'input[[1]]'!\n")
  }
  if ( !is.logical(para.verbose) ) {
    stop("c_cellID:: argument in 'input[[3]]' must be logical!\n")
  }

  cellID <- 1:g_nrVars(pI)
  for ( i in seq_along(para.names) ) {
    cellID <- intersect(cellID, which(!is.na(match(as.character(codesOrig[[indexVar[i]]]), para.codes[i]))))
  }
  if ( length(cellID) != 1) {
    stop("c_cellID:: check argument 'input' -> 0 or > 1 cells identified!\n")
  }
  return(cellID)
})

setMethod("c_finalize", signature=c("sdcProblem", "list"), definition=function(object, input) {
  Freq <- NULL
  time.start <- proc.time()

  pI <- g_problemInstance(object)
  dI <- g_dimInfo(object)
  levelObj <- g_dim_info(dI)
  strInfo <- g_str_info(dI)

  sdcStatus <- g_sdcStatus(pI)

  nrNonDuplicatedCells <- g_nrVars(pI)
  nrPrimSupps <- length(which(sdcStatus == 'u'))
  nrSecondSupps <- length(which(sdcStatus == 'x'))
  nrPublishableCells <- length(which(sdcStatus %in% c('z','s')))

  # merge codes and labels
  codesOriginal <- list()
  strIDs <- g_strID(pI)
  for ( i in seq_along(levelObj) ) {
    codesDefault <- mySplit(strIDs, strInfo[[i]][1]:strInfo[[i]][2])
    codesOriginal[[i]] <- c_match_orig_codes(object=levelObj[[i]], input=codesDefault)
  }
  out <- as.data.table(codesOriginal)
  setnames(out, g_varname(dI))
  out[,Freq:=g_freq(pI)]

  numVars <- g_numVars(pI)
  if ( !is.null(numVars) ) {
    data.obj <- g_dataObj(object)
    nV <- as.data.table(numVars)
    setnames(nV, colnames(g_raw_data(data.obj))[g_numvar_ind(data.obj)])
    out <- cbind(out, nV)
  }
  out[,sdcStatus:=g_sdcStatus(pI)]
  out[sdcStatus=="z", sdcStatus:="s"]

  # add duplicates
  hasDups <- sapply(1:length(levelObj), function(x) {
    g_has_dups(levelObj[[x]])
  })
  if ( any(hasDups) ) {
    for ( i in which(hasDups==TRUE) ) {
      dups <- g_dups(levelObj[[i]])
      dupsUp <- g_dups_up(levelObj[[i]])

      runInd <- TRUE
      while ( runInd ) {
        for ( j in 1:length(dups) ) {
          sub <- subset(out, out[[i]]==dupsUp[j])
          if ( nrow(sub) > 0 ) {
            sub[[i]] <- dups[j]
            out <- rbind(out, sub)
            dups <- dups[-j]
            dupsUp <- dupsUp[-j]
          }
        }
        if ( length(dups) == 0 )
          runInd <- FALSE
      }
    }
  }
  safeObj <- new("safeObj",
    finalData=out,
    dimInfo=dI,
    nrNonDuplicatedCells=nrNonDuplicatedCells,
    nrPrimSupps=nrPrimSupps,
    nrSecondSupps=nrSecondSupps,
    nrPublishableCells=nrPublishableCells,
    suppMethod=input$method,
    elapsedTime=g_elapsedTime(object) + (proc.time()-time.start)[3]
  )
  return(safeObj)
})

setMethod("c_ghmiter_diag_obj", signature=c("sdcProblem", "list"), definition=function(object, input) {
  cellToProtect <- input[[1]]
  indexList <- input[[2]]
  allowZeros <- input[[3]]

  if ( length(cellToProtect) != 1 ) {
    stop("makeDiagObj:: length of 'cellToProtect' must equal 1!\n")
  }
  pI <- g_problemInstance(object)
  nrVars <- g_nrVars(pI)
  if ( !cellToProtect %in% 1:nrVars ) {
    stop("makeDiagObj:: 'cellToProtect' must be >= 1 and <= ",nrVars,"!\n")
  }
  if ( !cellToProtect %in% 1:nrVars ) {
    stop("makeDiagObj:: 'cellToProtect' must be >= 1 and <= ",nrVars,"!\n")
  }
  sdcStatus <- g_sdcStatus(pI)

  if ( allowZeros == TRUE ) {
    indZeros <- which(sdcStatus=="z" & g_freq(pI)==0)
    if ( length(indZeros ) > 0 ) {
      sdcStatus[indZeros] <- "s"
    }
  }

  indToProtect <- sapply(indexList, function(x) { x[cellToProtect] } )
  diagIndices <- NULL
  for ( i in seq_along(indexList) ) {
    if ( length(unique(indexList[[i]])) > 1  ) {
      if ( is.null(diagIndices) ) {
        diagIndices <- which(indexList[[i]] != indexList[[i]][cellToProtect])
      } else {
        diagIndices <- intersect(diagIndices, which(indexList[[i]] != indexList[[i]][cellToProtect]))
      }
    }
  }
  diagIndices <- diagIndices[which(sdcStatus[diagIndices]!="z")]

  if ( length(diagIndices) == 0 )
    diagIndices <- NULL

  supp <- list()
  supp$cellToProtect <- cellToProtect
  supp$indToProtect <- indToProtect
  supp$diagIndices <- diagIndices
  return(supp)
})

setMethod("c_ghmiter_calc_info", signature=c("sdcProblem", "list"), definition=function(object, input) {
  # calculate some info on a given quader (normalization,...)
  calcQInfo <- function(g, d) {
    if ( length(g) != length(d) ) {
      stop("calcQInfo: 'g' and 'd' must have equal length!\n")
    }
    numberIndexVars <- length(g)

    quader <- expand(lapply(1:numberIndexVars, function(x) { c(g[x],d[x]) }  ), vector=FALSE)

    quader <- matrix(unlist(quader), length(quader[[1]]), length(quader))
    quader <- quader[!duplicated(quader),,drop=FALSE]
    quader <- lapply(1:ncol(quader), function(x) quader[,x])

    ### normquader
    normQ <- list()
    for ( i in seq_along(quader) ) {
      normQ[[i]] <- rep(0, length(quader[[1]]))
      normQ[[i]][which(quader[[i]]==g[i])] <- 1
    }

    ### g|u indication?
    indexing <- rep("g", length(quader[[1]]))
    indexing[which(apply(matrix(unlist(normQ),length(quader[[1]]),numberIndexVars),1,sum) %%2 != 0)] <- "u"
    return(list(quader=quader, normQ=normQ, indexing=indexing))
  }

  diagObj <- input[[1]]
  indexList <- input[[2]]
  protectionLevel <- input[[3]]
  allowZeros <- input[[4]]

  # TODO: error checking
  pI <- g_problemInstance(object)
  freqs <- g_freq(pI)
  sdcStatus <- g_sdcStatus(pI)

  if ( allowZeros == TRUE ) {
    indZeros <- which(sdcStatus=="z" & g_freq(pI)==0)
    if ( length(indZeros ) > 0 ) {
      sdcStatus[indZeros] <- "s"
    }
  }

  ### relevant Indices ### TODO FIXME CHECK!!!!
  relevantIndices <- which(sapply(indexList, function(x) { length(unique(x)) } ) > 1)

  resultObj <- list()
  # FIXME: What is with 1-dimensional data?
  limit <- length(diagObj$diagIndices)
  for ( z in 1:limit ) {
    g <- diagObj$indToProtect
    d <- sapply(indexList, function(x) { x[diagObj$diagIndices[z]] } )
    qInfo <- calcQInfo(g, d)

    # 2) position (indices==qPosition) of current quader in subTabObj
    valsQ <- pasteStrVec(unlist(qInfo$quader), length(qInfo$quader), coll=NULL)
    qPosition <- match(valsQ, g_strID(pI))
    suppStatus <- sdcStatus[qPosition]
    if ( !any(suppStatus == "z") ) {
      # 3) calculate various information about the selected quader (infoQuader)
      # 3.1) how many values would need to be suppressed for this quader
      indNonSupp <- which(sdcStatus[qPosition] == "s" & sdcStatus[qPosition] != "u")
      nrAdditionalSupps <- length(indNonSupp)

      # 3.2) whats the amount of information which needs to be suppressed?
      sumAdditionalSuppsFreq <- sum(freqs[qPosition[indNonSupp]])

      # 3.3) does the quader contains other single cells except for
      # the primary suppressed value (diaObj$cellToPretect) to check?
      # subIndices = current quader without primary suppressed cell to check
      indSingleItems <- setdiff(which(freqs[qPosition]==1),1)
      singleItems <- NULL
      indikatorSingleItems <- FALSE
      if( length(indSingleItems) >= 1 ) {
        indikatorSingleItems <- TRUE
        singleItems <- indSingleItems
      }

      # 3.5) is the quader protected enough? (protectionLevel)
      # we need to check for interval-protection only if protectionLevel > 0
      schutzInd <- TRUE
      schutz <- protectionLevel
      # FIXME: S|P (what to do with "x" that are temporarily "u"?)
      if( protectionLevel > 0 ) {
        if ( !all(qInfo$indexing =="u") ) {
          range <- min(freqs[qPosition[which(qInfo$indexing =="u")]], na.rm=TRUE) + min(freqs[qPosition[which(qInfo$indexing =="g")]], na.rm=TRUE)
          X <- freqs[diagObj$cellToProtect]
          if( X == 0 ) {
            tmpInd <- which(sdcStatus[qPosition] != "u" & freqs[qPosition] != 0)

            if( length(tmpInd) > 0 ) {
              # TODO: this needs testing !!! (page 60, repsilber)
              if( range <= min(freqs[tmpInd]) ) {
                schutzInd <- FALSE
                protectionLevel <- 0
              }
            }
          }
          else {
            schutz <- (100*range) / X
            if ( schutz < protectionLevel )
              schutzInd <- FALSE
          }
        }
      }

      # 4) return results
      # in this case, the cell is already protected, so we can stop!

      # allowZeros==TRUE: we have not found patterns without zeros
      # so we do not care for an 'optimal' solution that does not exist anyway
      if ( allowZeros == TRUE ) {
        if ( length(resultObj) == 100 ) {
          return(resultObj)
          break
        }
      } else {
        if( nrAdditionalSupps == 0 & schutzInd == TRUE & indikatorSingleItems == FALSE  ) {
          return(erg = NULL)
          break
        }
      }

      resultObj[[length(resultObj)+1]] <- list(
        quaderStrID = valsQ,
        indexing = qInfo$indexing,
        qPosition = qPosition,
        nrAdditionalSupps=nrAdditionalSupps,
        sumAdditionalSuppsFreq = sumAdditionalSuppsFreq,
        indikatorSingleItems = indikatorSingleItems,
        singleItems = singleItems,
        schutz = schutz,
        schutzInd = schutzInd
      )
    }
  }
  return(resultObj)
})

setMethod("c_ghmiter_suppress_quader", signature=c("sdcProblem", "list"), definition=function(object, input) {
  pI <- g_problemInstance(object)
  sdcStatus <- g_sdcStatus(pI)

  suppIndex <- setdiff(input$qPosition, input$qPosition[which(sdcStatus[input$qPosition]=="u")])
  s_sdcStatus(pI) <- list(index=suppIndex, vals=rep("x", length(suppIndex)))
  s_problemInstance(object) <- pI
  return(object)
})

setMethod("c_ghmiter_select_quader", signature=c("sdcProblem", "list"), definition=function(object, input) {
  infoObj <- input[[1]]
  suppMethod <- input[[2]]$suppMethod
  verbose <- input[[2]]$verbose

  sdcStatus <- g_sdcStatus(g_problemInstance(object))
  relevantIndices <- as.numeric(unlist(lapply(infoObj, '[', 'schutzInd'))[1])

  # already protected
  if ( is.null(infoObj) ) {
    suppObj <- NULL
    return(suppObj)
  }

  # not protected yet
  # which elements of iqsInfo are NULL?
  nullElements <- which(unlist(lapply(lapply(infoObj, '[[', 'qPosition'), function(x) { length(x) } )) == 0)
  if ( length(nullElements) > 0 ) {
    infoObj <- infoObj[-nullElements]
  }

  # put iqs together so that we can choose the optimal suppression scheme
  qIndexNr <- 1:length(infoObj)
  nrAdditionalSupps <- as.numeric(unlist(lapply(infoObj, '[', 'nrAdditionalSupps')))
  sumAdditionalSuppsFreq <- as.numeric(unlist(lapply(infoObj, '[', 'sumAdditionalSuppsFreq')))
  indikatorSingleItems <- as.logical(unlist(lapply(infoObj, '[', 'indikatorSingleItems')))
  schutz <- as.numeric(unlist(lapply(infoObj, '[', 'schutz')))
  schutzInd <- as.logical(unlist(lapply(infoObj, '[', 'schutzInd')))
  schutzInd <- as.logical(unlist(lapply(infoObj, '[', 'schutzInd')))

  possQuaders <- data.frame(qIndexNr, nrAdditionalSupps, sumAdditionalSuppsFreq, indikatorSingleItems, schutz, schutzInd)

  # are there any suppression schemes satisfying the necessary interval protection?
  indexIntervallOk <- FALSE
  if ( any(possQuaders$schutzInd==TRUE) )
    indexIntervallOk <- TRUE

  # do suppression schemes exist that do not contain single values?
  # these are preferred suppression schemes.
  existNonSingles <- FALSE
  if ( any(possQuaders$indikatorSingleItems == FALSE) ) {
    existNonSingles <- TRUE
    possQuaders <- possQuaders[possQuaders$indikatorSingleItems==FALSE,,drop=FALSE]
  }

  if( indexIntervallOk ) {
    if ( min(possQuaders$nrAdditionalSupps) > 0 & verbose == TRUE) {
      cat("# additional secondary Supps:", min(possQuaders$nrAdditionalSupps)," ")
    }
    if ( suppMethod == "minSupps" ) {
      possQuaders <- possQuaders[which(possQuaders$nrAdditionalSupps == min(possQuaders$nrAdditionalSupps)),,drop=FALSE]
      possQuaders <- possQuaders[which(possQuaders$sumAdditionalSuppsFreq == min(possQuaders$sumAdditionalSuppsFreq)),,drop=FALSE]
    }
    if ( suppMethod == "minSum" ) {
      possQuaders <- possQuaders[which(possQuaders$sumAdditionalSuppsFreq == min(possQuaders$sumAdditionalSuppsFreq)),,drop=FALSE]
      possQuaders <- possQuaders[which(possQuaders$nrAdditionalSupps == min(possQuaders$nrAdditionalSupps)),,drop=FALSE]
    }
    if ( suppMethod == "minSumLogs" ) {
      possQuaders <- possQuaders[which(possQuaders$sumAdditionalSuppsFreq == min(log(1+possQuaders$sumAdditionalSuppsFreq))),,drop=FALSE]
      possQuaders <- possQuaders[which(possQuaders$nrAdditionalSupps == min(possQuaders$nrAdditionalSupps)),,drop=FALSE]
    }

    # finally choose the suppression scheme
    possQuaders <- possQuaders[1,]
    suppObj <- infoObj[[possQuaders$qIndexNr]]
  }
  # problem: no suppression scheme is satisfying the
  # required interval protection
  else {
    # all cells in this subtable are already suppressed
    # -> everything is ok
    if( all(sdcStatus=="u" | sdcStatus == "x") )
      suppObj <- NULL

    # no suppression scheme satisfies the required interval protection
    # the suppression pattern with the max. protection level is selected
    else {
      possQuaders <- possQuaders[which(possQuaders$schutz == max(possQuaders$schutz)),]
      possQuaders <- possQuaders[which(possQuaders$nrAdditionalSupps == min(possQuaders$nrAdditionalSupps)),]
      possQuaders <- possQuaders[which(possQuaders$sumAdditionalSuppsFreq == min(possQuaders$sumAdditionalSuppsFreq)),]
      possQuaders <- possQuaders[1,]
      suppObj <- infoObj[[possQuaders$qIndexNr]]
    }
  }
  return(suppObj)
})

setMethod("c_ghmiter_supp_additional", signature=c("sdcProblem", "list"), definition=function(object, input) {
  diagObj <- input[[1]]
  infoObj <- input[[2]]
  suppObj <- input[[3]]
  suppMethod <- input[[4]]$suppMethod
  verbose <- input[[4]]$suppMethod

  ### Task: find quader (from) infoObj with following restrictions
  freqs <- g_freq(g_problemInstance(object))

  # - must not be suppObj itself
  cellToProtect <- diagObj$cellToProtect
  suppIndicesOrig <- suppObj$qPosition

  # the additional quader must non contain these indices
  # the singletons in the original suppressed pattern
  prohibitedIndices <- setdiff(suppIndicesOrig[which(freqs[suppIndicesOrig] == 1)],cellToProtect)

  # the possible indices
  possIndices <- lapply(infoObj, function(x) { x$qPosition } )

  # do the indices of the possible patterns contain any of the prohibited cells
  res <- sapply(1:length(possIndices), function(x) { any(prohibitedIndices %in% possIndices[[x]]) } )
  if ( all(res == TRUE ) ) {
    #warning("no additional cube could be found!\n")
    infoObj <- NULL
  } else {
    ind <- which(res==TRUE)
    infoObj <- infoObj[-ind]
  }

  suppObjNew <- c_ghmiter_select_quader(object, input=list(infoObj, input[[4]]))
  if ( !is.null(suppObjNew) ) {
    object <- c_ghmiter_suppress_quader(object, input=suppObjNew)
  }
  return(object)
})

setMethod("c_contributing_indices", signature=c("sdcProblem", "list"), definition=function(object, input) {
  strID <- input[[1]]
  dataObj <- g_dataObj(object)
  dimInfoObj <- g_dimInfo(object)
  dimInfo <- g_dim_info(dimInfoObj)
  pI <- g_problemInstance(object)

  if ( !strID %in% g_strID(pI) ) {
    stop("c_contributing_indices:: strID not found in the current problem!\n")
  }
  dims <- lapply(dimInfo, function(x) {
    g_dims(x)
  })
  indexVec <- which(g_str_id(dimInfoObj)==strID)
  # some (sub)totals need to be considered
  if( length(indexVec) == 0 ) {
    levInfo <- list()
    for ( z in 1:length(dimInfo) ) {
      subLevel <- substr(strID, g_str_info(dimInfoObj)[[z]][1], g_str_info(dimInfoObj)[[z]][2])
      if ( sum(as.numeric(subLevel)) == 0 ) {
        levInfo[[z]] <- sort(unique(unlist(dims[[z]])))
      } else {
        orderInd <- unlist(lapply(dims[[z]], function(x) { match(subLevel, x)}))
        if( min(orderInd, na.rm=TRUE) == 1 ) {
          levInfo[[z]] <- dims[[z]][[which(orderInd==1)]]
        } else {
          levInfo[[z]] <- subLevel
        }
      }
    }
    cellIndex <- pasteStrVec(unlist(expand.grid(levInfo)), length(levInfo))
    indexVec <- which(g_str_id(dimInfoObj) %in% cellIndex)
  }
  return(indexVec)
})

setMethod("c_reduce_problem", signature=c("sdcProblem", "list"), definition=function(object, input) {
  x <- object
  y <- input[[1]]

  pI <- g_problemInstance(x)
  dimInfo <- g_dimInfo(x)
  strInfo <- strInfoOrig <- g_str_info(dimInfo)

  if ( length(y) < 1 ) {
    stop("c_reduce_problem:: length of argument 'y' < 1!\n")
  }
  if ( !all(y %in% 1:g_nrVars(pI)) ) {
    stop("c_reduce_problem:: elements of indices y does not match with problem size!\n")
  }

  newDims <- lapply(1:length(strInfo), function(x) {
    substr(g_strID(pI)[y], strInfo[[x]][1], strInfo[[x]][2])
  })
  newDims2 <- lapply(1:length(newDims), function(x) {
    sort(unique(newDims[[x]]))
  })
  newDimsOrigCodes <- lapply(1:length(newDims), function(k) {
    c_match_orig_codes(object=dimInfo@dimInfo[[k]], input=newDims2[[k]])
  })

  lenNewDims <- sapply(newDims2, length)-1
  codesNew <- lapply(1:length(newDims), function(x) {
    c("@", rep("@@", lenNewDims[x]))
  })

  dimInfoOld <- lapply(1:length(newDims2), function(x) {
    init.dimVar(input=list(input=data.frame(codesNew[[x]], newDims2[[x]]), vName=paste('V',x,sep="")) )
  })
  dimInfoNew <- lapply(1:length(newDims2), function(x) {
    init.dimVar(input=list(input=data.frame(codesNew[[x]], newDimsOrigCodes[[x]]), vName=paste('V',x,sep="")) )
  })

  new.codes <- lapply(1:length(newDims), function(x) {
    dimInfoOld[[x]]@codesDefault[match(newDims[[x]], dimInfoOld[[x]]@codesOriginal)]
  })
  pI@strID <- pasteStrVec(unlist(new.codes), length(newDims))
  pI@Freq <- g_freq(pI)[y]
  if ( !is.null(g_w(pI)) ) {
    pI@w <- g_w(pI)[y]
  }
  numVars <- as.list(g_numVars(pI))
  if ( length(numVars) > 0 ) {
    for ( j in 1:length(numVars) ) {
      pI@numVars[[j]] <- numVars[[j]][y]
    }
  }
  pI@lb <- g_lb(pI)[y]
  pI@ub <- g_ub(pI)[y]
  pI@LPL <- g_LPL(pI)[y]
  pI@UPL <- g_UPL(pI)[y]
  pI@SPL <- g_SPL(pI)[y]
  pI@sdcStatus <- g_sdcStatus(pI)[y]
  x@dimInfo@dimInfo <- dimInfoNew

  # strInfo
  info <- c(0, cumsum(sapply(1:length(codesNew), function(x) { sum(sapply(table(codesNew[[x]]), nchar)) } )))
  for ( i in 2:length(info) ) {
    strInfo[[i-1]] <- c(info[i-1]+1, info[i] )
  }
  x@dimInfo@strInfo <- strInfo

  s_problemInstance(x) <- pI
  validObject(x)
  return(x)
})

setMethod("c_gen_structcuts", signature=c("sdcProblem", "list"), definition=function(object, input) {
  pI <- g_problemInstance(object)
  dimInfoObj <- g_dimInfo(object)
  partition <- c_make_partitions(input=list(objectA=pI, objectB=dimInfoObj))

  dimInfo <- g_dim_info(dimInfoObj)
  nrLevels <- length(dimInfo)
  nrVars <- g_nrVars(pI)
  primSupps <- g_primSupps(pI)
  strIDs <- g_strID(pI)
  indices <- partition$indices
  weights <- g_weight(pI)
  requiredCuts <- init.cutList(type='empty', input=list(nrCols=nrVars))
  strInfo <- g_str_info(dimInfoObj)
  x <- rep(0, nrVars)

  for ( z in seq_along(primSupps) ) {
    pSupp <- primSupps[z]
    currentPrimSupp <- strIDs[pSupp]
    matchInd <- unlist(lapply(1:length(indices), function(x) {
      lapply(1:length(indices[[x]]), function(y) {
        if ( !all(is.na(match(pSupp, indices[[x]][[y]]))) ) {c(x,y)}
      })
    }))
    if ( any(is.na(matchInd)) ) {
      stop('elements of matchInd must not be NA!\n')
    }

    splitMatchInd <- split(matchInd, rep(1:(length(matchInd)/2), each=2))

    for ( u in 1:length(splitMatchInd) ) {
      matchInd <- splitMatchInd[[u]]
      nrPow <- nrLevels - length(which(as.numeric(unlist(strsplit(partition$groups[[matchInd[1]]],"-")))==1))
      v1 <- v2 <- x
      index <- indices[[matchInd[1]]][[matchInd[2]]]
      v1[index] <- weights[index]
      v2[index] <- 1
      lim <- sum(sort(weights[index])[1:(2^nrPow)])
      if ( any(v1 != 0) ) {
        s_add_complete_constraint(requiredCuts) <- list(init.cutList(type='singleCut', input=list(vals=v1, dir=">=", rhs=lim)))
      }
      if ( any(v2 != 0) ) {
        s_add_complete_constraint(requiredCuts) <- list(init.cutList(type='singleCut', input=list(vals=v2, dir=">=", rhs=(2^nrPow))))
      }

      ### Todo: at least 2 suppressions in each dimension
      ### there is some error here! -> TODO: CHECK: FIXME!
      #for ( i in 1:length(dimInfo) ) {
      # lO <- dimInfo[[i]]
      # splitList <- lapply(strInfo[i], function(k) { seq(k[1], k[2]) } )
      # subStringToFix <- mySplitIndicesList(currentPrimSupp, splitList)
      # f <- mySplitIndicesList(strIDs, splitList)
      #
      # index <- which(f == subStringToFix)
      # v3 <- x
      # v4 <- x
      # v3[index] <- weights[index]
      # v4[index] <- 1
      # lim <- sum(sort(weights[index])[1:2])
      # if ( any(v3 != 0) ) {
      #   if ( !is.na(lim) ) {
      #     s_add_complete_constraint(requiredCuts) <- list(init.cutList(type='singleCut', input=list(vals=v3, dir=">=", rhs=lim)))
      #   }
      # }
      # if ( any(v4 != 0) ) {
      #   s_add_complete_constraint(requiredCuts) <- list(init.cutList(type='singleCut', input=list(vals=v4, dir=">=", rhs=2)))
      # }
      #}
    }
  }

  dupRows <- g_duplicated_rows(g_constraints(requiredCuts))
  if ( length(dupRows) > 0 ) {
    s_remove_complete_constraint(requiredCuts) <- list(dupRows)
  }
  return(requiredCuts)
})
