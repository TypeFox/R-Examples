#' @aliases get.problemInstance,problemInstance,character-method
#' @rdname get.problemInstance-method
setMethod(f='get.problemInstance', signature=c('problemInstance', 'character'),
  definition=function(object, type) {
    if ( !type %in% c('strID', 'nrVars', 'freq', 'w', 'numVars', 'sdcStatus',
        'lb', 'ub', 'LPL', 'SPL', 'UPL', 'primSupps', 'secondSupps',
        'forcedCells', 'hasPrimSupps', 'hasSecondSupps', 'hasForcedCells',
        'weight', 'suppPattern') ) {
      stop("get.problemInstance:: argument 'type' is not valid!\n")
    }

    if ( type == 'strID' ) {
      return(g_strID(object))
    }
    if ( type == 'nrVars' ) {
      return(g_nrVars(object))
    }
    if ( type == 'freq' ) {
      return(g_freq(object))
    }
    if ( type == 'w' ) {
      return(g_w(object))
    }
    if ( type == 'numVars' ) {
      return(g_numVars(object))
    }
    if ( type == 'sdcStatus' ) {
      return(g_sdcStatus(object))
    }
    if ( type == 'lb' ) {
      return(g_lb(object))
    }
    if ( type == 'ub' ) {
      return(g_ub(object))
    }
    if ( type == 'LPL' ) {
      return(g_LPL(object))
    }
    if ( type == 'UPL' ) {
      return(g_UPL(object))
    }
    if ( type == 'SPL' ) {
      return(g_SPL(object))
    }
    if ( type == 'primSupps' ) {
      return(g_primSupps(object))
    }
    if ( type == 'secondSupps' ) {
      return(g_secondSupps(object))
    }
    if ( type == 'forcedCells' ) {
      return(g_forcedCells(object))
    }
    if ( type == 'hasPrimSupps' ) {
      return(g_hasPrimSupps(object))
    }
    if ( type == 'hasSecondSupps' ) {
      return(g_hasSecondSupps(object))
    }
    if ( type == 'hasForcedCells' ) {
      return(g_hasForcedCells(object))
    }
    if ( type == 'weight' ) {
      return(g_weight(object))
    }
    if ( type == 'suppPattern' ) {
      return(g_suppPattern(object))
    }
  }
)

#' @aliases set.problemInstance,problemInstance,character,list-method
#' @rdname set.problemInstance-method
setMethod(f='set.problemInstance', signature=c('problemInstance', 'character', 'list'),
  definition=function(object, type, input) {
    index <- input[[1]]
    values <- input[[2]]
    if ( !type %in% c('lb', 'ub', 'LPL', 'UPL', 'SPL', 'sdcStatus') ) {
      stop("set.problemInstance:: check argument 'type'!\n" )
    }
    if ( !is.null(index) & length(values) != length(index) ) {
      stop("set.problemInstance:: arguments 'values' and 'index' differ in length!\n")
    }
    if ( !all(index %in% 1:g_nrVars(object)) ) {
      stop("set.problemInstance:: argument 'index' does not fit to given problem!\n")
    }
    if ( type == 'lb' ) {
      s_lb(object) <- list(index=index, vals=values)
    }
    if ( type == 'ub' ) {
      s_ub(object) <- list(index=index, vals=values)
    }
    if ( type == 'LPL' ) {
      s_LPL(object) <- list(index=index, vals=values)
    }
    if ( type == 'UPL' ) {
      s_UPL(object) <- list(index=index, vals=values)
    }
    if ( type == 'SPL' ) {
      s_SPL(object) <- list(index=index, vals=values)
    }
    if ( type == 'sdcStatus' ) {
      s_sdcStatus(object) <- list(index=index, vals=values)
    }
    #validObject(object)
    object
  }
)

#' @aliases calc.problemInstance,problemInstance,character,list-method
#' @rdname calc.problemInstance-method
setMethod(f='calc.problemInstance', signature=c('problemInstance', 'character','list'),
  definition=function(object, type, input) {
    if ( !type %in% c('makeMasterProblem', 'isProtectedSolution') ) {
      stop("calc.problemInstance:: check argument 'type'!\n" )
    }
    if ( type == 'makeMasterProblem' ) {
      return(c_make_masterproblem(object, input))
    }

    if ( type == 'isProtectedSolution' ) {
      return(c_is_protected_solution(object, input))
    }
  }
)

setMethod("g_sdcStatus", signature="problemInstance", definition=function(object) {
  object@sdcStatus
})

setMethod("g_primSupps", signature="problemInstance", definition=function(object) {
  which(g_sdcStatus(object)=="u")
})

setMethod("g_secondSupps", signature="problemInstance", definition=function(object) {
  which(g_sdcStatus(object)=="x")
})

setMethod("g_forcedCells", signature="problemInstance", definition=function(object) {
  which(g_sdcStatus(object)=="z")
})

setMethod("g_type", signature="problemInstance", definition=function(object) {
  object@type
})

setMethod("g_freq", signature="problemInstance", definition=function(object) {
  object@Freq
})

setMethod("g_strID", signature="problemInstance", definition=function(object) {
  object@strID
})

setMethod("g_UPL", signature="problemInstance", definition=function(object) {
  object@UPL
})

setMethod("g_LPL", signature="problemInstance", definition=function(object) {
  object@LPL
})

setMethod("g_SPL", signature="problemInstance", definition=function(object) {
  object@SPL
})

setMethod("g_nrVars", signature="problemInstance", definition=function(object) {
  length(g_strID(object))
})

setMethod("g_lb", signature="problemInstance", definition=function(object) {
  object@lb
})

setMethod("g_ub", signature="problemInstance", definition=function(object) {
  object@ub
})

setMethod("g_w", signature="problemInstance", definition=function(object) {
  object@w
})

setMethod("g_numVars", signature="problemInstance", definition=function(object) {
  object@numVars
})

setMethod("g_hasPrimSupps", signature="problemInstance", definition=function(object) {
  return(length(g_primSupps(object)) > 0)
})

setMethod("g_hasSecondSupps", signature="problemInstance", definition=function(object) {
  return(length(g_secondSupps(object)) > 0)
})

setMethod("g_hasForcedCells", signature="problemInstance", definition=function(object) {
  return(length(g_forcedCells(object)) > 0)
})

setMethod("g_weight", signature="problemInstance", definition=function(object) {
  wg <- g_w(object)
  if ( !is.null(wg) ) {
    return(wg)
  } else {
    return(g_freq(object))
  }
})

setMethod("g_suppPattern", signature="problemInstance", definition=function(object) {
  suppPattern <- rep(0, g_nrVars(object))
  if ( g_hasPrimSupps(object) ) {
    suppPattern[g_primSupps(object)] <- 1
  }
  secondSupps <- g_secondSupps(object)
  if ( length(secondSupps) > 0 ) {
    suppPattern[secondSupps] <- 1
  }
  return(suppPattern)
})

setReplaceMethod("s_sdcStatus", signature=c("problemInstance", "list"), definition=function(object, value) {
  sdcStatus <- g_sdcStatus(object)
  index <- value$index
  values <- value$vals
  if ( length(values) > length(sdcStatus) ) {
    stop("s_sdcStatus:: length of 'sdcVec' must be <=",length(sdcStatus),"\n")
  }
  if ( is.null(index) ) {
    indexVec <- 1:length(sdcStatus)
    if ( length(values) != length(sdcStatus) ) {
      stop("s_sdcStatus:: length of 'values' must be ==",length(sdcStatus), "if 'index' is NULL!\n")
    }
    object@sdcStatus[indexVec] <- values
  }
  if ( !is.null(index) ) {
    if ( !all(index %in% 1:length(sdcStatus)) ) {
      stop("s_sdcStatus:: elements of 'index' must be in 1:",length(sdcStatus),"!\n")
    }
    object@sdcStatus[index] <- values
  }
  validObject(object)
  object
})

setReplaceMethod("s_lb", signature=c("problemInstance", "list"), definition=function(object, value) {
  object@lb[value$index] <- value$vals
  validObject(object)
  object
})

setReplaceMethod("s_ub", signature=c("problemInstance", "list"), definition=function(object, value) {
  object@ub[value$index] <- value$vals
  validObject(object)
  object
})

setReplaceMethod("s_LPL", signature=c("problemInstance", "list"), definition=function(object, value) {
  object@LPL[value$index] <- value$vals
  validObject(object)
  object
})

setReplaceMethod("s_UPL", signature=c("problemInstance", "list"), definition=function(object, value) {
  object@UPL[value$index] <- value$vals
  validObject(object)
  object
})

setReplaceMethod("s_SPL", signature=c("problemInstance", "list"), definition=function(object, value) {
  object@SPL[value$index] <- value$vals
  validObject(object)
  object
})


setMethod("c_make_masterproblem", signature=c("problemInstance", "list"), definition=function(object, input) {
  mProb <- NULL
  if ( g_hasPrimSupps(object) ) {
    objective <- g_weight(object)
    primSupps <- g_primSupps(object)
    nrVars <- g_nrVars(object)

    M <- init.simpleTriplet(type='simpleTriplet', input=list(mat=matrix(0, nrow=0, ncol=nrVars)))
    direction <- rep("==", g_nr_rows(M))
    rhs <- rep(1, g_nr_rows(M))

    # cells with sdcStatus=="z" must be published
    if ( g_hasForcedCells(object) ) {
      forcedCells <- g_forcedCells(object)
      if ( length(forcedCells) > 0 ) {
        for ( i in seq_along(forcedCells) ) {
          M <- c_add_row(M, input=list(index=forcedCells[i], values=1))
        }
        direction <- c(direction, rep("==", length(forcedCells)))
        rhs <- c(rhs, rep(0, length(forcedCells)))
      }
    }
    types <- rep("C", nrVars)
    boundsLower <- list(ind=1:nrVars, val=rep(0, nrVars))
    boundsUpper <- list(ind=1:nrVars, val=rep(1, nrVars))

    if ( length(primSupps) > 0 ) {
      boundsLower$val[primSupps] <- 1
    }
    mProb <- new("linProb",
      objective=objective,
      constraints=M,
      direction=direction,
      rhs=rhs,
      boundsLower=boundsLower,
      boundsUpper=boundsUpper,
      types=types)
  }
  return(mProb)
})

setMethod("c_is_protected_solution", signature=c("problemInstance", "list"), definition=function(object, input) {
  input1 <- input$input1 # ~ limitsDown
  input2 <- input$input2 # ~ limitsUp
  primSupps <- g_primSupps(object)
  if ( length(input1) != length(input2) ) {
    stop("c_is_protected_solution:: parameters 'input1 (~limitsDown)' and 'input2 (~limitsUp)' differ in length!\n")
  }
  if ( length(input1) != length(primSupps) ) {
    stop("c_is_protected_solution:: parameter 'limits.x' and length of primary suppressed cells differ!\n")
  }

  protected <- TRUE
  weights <- g_weight(object)[primSupps]

  limits <- list()
  limits$LPL <- g_LPL(object)[primSupps]
  limits$UPL <- g_UPL(object)[primSupps]
  limits$SPL <- g_SPL(object)[primSupps]

  if ( any(weights - input1 < limits[[1]]) == TRUE ) {
    protected <- FALSE
  }

  if ( any(input2 - weights  < limits[[2]]) == TRUE ) {
    protected <- FALSE
  }

  if ( any(input2 - input1  < limits[[3]]) == TRUE ) {
    protected <- FALSE
  }
  return(protected)
})

