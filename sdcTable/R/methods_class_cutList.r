#' @aliases get.cutList,cutList,character-method
#' @rdname get.cutList-method
setMethod(f="get.cutList", signature=c("cutList", "character"),
  definition=function(object, type) {
    if ( !type %in% c("constraints", "direction", "rhs", "nrConstraints") ) {
      stop("get.cutList:: argument 'type' is not valid!\n")
    }
    if ( type == "constraints" ) {
      return(g_constraints(object))
    }
    if ( type == "direction" ) {
      return(g_direction(object))
    }
    if ( type == "rhs" ) {
      return(g_rhs(object))
    }
    if ( type == "nrConstraints" ) {
      return(g_nr_constraints(object))
    }
  }
)

#' @aliases set.cutList,cutList,character,list-method
#' @rdname set.cutList-method
setMethod(f='set.cutList', signature=c("cutList", "character", "list"),
  definition=function(object, type, input) {
    if ( !type %in% c("addCompleteConstraint", "removeCompleteConstraint") ) {
      stop("set.cutList:: argument 'type' is not valid!\n")
    }
    if ( type == "addCompleteConstraint" ) {
      s_add_complete_constraint(object) <- input[[1]]
    }
    if ( type == "removeCompleteConstraint" ) {
      s_remove_complete_constraint(object) <- input[[1]]
    }
    validObject(object)
    return(object)
  }
)

#' @aliases calc.cutList,cutList,character,list-method
#' @rdname calc.cutList-method
setMethod(f="calc.cutList", signature=c("cutList", "character", "list"),
  definition=function(object, type, input) {
    if ( !type %in% c("strengthen", "checkViolation", "bindTogether") ) {
      stop("calc.cutList:: argument 'type' is not valid!\n")
    }
    if ( type == "strengthen" ) {
      return(c_strengthen(object))
    }
    if ( type == "checkViolation" ) {
      return(c_check_violation(object, input))
    }
    if ( type == "bindTogether" ) {
      return(c_bind_together(object, input))
    }
  }
)

#' @aliases init.cutList,character,list-method
#' @rdname init.cutList-method
setMethod(f='init.cutList', signature=c('character', 'list'),
  definition=function(type, input) {
    if ( !type %in% c('empty', 'singleCut', 'multipleCuts') ) {
      stop("initialize.cutList:: argument 'type' is not valid!\n")
    }

    x <- new("cutList")
    if ( type == 'empty' ) {
      x@con <- init.simpleTriplet(type='simpleTriplet', input=list(mat=matrix(0, nrow=0, ncol=input$nrCols)))
    }

    if ( type == 'singleCut' ) {
      x@con <- init.simpleTriplet(type='simpleTriplet', input=list(mat=matrix(input$vals, nrow=1, ncol=length(input$vals))))
      x@direction <- input$dir
      x@rhs <- input$rhs
    }

    if ( type == 'multipleCuts' ) {
      x@con <- input$mat
      x@direction <- input$dir
      x@rhs <- input$rhs
    }
    validObject(x)
    x
  }
)

setMethod(f="g_constraints", signature=c("cutList"), definition=function(object) {
  return(object@con)
})

setMethod(f="g_direction", signature=c("cutList"), definition=function(object) {
  return(object@direction)
})

setMethod(f="g_rhs", signature=c("cutList"), definition=function(object) {
  return(object@rhs)
})

setMethod(f="g_nr_constraints", signature=c("cutList"), definition=function(object) {
  return(length(g_rhs(object)))
})

setReplaceMethod(f="s_add_complete_constraint", signature=c("cutList", "list"), definition=function(object, value) {
  input <- value[[1]]
  if ( g_nr_cols(g_constraints(object)) != g_nr_cols(g_constraints(input)) ) {
    stop("s_add_complete_constraint:: nrCols of 'object' and 'input' differ!\n")
  }
  if ( g_nr_constraints(input) > 0 ) {
    con <- g_constraints(input)
    for ( k in 1:g_nr_rows(con) ) {
      x <- g_row(con, input=list(k))
      object@con <- c_add_row(g_constraints(object), input=list(index=g_col_ind(x), values=g_values(x)))
    }
    object@direction <- c(g_direction(object), g_direction(input))
    object@rhs <- c(g_rhs(object), g_rhs(input))
  }
  return(object)
})

setReplaceMethod(f="s_remove_complete_constraint", signature=c("cutList", "list"), definition=function(object, value) {
  input <- value[[1]]
  if ( !all(input %in% 1:length(g_rhs(object))) ) {
    stop("elements of argument 'input' must be >=1 and <=",length(g_rhs(object)),"!\n")
  }
  object@con <- c_remove_row(g_constraints(object), input=list(input))
  object@direction <- g_direction(object)[-input]
  object@rhs <- g_rhs(object)[-input]
  return(object)
})

setMethod(f="c_strengthen", signature=c("cutList"), definition=function(object) {
  con <- g_constraints(object)
  nrRows <- g_nr_rows(con)
  nrCols <- g_nr_cols(con)
  rhs <- g_rhs(object)
  vals <- lapply(1:nrRows, function(x) {
    g_values(g_row(con, input=list(x)))
  })
  vals <- lapply(1:nrRows, function(x) {
    sapply(vals[[x]], function(y) { min(y, rhs[x]) } )
  })
  colInds <- lapply(1:nrRows, function(x) {
    g_col_ind(g_row(con, input=list(x)))
  })

  st <- init.simpleTriplet(type='simpleTriplet', input=list(mat=matrix(0, 0, nrCols)))
  lapply(1:nrRows, function(x) {
    st <<- c_add_row(st, input=list(index=colInds[[x]], values=vals[[x]]))
  })
  object@con <- st
  validObject(object)
  return(object)
})

setMethod(f="c_check_violation", signature=c("cutList"), definition=function(object, input) {
  sol <- input[[1]]
  w <- input[[2]]
  con <- g_constraints(object)
  if ( g_nr_cols(con) != length(sol) ) {
    stop("c_check_violation:: length of argument 'sol' and number of columns of 'someCuts' differ!\n")
  }
  if ( length(sol) != length(w) ) {
    stop("c_check_violation:: length of argument 'sol' 'w' differ!\n")
  }
  res <- rep(NA, g_nr_rows(con))
  rhs <- g_rhs(object)
  dir <- g_direction(object)
  for ( z in 1:g_nr_rows(con) ) {
    colInd <- g_col_ind(g_row(con, input=list(z)))
    expr <- paste(sum(sol[colInd]*w[colInd]), dir[z], rhs[z])
    res[z] <- eval(parse(text=expr))
  }
  return(any(res == FALSE))
})

setMethod(f="c_bind_together", signature=c("cutList"), definition=function(object, input) {
  object1 <- object
  object2 <- input[[1]]
  x <- new("cutList")
  x@con <- c_bind(object=g_constraints(object1), input=list(g_constraints(object2), bindRow=TRUE))
  x@direction <- c(g_direction(object1), g_direction(object2))
  x@rhs <- c(g_rhs(object1), g_rhs(object2))
  validObject(x)
  return(x)
})

