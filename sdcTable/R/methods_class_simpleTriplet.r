######################
### global methods ###
######################
setMethod(f="as.matrix", signature="simpleTriplet",
  definition=function(x, ...) {
    M <- matrix(0, nrow=g_nr_rows(x), ncol=g_nr_cols(x))
    i.x <- g_row_ind(x)
    j.x <- g_col_ind(x)
    v.x <- g_values(x)
    for ( i in 1:g_nr_cells(x) ) {
      M[i.x[i], j.x[i]] <- v.x[i]
    }
    return(M)
  }
)

##############################################
### methods only for class 'simpleTriplet' ###
##############################################
#' @aliases get.simpleTriplet,simpleTriplet,character,list-method
#' @rdname get.simpleTriplet-method
setMethod(f="get.simpleTriplet", signature=c("simpleTriplet", "character", 'list'),
  definition=function(object, type, input) {
    if ( !type %in% c("rowInd", "colInd", "values", "nrRows",
        "nrCols", "nrCells", "duplicatedRows", "transpose",
        "getRow", "getCol") ) {
      stop("get.simpleTriplet:: argument 'type' is not valid!\n")
    }
    if ( type == "rowInd" ) {
      return(g_row_ind(object))
    }
    if ( type == "colInd" ) {
      return(g_col_ind(object))
    }
    if ( type == "values" ) {
      return(g_values(object))
    }
    if ( type == "nrRows" ) {
      return(g_nr_rows(object))
    }
    if ( type == "nrCols" ) {
      return(g_nr_cols(object))
    }
    if ( type == "nrCells" ) {
      return(g_nr_cells(object))
    }
    if ( type == "duplicatedRows" ) {
      return(g_duplicated_rows(object))
    }
    if ( type == "transpose" ) {
      return(g_transpose(object))
    }
    if ( type == "getRow" ) {
      return(g_row(object, input))
    }
    if ( type == "getCol" ) {
      return(g_col(object, input))
    }
  }
)

#' @aliases calc.simpleTriplet,simpleTriplet,character,list-method
#' @rdname calc.simpleTriplet-method
setMethod(f="calc.simpleTriplet", signature=c("simpleTriplet", "character", "list"),
  definition=function(object, type, input) {
    if ( !type %in% c("removeRow", "removeCol", "addRow", "addCol",
      "modifyRow", "modifyCol", "modifyCell", "bind") ) {
      stop("calc.simpleTriplet:: check argument 'type'!\n")
    }
    if ( type == "removeRow" ) {
      return(c_remove_row(object, input))
    }
    if ( type == "removeCol" ) {
      return(c_remove_col(object, input))
    }
    if ( type == "addRow" ) {
      return(c_add_row(object, input))
    }
    if ( type == "addCol" ) {
      return(c_add_col(object, input))
    }
    if ( type == "modifyRow" ) {
      return(c_modify_row(object, input))
    }
    if ( type == 'modifyCol' ) {
      return(c_modify_col(object, input))
    }
    if ( type == "modifyCell" ) {
      return(c_modify_cell(object, input))
    }
    if ( type == "bind" ) {
      return(c_bind(object, input))
    }
  }
)

#' @aliases init.simpleTriplet,character,list-method
#' @rdname init.simpleTriplet-method
setMethod(f='init.simpleTriplet', signature=c('character', 'list'),
  definition=function(type, input) {
    if ( !type %in% c('simpleTriplet', 'simpleTripletDiag') ) {
      stop("init.simpleTriplet:: check argument 'type'!\n")
    }

    if ( type == 'simpleTriplet' ) {
      matA <- input$mat
      dims <- dim(matA)
      v <- as.vector(t(matA))
      ind <- v!=0
      i <- rep(1:dims[1], each=dims[2])[ind==TRUE]
      j <- rep(1:dims[2], length=dims[1]*dims[2])[ind]
      v <- v[ind]
      out <- new("simpleTriplet",
        i=i,
        j=j,
        v=v,
        nrRows=dims[1],
        nrCols=dims[2]
      )
    }

    if ( type == 'simpleTripletDiag' ) {
      nrRows <- input$nrRows
      negative <- input$negative
      i <- j <- 1:nrRows
      if ( negative ) {
        v <- rep(-1, nrRows)
      } else {
        v <- rep(1, nrRows)
      }
      out <- new("simpleTriplet",
        i=i,
        j=j,
        v=v,
        nrRows=nrRows,
        nrCols=nrRows
      )
    }
    validObject(out)
    return(out)
  }
)

# get-methods
setMethod("g_row_ind", signature=c("simpleTriplet"), definition=function(object) {
  return(object@i)
})

setMethod("g_col_ind", signature=c("simpleTriplet"), definition=function(object) {
  return(object@j)
})

setMethod("g_values", signature=c("simpleTriplet"), definition=function(object) {
  return(object@v)
})

setMethod("g_nr_rows", signature=c("simpleTriplet"), definition=function(object) {
  return(object@nrRows)
})

setMethod("g_nr_cols", signature=c("simpleTriplet"), definition=function(object) {
  return(object@nrCols)
})

setMethod("g_nr_cells", signature=c("simpleTriplet"), definition=function(object) {
  return(length(g_values(object)))
})

setMethod("g_duplicated_rows", signature=c("simpleTriplet"), definition=function(object) {
  i <- g_row_ind(object)
  j <- g_col_ind(object)
  v <- g_values(object)
  len <- g_nr_rows(object)
  o <- order(i, j)
  y <- split(paste(j[o], v[o], sep = "\r"), i[o])
  tmp <- character(len)
  names(tmp) <- seq_along(tmp)
  tmp[names(y)] <- sapply(y, paste, collapse = "\r")
  dupRows <- which(duplicated(tmp))
  if ( length(dupRows) == 0 ) {
    dupRows <- NULL
  }
  return(dupRows)
})

setMethod("g_transpose", signature=c("simpleTriplet"), definition=function(object) {
  out <- init.simpleTriplet(type='simpleTriplet', input=list(mat=matrix(0, nrow=g_nr_cols(object), ncol=0)))
  for ( i in 1:g_nr_rows(object) ) {
    r <- g_row(object, input=list(i))
    out <- c_add_col(out, input=list(index=g_col_ind(r), values=g_values(r)))
  }
  return(out)
})

setMethod("g_row", signature=c("simpleTriplet", "list"), definition=function(object, input) {
  # if somebody specifies a vector of length > 1, the row with the first index is returned
  index <- input[[1]][1]
  if ( !index %in% 1:g_nr_rows(object) ) {
    stop("g_row:: parameter 'index' must be >=1 and <=",g_nr_rows(object),"!\n")
  }

  out <- NULL
  indI <- which(g_row_ind(object) == index)
  if ( length(indI) > 0 ) {
    out <- new("simpleTriplet",
      i=rep(1, length(indI)),
      j=g_col_ind(object)[indI],
      v=g_values(object)[indI],
      nrRows=1,
      nrCols=g_nr_cols(object)
    )
  }
  return(out)
})
setMethod("g_col", signature=c("simpleTriplet", "list"), definition=function(object, input) {
  # if somebody specifies a vector of length > 1, the column with the first index is returned
  index <- input[[1]][1]
  if ( !index %in% 1:g_nr_cols(object) ) {
    stop("g_col:: parameter 'index' must be >=1 and <=",g_nr_cols(object),"!\n")
  }
  out <- NULL
  indJ <- which(g_col_ind(object) == index)
  if ( length(indJ) > 0 ) {
    out <- new("simpleTriplet",
      i=g_row_ind(object)[indJ],
      j=rep(1, length(indJ)),
      v=g_values(object)[indJ],
      nrRows=g_nr_rows(object),
      nrCols=1
    )
  }
  return(out)
})

setMethod("c_remove_row", signature=c("simpleTriplet"), definition=function(object, input) {
  index <- input[[1]]
  if ( !all(index %in% 1:g_nr_rows(object)) ) {
    stop("c_remove_row:: check dimensions of parameter 'index'!\n")
  }
  ind <- which(g_row_ind(object) %in% index)
  if ( length(ind) > 0 ) {
    object@i <- g_row_ind(object)[-ind]
    object@j <- g_col_ind(object)[-ind]
    object@v <- g_values(object)[-ind]
  }
  object@nrRows <- g_nr_rows(object)-length(index)
  object@i <- rep(1:length(unique(g_row_ind(object))), table(g_row_ind(object)))
  validObject(object)
  return(object)
})

setMethod("c_remove_col", signature=c("simpleTriplet"), definition=function(object, input) {
  index <- input[[1]]
  if ( !all(index %in% 1:g_nr_cols(object)) ) {
    stop("c_remove_col:: check dimensions of parameter 'index'!\n")
  }
  ind <- which(g_col_ind(object) %in% index)
  if ( length(ind) > 0 ) {
    object@i <- g_row_ind(object)[-ind]
    object@j <- g_col_ind(object)[-ind]
    object@v <- g_values(object)[-ind]
  }
  object@nrCols <- g_nr_cols(object)-length(index)
  object@j <- rep(1:length(unique(g_col_ind(object))), table(g_col_ind(object)))
  validObject(object)
  return(object)
})

setMethod("c_add_row", signature=c("simpleTriplet"), definition=function(object, input) {
  index <- input[[1]]
  values <- input[[2]]
  if ( !all(index %in% 1:g_nr_cols(object)) ) {
    stop("c_add_row:: check dimensions of parameter 'index'!\n")
  }
  if ( length(index) != length(values) ) {
    stop("c_add_row:: dimensions of 'index' and 'values' do not match!\n")
  }
  rowInd <- g_nr_rows(object)+1
  object@nrRows <- rowInd
  ind <- which(values != 0)
  nrAddedCells <- length(ind)
  if ( nrAddedCells > 0 ) {
    object@i <- c(g_row_ind(object), rep(rowInd, nrAddedCells))
    object@j <- c(g_col_ind(object), index[ind])
    object@v <- c(g_values(object), values[ind])
  }
  validObject(object)
  return(object)
})

setMethod("c_add_col", signature=c("simpleTriplet"), definition=function(object, input) {
  index <- input[[1]]
  values <- input[[2]]
  if ( !all(index %in% 1:g_nr_rows(object)) ) {
    stop("c_add_col:: check dimensions of parameter 'index'!\n")
  }
  if ( length(index) != length(values) ) {
    stop("c_add_col:: dimensions of 'index' and 'values' do not match!\n")
  }
  colInd <- g_nr_cols(object)+1
  object@nrCols <- colInd
  ind <- which(values != 0)
  nrAddedCells <- length(ind)
  if ( nrAddedCells > 0 ) {
    object@i <- c(g_row_ind(object), index[ind])
    object@j <- c(g_col_ind(object), rep(colInd, nrAddedCells))
    object@v <- c(g_values(object), values[ind])
  }
  validObject(object)
  return(object)
})

setMethod("c_modify_row", signature=c("simpleTriplet"), definition=function(object, input) {
  rowInd <- input[[1]]
  colInd <- input[[2]]
  values <- input[[3]]
  if ( length(rowInd) != 1 ) {
    stop("c_modify_row:: length of parameter 'rowInd' must equal 1!\n")
  }
  if ( !rowInd %in% 1:g_nr_rows(object) ) {
    stop("c_modify_row:: check dimensions of parameter 'rowInd'!\n")
  }
  if ( !all(colInd %in% 1:g_nr_cols(object)) ) {
    stop("c_modify_row:: check dimensions of parameter 'colInd'!\n")
  }
  if ( length(colInd) != length(values) ) {
    stop("c_modify_row:: dimensions of 'colInd' and 'values' do not match!\n")
  }
  ind <- which(g_row_ind(object) %in% rowInd)
  if ( length(ind) == 0 ) {
    stop("c_modify_row:: no row to modify!\n")
  }
  for ( j in seq_along(colInd) ) {
    object <- c_modify_cell(object, input=list(rowInd, colInd[j], values[j]))
  }
  validObject(object)
  return(object)
})

setMethod("c_modify_col", signature=c("simpleTriplet"), definition=function(object, input) {
  rowInd <- input[[1]]
  colInd <- input[[2]]
  values <- input[[3]]

  if ( length(colInd) != 1 ) {
    stop("c_modify_col:: length of parameter 'colInd' must equal 1!\n")
  }
  if ( !all(rowInd %in% 1:g_nr_rows(object)) ) {
    stop("c_modify_col:: check dimensions of parameter 'rowInd'!\n")
  }
  if ( !colInd %in% 1:g_nr_cols(object) ) {
    stop("c_modify_col:: check dimensions of parameter 'colInd'!\n")
  }
  if ( length(rowInd) != length(values) ) {
    stop("c_modify_col:: dimensions of 'rowInd' and 'values' do not match!\n")
  }
  ind <- which(g_col_ind(object) %in% colInd)
  if ( length(ind) == 0 ) {
    stop("c_modify_col:: no column to modify!\n")
  }
  for ( i in seq_along(rowInd) ) {
    object <- c_modify_cell(object, input=list(rowInd[i], colInd, values[i]))
  }
  validObject(object)
  return(object)
})

setMethod("c_modify_cell", signature=c("simpleTriplet"), definition=function(object, input) {
  rowInd <- input[[1]]
  colInd <- input[[2]]
  values <- input[[3]]

  if ( any(length(values), length(rowInd), length(colInd) != 1) ) {
    stop("c_modify_cell:: length of all arguments 'rowInd', 'colInd', 'values' must equal 1!\n")
  }
  if ( !all(colInd %in% 1:g_nr_cols(object)) ) {
    stop("c_modify_cell:: check dimensions of parameter 'colInd'!\n")
  }
  if ( !all(rowInd %in% 1:g_nr_rows(object)) ) {
    stop("c_modify_cell:: check dimensions of parameter 'rowInd'!\n")
  }
  ind <- which(g_row_ind(object) == rowInd & g_col_ind(object) == colInd)
  if ( length(ind) == 1 & values != 0 ) {
    object@v[ind] <- values
  } else if ( length(ind)==1 & values == 0 ) {
    object@i <- g_row_ind(object)[-ind]
    object@j <- g_col_ind(object)[-ind]
    object@v <- g_values(object)[-ind]
  } else if ( length(ind)==0 & values!=0 ) {
    object@i <- c(g_row_ind(object), rowInd)
    object@j <- c(g_col_ind(object), colInd)
    object@v <- c(g_values(object), values)
  }
  validObject(object)
  return(object)
})

setMethod("c_bind", signature=c("simpleTriplet"), definition=function(object, input) {
  object1 <- object
  object2 <- input[[1]]
  bindRow <- input[[2]]

  if ( bindRow == TRUE ) {
    # "rbind"
    if ( g_nr_cols(object1) != g_nr_cols(object2) ) {
      stop("c_bind:: nr of columns of 'object1' and 'object2' differ!\n")
    }
    out <- init.simpleTriplet(type='simpleTriplet',
      input=list(mat=matrix(0, nrow=g_nr_rows(object1) + g_nr_rows(object2), ncol=g_nr_cols(object1))))
    object2@i <- object2@i + g_nr_rows(object1)
  } else {
    # "cbind"
    if ( g_nr_rows(object1) != g_nr_rows(object2) ) {
      stop("c_bind:: nr of rows of 'object1' and 'object2' differ!\n")
    }
    out <- init.simpleTriplet(type='simpleTriplet',
      input=list(mat=matrix(0, nrow=g_nr_rows(object1), ncol=g_nr_cols(object1)+g_nr_cols(object2))))
    object2@j <- object2@j + g_nr_cols(object1)
  }
  out@i <- c(object1@i, object2@i)
  out@j <- c(object1@j, object2@j)
  out@v <- c(object1@v, object2@v)
  validObject(out)
  return(out)
})

