#' Cross Tabulation and Table Creation
#' 
#' Upgrades table to a generic function and implements a method 
#' for ff vectors which works for ff factors.
#' For other arguments passed on to table, uses \code{\link[base]{table}}\cr
#'
#' table.ff uses the cross-classifying factors to build a contingency table of 
#' the counts at each combination of factor levels.\cr
#' If \code{...} does not contain factors, \code{unique.ff} will add a levels 
#' attribute to the non-factors.
#' 
#' @seealso \code{\link[base]{table}}
#'
#' @param ... \code{ff} factors or \code{ff} integers
#' @param exclude see \code{\link[base]{table}}
#' @param useNA see \code{\link[base]{table}}
#' @param dnn see \code{\link[base]{table}}
#' @param deparse.level see \code{\link[base]{table}}
#' @usage table(..., exclude = if (useNA == "no") c(NA, NaN), 
#' useNA = c("no", "ifany", "always"), dnn = list.names(...), deparse.level = 1)
#' @return \code{\link[base]{table}} object
#' @export table
table <- function( ...
                 , exclude = if (useNA == "no") c(NA, NaN)
                 , useNA = c("no","ifany", "always")
                 , dnn = list.names(...)
                 , deparse.level = 1
){
  UseMethod("table")
}

#' @export 
table.default <- base::table

#' @export
#' @export table.ff
#' @usage table(..., exclude = if (useNA == "no") c(NA, NaN), 
#' useNA = c("no", "ifany", "always"), dnn = list.names(...), deparse.level = 1)
#' @rdname table
#' @aliases table
#' 
table.ff <- function( ...
                     , exclude = if (useNA == "no") c(NA, NaN)
                     , useNA = c("no","ifany", "always")
                     , dnn = list.names(...)
                     , deparse.level = 1
){
  ###
  args <- list(...)
  tab <- NULL
  useNA <- match.arg(useNA)
  
  dat <- do.call(ffdf, args) # create a ffdf  for estimating good chunking size 
                             #and checking if ... have equal length
  colnames(dat) <- names(args)
  ### Cover non-factors like integers by adding a levels attribute
  if(sum(!vmode(dat) %in% c("byte", "short", "integer")) > 0){  	
    stop(sprintf("Only vmodes integer currently allowed - are you sure ... contains only factors or integers?"))
  }
  nonfactors <- sapply(colnames(dat), FUN=function(column, dat) !is.factor.ff(dat[[column]]), dat=dat)
  nonfactors <- names(nonfactors)[nonfactors == TRUE]
  if(length(nonfactors) > 0){
    for(column in nonfactors){
      dat[[column]] <- as.character.ff(dat[[column]])
    } 
  }
  
  for (i in chunk(dat)){
    Log$chunk(i)
    factors <- unname(as.list(dat[i,, drop=FALSE]))
    factors$exclude <- exclude
    factors$useNA <- useNA
    factors$deparse.level <- deparse.level
    
    ttab <- do.call(table,factors)
    tab <- if (is.null(tab)){ 
      ttab
    }
    else { tab + ttab
    }
    #names(dimnames(tab)) <- names(dimnames(ttab))
  }
  return(tab)	
}

#borrowed from table
list.names <- function(...) {
  l <- as.list(substitute(list(...)))[-1L]
  nm <- names(l)
  fixup <- if (is.null(nm)) 
    seq_along(l)
  else nm == ""
  
  dep <- sapply(l[fixup], function(x) if (is.symbol(x)) as.character(x) else "")
  if (is.null(nm)) 
    dep
  else {
    nm[fixup] <- dep
    nm
  }
}

# setGeneric( "table"
# , signature="..."
# )

# setMethod( "table"
# , "ff"
# , table.ff
# )

# setMethod( "table"
# , "ff_vector"
# , table.ff
# )
