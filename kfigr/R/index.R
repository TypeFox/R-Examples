#' Chunk Indexing
#'
#' (Internal) index a chunk. Should not be called by the user directly.
#'
#' @seealso \code{\link{figr}}, \code{\link{anchors}}
#'
#' @param label The chunk label.
#' @param type The type of chunk to be indexed.
#' @return The rank of the indexed chunk.
#'
#' @importFrom stats setNames
index <- function(label, type){
  if (!type %in% get("types", envir=anchorenv)){  # check or define the figr type
    assign("types", c(get("types", envir=anchorenv), type), 
	       envir=anchorenv)
    assign(type, character(0), envir=anchorenv)
  }
  if (!label %in% get(type, envir=anchorenv)){  # index the chunk if needed
    assign(type, c(get(type, envir=anchorenv), label), 
	       envir=anchorenv)
    assign("index", c(get("index", envir=anchorenv), setNames(list(type), label)),
	       envir=anchorenv)
  }
  # get the figr index
  number <- grep(label, anchors(type))
  # update history
  updatehist <- function(h, caller){
    caller=deparse(caller)
	if(caller == "index(label = options$label, type = options$anchor)")
      called.by <- 'hook_anchor'
    else if(caller == "index(label = label, type = type)")
      called.by <- "figr"
    # seems to be the only way to append a vector as a single element
    h[[length(h) + 1L]] <- data.frame(label=label, type=type, 
                                      number=number, called.by=called.by,
                                      stringsAsFactors=FALSE)
    return(h)
  }
  assign("history", updatehist(get("history", envir=anchorenv), match.call()),
         envir=anchorenv)  
  return(number)
}
