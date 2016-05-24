### Linking Routines ###

##' A utility for creating linking functions that operate in both
##' directions (full duplex).
##'
##' The generated linker function takes two arguments:
##' \code{from_selection} and \code{new_selection}. If
##' \code{new_selection} is specified, \code{new_selection} is mapped
##' from \code{to_data} to \code{from_data}. Otherwise,
##' \code{from_selection} is mapped from \code{from_data} to
##' \code{to_data}.
##' 
##' @title Duplex linking
##' @param delegate The linking function that performs the mapping,
##' such as \code{\link{match_any_linker}}.
##' @param from_data A \code{data.frame} of keys
##' @param to_data A \code{data.frame} of keys
##' @return A two-way linking function as described in the details.
##' @author Michael Lawrence
##' @export
duplex_data_linker <- function(delegate, from_data, to_data = from_data) {
  function(from_selection, new_selection) {
    if (!missing(new_selection))
      delegate(new_selection, to_data, from_data)
    else delegate(from_selection, from_data, to_data)
  }
}

generate_key <- function(data) {
  if (ncol(data) == 1L)
    data[[1]]
  else do.call(paste, c(as.list(data), sep = "\r"))
}

##' Linking functions return a logical vector, with the \code{TRUE}
##' elements indicating rows in the data that are linked. 
##'
##' The \code{match_any_linker} function links rows in
##' \code{from_data} to rows in \code{to_data} that share the same
##' key.
##'
##' By convention, a key is defined as the combination of the values
##' in every column of \code{from_data} and \code{to_data}. Thus,
##' \code{from_data} and \code{to_data} should contain only the
##' columns necessary for key generation. They should not be an
##' entire dataset. 
##' 
##' @title match_any_linker
##' @param from_data A \code{data.frame}-like object containing the
##' keys for linking the corresponding rows to rows in \code{to_data}
##' @param to_data A \code{data.frame}-like object containing the keys
##' that will be matched against the keys in \code{from_data}
##' @return a logical vector, indicating which \code{from_data} rows are linked
##' @author Michael Lawrence
##' @export
match_any_linker <- function(from_data, to_data = from_data)
{
  duplex_data_linker(function(from_selection, from_data, to_data) {
    from_logical <- as.logical(from_selection)
    generate_key(to_data) %in% generate_key(from_data)[from_logical]
  }, from_data, to_data)
}

## requires that all records in from_data that map to to_data are selected
match_all_linker <- function(from_data, to_data = from_data) {
  match_weight <- match_weight_linker(from_data, to_data)
  duplex_data_linker(function(from_selection, from_data, to_data) {
    match_weight(from_selection, from_data, to_data) == 1.0
  }, from_data, to_data)
}

## weights by selected prsoportion of records in to_data that map to from_data
match_weight_linker <- function(from_data, to_data = from_data) {
  duplex_data_linker(function(from_selection, from_data, to_data) {
    m <- match(generate_key(from_data), generate_key(to_data))
    tabulate(m[as.logical(from_selection)]) / tabulate(m)
  }, from_data, to_data)
}

## FANCIER STUFF IN PROGRESS

distance_cutoff_linker <- function(from_data, to_data = from_data,
                                   cutoff = 0, method)
{
  
}

distance_weight_linker <- function(from_data, to_data = from_data, method)
{    
  duplex_data_linker(function(from_selection, from_data, to_data) {
    
  })
}

## would help from the 'graph' package
## to do enough cool stuff (p-neighborhood, connected components)
## maybe separate package?
graph_linker <- function(from_data, to_data = from_data, graph, p = 1L) { }


