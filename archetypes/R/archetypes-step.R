#' @include archetypes-class.R
{}


#' Run archetypes algorithm repeatedly
#'
#' @param ... Passed to the specific archetype function.
#' @param k A vector of integers passed in turn to the k argument of
#'   \code{\link{archetypes}}.
#' @param nrep For each value of \code{k} run \code{\link{archetypes}}
#'   \code{nrep} times.
#' @param method Archetypes function to use, typically
#'   \code{\link{archetypes}}, \code{\link{weightedArchetypes}} or
#'   \code{\link{robustArchetypes}},
#' @param verbose Show progress during exection.
#'
#' @return A list with \code{k} elements and class attribute
#'   \code{stepArchetypes}. Each element is a list of class
#'   \code{repArchetypes} with \code{nrep} elements; only for internal
#'   usage.
#'
#' @seealso \code{\link{archetypes}}
#'
#' @examples
#'   \dontrun{
#'   data(skel)
#'   skel2 <- subset(skel, select=-Gender)
#'   as <- stepArchetypes(skel2, k=1:5, verbose=FALSE)
#'
#'   ## Residual sum of squares curve:
#'   screeplot(as)
#'
#'   ## Select three archetypes and from that the best
#'   ## recurrence:
#'   a3 <- bestModel(as[[3]])
#'   }
#'
#' @export
stepArchetypes <- function(..., k, nrep = 3, method = archetypes, verbose = TRUE) {

  mycall <- match.call()
  as <- list()

  for ( i in 1:length(k) ) {
    as[[i]] <- list()
    class(as[[i]]) <- 'repArchetypes'

    for ( j in seq_len(nrep) ) {
      if ( verbose )
        cat('\n*** k=', k[i], ', rep=', j, ':\n', sep='')

      as[[i]][[j]] <- method(..., k=k[i])
    }
  }

  return(structure(as, class='stepArchetypes', call=mycall))
}


#' @import methods
setOldClass('repArchetypes')
setOldClass('stepArchetypes')



#' Extract method
#'
#' An extraction on a \code{stepArchetypes} object returns again a
#' \code{stepArchetypes} object.
#'
#' @param x A \code{stepArchetypes} object.
#' @param i The indizes to extract.
#' @return A \code{stepArchetypes} object containing only the parts
#'   defined in \code{i}.
#' @method [ stepArchetypes
#' @rdname extract
#'
#' @S3method "[" stepArchetypes
`[.stepArchetypes` <- function(x, i) {
  y <- unclass(x)[i]
  attributes(y) <- attributes(x)

  return(y)
}



#' @S3method print stepArchetypes
print.stepArchetypes <- function(x, ...) {
  cat('StepArchetypes object\n\n')
  cat(deparse(attr(x, 'call')), '\n')
}



#' Summary method for stepArchetypes object
#'
#' @param object A \code{stepArchetypes} object.
#' @param ... Ignored.
#' @return Undefined.
#'
#' @method summary stepArchetypes
#' @rdname summary
#'
#' @S3method summary stepArchetypes
summary.stepArchetypes <- function(object, ...) {
  print(object)

  ps <- nparameters(object)

  for ( i in seq_along(object) ) {
    cat('\nk=', ps[i], ':\n', sep='')
    print(object[[i]], full=FALSE)
  }
}



#' @rdname parameters
#' @aliases parameters,stepArchetypes-method
#' @importFrom modeltools parameters
#' @import methods
#' @exportMethod parameters
setMethod('parameters', signature = c(object = 'stepArchetypes'),
function(object, ...) {
  lapply(object, parameters)
})



#' @rdname nparameters
#' @method nparameters stepArchetypes
#'
#' @S3method nparameters stepArchetypes
nparameters.stepArchetypes <- function(object, ...) {
  return(sapply(object, nparameters))
}



#' @rdname rss
#' @method rss stepArchetypes
#'
#' @S3method rss stepArchetypes
rss.stepArchetypes <- function(object, ...) {
  ret <- t(sapply(object, rss))
  rownames(ret) <- paste('k', nparameters(object), sep='')
  return(ret)
}



#' Return best model
#'
#' @param object An \code{archetypes} object.
#' @param ... Ignored
#'
#' @rdname bestModel
#' @method bestModel stepArchetypes
#'
#' @S3method bestModel stepArchetypes
bestModel.stepArchetypes <- function(object, ...) {
  zsmin <- lapply(object, bestModel)

  if ( length(zsmin) == 1 )
    return(zsmin[[1]])
  else
    return(zsmin)
}



#' @S3method print repArchetypes
print.repArchetypes <- function(x, ...) {
  for ( i in seq_along(x) )
    print(x[[i]], ...)

  invisible(x)
}



#' @rdname parameters
#' @aliases parameters,repArchetypes-method
#' @importFrom modeltools parameters
#' @exportMethod parameters
setMethod('parameters', signature = signature(object = 'repArchetypes'),
function(object, ...) {
  lapply(object, parameters)
})

 

#' @rdname rss
#' @method rss repArchetypes
#'
#' @S3method rss repArchetypes
rss.repArchetypes <- function(object, ...) {
  ret <- sapply(object, rss)
  names(ret) <- paste('r', seq_along(ret), sep='')

  return(ret)
}



#' @rdname nparameters
#' @method nparameters repArchetypes
#'
#' @S3method nparameters repArchetypes
nparameters.repArchetypes <- function(object, ...) {
  nparameters(object[[1]])
}



#' @rdname bestModel
#' @method bestModel repArchetypes
#'
#' @S3method bestModel repArchetypes
bestModel.repArchetypes <- function(object, ...) {
  m <- which.min(rss(object))

  if ( length(m) == 0 )
    return(object[[1]])
  else
    return(object[[m]])
}


