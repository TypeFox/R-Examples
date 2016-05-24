.checkTreeMen <- function(object) {
  .check <- function(i, invlds) {
    if(class(object@treelist[[i]])[1] != "TreeMan") {
      invlds <<- c(i, invlds)
    }
    NULL
  }
  invlds <- NULL
  sapply(1:length(object@treelist), .check, invlds=invlds)
  if(length(invlds) > 0) {
    for(i in invlds) {
      cat("[", i, "] in treelist is not a TreeMan object\n", sep="")
    }
    return(FALSE)
  }
  TRUE
}

#' @name TreeMen-class
#' @title TreeMen-class
#' @param x \code{TreeMen} object
#' @param i tree index (integer or character)
#' @param object \code{TreeMen} object
#' @param max.level \code{str()} maximum level
#' @param ... additional tree objects
#' @description S4 class for multiple phylogenetic trees
#' @slot treelist list of \code{TreeMan} objects
#' @slot ntips sum of tips per tree
#' @slot ntrees total number of trees
#' @exportClass TreeMen
#' @seealso 
#' \code{\link{cTrees}}
setClass('TreeMen', representation=representation(
  treelist='list',       # list of TreeMan objects
  ntips='numeric',       # sum of tips per tree
  ntrees='numeric'),     # number of trees in object
  validity=.checkTreeMen)

# concatenate methods
.cMenToMen <- function(treemen_1, treemen_2) {
  treelist <- c(treemen_1@treelist, treemen_2@treelist)
  treemen_1@treelist <- treelist
  treemen_1@ntips <- treemen_1@ntips + treemen_2@ntips
  treemen_1@ntrees <- treemen_1@ntrees + treemen_2@ntrees
  treemen_1
}

.cMenToMan <- function(treemen, treeman) {
  treelist <- c(treemen@treelist, treeman)
  treemen@treelist <- treelist
  treemen@ntips <- treeman@ntips + treemen@ntips
  treemen@ntrees <- treemen@ntrees + 1
  treemen
}

.cMenToAny <- function(treemen, treeobj) {
  if(class(treeobj)[1] == "TreeMan") {
    treemen <- .cMenToMan(treemen, treeobj)
  } else if (class(treeobj)[1] == "TreeMen") {
    treemen <- .cMenToMen(treemen, treeobj)
  }
  treemen
}

.cTreeObjs <- function(treemen, treeobj, ...) {
  if(nargs() > 2) {
    treemen <- .cMenToAny(treemen, treeobj)
    treemen <- .cTreeObjs(treemen, ...)
  } else {
    treemen <- .cMenToAny(treemen, treeobj)
  }
  treemen
}

#' @title cTrees
#' @description Return \code{TreeMen} of concatenated trees.
#' @details Concatenate trees into single \code{TreeMen} object.
#' @param x \code{TreeMan} or \code{TreeMen} objects
#' @param ... more \code{TreeMan} or \code{TreeMen} objects
#' @seealso 
#' \code{\link{TreeMen-class}}, \code{\link{TreeMan-class}}, \code{\link{as-TreeMen}}
#' @examples 
#' library(treeman)
#' trees <- cTrees(randTree(10), randTree(10))
#TODO: fix this for when just two trees cTrees(tree, tree)
#' @export
setGeneric("cTrees", signature=c("x"),
           function(x, ...) {
             standardGeneric("cTrees")
           })
#' @rdname TreeMan-class
#' @aliases TreeMan-method
#' @exportMethod cTrees
setMethod("cTrees", c("TreeMan"),
          function(x, ...) {
            x <- list(x)
            x <- as(x, 'TreeMen')
            x <- .cTreeObjs(x, ...)
            x
            })
#' @rdname TreeMen-class
#' @aliases TreeMen-method
#' @exportMethod cTrees
setMethod('cTrees', c('TreeMen'),
          function(x, ...) {
            x <- .cTreeObjs(x, ...)
            x
          })

# Accessor methods
#' @rdname TreeMen-class
#' @aliases TreeMen-method
#' @exportMethod [[
setMethod('[[', c('TreeMen', 'ANY'),
          function(x, i) {
            if(!i %in% 1:length(x@treelist)) {
              stop(paste0(i, ' not in tree'))
            }
            x@treelist[[i]]
          })
#' @rdname TreeMen-class
#' @aliases TreeMen-method
#' @exportMethod [
setMethod('[', c('TreeMen', 'character'),
          function(x, i) {
            if(!i %in% slotNames(x)) {
              stop(paste0(i, '  not in tree'))
            }
            slot(x, i)
          })

#' @name as-TreeMen
#' @title Convert list to a TreeMen
#' @description Return a \code{TreeMen} object from a list of \code{TreeMans}
#' @seealso 
#' \code{\link{TreeMen-class}}
#' @examples 
#' library(treeman)
#' trees <- list('tree_1'=randTree(10), 'tree_2'=randTree(10))
#' trees <- as(trees, 'TreeMen')
# Conversion method
setAs(from="list", to="TreeMen", def=function(from, to) {
  ntips <- sum(unlist(lapply(from, function(tree) tree@ntips)))
  ntrees <- length(from)
  new(to, treelist=from, ntips=ntips, ntrees=ntrees)
  })

# display methods
#' @rdname TreeMen-class
#' @aliases TreeMen-method
#' @exportMethod as.character
setMethod('as.character', c('x'='TreeMen'),
          function(x) {
            paste0('TreeMen Object of [', x@ntrees,'] trees')
          })
#' @rdname TreeMen-class
#' @aliases TreeMen-method
#' @exportMethod show
setMethod('show', 'TreeMen',
          function(object){
            msg <- as.character(object)
            cat(msg)
          })
#' @rdname TreeMen-class
#' @aliases TreeMen-method
#' @exportMethod str
setMethod('str', c('object'='TreeMen'),
          function(object, max.level=2L, ...) {
            if(is.na(max.level)) {
              stop('max.level must be numeric')
            }
            str@default(object, max.level=max.level, ...)
          })
#' @rdname TreeMen-class
#' @aliases TreeMen-method
#' @exportMethod print
setMethod('print', c('x'='TreeMen'),
          function(x){
            msg <- 'Trees (TreeMen Object):\n'
            msg <- paste0(msg, '  + ', x@ntrees, ' trees\n')
            msg <- paste0(msg, '  + ', x@ntips, ' tips\n')
            cat(msg)
          })