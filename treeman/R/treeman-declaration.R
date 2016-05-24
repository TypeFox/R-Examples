# roxygen imports
#' @import methods
#' @importFrom graphics lines plot.default text
#' @importFrom utils combn write.table
#' @importFrom stats runif
.checkTreeMan <- function(object) {
  .check <- function(nd) {
    # must have id
    if(!'id' %in% names(nd)) {
      return(FALSE)
    }
    # must have either prid/ptid or both
    if(!('ptid' %in% names(nd) | 'prid' %in% names(nd))){
      return(FALSE)
    }
    # must have pd if spns
    if(length(nd[['spn']]) > 0 & is.null(nd[['pd']])) {
      return(FALSE)
    }
    # kids must be all unique
    if(any(duplicated(nd[['kids']]))) {
      return(FALSE)
    }
    test_1 <- nd[['id']] %in% nds
    test_2 <- is.null(nd[['prid']]) || (nd[['prid']] %in% nds)
    test_3 <- is.null(nd[['prid']]) || all(nd[['ptid']] %in% nds)
    if(test_1 & test_2 & test_3) {
      return(TRUE)
    }
    FALSE
  }
  nds <- names(object@ndlst)
  nd_checks <- unlist(lapply(object@ndlst, .check))
  if(!all(nd_checks)) {
    msg <- 'These nodes are invalid:\n'
    bad <- which(!nd_checks)
    for(i in bad[-length(bad)]) {
      msg <- paste0(msg, nds[i], ', ')
    }
    msg <- paste0(msg, nds[bad[length(bad)]], '\n\n')
    msg <- paste0(msg, 'They may be pointing to non-existent nodes in tree, their ID may not be a named element in `@ndlst` or they may have missing expected node slots.')
    cat(msg)
    return(FALSE)
  }
  TRUE
}

#' @name TreeMan-class
#' @title TreeMan-class
#' @description S4 class for representing phylogenetic trees as a list of nodes.
#' @param x \code{TreeMan} object
#' @param i node ID or slot name
#' @param object \code{TreeMan} object
#' @param max.level \code{str()} maximum number of levels to show
#' @param ... additional tree objects
#' @slot ndlst list of nodes
#' @slot nds vector of node ids that are internal nodes
#' @slot nnds numeric of number of internal nodes in tree
#' @slot tips vector of node ids that are tips
#' @slot ntips numeric of number of internal nodes in tree
#' @slot all vector of all node ids
#' @slot nall numeric of number of all nodes in tree
#' @slot age numeric of max root to tip distance
#' @slot pd numeric of total branch length of tree
#' @slot ext vector of node ids of all tips with 0 age
#' @slot exc vector of node ids of all tips with age > 0
#' @slot wspn logical, do nodes have spans
#' @slot ultr logical, do all tips end at 0
#' @slot ply logical, is tree bifurcating
#' @slot tol numeric of tolerance for determining extant
#' @slot root character of node id of root, if no root then empty character
#' @details
#' A \code{TreeMan} object holds a list of nodes. The idea of the \code{TreeMan}
#' class is to make adding and removing nodes as similar as possible to adding
#' and removing elements in a list. Note that internal nodes and tips are
#' both considered nodes. Trees can be polytomous but not unrooted.
#' 
#' 
#' Each node within the \code{TreeMan} \code{ndlst} contains the following data slots:
#' \itemize{
#'    \item \code{id}: character string for the node ID
#'    \item \code{txnym}: name of taxonomic clade (optional)
#'    \item \code{spn}: length of the preceding branch
#'    \item \code{prid}: IDs of the preceding nodes to the root
#'    \item \code{ptid}: IDs of the immediately connecting nodes
#'    \item \code{kids}: descending tip IDs
#'    \item \code{pd}: phylogenetic diversity represented by node
#'    \item \code{prdst}: pre distance(distance to root if rooted or
#'    most distal tip if unrooted)
#' }
#' These data slots are updated whenever a node is modified, added or removed.
#' 
#' See below in 'Examples' for these methods in use.
#' @seealso
#' \code{\link{randTree}}, \code{\link{Node-class}}, \code{\link{viz}}
#' @examples
#' library(treeman)
#' # Generate random tree
#' tree <- randTree(10)
#' # Print to get basic stats
#' print(tree)
#' # Currently available methods
#' tree['tips']  # return all tips IDs
#' tree['nds']  # return all internal node IDs
#' tree['ntips']  # count all tips
#' tree['nnds']  # count all internal nodes
#' tree['root']  # identify root node
#' tree[['t1']]  # return t1 node object
#' tree['pd']  # return phylogenetic diversity
#' tree['age']  # return age of tree
#' tree['ultr']  # is ultrametric?
#' tree['ply']  # is polytomous?
#' tree['ext']  # return all extant tip IDs
#' tree['exc']  # return all extinct tip IDs
#' tree <- setTol(tree, 10)  # reset tolerance, default 1e-8
#' # now tol is higher more tips will be classed as extant
#' tree['ext']
#' # Because all nodes are lists with metadata we can readily
#' #  get specific information on nodes of interest
#' nd <- tree[['n2']]
#' print(nd)
#' # And then use the same syntax for the tree
#' nd['age']  # .... nkids, pd, etc.
#' @exportClass TreeMan
setClass('TreeMan', representation=representation(
  ndlst='list',         # list of node lists
  nds='vector',          # vector of node ids that are internal nodes
  nnds='numeric',        # numeric of number of internal nodes in tree
  tips='vector',         # vector of node ids that are tips
  ntips='numeric',       # numeric of number of internal nodes in tree
  all='vector',          # vector of all Node ids
  nall='numeric',        # numeric of number of all nodes in tree
  age='numeric',         # numeric of max root to tip distance
  pd='numeric',          # numeric of total branch length of tree
  ext='vector',          # vector of node ids of all tips with 0 age
  exc='vector',          # vector of node ids of all tips with age > 0
  wspn='logical',        # logical, do nodes have spans
  ultr='logical',        # logical, do all tips end at 0
  ply='logical',         # logical, is tree bifurcating
  tol='numeric',         # numeric of tolerance for determining extant
  root='character'),     # character of node id of root, if no root then empty character
  prototype=prototype(tol=1e-8), validity=.checkTreeMan)

# Accessor methods
#' @rdname TreeMan-class
#' @aliases TreeMan-method
#' @exportMethod [[
setMethod('[[', c('TreeMan', 'character'),
          function(x, i) {
            if(!i %in% names(x@ndlst)) {
              srch_trm <- gsub(' ', '_', i)  # usual mistake
              pssbls <- which(agrepl(srch_trm, names(x@ndlst), ignore.case=TRUE,
                                     max.distance=0.25))
              pssbls <- names(x@ndlst)[pssbls]
              if(length(pssbls) > 0 & length(pssbls) < 50) {
                msg <- paste0("Can't find [", i, "]. Did you mean ....\n")
                for(p in pssbls) {
                  msg <- paste0(msg, '"', p, '"\n')
                }
                msg <- paste0(msg, "?\n")
              } else {
                msg <- paste0("Can't find [", i, "] in tree.")
              }
              stop(msg)
            }
            .newNd(x, i)
          })
#' @rdname TreeMan-class
#' @aliases TreeMan-method
#' @exportMethod [
setMethod('[', c('TreeMan', 'character'),
          function(x, i) {
            slt_nms <- slotNames(x)
            slt_nms <- slt_nms[slt_nms != 'ndlst']
            if(!i %in% slt_nms) {
              slt_nms <- paste0(slt_nms, collapse=', ')
              stop(paste0('`', i, '` not a tree slot. Available slots: ', slt_nms))
            }
            slot(x, i)
          })

# display methods
#' @rdname TreeMan-class
#' @aliases TreeMan-method
#' @exportMethod as.character
setMethod('as.character', c('x'='TreeMan'),
          function(x) {
            paste0('TreeMan Object of [', length(x@tips),'] tips')
          })
#' @rdname TreeMan-class
#' @aliases TreeMan-method
#' @exportMethod show
setMethod('show', 'TreeMan',
          function(object){
            msg <- as.character(object)
            cat(msg)
          })
#' @rdname TreeMan-class
#' @aliases TreeMan-method
#' @exportMethod str
setMethod('str', c('object'='TreeMan'),
          function(object, max.level=2L, ...) {
            if(is.na(max.level)) {
              stop('max.level must be numeric')
            }
            str@default(object, max.level=max.level, ...)
          })
#' @rdname TreeMan-class
#' @aliases TreeMan-method
#' @exportMethod print
setMethod('print', c('x'='TreeMan'),
          function(x){
            msg <- 'Tree (TreeMan Object):\n'
            msg <- paste0(msg, '  + ', x@ntips, ' tips\n')
            msg <- paste0(msg, '  + ', x@nnds, ' internal nodes\n')
            if(x@ply) {
              msg <- paste0(msg, '  + Polytomous\n')
            } else {
              msg <- paste0(msg, '  + Binary\n')
            }
            if(length(x@root) == 0) {
              if(!x@wspn) {
                msg <- paste0(msg, '  + Unrooted and without node spans\n')
              } else {
                msg <- paste0(msg, '  + Unrooted, with node spans\n')
                msg <- paste0(msg, '  + PD ', signif(x@pd, 3), '\n')
              }
            } else {
              if(x@wspn) {
                msg <- paste0(msg, '  + Age ', signif(x@age, 3), '\n')
                msg <- paste0(msg, '  + PD ', signif(x@pd, 3), '\n')
                if(x@ultr) {
                  msg <- paste0(msg, '  + Ultrametric (all tips are extant)\n')
                } else {
                  msg <- paste0(msg, '  + Not ultrametric (with extinct tips)\n')
                }
              } else {
                msg <- paste0(msg, '  + Without node spans\n')
              }
              msg <- paste0(msg, '  + Root node is \"', x@root, '\"\n')
            }
            cat(msg)
          })

#' Method viz
#' @name viz
#' @param tree \code{TreeMan} object
#' @param taxonyms Boolean, show taxonyms rather than IDs?
#' @rdname TreeMan-method
#' @description Crude plot of \code{TreeMan} objects
#' @seealso 
#' \code{\link{TreeMan-class}}
#' @exportMethod viz
setGeneric("viz", signature=c("tree", "taxonyms"),
           function(tree, taxonyms=FALSE) {
             standardGeneric("viz")
           })
#' @rdname TreeMan-method
#' @aliases viz
#' @exportMethod viz
setMethod('viz', 'TreeMan',
          function(tree, taxonyms){
            get_pnts <- function(nd, y, pnts) {
              pstids <- nd[['ptid']]
              low_y_diff <- -nd[['pd']]/2
              high_y_diff <- nd[['pd']]/2
              y_diffs <- seq(from=low_y_diff, to=high_y_diff,
                             length.out=length(pstids))
              counter <- 1
              for(pstid in pstids) {
                pstnd <- tree@ndlst[[pstid]]
                pstnd_x <- pstnd[['prdst']]
                pstnd_y <- y + y_diffs[counter]
                pnts <- rbind(pnts, data.frame(nd=pstid,
                                               x=pstnd_x, y=pstnd_y))
                pnts <- get_pnts(pstnd, pstnd_y, pnts)
                counter <- counter + 1
              }
              pnts
            }
            if(!tree@wspn) {
              # TODO: switch to setNdsspn
              for(i in 1:length(tree@ndlst)) {
                tree@ndlst[[i]][['spn']] <- 1
                tree@ndlst[[i]][['pd']] <- length(tree@ndlst[[i]][['kids']])
                prids <- getNdPrid(tree, tree@ndlst[[i]][['id']])
                tree@ndlst[[i]][['prdst']] <- length(prids)
              }
              tree@pd <- length(tree@ndlst) - 1
            }
            # start with root node
            # TODO: handle unrooted tree
            pnts <- data.frame(nd=tree@root, x=0, y=tree@pd, stringsAsFactors=FALSE)
            root_nd <- tree@ndlst[[tree@root]]
            pnts <- get_pnts(root_nd, y=tree@pd, pnts=pnts)
            # add 10% to min y limit for node label
            min_y <- abs(min(pnts$y))
            min_y <- min_y + (min_y*.1)
            min_y <- ifelse(min(pnts$y) > 0, min_y, -1*min_y)
            y_lmts <- c(min_y, max(pnts$y))
            plot.default(x=pnts$x, y=pnts$y, col='black', pch=19, yaxt='n', ylab='',
                         xlab='', bty='n', ylim=y_lmts)
            if(taxonyms) {
              text(x=pnts$x, y=pnts$y,
                   labels=sapply(pnts$nd, function(n) tree@ndlst[[n]][['taxonym']]),
                   pos=1)
            } else {
              text(x=pnts$x, y=pnts$y, labels=pnts$nd, pos=1)
            }
            # draw lines
            for(i in 2:nrow (pnts)) {
              prend <- tree@ndlst[[pnts$nd[i]]][['prid']][1]
              ind <- c(i, which(pnts$nd == prend))
              lines(x=pnts$x[ind], y=pnts$y[ind])
            }
          })