.newNd <- function(tree, nd) {
  nd <- tree@ndlst[[nd]]
  if(is.null(nd[['spn']]) | !tree@wspn) {
    spn <- pd <- prdst <- numeric()
  } else {
    spn <- nd[['spn']]
    pd <- nd[['pd']]
    prdst <- nd[['prdst']]
  }
  if(length(tree@age) > 0) {
    age <- tree@age - nd[['prdst']]
  } else {
    age <- numeric()
  }
  if(is.null(nd[['txnym']])) {
    txnym <- vector()
  } else {
    txnym <- nd[['txnym']]
  }
  new('Node', id=nd[['id']], spn=spn, prid=as.character(nd[['prid']][1]),
     ptid=as.character(nd[['ptid']]), kids=as.character(nd[['kids']]),
     nkids=length(as.character(nd[['kids']])), pd=pd, txnym=txnym,
     prdst=prdst, root=tree@root == nd[['id']],
     age=age, tip=length(nd[['ptid']]) == 0)
}

#' @name Node-class
#' @param x \code{Node} object
#' @param object \code{Node} object
#' @param i slot name
#' @title Node-class
#' @description The \code{Node} is an S4 class used for displaying node information.
#' It is only generated when a user implements the \code{[[]]} on a tree.
#' @slot id unique ID for node in tree['ndlst']
#' @slot spn length of preceding branch
#' @slot prid parent node ID
#' @slot ptid child node ID
#' @slot kids descending tip IDs
#' @slot nkids number of descending tip IDs
#' @slot txnym list of associated taxonyms
#' @slot pd total branch length represented by node
#' @slot prdst total branch length of connected prids
#' @slot age age of node in tree
#' @slot root T/F root node?
#' @slot tip T/F tip node?
#' @exportClass Node
#' @seealso 
#' \code{\link{cTrees}}
setClass ('Node', representation=representation (
  id='character',        # unique ID for node in tree@nodelist
  spn='numeric',        # length of preceding branch
  prid='character',      # parent node ID
  ptid='vector',         # child node IDs
  kids='vector',         # descending tip IDs
  nkids='numeric',       # number of descending tips
  txnym="vector",        # list of associated taxonyms
  pd='numeric',          # total branch length represented by node
  prdst='numeric',       # total branch length of connected pres
  age='numeric',         # age of node in tree
  root='logical',        # T/F root node?
  tip='logical')         # T/F tip node?
)

#' @rdname Node-class
#' @aliases Node-method
#' @exportMethod as.character
setMethod ('as.character', c('x'='Node'),
           function(x) {
             x@id
           })
#' @rdname Node-class
#' @aliases Node-method
#' @exportMethod show
setMethod ('show', 'Node',
           function(object){
             print (object)
           })
#' @rdname Node-class
#' @aliases Node-method
#' @exportMethod print
setMethod ('print', c('x'='Node'),
           function(x){
             if(x@root) {
               msg <- paste0('Node (root node):\n')
             } else if (x@tip){
               msg <- paste0('Node (tip node):\n')
             } else {
               msg <- paste0('Node (internal node):\n')
             }
             msg <- paste0(msg, '  + ID: \"', x@id, '\"\n')
             if(length(x@txnym) > 0) {
               msg <- paste0(msg, '  + txnym: \"', paste0(x@txnym, collapse='\", \"'), '\"\n')
             }
             if(!x@root) {
               msg <- paste0(msg, '  + prid: \"', x@prid, '\"\n')
             }
             if(!x@tip) {
               msg <- paste0(msg, '  + ptid: \"', paste0(x@ptid, collapse='\", \"'), '\"\n')
               msg <- paste0(msg, '  + nkids: ', length(x@kids), '\n')
             }
             if(length(x@spn) > 0) {
               if(!x@root) {
                 msg <- paste0(msg, '  + spn: ', signif(x@spn, 2), '\n')
               }
               if(length(x@age) > 0) {
                 msg <- paste0(msg, '  + age: ', signif(x@age, 2), '\n')
               } else {
                 msg <- paste0(msg, '  + predist: ', signif(x@prdst, 2), '\n') 
               }
               msg <- paste0(msg, '  + pd: ', signif(x@pd, 2), '\n')
             }
             cat (msg)
           })
#' @rdname Node-class
#' @aliases Node-method
#' @exportMethod [
setMethod('[', c('Node', 'character'),
          function(x, i) {
            if(!i %in% slotNames(x)) {
              stop(paste0(i, '  not in node'))
            }
            slot(x, i)
          })