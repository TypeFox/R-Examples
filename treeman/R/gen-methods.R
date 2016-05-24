#' @name randTree
#' @title Generate a random tree
#' @description Returns a random \code{TreeMan} tree with \code{n}
#' tips.
#' @details Equivalent to \code{ape}'s \code{rtree()} but returns a
#' \code{TreeMan} tree. Tree is always rooted and bifurcating.
#' @param n number of tips, integer, must be greater than 1
#' @seealso
#' \code{\link{TreeMan-class}}
#' @export
#' @examples
#' library(treeman)
#' tree <- randTree(5)

randTree <- function(n) {
  .nd <- function(n_left, id, spn, pre, predist, ndlst) {
    post <- kids <- c()
    pd <- 0
    nd <- list('id'=id,
                  'spn'=spn,
                  'prid'=pre,
                  'ptid'=post,
                  'kids'=kids,
                  'pd'=pd,
                  'prdst'=predist)
    ndlst[[id]] <- nd
    # if there are enough ns left to have kids
    n_left <- n_left - 1
    if(n_left > 1) {
      nls <- n_left
      if(n_left == 2) {
        nls <- c(1, 1)
      } else {
        # split must be binary
        splits <- seq(from=1, to=((n_left)-1)/2, by=2)
        nls[2] <- sample(splits, 1)
        nls[1] <- n_left - nls[2]
      }
      for(nl in nls) {
        if(nl == 1) {
          ntips <<- ntips + 1
          new_id <- paste0('t', ntips)
          kids <- c(new_id, kids)
        } else {
          nnds <<- nnds + 1
          new_id <- paste0('n', nnds)
        }
        post <- c(post, new_id)
        new_spn <- runif(min=0, max=1, n=1)
        new_pre <- c(id, pre)
        new_predist <- predist + new_spn
        ndlst <- .nd(nl, new_id, new_spn,
                           new_pre, new_predist, ndlst)
        kids <- c(kids, ndlst[[new_id]][['kids']])
        pd <- pd + new_spn + ndlst[[new_id]][['pd']]
      }
    }
    ndlst[[id]][['kids']] <- kids
    ndlst[[id]][['ptid']] <- post
    ndlst[[id]][['pd']] <- pd
    ndlst
  }
  if(n < 2) {
    stop('`n` must be greater than 1')
  }
  # init empty ndlst
  ntips <- 0
  nnds <- 1
  ndlst <- list()
  n_left <-(n - 1) + n
  # generate root node
  id <- paste0('n', 1)
  predist <- spn <- 0
  pre <- c()
  # gen nodelist
  ndlst <- .nd(n_left, id, spn, pre, predist, ndlst)
  # init new tree object
  tree <- new('TreeMan', ndlst=ndlst, root='n1')
  .updateTreeSlts(tree)
}

blncdTree <- function(...) {
  cat('This function is in progress.... \n')
}

imblncdTree <- function(...) {
  cat('This function is in progress.... \n')
}

#@name blncdTree
# blncdTree <- function(n) {
#   mkNd <- function(id, prid, ptid, pd) {
#     nd <- list('id'=id,
#                   'spn'=1,
#                   'prid'=pre,
#                   'ptid'=post,
#                   'kids'=kids,
#                   'pd'=pd,
#                   'prdst'=prdst)
#     nd
#   }
#   prdst <- n/2
#   is <- seq(from=2, to=n, by=2)
#   cntr <- 0
#   ndlst <- list()
#   for(i in is) {
#     id <- paste0()
#     ndlst <- mkNd
#   }
#   
#   
# }