#' @name writeTree
#' @title Write a Newick tree
#' @description Creates a Newick tree from a \code{TreeMan} object.
#' @details The \code{ndLabels} argument can be used to add a user defined node label in
#' the Newick tree. It should take only 1 argument, \code{nd}, the node represented as a list.
#' It should only return a single character value that can be added to a newick string.
#' @param tree \code{TreeMan} object
#' @param file file path
#' @param ndLabels node label function
#' @seealso
#' \code{\link{readTree}}, \code{\link{randTree}}, \url{https://en.wikipedia.org/wiki/Newick_format}
#' @export
#' @examples
#' library(treeman)
#' tree <- randTree(10)
#' ndLabels <- function(n) {
#' paste0(n[['id']], '_ndlabel')
#' }
#' writeTree(tree, file='example.tre', ndLabels)
#' file.remove('example.tre')
# TODO: test this with unrooted trees, adapt for TreeMen
writeTree <- function(tree, file, ndLabels=function(nd){
  return(NULL)
  }) {
  tipBytip <- function(i) {
    ids <- c(ndlst[[prid]][['kids']], prid,
             ndlst[[prid]][['prid']])
    id <<- ids[!ids %in% deja_vues][1]
    deja_vues[i] <<- id
    spn <- ndlst[[id]][['spn']]
    if(id %in% tids) {
      dpth <- which(ndlst[[id]][['prid']] == prid) - 1
      prid <<- ndlst[[id]][['prid']][[1]]
      tpstr <- paste0(id, ':', spn)
      if(dpth > 0) {
        brckts <- paste0(rep('(', dpth), collapse='')
        trstr <<- paste0(trstr, ',', brckts, tpstr)
      } else {
        trstr <<- paste0(trstr, ',', tpstr)
      }
    } else {
      prid <<- ndlst[[id]][['prid']][[1]]
      ndlbl <- ndLabels(ndlst[[id]])
      trstr <<- paste0(trstr, ')', ndlbl,':', spn)
    }
    NULL
  }
  # start with first tip
  # loop through tree structure adding tip by tip to string
  # unpack
  ndlst <- tree@ndlst
  tids <- tree@tips
  nids <- tree@nds
  rid <- tree@root
  # add first tip
  id <- tids[1]
  trstr <-  ''
  deja_vues <- rep(NA, length(ndlst))
  deja_vues[1] <- id
  spn <- ndlst[[id]][['spn']]
  dpth <- length(ndlst[[id]][['prid']])
  prid <- ndlst[[id]][['prid']][[1]]
  tpstr <- paste0(id, ':', spn)
  trstr <- paste0(rep('(', dpth), collapse='')
  trstr <- paste0(trstr, tpstr)
  # loop through nodes
  plyr::m_ply(2:(length(ndlst) - 1), .fun=tipBytip)
  trstr <- paste0(trstr, ');')
  write.table(x=trstr, file=file, quote=FALSE, row.names=FALSE,
              col.names=FALSE)
}

#' @name readTree
#' @title Read a Newick tree
#' @description Return a \code{TreeMan} or \code{TreeMen} object from a Newick treefile
#' @details Read a single or multiple trees from a file, or a text string. Parallelizable.
#' @param file file path
#' @param text Newick character string
#' @param ... \code{plyr} arguments
#' @seealso
#' \code{\link{writeTree}}, \code{\link{randTree}}, \url{https://en.wikipedia.org/wiki/Newick_format}
#' @export
#' @examples
#' library(treeman)
#' tree <- readTree(text="((A:1.0,B:1.0):1.0,(C:1.0,D:1.0):1.0);")
readTree <- function(file=NULL, text=NULL, ...) {
  if(!is.null(file)) {
    trstr <- scan(file, what="raw", quiet=TRUE)
  } else {
    trstr <- text
  }
  if(length(trstr) > 1) {
    trees <- plyr::mlply(trstr, .fun=.readTree, ...)
    tree <- as(trees, 'TreeMen')
  } else {
    tree <- .readTree(trstr)
  }
  tree
}

.readTree <- function(trstr) {
  cuts <- gregexpr("(\\(|\\)|,|;)", trstr)[[1]]
  cuts <- c(cuts[1], cuts[2:length(cuts)] - cuts[1:(length(cuts)-1)])
  rdrenv <- .getRdrEnv(trstr)
  l_data <- data.frame(end_pos=cuts, stringsAsFactors=FALSE)
  plyr::m_ply(l_data, .mkNdLst, rdrenv=rdrenv)
  .addRoot(rdrenv)
  if(length(rdrenv[['root']]) == 0) {
    ndlst <- .globalUpdateKids(rdrenv[['ndlst']])
  } else {
    ndlst <- .globalUpdateAll(rdrenv[['ndlst']])
  }
  tree <- new('TreeMan', ndlst=ndlst, root=rdrenv[['root']])
  .updateTreeSlts(tree)
}

# set-up reader env
.getRdrEnv <- function(trstr) {
  rdrenv <- new.env()
  rdrenv$wspn <- grepl(':', trstr)
  rdrenv$trstr <- trstr
  rdrenv$ndlst <- list()
  rdrenv$prnds <- list()
  rdrenv$cntr <- 0L
  rdrenv$i <- 0L
  rdrenv$nxt_is_intrnl <- FALSE
  rdrenv
}

# extract ID and span from ndstr
.getIDandSpn <- function(ndstr, nints, wspn) {
  nd <- .mkNd(id='', wspn)
  ndstr <- gsub("(\\(|\\)|\\;|,)", "", ndstr)
  ndstr <- strsplit(ndstr, ":")[[1]]
  if(length(ndstr) > 1) {
    nd[['spn']] <- as.numeric(ndstr[2])
  }
  if(length(ndstr) == 0 || ndstr[1] == "") {
    nd[['id']] <- paste0("n", nints)
  } else {
    nd[['id']] <- ndstr[1]
  }
  nd
}

.mkNd <- function(id, wspn) {
  nd <- list()
  nd[['id']] <- id
  if(wspn) {
    nd[['spn']] <- nd[['prdst']] <- nd[['pd']] <- 0
  }
  nd
}

.getPrid <- function(prnds) {
  unlist(lapply(prnds[length(prnds):1], function(nd) nd[['id']]))
}

# cut trstr and generate nodelist
.mkNdLst <- function(end_pos, rdrenv) {
  ndstr <- substr(rdrenv$trstr, 1, end_pos)
  if(grepl("^\\(", ndstr)) {
    rdrenv$cntr <- rdrenv$cntr + 1
    rdrenv$i <- rdrenv$i + 1
    rdrenv$prnds[[rdrenv$i]] <- .mkNd(id=paste0("n", rdrenv$cntr),
                                      wspn=rdrenv$wspn)
  } else {
    nd <- .getIDandSpn(ndstr, rdrenv$cntr, rdrenv$wspn)
    if(rdrenv$nxt_is_intrnl) {
      # TODO: utilise nd$id, e.g. node labels are taxonym or support
      nd$id <- rdrenv$prnds[[rdrenv$i]][['id']]
      nd$ptid <- rdrenv$prnds[[rdrenv$i]][['ptid']]
      rdrenv$prnds <- rdrenv$prnds[-rdrenv$i]
      rdrenv$i <- rdrenv$i - 1
      nd$prid <- .getPrid(rdrenv$prnds)
      rdrenv$ndlst[[nd[['id']]]] <- nd
      rdrenv$nxt_is_intrnl <- FALSE
    } else {
      nd[['prid']] <- .getPrid(rdrenv$prnds)
      rdrenv$ndlst[[nd[['id']]]] <- nd
    }
    if(length(rdrenv$prnds) > 0) {
      rdrenv$prnds[[rdrenv$i]][['ptid']] <-
        c(rdrenv$prnds[[rdrenv$i]][['ptid']], nd[['id']])
    }
    if(grepl("\\)$", ndstr)) {
      rdrenv$nxt_is_intrnl <- TRUE
    }
  }
  rdrenv$trstr <- substr(rdrenv$trstr, end_pos+1, nchar(rdrenv$trstr))
  NULL
}

.addRoot <- function(rdrenv) {
  root_i <- which(unlist(lapply(rdrenv$ndlst, function(n) length(n[['prid']]) == 0)))
  if(length(root_i) > 0) {
    rdrenv$ndlst[[root_i]][['prid']] <- NULL
    rdrenv$ndlst[[root_i]][['spn']] <- 0
    rdrenv$root <- names(rdrenv$ndlst)[root_i]
  } else {
    rdrenv$root <- character()
  }
  NULL
}