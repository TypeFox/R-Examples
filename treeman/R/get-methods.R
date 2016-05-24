#' @name getNdSstr
#' @title Get sister id
#' @description Returns the id of the sister(s) of node id given.
#' @details An error is raised if there is no sister (e.g. for the root).
#'  There can be more than one sister if tree is polytomous.
#' @param tree \code{TreeMan} object
#' @param id node id
#' @seealso
#' \code{\link{getNdsSstr}}, 
#' \url{https://github.com/DomBennett/treeman/wiki/get-methods}
#' @export
#' @examples
#' library(treeman)
#' tree <- randTree(10)
#' getNdSstr(tree, id='t1')
getNdSstr <- function(tree, id) {
  prid <- tree@ndlst[[id]][['prid']][[1]]
  ptids <- tree@ndlst[[prid]][['ptid']]
  ptids[ptids != id]
}

#' @name getNdsSstr
#' @title Get sister id
#' @description Returns the ids of the sister(s) of nd ids given.
#' @details An error is raised if there is no sister (e.g. for the root).
#'  There can be more than one sister if tree is polytomous. Parallelizable.
#' @param tree \code{TreeMan} object
#' @param ids nd ids
#' @param ... \code{plyr} arguments
#' @seealso
#' \code{\link{getNdSstr}}, 
#' \url{https://github.com/DomBennett/treeman/wiki/get-methods}
#' @export
#' @examples
#' library(treeman)
#' tree <- randTree(10)
#' getNdsSstr(tree, ids=tree['tips'])
getNdsSstr <- function(tree, ids, ...) {
  l_data <- data.frame(id=ids, stringsAsFactors=FALSE)
  res <- plyr::mdply(.data=l_data, .fun=getNdSstr, tree=tree, ...)
  res[ ,2]
}

#' @name getTxnyms
#' @title Get node id(s) for txonyms
#' @description Returns the node ids of nodes with given taxonyms.
#' @details Returns a \code{list}, parallelizable.
#' @param tree \code{TreeMan} object
#' @param txnyms vector of taxonomic names
#' @param ... \code{plyr} arguments
#' @seealso
#' \code{\link{getNdLng}}, 
#' \code{\link{getNdsLng}}, 
#' \url{https://github.com/DomBennett/treeman/wiki/get-methods}
#' @export
#' @examples
#' library(treeman)
#' data(mammals)
#' homo_ids <- getTxnyms(mammals, txnyms='Homo')
getTxnyms <- function(tree, txnyms, ...) {
  # get node id(s) for taxonyms
  .get <- function(id, txnym, ...) {
    for(t in txnyms) {
      if(t %in% txnym) {
        res[[t]] <<- c(res[[t]], id)
      }
    }
  }
  res <- list()
  plyr::m_ply(tree@ndlst, .fun=.get)
  res
}

#' @name getOtgrp
#' @title Get outgroup
#' @description Return the outgroup based on a tree and a vector of IDs.
#' @details Returns a id, character. If there are multiple possible outgroups, returns NULL.
#' @param tree \code{TreeMan} object
#' @param ids vector of node ids
#' @seealso
#' \url{https://github.com/DomBennett/treeman/wiki/get-methods}
#' @export
#' @examples
#' library(treeman)
#' data(mammals)
#' # orangutan is an outgroup wrt humans and chimps
#' getOtgrp(mammals, ids=c('Homo_sapiens', 'Pan_troglodytes', 'Pongo_pygmaeus'))
getOtgrp <- function(tree, ids) {
  .cntr <- function(id) {
    kids <- tree@ndlst[[id]][['kids']]
    sum(ids %in% kids)
  }
  prnt <- getPrnt(tree, ids)
  ptids <- tree@ndlst[[prnt]][['ptid']]
  cnts <- sapply(ptids, .cntr)
  outnd <- names(cnts)[which.min(cnts)]
  kids <- tree@ndlst[[outnd]][['kids']]
  if(length(kids) == 0) {
    return(outnd)
  }
  outgroup <- ids[ids %in% kids]
  if(length(outgroup) > 1) {
    return(NULL)
  }
  outgroup
}

#' @name getNdSlt
#' @title Get a node slot
#' @description Returns the value of named slot.
#' @details Returned object depends on name, either character, vector or numeric.
#' @param tree \code{TreeMan} object
#' @param slt_nm slot name
#' @param id node id
#' @seealso
#' \code{\link{getNdsSlt}}, 
#' \url{https://github.com/DomBennett/treeman/wiki/get-methods}
#' @export
#' @examples
#' library(treeman)
#' tree <- randTree(10)
#' getNdSlt(tree, slt_nm='spn', id='t1')  # return span of t1
getNdSlt <- function(tree, slt_nm, id) {
  tree@ndlst[[id]][[slt_nm]]
}

#' @name getNdsSlt
#' @title Get a node slot for multiple nodes
#' @description Returns the value of named slot.
#' @details Returned object depends on name, either character, vector or numeric. Parallelizable.
#' @param tree \code{TreeMan} object
#' @param slt_nm slot name
#' @param ids vector of node ids
#' @param ... \code{plyr} arguments
#' @seealso
#' \code{\link{getNdSlt}}, 
#' \url{https://github.com/DomBennett/treeman/wiki/get-methods}
#' @export
#' @examples
#' library(treeman)
#' tree <- randTree(10)
#' getNdsSlt(tree, slt_nm='spn', ids=tree['tips'])  # return spans of all tips
getNdsSlt <- function(tree, slt_nm, ids, ...) {
  .get <- function(i) {
    getNdSlt(tree, slt_nm, ids[i])
  }
  l_data <- data.frame(i=1:length(ids), stringsAsFactors=FALSE)
  res <- plyr::mdply(.data=l_data, .fun=.get, ...)
  res[ ,2]
}

#' @name getNdKids
#' @title Get children IDs
#' @description Return the node ids of all tips that descend from node.
#' @details Returns a vector
#' @param tree \code{TreeMan} object
#' @param id node id
#' @seealso
#' \code{\link{getNdsKids}}, 
#' \url{https://github.com/DomBennett/treeman/wiki/get-methods}
#' @export
#' @examples
#' library(treeman)
#' tree <- randTree(10)
#' # everyone descends from root
#' getNdKids(tree, id=tree['root'])
getNdKids <- function(tree, id) {
  nd <- tree@ndlst[[id]]
  nd[['kids']]
}

#' @name getNdsKids
#' @title Get children IDs for multiple nodes
#' @description Return the node ids of all tips that descend from each node in \code{ids}.
#' @details Returns a list, parallelizable.
#' @param tree \code{TreeMan} object
#' @param ids vector of node ids
#' @param ... \code{plyr} arguments
#' @seealso
#' \code{\link{getNdKids}}, 
#' \url{https://github.com/DomBennett/treeman/wiki/get-methods}
#' @export
#' @examples
#' library(treeman)
#' tree <- randTree(10)
#' getNdsKids(tree, id=tree['nds'])

getNdsKids <- function(tree, ids, ...) {
  l_data <- data.frame(id=ids, stringsAsFactors=FALSE)
  res <- plyr::mlply(.data=l_data, .fun=getNdKids, tree=tree, ...)
  names(res) <- ids
  res[1:length(res)]
}

#' @name getNdAge
#' @title Get age
#' @description Return the root to tip distance for \code{id}.
#' @details Returns a numeric.
#' @param tree \code{TreeMan} object
#' @param id node id
#' @seealso
#' \code{\link{getNdsAge}}, 
#' \code{\link{getSpnAge}}, 
#' \code{\link{getSpnsAge}}, 
#' \code{\link{getPrnt}}, 
#' \url{https://github.com/DomBennett/treeman/wiki/get-methods}
#' @export
#' @examples
#' library(treeman)
#' data(mammals)
#' # when did apes emerge?
#' # get parent id for all apes
#' prnt_id <- getPrnt(mammals, ids=c('Homo_sapiens', 'Hylobates_concolor'))
#' getNdAge(mammals, id=prnt_id)
#TODO: how to effectively handle unrooted trees, age has no meaning
getNdAge <- function(tree, id) {
  nd <- tree@ndlst[[id]]
  age <- tree@age - nd[['prdst']]
  age
}

#' @name getNdsAge
#' @title Get ages for multiple nodes
#' @description Return the root to tip distances for \code{ids}.
#' @details Returns a vector, parallelizable.
#' @param tree \code{TreeMan} object
#' @param ids vector of node ids
#' @param ... \code{plyr} arguments
#' @seealso
#' \code{\link{getNdAge}}, 
#' \code{\link{getSpnAge}}, 
#' \code{\link{getSpnsAge}}, 
#' \url{https://github.com/DomBennett/treeman/wiki/get-methods}
#' @export
#' @examples
#' library(treeman)
#' tree <- randTree(10)
#' getNdsAge(tree, ids=tree['nds'])
getNdsAge <- function(tree, ids, ...) {
  l_data <- data.frame(id=ids, stringsAsFactors=FALSE)
  res <- plyr::mdply(.data=l_data, .fun=getNdAge, tree=tree, ...)
  res[ ,2]
}

#' @name getSpnAge
#' @title Get age range
#' @description Return start and end root to tip distances for \code{id}.
#' @details Returns a dataframe.
#' @param tree \code{TreeMan} object
#' @param id node id
#' @seealso
#' \code{\link{getNdAge}}, 
#' \code{\link{getNdsAge}}, 
#' \code{\link{getSpnsAge}}, 
#' \url{https://github.com/DomBennett/treeman/wiki/get-methods}
#' @export
#' @examples
#' library(treeman)
#' data(mammals)
#' getSpnAge(mammals, id='Homo_sapiens')
getSpnAge <- function(tree, id) {
  start <- getNdAge(tree, tree@ndlst[[id]][['prid']][1])
  end <- getNdAge(tree, id)
  data.frame(spn=id, start, end)
}

#' @name getSpnsAge
#' @title Get age ranges for multiple nodes
#' @description Return start and end root to tip distances for \code{ids}.
#' @details Returns a dataframe, parallelizable.
#' @param tree \code{TreeMan} object
#' @param ids vector of node ids
#' @param ... \code{plyr} arguments
#' @seealso
#' \code{\link{getNdAge}}, 
#' \code{\link{getNdsAge}}, 
#' \code{\link{getSpnAge}}, 
#' \url{https://github.com/DomBennett/treeman/wiki/get-methods}
#' @export
#' @examples
#' library(treeman)
#' tree <- randTree(10)
#' # all nodes but root
#' ids <- tree['nds'][tree['nds'] != tree['root']]
#' getSpnsAge(tree, ids=ids)
getSpnsAge <- function(tree, ids, ...) {
  l_data <- data.frame(id=ids, stringsAsFactors=FALSE)
  res <- plyr::mdply(.data=l_data, .fun=getSpnAge, tree=tree, ...)
  res <- res[ ,colnames(res) != 'id']
  res
}

#' @name getPrnt
#' @title Get parent
#' @description Return parental (most recent common ancestor) node id for \code{ids}.
#' @details Returns a character.
#' @param tree \code{TreeMan} object
#' @param ids vector of node ids
#' @seealso
#' \code{\link{getSubtree}}, 
#' \url{https://github.com/DomBennett/treeman/wiki/get-methods}
#' @export
#' @examples
#' library(treeman)
#' data(mammals)
#' # choosing ids from the two main branches of apes allows to find the parent for all apes
#' ape_id <- getPrnt(mammals, ids=c('Homo_sapiens', 'Hylobates_concolor'))
getPrnt <- function(tree, ids) {
  prids <- getNdsPrid(tree, ids)
  rf <- prids[[1]]
  mn_rnk <- 0
  for(n in prids[-1]) {
    rnk <- min(match(n, rf), na.rm=TRUE)
    if(rnk > mn_rnk) mn_rnk <- rnk
  }
  rf[mn_rnk]
}

#' @name getPath
#' @title Get path between nodes
#' @description Return node ids for connecting \code{from} to \code{to}.
#' @details Returns a vector, first id is \code{from} to \code{to}.
#' @param tree \code{TreeMan} object
#' @param from starting node id
#' @param to ending node id
#' @seealso
#' \url{https://github.com/DomBennett/treeman/wiki/get-methods}
#' @export
#' @examples
#' library(treeman)
#' data(mammals)
#' # what's the phylogenetic distance from humans to gorillas?
#' ape_id <- getPrnt(mammals, ids=c('Homo_sapiens', 'Hylobates_concolor'))
#' pth <- getPath(mammals, from='Homo_sapiens', to='Gorilla_gorilla')
#' sum(getNdsSlt(mammals, ids=pth, slt_nm='spn'))
getPath <- function(tree, from, to) {
  pre_1 <- c(from, getNdPrid(tree, from))
  pre_2 <- c(to, getNdPrid(tree, to))
  parent <- pre_1[which(pre_1 %in% pre_2)[1]]
  path_1 <- pre_1[!pre_1 %in% pre_2]
  path_2 <- pre_2[!pre_2 %in% pre_1]
  path_2 <- path_2[length(path_2):1]
  c(path_1, parent, path_2)
}

#' @name getNdPrid
#' @title Get pre-nodes to root
#' @description Return node ids for connecting \code{id} to root.
#' @details Returns a vector.
#' @param tree \code{TreeMan} object
#' @param id node id
#' @seealso
#' \code{\link{getNdsPrid}}, 
#' \code{\link{getNdPtid}}, 
#' \code{\link{getNdsPtid}}, 
#' \url{https://github.com/DomBennett/treeman/wiki/get-methods}
#' @export
#' @examples
#' library(treeman)
#' tree <- randTree(10)
#' # get all nodes to root
#' getNdPrid(tree, id='t1')
getNdPrid <- function(tree, id) {
  tree@ndlst[[id]][['prid']]
}

#' @name getNdsPrid
#' @title Get pre-nodes to root for multiple nodes
#' @description Return node ids for connecting \code{ids} to root.
#' @details Returns a list, parallizable.
#' @param tree \code{TreeMan} object
#' @param ids vector of node ids
#' @param ... \code{plyr} arguments
#' @seealso
#' \code{\link{getNdPrid}}, 
#' \code{\link{getNdPtid}}, 
#' \code{\link{getNdsPtid}}, 
#' \url{https://github.com/DomBennett/treeman/wiki/get-methods}
#' @export
#' @examples
#' library(treeman)
#' tree <- randTree(10)
#' # get all nodes to root
#' getNdsPrid(tree, ids=tree['tips'])
getNdsPrid <- function(tree, ids, ...) {
  l_data <- data.frame(id=ids, stringsAsFactors=FALSE)
  res <- plyr::mlply(.data=l_data, .fun=getNdPrid, tree=tree, ...)
  names(res) <- ids
  res[1:length(res)]
}

#' @name getNdPtid
#' @title Get post-nodes to tips
#' @description Return node ids for connecting \code{id} to kids.
#' @details Returns a vector.
#' @param tree \code{TreeMan} object
#' @param id node id
#' @seealso
#' \code{\link{getNdsPtid}}, 
#' \code{\link{getNdPrid}}, 
#' \code{\link{getNdsPrid}}, 
#' \url{https://github.com/DomBennett/treeman/wiki/get-methods}
#' @export
#' @examples
#' library(treeman)
#' tree <- randTree(10)
#' # get all nodes from root to tip
#' getNdPtid(tree, id='n1')
# reduce dependence on the recursive, by getting prenodes
# tip ids to id
getNdPtid <- function(tree, id) {
  .get <- function(id) {
    tmp <- c(id, getNdPrid(tree, id))
    index <- seq(from=1, to=(which(tmp %in% pstids)[1]-1), by=1)
    pstids <<- c(tmp[index], pstids)
    NULL
  }
  pstids <- id
  l_data <- data.frame(id=tree@ndlst[[id]][['kids']],
                       stringsAsFactors=FALSE)
  plyr::m_ply(.data=l_data, .fun=.get)
  pstids
}

#' @name getNdsPtid
#' @title Get post-nodes to tips for multiple nodes
#' @description Return node ids for connecting \code{ids} to kids.
#' @details Returns a list, parallizable.
#' @param tree \code{TreeMan} object
#' @param ids vector of node ids
#' @param ... \code{plyr} arguments
#' @seealso
#' \code{\link{getNdPtid}}, 
#' \code{\link{getNdPrid}}, 
#' \code{\link{getNdsPrid}}, 
#' \url{https://github.com/DomBennett/treeman/wiki/get-methods}
#' @export
#' @examples
#' library(treeman)
#' tree <- randTree(10)
#' # get all nodes to tip for all nodes
#' getNdsPtid(tree, ids=tree['nds'])
getNdsPtid <- function(tree, ids, ...) {
  l_data <- data.frame(id=ids, stringsAsFactors=FALSE)
  res <- plyr::mlply(.data=l_data, .fun=getNdPtid, tree=tree, ...)
  names(res) <- ids
  res[1:length(res)]
}

#' @name getNdLng
#' @title Get lineage
#' @description Return unique taxonyms for connecting \code{id} to root.
#' @details Returns a vector.
#' @param tree \code{TreeMan} object
#' @param id node id
#' @seealso
#' \code{\link{getNdsLng}}, 
#' \url{https://github.com/DomBennett/treeman/wiki/get-methods}
#' @export
#' @examples
#' library(treeman)
#' data(mammals)
#' # return human lineage
#' getNdLng(mammals, id='Homo_sapiens')
getNdLng <- function(tree, id) {
  .get <- function(txnym, ...) {
    lng <<- c(txnym, lng)
  }
  prids <- c(id, getNdPrid(tree, id))
  lng <- NULL
  plyr::m_ply(tree@ndlst[prids], .fun=.get)
  unique(lng)
}

#' @name getNdsLng
#' @title Get lineage for multiple nodes
#' @description Return unique taxonyms for connecting \code{ids} to root.
#' @details Returns a list, parallelizable.
#' @param tree \code{TreeMan} object
#' @param ids vector of node ids
#' @param ... \code{plyr} arguments
#' @seealso
#' \code{\link{getNdLng}}, 
#' \url{https://github.com/DomBennett/treeman/wiki/get-methods}
#' @export
#' @examples
#' library(treeman)
#' data(mammals)
#' # return human and gorilla lineages
#' getNdLng(mammals, id=c('Homo_sapiens', 'Gorilla_gorilla'))
getNdsLng <- function(tree, ids, ...) {
  l_data <- data.frame(id=ids, stringsAsFactors=FALSE)
  plyr::mlply(.data=l_data, .fun=getNdLng, tree=tree, ...)
}

#' @name getSubtree
#' @title Get subtree
#' @description Return tree descending from \code{id}.
#' @details Returns a \code{TreeMan}, parallelizable.
#' @param tree \code{TreeMan} object
#' @param id node id
#' @seealso
#' \code{\link{getPrnt}}, 
#' \url{https://github.com/DomBennett/treeman/wiki/get-methods}
#' @export
#' @examples
#' library(treeman)
#' data(mammals)
#' # get tree of apes
#' ape_id <- getPrnt(mammals, ids=c('Homo_sapiens', 'Hylobates_concolor'))
#' apes <- getSubtree(mammals, id=ape_id)
getSubtree <- function(tree, id) {
  .prdst <- function(nd) {
    nd[['prdst']] <- nd[['prdst']] - nd_prdst
    nd[['prid']] <- nd[['prid']][nd[['prid']] %in% nids]
    nd
  }
  pstids <- getNdPtid(tree, id)
  ndlst <- tree@ndlst[pstids]
  nids <- names(ndlst)
  nd_prdst <- ndlst[[id]][['prdst']]
  ndlst <- plyr::llply(.data=ndlst, .fun=.prdst)
  ndlst[[id]][['prid']] <- NULL
  ndlst[[id]][['spn']] <- 0
  new_tree <- new('TreeMan', ndlst=ndlst, root=id)
  .updateTreeSlts(new_tree)
}