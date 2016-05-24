
# TODO: need to rethink this, make it more logical
# -- update tips or nodes
# -- update upstream or downstream or all
# -- update prdst, pd, kids or all

.updateNdsSlt <- function(ndlst, nids, updater) {
  # update nids using updater function
  ndlst[nids] <- plyr::llply(ndlst[nids], .fun=updater)
  ndlst
}

.dwndateNd <- function(ndlst, nid, rid) {
  .add <- function(nd) {
    nd[['pd']] <- nd[['pd']] - nd_spn
    nd[['kids']] <- nd[['kids']][nd[['kids']] != nid]
    nd
  }
  nd_spn <- ndlst[[nid]][['spn']]
  prids <- ndlst[[nid]][['prid']]
  ndlst[prids] <- plyr::llply(ndlst[prids], .fun=.add)
  ndlst
}

.updateNd <- function(ndlst, nid, rid) {
  .add <- function(nd) {
    nd[['pd']] <- nd[['pd']] + nd_spn
    prdst <<- nd[['spn']] + prdst
    nd
  }
  nd_spn <- ndlst[[nid]][['spn']]
  prdst <- nd_spn
  prids <- ndlst[[nid]][['prid']]
  ndlst[prids] <- plyr::llply(ndlst[prids], .add)
  ndlst[[nid]][['prdst']] <- prdst
  ndlst
}

.dwndateTip <- function(ndlst, tid, rid) {
  .add <- function(nd) {
    kids <- nd[['kids']]
    nd[['kids']] <- kids[kids != tid]
    nd[['pd']] <- nd[['pd']] - tp_spn
    nd
  }
  tp_spn <- ndlst[[tid]][['spn']]
  prids <- ndlst[[tid]][['prid']]
  ndlst[prids] <- plyr::llply(ndlst[prids], .fun=.add)
  ndlst
}

.updateTip <- function(ndlst, tid, rid) {
  .add <- function(nd) {
    kids <- nd[['kids']]
    nd[['kids']] <- c(kids, tid)
    nd[['pd']] <- nd[['pd']] + tp_spn
    prdst <<- nd[['spn']] + prdst
    nd
  }
  tp_spn <- ndlst[[tid]][['spn']]
  prids <- ndlst[[tid]][['prid']]
  prdst <- tp_spn
  ndlst[prids] <- plyr::llply(ndlst[prids], .add)
  ndlst[[tid]][['prdst']] <- prdst
  ndlst
}

.globalUpdateAll <- function(ndlst, just_spn_data=FALSE) {
  tip <- function(tid) {
    ndlst <- .updateTip(ndlst, tid, rid)
    ndlst <<- ndlst
  }
  nd <- function(nid) {
    ndlst <- .updateNd(ndlst, nid, rid)
    ndlst <<- ndlst
  }
  wo_prnds <- sapply(ndlst, function(n) length(n[['prid']]) == 0)
  if(!just_spn_data) {
    wo_pstnds <- sapply(ndlst, function(n) length(n[['ptid']]) == 0)
    nids <- names(ndlst)[(!wo_pstnds) & (!wo_prnds)]
    tids <- names(ndlst)[wo_pstnds]
    rid <- names(ndlst)[wo_prnds]
    l_data <- data.frame(tid=tids, stringsAsFactors=FALSE)
    plyr::m_ply(.data=l_data, .fun=tip)
  } else {
    # just run updateNd for all nodes if just spn data needs updating
    nids <- names(ndlst)[!wo_prnds]
    rid <- names(ndlst)[wo_prnds]
  }
  l_data <- data.frame(nid=nids, stringsAsFactors=FALSE)
  plyr::m_ply(.data=l_data, .fun=nd)
  ndlst
}

.updateKids <- function(ndlst, tid, rid) {
  .add <- function(nd) {
    kids <- nd[['kids']]
    nd[['kids']] <- c(kids, tid)
    nd
  }
  prids <- ndlst[[tid]][['prid']]
  ndlst[prids] <- plyr::llply(ndlst[prids], .fun=.add)
  ndlst
}

.dwndateKids <- function(ndlst, tid, rid) {
  .add <- function(nd) {
    kids <- nd[['kids']]
    nd[['kids']] <- kids[kids != tid]
    nd
  }
  prids <- ndlst[[tid]][['prid']]
  ndlst[prids] <- plyr::llply(ndlst[prids], .fun=.add)
  ndlst
}

.globalUpdateKids <- function(ndlst) {
  tip <- function(tid) {
    ndlst <- .updateKids(ndlst, tid, rid)
    ndlst <<- ndlst
  }
  wo_pstnds <- sapply(ndlst, function(n) length(n[['ptid']]) == 0)
  w_prnds <- sapply(ndlst, function(n) length(n[['prid']]) == 0)
  tids <- names(ndlst)[wo_pstnds]
  rid <- names(ndlst)[w_prnds]
  l_data <- data.frame(tid=tids, stringsAsFactors=FALSE)
  plyr::m_ply(.data=l_data, .fun=tip)
  ndlst
}

.updateTreeSlts <- function(tree) {
  wo_pstndes <- sapply(tree@ndlst,
                       function(n) length(n[['ptid']]) == 0)
  tree@tips <- sort(names(wo_pstndes)[wo_pstndes])
  tree@ntips <- length(tree@tips)
  tree@nds <- sort(names(wo_pstndes)[!wo_pstndes])
  tree@nnds <- length(tree@nds)
  tree@all <- c(tree@tips, tree@nds)
  tree@nall <- length(tree@all)
  if(length(tree@root) > 0) {
    wspn <- names(tree@ndlst)[names(tree@ndlst) != tree@root]
  } else {
    wspn <- names(tree@ndlst)
  }
  tree@wspn <- all(sapply(tree@ndlst[wspn], function(n) !is.null(n[['spn']])))
  if(tree@wspn) {
    if(length(tree@root) > 0) {
      tree@age <- max(sapply(tree@ndlst[wspn], function(n) n[['prdst']]))
      extant_is <- unlist(sapply(tree@tips, function(i) {
        (tree@age - tree@ndlst[[i]][['prdst']]) <= tree@tol}))
      tree@ext <- names(extant_is)[extant_is]
      tree@exc <- tree@tips[!tree@tips %in% tree@ext]
      tree@ultr <- all(tree@tips %in% tree@ext)
    } else {
      tree@ext <- tree@exc <- vector()
      tree@ultr <- FALSE
      tree@age <- numeric()
    }
    tree@pd <- sum(sapply(tree@ndlst, function(n) n[['spn']]))
  } else {
    tree@age <- tree@pd <- numeric()
    tree@ext <- tree@ext <- vector()
    tree@ultr <- logical()
  }
  tree@ply <- any(sapply(tree@ndlst, function(n) length(n[['ptid']]) > 2))
  initialize(tree)
}