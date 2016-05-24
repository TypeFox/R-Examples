# TODO: mergeTree, collapseNode, removeNode, reRoot

#' @name rmTip
#' @title Remove tip from a tree
#' @description Returns a tree with a tip ID remove
#' @details Removes a tip in a tree. Set drp_intrnl to FALSE to convert
#' internal nodes into new tips.
#' @param tree \code{TreeMan} object
#' @param tid tip ID
#' @param drp_intrnl Boolean, drop internal branches, default FALSE
#' @seealso
#' \code{\link{addTip}}, 
#' \url{https://github.com/DomBennett/treeman/wiki/manip-methods}
#' @export
#' @examples
#' library(treeman)
#' tree <- randTree(10)
#' tree <- rmTip(tree, 't1')
rmTip <- function(tree, tid, drp_intrnl=TRUE) {
  updater <- function(nd) {
    nd[['prid']] <- nd[['prid']][nd[['prid']] != tid]
    nd[['prid']] <- nd[['prid']][nd[['prid']] != prid]
    nd
  }
  # unpack
  ndlst <- tree@ndlst
  rid <- tree@root
  sids <- getNdSstr(tree, tid)
  # get prid
  prid <- ndlst[[tid]][['prid']][[1]]
  # remove tid
  ndlst <- .dwndateNd(ndlst, nid=tid, rid=rid)
  ndlst <- ndlst[names(ndlst) != tid]
  ndlst[[prid]][['ptid']] <-
    ndlst[[prid]][['ptid']][ndlst[[prid]][['ptid']] != tid]
  # remove prnd if specified and not polytomous
  if(drp_intrnl & length(sids) == 1) {
    nids <- unlist(getNdsPtid(tree, sids))
    ptid <- ndlst[[prid]][['ptid']][[1]]
    ndlst[[ptid]][['spn']] <- ndlst[[prid]][['spn']] +
      ndlst[[ptid]][['spn']]
    if(prid != rid) {
      gprid <- ndlst[[prid]][['prid']][[1]]
      ndlst[[gprid]][['ptid']] <-
        ndlst[[gprid]][['ptid']][ndlst[[gprid]][['ptid']] != prid]
      ndlst[[gprid]][['ptid']] <- c(ndlst[[gprid]][['ptid']], ptid)
    } else {
      # if prid to be dropped is root, set root to ptid
      tree@root <- ptid
    }
    ndlst <- ndlst[names(ndlst) != prid]
    ndlst <- .updateNdsSlt(ndlst, nids, updater)
  } else {
    prid <- NULL
  }
  tree@ndlst <- ndlst
  tree <- .updateTreeSlts(tree)
}

#' @name addTip
#' @title Add tip to a tree
#' @description Returns a tree with a tip ID added
#' @details User must provide a new tip ID, the ID of a node
#' which will become the new tip's sister, a start time point to
#' specify when the new branch will start in time and, an end time point
#' (0 for extant tips).
#' @param tree \code{TreeMan} object
#' @param tid tip ID
#' @param sid ID of node that will become new tip's sister
#' @param start start time
#' @param end end time
#' @param pid parent ID (default is 'p_' + tid)
#' @seealso
#' \code{\link{rmTip}}, 
#' \url{https://github.com/DomBennett/treeman/wiki/manip-methods}
#' @export
#' @examples
#' library(treeman)
#' tree <- randTree(10)
#' # add a new tip to the branch preceding t1
#' # calculate the span and find a point in that time frame for start
#' t1_spn <- getSpnAge(tree, 't1')
#' start <- runif(max=t1_spn[1, 'start'], min=t1_spn[1, 'end'], n=1)
#' end <- runif(max=start, min=0, n=1)
#' tree <- addTip(tree, tid='t11', sid='t1', start=start, end=end)
addTip <- function(tree, tid, sid, start, end,
                   pid=paste0("p_", tid)) {
  # terminology
  # snd, sid -- old sister node and id
  # tnd, tid -- new tip node and id
  # pnd, pid -- new parent node and id
  # gpnd, gpid -- grand parent (prid of old sister)
  updater <- function(nd) {
    # node operation
    i <- which(nd[['prid']] == gpid) - 1
    mxi <- length(nd[['prid']])
    nd[['prid']] <- c(nd[['prid']][0:i], pid, 
                      nd[['prid']][(i+1):mxi])
    nd
  }
  # unpack
  ndlst <- tree@ndlst
  # get key data from tree
  nids <- getNdPtid(tree, sid)
  age <- getNdAge(tree, sid)
  # init new nodes
  tnd <- list('id'=tid)
  snd <- ndlst[[sid]]
  gpid <- snd[['prid']][[1]]
  gpnd <- ndlst[[gpid]]
  pnd <- list('id'=pid, 'kids'=sid)
  # update spans
  tnd[['spn']] <- start - end
  pnd[['spn']] <- snd[['spn']] - (start - age)
  snd[['spn']] <- start - age
  # update ptid
  gpnd[['ptid']] <- gpnd[['ptid']][!gpnd[['ptid']] %in% snd[['id']]]
  gpnd[['ptid']] <- c(gpnd[['ptid']], pnd[['id']])
  pnd[['ptid']] <- c(tid, sid)
  # set prid
  tnd[['prid']] <- snd[['prid']]
  pnd[['prid']] <- snd[['prid']]
  # set prdst
  pnd[['prdst']] <- snd[['prdst']] - snd[['spn']]
  tnd[['prdst']] <- pnd[['prdst']] + tnd[['spn']]
  # set kids
  if(is.null(snd[['kids']])) {
    pnd[['kids']] <- sid
  } else {
    pnd[['kids']] <- snd[['kids']]
  }
  # set pd
  pnd[['pd']] <- snd[['pd']]
  tnd[['pd']] <- 0
  # add to ndlst
  ndlst[[tid]] <- tnd
  ndlst[[pid]] <- pnd
  ndlst[[sid]] <- snd
  ndlst[[gpid]] <- gpnd
  # update upstream prids from sid onwards
  ndlst <- .updateNdsSlt(ndlst, c(tid, nids), updater)
  # update downstream using updateTip
  tree@ndlst <- .updateTip(ndlst, tid=tid, rid=tree@root)
  .updateTreeSlts(tree)
}

#' @name pinTips
#' @title Pin tips to a tree
#' @description Returns a tree with new tips added based on given lineages and time points
#' @details User must provide a vector of new tip IDs, a list of the ranked lineages
#' for these IDs (in ascending order) and a vector of end time points for each new ID
#' (0s for extant tips). The function expects the given tree to be taxonomically informed;
#' the \code{txnym} slot for every node should have a taxonomic label. The function takes
#' the lineage and tries to randomly add the new tip at the lowest point in the taxonomic rank
#' before the end time point. Parallelizable.
#' @param tree \code{TreeMan} object
#' @param tids new tip ids
#' @param lngs list of vectors of the lineages of each tid
#' @param ends end time points for each tid
#' @param ... \code{plyr} arguments
#' @seealso
#' \url{https://github.com/DomBennett/treeman/wiki/manip-methods}
#' @export
#' @examples
#' # see https://github.com/DomBennett/treeman/wiki/Pinning-tips for a detailed example
pinTips <- function(tree, tids, lngs, ends, ...) {
  .pin <- function(i) {
    # unpack
    tid <- tids[i]
    lng <- lngs[[i]]
    end <- ends[i]
    for(j in length(lng):1) {
      spns <- names(txnyms)[which(txnyms %in% lng[j])]
      if(length(spns) == 0) {
        next
      }
      spns <- c(spns, unlist(sapply(spns, function(n) tree@ndlst[[n]][['ptid']])))
      spns <- spns[spns != tree@root]
      rngs <- getSpnsAge(tree, ids=spns)
      bool <- rngs[ ,'start'] > end
      if(any(bool)) {
        rngs <- rngs[bool, ]
        rngs[rngs[ ,'end'] <= end, "end"] <- end
        # pinning is based on branch length
        prbs <- rngs$start - rngs$end
        e <- as.vector(sample(rngs$spn, prob=prbs, size=1))
        e_i <- which(rngs$spn == e)
        start <- runif(min=rngs$end[e_i], max=rngs$start[e_i], n=1)
        if(j != length(lng)) {
          tip_txnym <- lng[j+1]
        } else {
          tip_txnym <- lng[j]
        }
        pid <- paste0('p_', tid, sep='')
        tree <- addTip(tree, tid=tid, sid=e, start=start, end=end,
                        pid=pid)
        tree@ndlst[[tid]][['txnym']] <- tip_txnym
        tree@ndlst[[pid]][['txnym']] <- lng[j]
        # add to txnyms list
        txnyms[[tid]] <<- tip_txnym
        # push out
        tree <<- tree
        break
      }
    }
  }
  .getTxnyms <- function(txnym, ...) {
    txnym
  }
  txnyms <- plyr::mlply(tree@ndlst, .fun=.getTxnyms)
  txnyms <- txnyms[1:length(txnyms)]
  names(txnyms) <- names(tree@ndlst)
  plyr::m_ply(1:length(tids), .pin)
  tree
}
