extractCPT <- function(x, graph, smooth=0){
  if (!inherits(x, c("data.frame","table", "xtabs")))
    stop("'x' must be one of dataframe, table, xtabs")
  if (!inherits(graph, "graphNEL"))
    stop("'graph' must be a graphNEL object")
  if (!is.DAG(graph))
    stop("'graph' must be a DAG")

  V   <- graph::nodes(graph)
  vpa <- vpar(graph)[V]

  if (class(x)[1]=="data.frame"){
      ### cat("extractCPT - data.frame\n")
      ans <- lapply(vpa, function(ss){
          ##cat(sprintf("---- %s ----\n", toString( ss )))
          zzz <- xtabs(~., data=x[, ss, drop=FALSE])
          ##print(zzz)
          zzz
      })
  } else {
    ans <- lapply(vpa, function(ss){tableMargin(x, ss)})
  }

  #.ans <<- ans
  ans <- lapply(ans, as.parray, normalize="first", smooth=smooth)
  chk <- unlist(lapply(ans, function(zz) any(is.na(zz))))
  nnn <- names(chk)[which(chk)]
  if (length(nnn)>0){
      cat(sprintf("NAs found in conditional probability table(s) for nodes: %s\n", toString(nnn)))
      cat(sprintf("  ... consider using the smooth argument\n"))
  }

  class(ans) <- c("extractCPT","list")
  ans
}

extractPOT <- function(x, graph, smooth=0){
  if (!inherits(x, c("data.frame","table")))
    stop("'x' must be a dataframe or a table")
  if (!inherits(graph, "graphNEL"))
    stop("'graph' must be a graphNEL object")
  if (!is.TUG(graph))
    stop("'graph' must be a triangulated undirected graph")

  .rip  <- rip( graph )

  if (class(x)[1]=="data.frame"){
    ans <- .extractPOT_dataframe(x, .rip$cliques, .rip$sep, smooth=smooth)
  } else {
    ans <- .extractPOT_table(x, .rip$cliques, .rip$sep, smooth=smooth)
  }

  attr(ans, "rip")     <- .rip

  ## FIXME: Fix of ug2dag  ##dg 	  <- .ug2dag(graph)
  dg 	    <- ug2dag(graph)
  cptlist <- extractCPT(x, dg, smooth=smooth)
  attr(ans, "dag")     <- dg        ## Needed to save network in Hugin format
  attr(ans, "cptlist") <- cptlist   ## Needed to save network in Hugin format

  class(ans) <- c("extractPOT","list")
  ans
}

.extractPOT_table <- function(x, cliq, seps=NULL, smooth=0){
  ans <- vector("list", length(cliq))
  for ( i  in seq_along(cliq)){
    cq    <- cliq[[ i ]]
    sp    <- seps[[ i ]]
    t.cq  <- tableMargin(x, cq) + smooth
    names(dimnames(t.cq)) <- cq
    if (!is.null(seps) && length(sp)>0){
      t.sp      <- tableMargin(t.cq, sp)
      ans[[ i ]] <- tableOp2(t.cq, t.sp, op=`/`)
    } else {
      ans[[ i ]] <- t.cq / sum(t.cq)
    }
  }
  ans
}

.extractPOT_dataframe <- function(x, cliq, seps=NULL, smooth=0){
  ans <- vector("list", length(cliq))
  for ( i  in seq_along(cliq)){
    cq   <- cliq[[ i ]]
    sp   <- seps[[ i ]]
    ## FIXME: Would say that the following two lines do the same...
    xxx  <- xtabs(~., data=x[ , cq, drop=FALSE])
    t.cq <- tableMargin(xxx, cq) + smooth
    names(dimnames(t.cq)) <- cq
    if (!is.null(seps) && length(sp)>0){
      t.sp       <- tableMargin(t.cq, sp)
      ans[[ i ]] <- tableOp2(t.cq, t.sp, op=`/`)
    } else {
      ans[[ i ]] <- t.cq / sum(t.cq)
    }
  }
  ans
}


## FIXME: gRbase (1.7-2) now has ug2dag
## FIXME: remove .ug2dag after in next release.
#' .ug2dag <- function(ug){
#'   m <- mcs(ug)
#'   if (length(m)==0)
#'     return(NULL)
#'   adjList <- graph::adj(ug, m)
#'   vparList <- vector("list",length(m))
#'   names(vparList) <- m

#'   ii <- 2
#'   vparList[[1]] <- m[1]
#'   for (ii in 2:length(m)){
#'     vparList[[ii]] <- c(m[ii],intersectPrim(adjList[[ii]], m[1:ii]))
#'   }

#'   dg <- dagList(vparList)
#'   dg
#' }
