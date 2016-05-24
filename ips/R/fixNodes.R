## package: ips
## called by: USER
## author: Christoph Heibl (at gmx.net)
## last update: 2014-07-30

fixNodes <- function(phy){
  
  if (!inherits(phy, "phylo"))                             
    stop("object 'phy' is not of class 'phylo'")
  unfixed <- phy
  s <- 999999
  tip <- 1:Ntip(phy)
  internal <- Ntip(phy) + (1:Nnode(phy))
  pool <- union(tip, internal)
  root <- min(internal)
  is.internal <- function(x, maxtip){x > maxtip}
  phy$edge <- phy$edge + s
  
  ## fix root node
  ## -------------
  phy$edge[phy$edge[, 1] == min(internal) + s, 1] <- root
  pool <- setdiff(pool, root)
  internal <- setdiff(internal, root)
  
  ## fix non-root nodes
  ## ------------------
  for ( i in 1:nrow(phy$edge) ){
    n <- phy$edge[i, 2]
    if ( n <= Ntip(phy) + s ){
      phy$edge[i, 2] <- tip[1]
      tip <- tip[-1]
    } else {
      phy$edge[phy$edge[, 1] == n, 1] <- internal[1]
      phy$edge[i, 2] <- internal[1]
      internal <- internal[-1]
    }
  }
  
  ## fix tip- and nodelabels
  ## -----------------------
  trans <- cbind(unfixed$edge[, 2], phy$edge[, 2])
  ttrans <- trans[trans[, 2] <= Ntip(phy), ]
  ntrans <- trans[trans[, 2] > Ntip(phy), ]
  phy$tip.label <- phy$tip.label[match(ttrans[, 1], ttrans[, 2])]
  phy$node.label[-1] <- phy$node.label[-1][match(ntrans[, 1], ntrans[, 2])]
  
#   ADD: FIX NON-CANONICAL NODELABELS
  
  # fix <node.label> element
  # ------------------------
#   id <- node.trans(source = unfixed, target = phy, index = TRUE)
#   sdtnames <- c("edge", "edge.length", "Nnode", "tip.label")
#   if ( !all(names(phy) %in% sdtnames) ){
#     nls <- which(!names(phy) %in% sdtnames)
#     for ( i in nls )
#       phy[[i]] <- phy[[i]][id]
#   }  
  
  phy
}