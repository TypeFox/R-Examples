addTip <- function(phy, tip, tax, insert = "crown", stem.edge = 0.5){
  
  if ( !inherits(phy, "phylo") ) stop("'phy' is not of class 'phylo'")
  #   if ( !is.ultrametric(phy) ) stop("'phy' must be ultrametric")
  
  tip <- gsub(" ", "_", tip)
  
  if ( tip %in% phy$tip.label ) stop("'tip' is already contained in 'phy'")
  insert <- match.arg(insert, c("crown", "stem", "randomly"))
  
  # number of tips, internal nodes and edges:
  nt <- Ntip(phy); ni <- Nnode(phy); ne <- Nedge(phy)
  
  it <- ifelse(insert == "randomly", insert, paste("at", insert))
  cat(" add", it, "of ")
  an <- whereToInsert(phy, tax, tip)
  
  if ( an <= nt ){
    
#     cat(" - add to one congeneric species")
    ## add tip to one congeneric
    ## -------------------------
    pretip <- an
    
    ## an: This is now the subtending (or stem) node
    ## of the node where the new tip is to be inserted
    an <- phy$edge[phy$edge[, 2] == pretip, 1]
    
    ## new internal node number
    newinternal <- descendants(phy, an, type = "i")
    if ( length(newinternal) == 0 ){
      newinternal <- an + 2
    } else {
      newinternal <- max(descendants(phy, an, type = "i")) + 2
    }
    
    
    # increase node number greater than 'pretip' by 1
    # to create a gap to insert new tip at number pretip + 1
    phy$edge[phy$edge > pretip] <- phy$edge[phy$edge > pretip] + 1
    an <- an + 1 # step ancestral node accordingly
    ## create splits above and below gap
    id <- which(phy$edge[, 1] == an & phy$edge[, 2] == pretip)
    upper <- 1:id; lower <- (id + 1):ne
    phy$edge[phy$edge >= newinternal] <- phy$edge[phy$edge >= newinternal] + 1
    phy$edge[id, 2] <- newinternal
    
    # add edges
    phy$edge <- rbind(phy$edge[upper, ], 
                      c(newinternal, pretip),
                      c(newinternal, pretip + 1),
                      phy$edge[lower, ])
    
    # add edge.lengths
    phy$edge.length <- c(phy$edge.length[head(upper, -1)],
                         phy$edge.length[id] * stem.edge,
                         rep(phy$edge.length[id] * (1 - stem.edge), 2),
                         phy$edge.length[lower])
    
    # add tip.label
    phy$tip.label <- c(phy$tip.label[1:pretip], 
                       tip,
                       phy$tip.label[(pretip + 1):nt])
    
    # ajust internal node number
    phy$Nnode <- ni + 1 # and not nt, which would be valid 
    # only for binary trees!
    
  } else {
    
    ## add tip to more than one congenerics or a higher rank
    ## -----------------------------------------------------
#     cat(" - add to more congeneric or higher rank relatives")
    
    if ( insert == "stem" ) an <- noi(phy, strip.spec(tip), 
                                      regex = TRUE, stem = TRUE)
    if ( insert == "randomly" ) {
      an <- descendants(phy, an, type = "i")
      an <- sample(an, 1)
    }
    
    # crown group age:
    if ( is.ultrametric(phy) ){
      ael <- branching.times(phy)[an - Ntip(phy)]
    } else {
      tip.heights <- tipHeights(extract.clade(phy, an))
      ael <- runif(1, min(tip.heights), max(tip.heights))
    }
    
    # number of tip to be inserted: the smallest tip number,
    # because the existing tip numberd will be increased by 1
    newtip <- min(descendants(phy, an, "t"))
    
    # increase number of internal nodes by 1
    phy$edge[phy$edge >= newtip] <- phy$edge[phy$edge >= newtip] + 1
    an <- an + 1
    
    # insert before id (after would be more difficult)
    id <- min(which(phy$edge[, 1] == an))
    upper <- 1:(id - 1); lower <- id:ne
    
    # add edge
    phy$edge <- rbind(phy$edge[upper, ], 
                      c(an, newtip),
                      phy$edge[lower, ])
    # add edge.length
    phy$edge.length <- c(phy$edge.length[upper], 
                         ael,
                         phy$edge.length[lower])
    # add tip.label
    phy$tip.label <- c(phy$tip.label[1:(newtip - 1)], 
                       tip,
                       phy$tip.label[newtip:nt])  
  }
  #   fixNodes(phy)
  phy
}