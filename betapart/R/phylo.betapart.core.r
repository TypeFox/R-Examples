phylo.betapart.core<-function(x, tree)
{

    if (!is.matrix(x)) {
        x <- as.matrix(x)
    }

    if(nrow(x)<2)
        stop("Computing dissimilairty requires at least 2 communities", call. = TRUE)

    if (!is.numeric(x))
        stop("The data in 'x' is not numeric.", call. = TRUE)

    xvals <- unique(as.vector(x))
    if (any(!is.element(xvals, c(0, 1))))
        stop("The community matrix contains values other than 0 and 1: data should be presence/absence.", call. = TRUE)

    if (class(tree)!="phylo")
        stop("### invalid tree's format: \"phylo\" format required ###\n\n", call. = TRUE)

    if(any(!(colnames(x)%in%tree$tip)))
        warning("At least one species in community matrix is not included in the tree" , call. = TRUE)

    ############ Phylogenetic Diversity Faith (adapted from pd function in picante library) #######################

    pdnew <- function (samp, tree, include.root = TRUE) {

      if (is.null(tree$edge.length)) {
        stop("Tree has no branch lengths, cannot compute pd")
        }
      species <- colnames(samp)
      tree  <- node.age(tree)
      PDout <- apply(samp,1, function(x) {
        present <- species[x > 0]
        treeabsent <- tree$tip.label[which(!(tree$tip.label %in%present))]
        if (length(present) == 0) {
          PD <- 0
          }
        else if (length(present) == 1) {
          if (!is.rooted(tree) || !include.root) {
            warning("Rooted tree and include.root=TRUE argument required to calculate PD of single-species sampunities. Single species sampunity assigned PD value of NA.")
            PD <- NA
            }
          else {
            PD <- tree$ages[which(tree$edge[, 2] ==
            which(tree$tip.label == present))]
            }
          }
        else if (length(treeabsent) == 0) {
          PD <- sum(tree$edge.length)
          }
        else {
          sub.tree <- drop.tip(tree, treeabsent)
          if (include.root) {
            if (!is.rooted(tree)) {
              stop("Rooted tree required to calculate PD with include.root=TRUE argument")
              }
            sub.tree <- node.age(sub.tree)
            sub.tree.depth <- max(sub.tree$ages)
            orig.tree.depth <- max(tree$ages[which(tree$edge[,2] %in% which(tree$tip.label %in% present))])
            PD <- sum(sub.tree$edge.length) + (orig.tree.depth - sub.tree.depth)
            }
          else {
            PD <- sum(sub.tree$edge.length)
            }
          }
        SR <- length(present)
        PDout <- c(PD,SR)
        } )
      PDout <- t(PDout)
      rownames(PDout) <- rownames(samp)
      colnames(PDout) <- c("PD","SR")
      return(PDout)
    } # end of function pdnew
    
    ############ Paired matrix to distance matrix conversion (utility function) #######################

    dist.mat <- function(com,pair) {

      ncom <- nrow(com)
      distmat <- matrix(nrow=ncom,ncol=ncom,0,dimnames=list(rownames(com),rownames(com)))
      st <- c(0,cumsum(seq(ncom-1,2)))+1
      end <- cumsum(seq(ncom-1,1))
      for (i in 1:(ncom-1)) distmat[i,(ncom:(seq(1,ncom)[i]))]=c(pair[end[i]:st[i]],0)
      distmat <- as.dist(t(distmat))
      return(distmat)

    } # end of function dist.mat

    # pariwise comparisons
    com=x
    combin  <-  combn(nrow(com),2) # table with all pairs

    labcomb <-  apply(combin,2,function(x) paste(rownames(com)[x],collapse="-"))

    pd <-  pdnew(com,tree)[,"PD"] # PD for each community of the community matrix
    com.tot.pair <- t(apply(combin,2,function(x) colSums(com[x,])>0))
    pd.tot.pair <- pdnew(com.tot.pair,tree)[,"PD"]  # PD of the two communities combined
    sum.pd.pair <- apply(combin,2,function(x) sum(pd[x])) # Sum of PD for each community, separetely
    com.tot.multi <- t(as.matrix(colSums(com)>0))
    pd.tot.multi <- as.numeric(pdnew(com.tot.multi,tree)[,"PD"])  # PD of all communities combined

    min.not.shared <- apply(pd.tot.pair-t(combn(pd,2)),1,min) # minimun (b,c)
    max.not.shared <- apply(pd.tot.pair-t(combn(pd,2)),1,max) # maximum (b,c)
    sum.not.shared <- 2*pd.tot.pair - sum.pd.pair  # b+c
    shared <- pd.tot.pair - sum.not.shared    # a

    # returning results of functional.betapart.core
    phylo.computations<-list( sumSi=sum(pd),St=pd.tot.multi, shared = dist.mat(com,shared), sum.not.shared = dist.mat(com,sum.not.shared),
    								max.not.shared = dist.mat(com,max.not.shared), min.not.shared = dist.mat(com,min.not.shared))
	class(phylo.computations) <- "phylo.betapart"
    return(phylo.computations)

} # end of function