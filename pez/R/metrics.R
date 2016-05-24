#' Phylogenetic and functional trait metrics within pez
#'
#' Using these functions, you can calculate any of the phylogenetic
#' metrics within pez, using \code{\link{comparative.comm}}
#' objects. While you can call each individually, using the
#' \code{\link{pez.shape}}, \code{\link{pez.evenness}},
#' \code{\link{pez.dispersion}}, and \code{\link{pez.dissimilarity}}
#' wrapper functions (and the more flexible
#' \code{\link{generic.metrics}} and null model functions) are probably
#' your best bet. Note that *all of these functions* take a common
#' first parameter: a \code{\link{comparative.comm}} object. There are
#' additional parameters that can be passed, which are described
#' below.
#'
#' \code{.pd} returns two metrics: Faith's PD (which does not take
#' into account abundance) and Faith's PD corrected for species
#' richness or total abundance (depending on
#' \code{abundance.weighted}). I am almost certain that I got the idea
#' for this from somewhere, but I can't find the reference: if you
#' published on this before 2012, please get in touch with me.
#'
#' @note Many (but not all) of these functions are fairly trivial
#' wrappers around functions in other packages. In the citations for
#' each metric, * indicates a function that's essentially written in
#' \code{\link{picante}}. The Pagel family of measures are also fairly
#' trivial wrapper around \code{\link{caper}} code, functional
#' dissimilarity \code{\link{FD}} code, \code{gamma} \code{\link{ape}}
#' code, and \code{colless}
#' \code{\link[apTreeshape:colless]{apTreeshape}} code. I can't demand
#' it, but I would be grateful if you would cite these authors when
#' using these wrappers.
#'
#' The \code{\link{pez.shape}}, \code{\link{pez.evenness}},
#' \code{\link{pez.dispersion}}, and \code{\link{pez.dissimilarity}}
#' wrapper functions go to some trouble to stop you calculating
#' metrics using inappropriate data (see their notes). These functions
#' give you access to the underlying code within \code{pez}; there is
#' nothing I can do to stop you calculating a metric that, in my
#' opinion, doesn't make any sense. You have been warned :D
#'
#' If you're a developer hoping to make your metric(s) work in this
#' framework, please use the argument naming convention for arguments
#' described in this help file, and use the \code{...} operator in
#' your definition. That way functions that don't need particular
#' arguments can co-exist peacefully with those that do. The first
#' argument to one of these functions should \emph{always} be a
#' \code{\link{comparative.comm}} object; there is no method dispatch
#' on any of these functions and I foresee future pain without this
#' rule.
#' @export
#' @param x \code{\link{comparative.comm}} object
#' @param dist distance matrix for use with calculations; could be
#' generated from traits, a square-root-transformed distance matrix
#' (see \code{\link{.sqrt.phy}} for creating a
#' \code{\link{comparative.comm}} object with a square-root
#' transformed phylogeny). Default: NULL (--> calculate distance
#' matrix from phylogeny)
#' @param abundance.weighted whether to include species' abundances in
#' metric calculation, often dictating whether you're calculating a
#' \code{\link{pez.shape}} or \code{\link{pez.evenness}}
#' metric. Default: FALSE
#' @param na.rm remove NAs in calculations (altering this can obscure
#' errors that are meaningful; I would advise leaving alone)
#' @param include.root include root in PD calculations (default is
#' TRUE, as in picante, but within \code{\link{pez.shape}} I specify
#' FALSE
#' @param method whether to calculate using phylogeny ("phy";
#' default) or trait data ("traits")
#' @param which.eigen which phylo-eigenvector to be used for PVR
#' metric
#' @param q the q parameter for \code{.scheiner}; default 0.0001
#' @param permute number of permutations of null randomisations
#' (mostly only applies to \code{\link[pez:pez.dispersion]{dispersion
#' metrics}})
#' @param null.model one of "taxa.labels", "richness", "frequency",
#' "sample.pool", "phylogeny.pool", "independentswap", or
#' "independentswap". These correspond to the null models available in
#' \code{\link{picante}}; only \code{d} does not use these null models
#' @param ... ignored
#' @importFrom apTreeshape colless tipsubtree
#' @references \code{colless} Colless D.H. (1982). Review of
#' phylogenetics: the theory and practice of phylogenetic
#' systematics. Systematic Zoology, 31, 100-104.
#' @export
#' @rdname pez.metrics
#' @name pez.metrics
#' @examples
#' data(laja)
#' data <- comparative.comm(invert.tree, river.sites)
#' .psv(data)
#' @export
.colless <- function(x, ...)
{
    if(!inherits(x, "comparative.comm"))  stop("'data' must be a comparative community ecology object")
    output <- numeric(nrow(x$comm))
    for(i in seq(nrow(x$comm)))
        output[i] <- colless(as.treeshape(drop_tip(x$phy, colnames(x$comm)[x$comm[i,]==0])))
    names(output) <- rownames(x$comm)
    return(output)
}
#' @importFrom picante pd evol.distinct
#' @references \code{eed,hed} (i.e., \emph{Eed, Hed}) Cadotte M.W.,
#' Davies T.J., Regetz J., Kembel S.W., Cleland E. & Oakley
#' T.H. (2010). Phylogenetic diversity metrics for ecological
#' communities: integrating species richness, abundance and
#' evolutionary history. Ecology Letters, 13, 96-105.
#' @rdname pez.metrics
#' @name pez.metrics
#' @export
.hed <- function(x, ...){
    #Argument handling
    if(!inherits(x, "comparative.comm"))  stop("'data' must be a comparative community ecology object")

    #Setup
    ed <- evol.distinct(x$phy, "fair.proportion")$w
    pd <- pd(x$comm, x$phy)$PD

    #Internal assemblage calc.
    ..hed <- function(ed, pd.comm){
        hed <- ed / pd.comm
        return(-sum(hed * log(hed)))
    }

    #Calculate, clean, and return
    output <- numeric(nrow(x$comm))
    names(output) <- rownames(x$comm)
    for(i in seq(nrow(x$comm)))
        output[i] <- ..hed(ed, pd[i])
    return(output)
}

#' @rdname pez.metrics
#' @name pez.metrics
#' @export
.eed <- function(x, na.rm=TRUE, ...) {
    if(!inherits(x, "comparative.comm"))  stop("'x' must be a comparative community ecology object")
    output <- .hed(x) / log(apply(x$comm, 1, function(x) sum(x != 0)))
    names(output) <- rownames(x)
    return(output)
}

#' @references \code{PSV,PSR,PSE} Helmus M.R., Bland T.J., Williams
#' C.K. & Ives A.R. (2007). Phylogenetic measures of
#' biodiversity. American Naturalist, 169, E68-E83.
#' @importFrom picante psv
#' @rdname pez.metrics
#' @name pez.metrics
#' @export
.psv <- function(x, ...){
    if(!inherits(x, "comparative.comm"))
        stop("'x' must be a comparative.comm object")
    psv(x$comm, x$phy)[,1]
}

#' @importFrom picante psd
#' @rdname pez.metrics
#' @name pez.metrics
#' @export
.psr <- function(x, ...){
    if(!inherits(x, "comparative.comm"))
        stop("'x' must be a comparative.comm object")
    psd(x$comm, x$phy)[,4]
}

#' @importFrom picante mpd
#' @importFrom ape cophenetic.phylo
#' @rdname pez.metrics
#' @name pez.metrics
#' @export
.mpd <- function(x, dist=NULL, abundance.weighted=FALSE, ...){
    if(!inherits(x, "comparative.comm"))
        stop("'x' must be a comparative.comm object")
    if(is.null(dist))
        dist <- cophenetic(x$phy)
    mpd(x$comm, as.matrix(dist), abundance.weighted=abundance.weighted)
}

#' @references \code{PD} Faith D.P. (1992). Conservation evaluation
#' and phylogenetic diversity. Biological Conservation, 61, 1-10.
#' @importFrom picante pd
#' @importFrom stats lm
#' @rdname pez.metrics
#' @name pez.metrics
#' @export
.pd <- function(x, include.root=TRUE, abundance.weighted=FALSE, ...){
    if(!inherits(x, "comparative.comm"))
        stop("'x' must be a comparative.comm object")
    if(abundance.weighted==FALSE)
        x$comm[x$comm > 1] <- 1
    pd <- pd(x$comm, x$phy, include.root)[,1]
    pd.ivs <- unname(resid(lm(pd ~ rowSums(x$comm))))
    return(cbind(pd, pd.ivs))
}

#' @importFrom picante mntd
#' @importFrom ape cophenetic.phylo
#' @rdname pez.metrics
#' @name pez.metrics
#' @export
.mntd <- function(x, dist=NULL, abundance.weighted=FALSE, ...){
    if(!inherits(x, "comparative.comm"))
        stop("'x' must be a comparative.comm object")
    if(is.null(dist))
        dist <- cophenetic(x$phy)
    return(mntd(x$comm, as.matrix(dist), abundance.weighted=abundance.weighted))
}

#' @references \code{gamma} Pybus O.G. & Harvey P.H. (2000) Testing
#' macro-evolutionary models using incomplete molecular
#' phylogenies. _Proceedings of the Royal Society of London. Series
#' B. Biological Sciences 267: 2267--2272.
#' @importFrom ape gammaStat
#' @importFrom apTreeshape as.treeshape
#' @rdname pez.metrics
#' @name pez.metrics
#' @export
.gamma <- function(x, ...){
    if(!inherits(x, "comparative.comm"))
        stop("'x' must be a comparative.comm object")
    
    ..gamma <- function(pa.vec,tree,nams){
        if(sum(pa.vec)<3){
            return(NA)
        } else {
            if(length(setdiff(tree$tip.label, nams)) != 0){
                tree <- drop_tip(setdiff(tree$tip.label, ))
            }
            return(gammaStat(drop_tip(tree,nams[pa.vec==0])))
        }
    }
    
    tree.shape <- as.treeshape(x$phy)
    nams <- tree.shape$names
    return(apply(x$comm, 1, ..gamma, x$phy, nams))
}

#' @importFrom ape cophenetic.phylo
#' @importFrom vegan taxondive
#' @references \code{taxon} Clarke K.R. & Warwick R.M. (1998). A
#' taxonomic distinctness index and its statistical
#' properties. J. Appl. Ecol., 35, 523-531.
#' @rdname pez.metrics
#' @name pez.metrics
#' @export
.taxon <- function(x, dist=NULL, abundance.weighted=FALSE, ...){
    if(!inherits(x, "comparative.comm"))
        stop("'x' must be a comparative.comm object")
    if(is.null(dist))
        dist <- cophenetic(x$phy)
    if(abundance.weighted==FALSE)
        x$comm[x$comm > 1] <- 1
    output <- taxondive(x$comm, dist)
    output <- with(output, data.frame(Delta=D, DeltaStar=Dstar, LambdaPlus=Lambda, DeltaPlus=Dplus, S.DeltaPlus=SDplus))
    return(output)
}

#' @references \code{eigen.sum} Diniz-Filho J.A.F., Cianciaruso M.V.,
#' Rangel T.F. & Bini L.M. (2011). Eigenvector estimation of
#' phylogenetic and functional diversity. Functional Ecology, 25,
#' 735-744.
#' @importFrom ape cophenetic.phylo
#' @importFrom stats var
#' @rdname pez.metrics
#' @name pez.metrics
#' @export
.eigen.sum <- function(x, dist=NULL, which.eigen=1, ...){
    if(!inherits(x, "comparative.comm"))
        stop("'x' must be a comparative.comm object")
    ..eigen.sum <- function(x, evc, vecnums) {
        if(sum(x>0)) {
            return(sum(apply(as.matrix(evc[x>0,vecnums]),2,var)))
        } else {
            return(NA)
        }
    }
    if(is.null(dist))
        dist <- cophenetic(x$phy)
    
    eigen <- -0.5 * dist
    l <- matrix(1/nrow(eigen), nrow=nrow(eigen), ncol=ncol(eigen))
    eigen <- (diag(nrow(eigen)) - l) %*% eigen %*% (diag(nrow(eigen)) - l)
    eigen <- eigen(eigen, symmetric=TRUE)$vectors
    return(apply(x$comm, 1, ..eigen.sum, eigen, which.eigen))
}

#' @importFrom FD dbFD
#' @importFrom ape cophenetic.phylo
#' @rdname pez.metrics
#' @name pez.metrics
#' @export
.dist.fd <- function(x, method="phy", abundance.weighted=FALSE, ...){
    x <- comparative.comm()
    if(!inherits(x, "comparative.comm"))
        stop("'x' must be a comparative.comm object")
    if(method == "phy")
        data <- as.dist(cophenetic(x$phy))
    if(method == "traits")
        data <- x$data
    if(is.matrix(method) | is.data.frame(method))
        data <- method
    if(class(method)=="dist")
        data <- method
    output <- capture.output(dbFD(data, x$comm, w.abun=abundance.weighted, messages=TRUE), file=NULL)
    coefs <- with(output, cbind(coefs, cbind(FRic, FEve, FDiv, FDis, RaoQ)))
    
    #Only bother getting CWMs if we have trait data
    if(method=="traits" | is.data.frame(method)){
        t <- output$dist.fd$CWM
        colnames(t) <- paste(colnames(t), "cmw", sep=".")
        coefs$dist.fd <- rbind(coef$dist.fd, t)
    }
}

#' @importFrom ape is.ultrametric as.phylo cophenetic.phylo
#' @rdname pez.metrics
#' @name pez.metrics
#' @export
.sqrt.phy <- function(x){
    if(!inherits(x, "comparative.comm"))
        stop("'x' must be a comparative.comm object")
    if(!is.ultrametric(x$phy))
        stop("Phylogeny is not ultrametric; cannot square-root (known 'bug', see help)")
    dist <- sqrt(cophenetic(x$phy))
    x$phy <- as.phylo(hclust(as.dist(dist)))
    return(x)
}

#' @references \code{entropy} Allen B., Kon M. & Bar-Yam Y. (2009). A
#' New Phylogenetic Diversity Measure Generalizing the Shannon Index
#' and Its Application to Phyllostomid Bats. The American Naturalist,
#' 174, 236-243.
#' @importFrom ape write.tree
#' @importFrom ade4 newick2phylog
#' @rdname pez.metrics
#' @name pez.metrics
#' @export
.phylo.entropy <- function(x, ...)
{
  #Assertions and argument handling
  if(!inherits(x, "comparative.comm"))  stop("'x' must be a comparative community ecology object")
  
  #Setup
  tree.phylog <- newick2phylog(write.tree(x$phy))
  species <- colnames(x$comm)
  hp.sites <- numeric(nrow(x$comm))
  
  ## beginning of the sites iteration
  for (i in seq(nrow(x$comm)))
  {
    ## species which occured at the site
    site.species <- species[which(x$comm[i,] > 0)]
    
    ## species which NOT occured at the site. They will be removed
    ## from the phylogenetic tree.
    other.species <- species[which(x$comm[i,] == 0)]
    
    ## proportions of occurrence of each species
    proportions <- x$comm[i,site.species] / sum(x$comm[i,])
    
    ## god, how I hate these tree conversions...
    ## also, this comparison is very ugly, it has to be a better
    ## way...
    ## TODO: Search for a better way to verify if an object is empty
    if (length(other.species) == 0) {
      partial.tree <- tree.phylog
    } else {
      if(length(site.species) == 1) {other.species<-c(other.species,site.species)}
      partial.tree <- drop.tip(x$phy, other.species)
      if (all(partial.tree$edge.length[1] == partial.tree$edge.length) | length(site.species) == 2)
      {
        hp.sites[i] <- abs(sum(proportions * log(proportions) * partial.tree$edge.length[1]))
        next
      }
      partial.tree <- newick2phylog(write.tree(drop.tip(x$phy, other.species)))
    }
    ## TODO: Some (I think) partial trees are not rooted. The results
    ## seem to be ok, but the paper states that Hp should be
    ## calculated from a rooted tree. Do not know why though. Check
    ## what can be done to keep partial trees rooted.
    if ("Root" %in% names(partial.tree$nodes))
    {
      partial.branches <-
        partial.tree$nodes[-c(length(partial.tree$nodes))]
    } else {
      partial.branches <- partial.tree$nodes
    }
    ## terminal branches sizes
    partial.leaves <- partial.tree$leaves
    
    ## first part of the calculations. Here we calculate the index for
    ## each terminal branch
    sum.leaves <- sum(partial.leaves *
                        proportions[names(partial.leaves)] *
                        log(proportions[names(partial.leaves)]))
    
    ## storing the first part of the calculation
    hp <- c(sum.leaves)
    
    ## initilizing the list that will hold the descending leaves for
    ## each branch
    descending.leaves <- list()
    
    ## determining the descending leaves for each branch
    for (j in names(partial.branches)) {
      if (all(partial.tree$parts[[j]] %in% names(partial.leaves))) {
        descending.leaves[[j]] <- partial.tree$parts[[j]]
      } else {
        branches <- partial.tree$parts[[j]][!partial.tree$parts[[j]]
                                            %in% names(partial.leaves)]
        leaves <- partial.tree$parts[[j]][partial.tree$parts[[j]] %in%
                                            names(partial.leaves)]
        for (k in branches) {
          leaves <- c(leaves, descending.leaves[[k]])
        }
        descending.leaves[[j]] <- leaves
      }
    }
    ## calculating the index for each internal branch
    for (j in names(partial.branches)) {
      sum.proportions.desc.leaves <-
        sum(proportions[descending.leaves[[j]]])
      hp <- c(hp, (partial.branches[[j]] * sum.proportions.desc.leaves
                   * log(sum.proportions.desc.leaves)))
    }
    ## putting it all together
    hp.sites[i] <- abs(sum(hp))
  }
  #Make the 0 values NAs
  hp.sites[hp.sites==0]<-NA
  names(hp.sites) <- rownames(x$comm)
  ## the end.
  return(hp.sites)
}

#' @importFrom ape extract.clade
#' @importFrom caper clade.matrix
#' @rdname pez.metrics
#' @name pez.metrics
#' @export
.aed <- function(x, ...){
    #Setup
    if(!inherits(x, "comparative.comm"))  stop("'x' must be a comparative community ecology object")
    uni.nodes <- x$phy$edge[,2][!x$phy$edge[,2] %in% seq_along(x$phy$tip.label)]

    #Internal AED for each assemblage
    ..aed <- function(assemblage, tree, uni.nodes, clade.matrix){
        #Nodal values
        node.values <- numeric(length(uni.nodes))
        for(i in seq_along(uni.nodes)){
            t <- extract.clade(tree, uni.nodes[i])
            t.abund <- assemblage[names(assemblage) %in% t$tip.label]
            node.values[i] <- (tree$edge.length[which(tree$edge[,2]==uni.nodes[i])]) / sum(t.abund)
        }
        
        #AED
        aed <- numeric(length(assemblage))
        for(i in seq_along(tree$tip.label)){
            sp <- tree$tip.label[i]
            nodes <- rownames(clade.matrix)[clade.matrix[,i] == 1]
            splength <- tree$edge.length[tree$edge[,2] == i]
            t <- assemblage[i]
            aed[i] <- sum(node.values[which(uni.nodes %in% nodes)]) + unname(ifelse(t==0,0,splength/t))
        }
        return(aed)
    }

    #Calculate, neaten, and return
    aed <- apply(x$comm, 1, ..aed, x$phy, uni.nodes, clade.matrix(x$phy)$clade.matrix[-seq_along(x$phy$tip.label),])
    rownames(aed) <- x$phy$tip.label
    colnames(aed) <- rownames(x$comm)
    aed[aed == Inf | aed == -Inf] <- NA
    return(aed)
}

#' @references \code{pae,aed,iac,haed,eaed} Cadotte M.W., Davies T.J.,
#' Regetz J., Kembel S.W., Cleland E. & Oakley
#' T.H. (2010). Phylogenetic diversity metrics for ecological
#' communities: integrating species richness, abundance and
#' evolutionary history. Ecology Letters, 13, 96-105.
#' @importFrom picante pd
#' @importFrom picante evol.distinct
#' @rdname pez.metrics
#' @name pez.metrics
#' @export
.haed <- function(x, ...){
    #Argument handling
    if(!inherits(x, "comparative.comm"))  stop("'x' must be a comparative community ecology object")

    #Setup
    ed <- evol.distinct(x$phy, "fair.proportion")$w
    pd <- pd(x$comm, x$phy)$PD
    aed <- .aed(x)

    #Internal assemblage calc.
    ..haed <- function(ed, pd.comm, aed.comm, assemblage){
        s.aed <- rep(aed.comm, assemblage) / pd.comm      
        haed <- -sum(s.aed * log(s.aed), na.rm=TRUE)
        return(haed)
    }

    #Calculate, clean, and return
    output <- numeric(nrow(x$comm))
    names(output) <- rownames(x$comm)
    for(i in seq(nrow(x$comm)))
        output[i] <- ..haed(ed, pd[i], aed[,i], x$comm[i,])
    return(output)
}

#' @importFrom ape cophenetic.phylo
#' @rdname pez.metrics
#' @name pez.metrics
#' @export
.simpson.phylogenetic <- function(x) {
    N.relative <- prop.table(x$comm, 2)
    dmat <- cophenetic.phylo(x$phy)
    out <- apply(N.relative, 1, function(n) sum((n %o% n)*dmat))
    return(out) 
}

#' @rdname pez.metrics
#' @name pez.metrics
#' @export
.iac <- function(x, na.rm=TRUE, ...) {
    #Assertions and argument handling
    if(!inherits(x, "comparative.comm"))  stop("'x' must be a comparative community ecology object")
    
    subtrees <- assemblage.phylogenies(x)
    .ancestors <- function(tree, no.root=TRUE){
        mat <- clade.matrix(tree)$clade.matrix
        if(no.root)
            mat <- mat[!apply(mat, 1, function(x) all(x==1)), ]
        #Stop self-references
        diag(mat) <- 0
        apply(mat, 2, function(x) which(x == 1))
    }
        
    .denom <-  function(tree) {
        # Count number of lineages originating at each internal node
        # (i.e. number of splits)
        nSplits <- table(tree$edge[,1])
        # For each tip, take the product of the number of splits across
        # all of its ancestral nodes
        res <- sapply(.ancestors(tree), function(x)
            prod(nSplits[as.character(x)]))
        return(res * 2)
    }

    # now for each subtree...
    denom <- lapply(subtrees, .denom)
    denom <- do.call("cbind", denom)
    nnodes <- sapply(subtrees, function(x) x$Nnode)

    # Calculate expected number of individuals under null hypothesis
    # of equal allocation to each lineage at each (node) split 
    expected <- rowSums(x$comm, na.rm=na.rm) / t(denom)

    # IAC: summed absolute difference between expected and observed
    # abundances, divided by number of nodes
    return(rowSums(abs(expected - x$comm), na.rm=na.rm) / nnodes)
}

#' @rdname pez.metrics
#' @name pez.metrics
#' @importFrom stats setNames
#' @export
.pae <- function(x, na.rm=TRUE, ...) {
    #Assertions and argument handling
    if(!inherits(x, "comparative.comm"))  stop("'x' must be a comparative community ecology object")
    
    subtrees <- assemblage.phylogenies(x)
    PD <- pd(x$comm, x$phy)
    tmp <- setNames(rep(0, ncol(x$comm)), colnames(x$comm))
    TL <- lapply(subtrees, function(tree) {
        #Get terminal edge length
        res <- x$phy$edge.length[x$phy$edge[,2] <= length(x$phy$tip.label)]
        tmp[match(x$phy$tip.label, names(tmp))] <- res
        tmp
    })
    TL <- do.call("cbind", TL)
    numer <- PD$PD + colSums(TL * (t(x$comm) - 1))
    denom <- PD$PD + (rowSums(x$comm, na.rm = na.rm) / rowSums(x$comm,
        na.rm=na.rm) - 1) * colSums(TL)
    res <- numer/denom
    names(res) <- rownames(x$comm)
    return(res)
}

#' @importFrom picante evol.distinct
#' @rdname pez.metrics
#' @name pez.metrics
#' @references \code{scheiner} Scheiner, S.M. (20120). A metric of
#' biodiversity that integrates abundance, phylogeny, and function.
#' Oikos, 121, 1191-1202.
#' @export
.scheiner <- function(x, q=0.0001, abundance.weighted = FALSE, ...){
    #Assertions and argument handling
    if(!inherits(x, "comparative.comm")) stop("'x' must be a comparative community ecology object")
    
    #Setup
    ed <- evol.distinct(x$phy, "fair.proportion")$w
    if(!abundance.weighted)
        x$comm[x$comm > 0] <- 1
    
    #Calculate scheiner; beware dividing by zero inadvertantly
    output <- numeric(nrow(x$comm))
    for(i in seq(nrow(x$comm))){
        if(q==1)
            output[i] <- exp(-1*sum(((x$comm[i,]*ed[i])/(sum(x$comm[i,])*ed[i])) * log((x$comm[i,]*ed[i])/(sum(x$comm[i,])*ed[i]))))
        else
            output[i] <- sum(((x$comm[i,]*ed[i])/(sum(x$comm[i,])*ed[i]))^q)^(1/(1-q))
    }
    return(output)
}

#' @importFrom picante pse
#' @rdname pez.metrics
#' @name pez.metrics
#' @export
.pse <- function(x, ...){
    if(!inherits(x, "comparative.comm"))
        stop("'x' must be a comparative.comm object")
    pse(x$comm, x$phy)[,1]
}

#' @references \code{rao} Webb C.O. (2000). Exploring the phylogenetic
#' structure of ecological communities: An example for rain forest
#' trees. American Naturalist, 156, 145-155.
#' @importFrom picante raoD
#' @rdname pez.metrics
#' @name pez.metrics
#' @export
.rao <- function(x, ...){
    if(!inherits(x, "comparative.comm"))
        stop("'x' must be a comparative.comm object")
    raoD(x$comm, x$phy)$Dkk
}

#' @references \code{lambda,delta,kappa} Mark Pagel (1999) Inferring
#' the historical patterns of biological evolution. Nature 6756(401):
#' 877--884.
#' @importFrom caper pgls comparative.data
#' @rdname pez.metrics
#' @name pez.metrics
#' @export
.lambda <- function(x, ...){
    if(!inherits(x, "comparative.comm"))
        stop("'x' must be a comparative.comm object")
    
    output <- numeric(length(sites(x)))
    for(i in seq(nrow(x$comm))){
        c.data <- x$comm[i,][x$comm[i,] > 0]
        if(length(unique(c.data)) > 1){
            c.data <- data.frame(comm=c.data, names=names(c.data))
            c.data <- comparative.data(phy=drop_tip(x$phy, setdiff(x$phy$tip.label, c.data$names)), data=c.data, names.col=names)
            model <- pgls(comm ~ 1, c.data, lambda="ML")
            output[i] <- summary(model)$param.CI$lambda$opt
        } else output[i] <- NA
    }
    return(output)
}

#' @importFrom caper pgls comparative.data
#' @rdname pez.metrics
#' @name pez.metrics
#' @export
.delta <- function(x, ...){
    if(!inherits(x, "comparative.comm"))
        stop("'x' must be a comparative.comm object")
    
    output <- numeric(length(sites(x)))
    for(i in seq(nrow(x$comm))){
        c.data <- x$comm[i,][x$comm[i,] > 0]
        if(length(unique(c.data)) > 1){
            c.data <- data.frame(comm=c.data, names=names(c.data))
            c.data <- comparative.data(phy=drop_tip(x$phy, setdiff(x$phy$tip.label, c.data$names)), data=c.data, names.col=names)
            model <- pgls(comm ~ 1, c.data, delta="ML")
            output[i] <- summary(model)$param.CI$delta$opt
        } else output[i] <- NA
    }
    return(output)
}

#' @importFrom caper pgls comparative.data
#' @rdname pez.metrics
#' @name pez.metrics
#' @export
.kappa <- function(x, ...){
    if(!inherits(x, "comparative.comm"))
        stop("'x' must be a comparative.comm object")
    
    output <- numeric(length(sites(x)))
    for(i in seq(nrow(x$comm))){
        c.data <- x$comm[i,][x$comm[i,] > 0]
        if(length(unique(c.data)) > 1){
            c.data <- data.frame(comm=c.data, names=names(c.data))
            c.data <- comparative.data(phy=drop_tip(x$phy, setdiff(x$phy$tip.label, c.data$names)), data=c.data, names.col=names)
            model <- pgls(comm ~ 1, c.data, kappa="ML")
            output[i] <- summary(model)$param.CI$kappa$opt
        } else output[i] <- NA
    }
    return(output)
}

#' @rdname pez.metrics
#' @name pez.metrics
#' @export
.eaed <- function(x, ...){
    if(!inherits(x, "comparative.comm"))
        stop("'x' must be a comparative.comm object")
    return(.haed(x)/log(rowSums(x$comm)))
}

#' @references \code{unifrac} Lozupone C.A. & Knight
#' R. (2005). UniFrac: a new phylogenetic method for comparing
#' microbial communities. Applied and Environmental Microbiology, 71,
#' 8228-8235.
#' @importFrom picante unifrac
#' @rdname pez.metrics
#' @name pez.metrics
#' @export
.unifrac <- function(x, ...){
    if(!inherits(x, "comparative.comm"))
        stop("'x' must be a comparative.comm object")
    return(unifrac(x$comm, x$phy))
}

#' @references \code{pcd} Ives A.R. & Helmus M.R. (2010). Phylogenetic
#' metrics of community similarity. The American Naturalist, 176,
#' E128-E142.
#' @importFrom picante ses.mpd
#' @rdname pez.metrics
#' @name pez.metrics
#' @export
.pcd <- function(x, permute=1000, ...){
    if(!inherits(x, "comparative.comm"))
        stop("'x' must be a comparative.comm object")
    return(pcd(x$comm, x$phy, reps=permute))
}

#' @references \code{comdist} C.O. Webb, D.D. Ackerly, and
#' S.W. Kembel. 2008. Phylocom: software for the analysis of
#' phylogenetic community structure and trait
#' evolution. Bioinformatics 18:2098-2100.
#' @importFrom picante comdist
#' @rdname pez.metrics
#' @name pez.metrics
#' @export
.comdist <- function(x, dist=NULL, abundance.weighted=FALSE, ...){
    if(!inherits(x, "comparative.comm"))
        stop("'x' must be a comparative.comm object")
    if(is.null(dist))
        dist <- cophenetic(x$phy)
    return(comdist(x$comm, as.matrix(dist), abundance.weighted=abundance.weighted))
}

#' @references \code{phylosor} Bryant J.A., Lamanna C., Morlon H.,
#' Kerkhoff A.J., Enquist B.J. & Green J.L. (2008). Microbes on
#' mountainsides: Contrasting elevational patterns of bacterial and
#' plant diversity. Proceedings of the National Academy of Sciences of
#' the United States of America, 105, 11505-11511.
#' @importFrom picante phylosor
#' @rdname pez.metrics
#' @name pez.metrics
#' @export
.phylosor <- function(x, dist=NULL, abundance.weighted=FALSE, ...){
    if(!inherits(x, "comparative.comm"))
        stop("'x' must be a comparative.comm object")
    if(is.null(dist))
        dist <- cophenetic(x$phy)
    output <- phylosor(x$comm, x$phy)
    output <- as.dist(1 - as.matrix(output))
    return(output)
}


#' @references \code{d} Fritz S.A. & Purvis A. (2010). Selectivity in
#' Mammalian Extinction Risk and Threat Types: a New Measure of
#' Phylogenetic Signal Strength in Binary Traits. Conservation
#' Biology, 24, 1042-1051.
#' @importFrom caper contrCalc VCV.array
#' @importFrom mvtnorm rmvnorm
#' @importFrom stats reorder
#' @rdname pez.metrics
#' @name pez.metrics
#' @export
.d <- function(x, permute=1000, ...) {
  #Checking
  if(! inherits(x, "comparative.comm"))  stop("'x' must be a comparative community ecology object")
  if (!is.numeric(permute)) (stop("'", permute, "' is not numeric."))
  x$comm[x$comm > 1] <- 1
  # check tree branch lengths
  el    <- x$phy$edge.length
  elTip <- x$phy$edge[,2] <= length(x$phy$tip.label)
  
  if(any(el[elTip] == 0)) 
    stop('Phylogeny contains pairs of tips on zero branch lengths, cannot currently simulate')
  if(any(el[! elTip] == 0)) 
    stop('Phylogeny contains zero length internal branches. Use di2multi.')


  #Internal D calculation
  ..d <- function(ds, vcv, permute, phy){
      dsSort <- sort(ds)
      
      ## Random Association model
      ds.ran <- replicate(permute, sample(ds))
      
      ## Brownian Threshold model random data
      ds.phy <- rmvnorm(permute, sigma=unclass(vcv)) # class of 'VCV.array' throws the method dispatch
      ds.phy <- as.data.frame(t(ds.phy))
      
                                        # turn those into rank values, then pull that rank's observed value
      ds.phy <- apply(ds.phy, 2, rank, ties="random")
      ds.phy <- apply(ds.phy, 2, function(x) as.numeric(dsSort[x]))
      
      ## Get change along edges
      ## insert observed and set dimnames for contrCalc
      ds.ran <- cbind(Obs=ds, ds.ran)
      ds.phy <- cbind(Obs=ds, ds.phy)
      dimnames(ds.ran) <- dimnames(ds.phy) <- list(x$phy$tip.label, c('Obs', paste('V',1:permute, sep='')))
      
      ## now run that through the contrast engine 
      ds.ran.cc <- contrCalc(vals=ds.ran, phy=phy, ref.var='V1', picMethod='phylo.d', crunch.brlen=0)
      ds.phy.cc <- contrCalc(vals=ds.phy, phy=phy, ref.var='V1', picMethod='phylo.d', crunch.brlen=0)
      
      ## get sums of change and distributions
      ransocc <- colSums(ds.ran.cc$contrMat)
      physocc <- colSums(ds.phy.cc$contrMat)
      # double check the observed, but only to six decimal places or you can get floating point errors
      if(round(ransocc[1], digits=6) != round(physocc[1], digits=6)) stop('Problem with character change calculation in phylo.d')
      obssocc <- ransocc[1]
      ransocc <- ransocc[-1]
      physocc <- physocc[-1]
      
      soccratio <- (obssocc - mean(physocc)) / (mean(ransocc) - mean(physocc))
      soccpval1 <- sum(ransocc < obssocc) / permute
      soccpval0 <- sum(physocc > obssocc) / permute
      
      return(c(soccratio, soccpval1, soccpval0))
  }

  ## being careful with the edge order - pre-reorder the phylogeny
  phy <- reorder(x$phy, 'pruningwise')
  
  vcv <- VCV.array(x$phy)
  vals <- matrix(ncol=3, nrow=nrow(x$comm))
  rownames(vals) <- rownames(x$comm)
  colnames(vals) <- c("D", "P(D=1)", "P(D=0)")
  for(i in seq(nrow(x$comm)))
      vals[i,] <- ..d(x$comm[i,], vcv, permute, x$phy)
  
  return(vals)
}

#' @references \code{sesmpd,sesmntd} Webb C.O. (2000). Exploring the
#' phylogenetic structure of ecological communities: An example for
#' rain forest trees. American Naturalist, 156, 145-155.
#' @importFrom picante ses.mpd
#' @rdname pez.metrics
#' @name pez.metrics
#' @importFrom picante ses.mpd
#' @export
.ses.mpd <- function(x, dist=NULL, null.model="taxa.labels", abundance.weighted=FALSE, permute=1000, ...){
    if(!inherits(x, "comparative.comm"))
        stop("'x' must be a comparative.comm object")
    if(is.null(dist))
        dist <- cophenetic(x$phy)
    return(ses.mpd(x$comm, dis=as.matrix(dist), null.model=null.model, abundance.weighted=abundance.weighted, runs=permute))
}

#' @importFrom picante ses.mntd
#' @rdname pez.metrics
#' @name pez.metrics
#' @importFrom picante ses.mntd
#' @export
.ses.mntd <- function(x, dist=NULL, null.model="taxa.labels", abundance.weighted=FALSE, permute=1000, ...){
    if(!inherits(x, "comparative.comm"))
        stop("'x' must be a comparative.comm object")
    if(is.null(dist))
        dist <- cophenetic(x$phy)
    return(ses.mntd(x$comm, dis=as.matrix(dist), null.model=null.model, abundance.weighted=abundance.weighted, runs=permute))
}

#' @references \code{innd,mipd} Ness J.H., Rollinson E.J. & Whitney
#' K.D. (2011). Phylogenetic distance can predict susceptibility to
#' attack by natural enemies. Oikos, 120, 1327-1334.
#' @importFrom picante ses.mpd
#' @rdname pez.metrics
#' @name pez.metrics
#' @export
.ses.mipd <- function(x, dist=NULL, null.model="taxa.labels", abundance.weighted=FALSE, permute=1000, ...){
    if(!inherits(x, "comparative.comm"))
        stop("'x' must be a comparative.comm object")
    if(is.null(dist))
        dist <- cophenetic(x$phy)
    return(ses.mpd(x$comm, dis=1/as.matrix(dist), null.model=null.model, abundance.weighted=abundance.weighted, runs=permute))
}

#' @importFrom picante ses.mntd
#' @rdname pez.metrics
#' @name pez.metrics
#' @export
.ses.innd <- function(x, dist=NULL, null.model="taxa.labels", abundance.weighted=FALSE, permute=1000, ...){
    if(!inherits(x, "comparative.comm"))
        stop("'x' must be a comparative.comm object")
    if(is.null(dist))
        dist <- cophenetic(x$phy)
    return(ses.mntd(x$comm, dis=1/as.matrix(dist), null.model=null.model, abundance.weighted=abundance.weighted, runs=permute))
}

#' @importFrom picante mpd
#' @rdname pez.metrics
#' @name pez.metrics
#' @export
.mipd <- function(x, dist=NULL, abundance.weighted=FALSE, ...){
    if(!inherits(x, "comparative.comm"))
        stop("'x' must be a comparative.comm object")
    if(is.null(dist))
        dist <- cophenetic(x$phy)
    return(mpd(x$comm, dis=1/as.matrix(dist), abundance.weighted=abundance.weighted))
}

#' @importFrom picante mntd
#' @rdname pez.metrics
#' @name pez.metrics
#' @export
.innd <- function(x, dist=NULL, abundance.weighted=FALSE, ...){
    if(!inherits(x, "comparative.comm"))
        stop("'x' must be a comparative.comm object")
    if(is.null(dist))
        dist <- cophenetic(x$phy)
    return(mntd(x$comm, dis=1/as.matrix(dist), abundance.weighted=abundance.weighted))
}
