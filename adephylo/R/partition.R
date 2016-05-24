##
## Functions to obtain partitions of tips from a tree.
## For instance to obtain dummy vectors used in the orthogram.
##





############
# listTips
############
listTips <- function(x){
    ## if(!require(phylobase)) stop("phylobase package is not installed")

    ## conversion from phylo, phylo4 and phylo4d
    x <- as(x, "phylo4")

    ## check phylo4 object
    if (is.character(checkval <- checkPhylo4(x))) stop(checkval)

    ## computations
    nodIdx <- nTips(x)+1
    nodIdx <- nodIdx:(nodIdx+nNodes(x)-1)
    res <- lapply(nodIdx, function(i) descendants(x, i))

    if(hasNodeLabels(x)) {names(res) <- nodeLabels(x)}

    return(res)
} # end listTips





###########
# treePart
###########
treePart <- function(x, result=c("dummy", "orthobasis")){
    ## if(!require(phylobase)) stop("phylobase package is not installed")

    ## conversion from phylo, phylo4 and phylo4d
    x <- as(x, "phylo4")
    result <- match.arg(result)

    ## check phylo4 object
    if (is.character(checkval <- checkPhylo4(x))) stop(checkval)

    n <- nTips(x) # number of tips
    HTU.idx <- (n+1):(n+nNodes(x)) # index of internal nodes (HTU)

    if(!hasNodeLabels(x)) { # node labels will be used after
        nodeLabels(x) <- as.character(HTU.idx)
    }

    ## function coding one dummy vector
    fDum <- function(vec){ # vec is a vector of tip numbers
        dum <- integer(n)
        dum[vec] <- 1
        return(dum)
    }

    ## main computations
    temp <- listTips(x)
    res <- data.frame(lapply(temp,fDum))
    row.names(res) <- tipLabels(x)
    res <- res[,-1, drop=FALSE]

    if(result=="dummy"){
        return(res) # res is a data.frame of dummy vectors
    }



    ## If orthobasis is required ##

    ## Find values 'w' for all nodes
    ##
    ## Notations:
    ## - n: an internal node (HTU)
    ## - Dn: the set of all internal nodes descending from 'n'
    ## - En: the set 'n U Dn' (that is, Dn plus n itself)
    ## - ndd(e): the number of direct descendants from a node 'e'
    ##
    ## Then the values 'w' are computed as:
    ##
    ## w(n) = sum_{e \in En} lgamma( ndd(e) + 1)
    ##

    listDDx <- listDD(x)

    nbOfDD <- sapply(listDDx, length) # nb of DD for each node
    names(nbOfDD) <- HTU.idx # used to match the results of Dn

    findAlldHTU <- function(node){ # find all HTU descending from a node
        res <- descendants(x, node, type="all") # tips and HTU
        res <- res[res > n] # only HTU (here, just node numbers are kept
        if(length(res)==0) return(NULL)
        return(res)
    }


    listAlldHTU <- lapply(HTU.idx, function(node) c(node,findAlldHTU(node))) # ='Dn': for each HTU, list all HTU descending from it

    w <- sapply(listAlldHTU, function(e) sum(lgamma(nbOfDD[as.character(e)]+1))) # w(n)
    ## from now on, 'w' stores the w(n) values.

    ## add dummy vectors for tips
    res <- cbind(diag(1, n), root=rep(1,n), res) # sorted from first tip to last node
    colnames(res) <- 1:(nTips(x) + nNodes(x))
    valDum <- c(rep(-1, n), w) # dummy vectors of tips are given a negative value
    ## note: valDum is the w values for all nodes, sorted from first tip to last node

    ## Discard dummy vectors with lowest valDum (value of dummy vectors, w).
    ## -> for each node, a dummy vector associated to its DD is removed
    ## this one is that with the lowest valDum.

    discardOneDum <- function(node, DDnode){ # node is a node label, not a node number
        if(length(DDnode)==1) return(NULL)
        val <- valDum[DDnode]
        toRemove <- which.min(val)
        keptDD <- DDnode[-toRemove]
        return(keptDD)
    } # end discardOneDum

    dumToKeep <- lapply(1:length(listDDx), function(i) discardOneDum(i, listDDx[[i]]))
    dumToKeep <- unlist(dumToKeep) # contains indices of kept dummy vectors

    res <- res[dumToKeep] # retained dummy vectors
    res <- res[,order(valDum[dumToKeep], decreasing=TRUE)] # reorder vectors by decreasing w

    ## orthonormalization
    res <- cbind(root=rep(1,n), res) # for centring: vectors will be orthogonal to 1_n
    res <- qr.Q(qr(res)) # Gram-Schmidt orthogonalization
    res <- res[,-1] # keep only centred vectors; orthogonal for identity
    res <- res * sqrt(n) # render vectors orthogonal for 1/n

    rownames(res) <- tipLabels(x)
    colnames(res) <- paste("V",1:ncol(res))

    return(as.data.frame(res))

} # end treePart



## EXAMPLE
##
## plot(x <- read.tree(te=newick.eg[[2]]))
## plot(y <- newick2phylog(newick.eg[[2]]), clabel.node=1)
##
##
