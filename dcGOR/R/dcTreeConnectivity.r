#' Function to calculate the sparse connectivity matrix between parents and children from a phylo-formatted phylogenetic tree
#'
#' \code{dcTreeConnectivity} is supposed to calculate the sparse connectivity matrix between parents and children from a phylo-formatted phylogenetic tree. The matrix has internal nodes (in rows) and tips plus internal nodes (in columns). For a row (an internal node; as a parent), the non-zeros indicate all its descendants/children.
#'
#' @param phy an object of class 'phylo'
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to TRUE for display
#' @return
#' a sparse matrix of \eqn{Nnode} X \eqn{Ntip+Nnode}, where \eqn{Ntip} and \eqn{Nnode} are the number of tips and internal nodes. A non-zero entry indicates a pair of a parent and its child.
#' @note
#' None
#' @export
#' @seealso \code{\link{dcTreeConnectivity}}
#' @include dcTreeConnectivity.r
#' @examples
#' # a newick tree
#' tree <- "(((t1:5,t2:5):2,(t3:4,t4:4):3):2,(t5:4,t6:4):6);"
#' phy <- ape::read.tree(text=tree)
#'
#' # connectivity matrix
#' res <- dcTreeConnectivity(phy)
#' dim(res)
#' # convert to a full Matrix
#' as.matrix(res)

dcTreeConnectivity <- function(phy, verbose=T)
{
    
    if (class(phy) != "phylo"){
        stop("The input 'phy' must belong to the class 'phylo'!")
    }
    
    ## Define the tree struct dimensions and indexing
    Ntip <- ape::Ntip(phy)
    Nnode <- ape::Nnode(phy)
    Ntot <- Ntip+Nnode
    
    ## Calculate most recent common ancestor (MRCA) for each pair of tips and nodes
    mrca_node <- ape::mrca(phy, full=T)
    ### exclude self-self
    diag(mrca_node) <- 0
    
    ## Calculate connectivity linking each ancestor to its all children
    ### A matrix of Nnode X (Ntip+Nnode)
    connectivity <- matrix(0, nrow=Nnode, ncol=Ntip+Nnode)
    rownames(connectivity) <- (Ntip+1):Ntot
    colnames(connectivity) <- 1:Ntot
    for (i in 1:Nnode) {
        node_tmp <- i+Ntip
        child <- which(mrca_node[node_tmp,]==node_tmp, arr.ind=T)
        connectivity[i,child] <- 1
    }
    
    if(verbose){
        message(sprintf("The connectivity for %d internal nodes X %d tips+internal nodes", Nnode, Ntot), appendLF=T)
    }
    
    res <- Matrix::Matrix(connectivity, sparse=T)
    
    invisible(res)
}
