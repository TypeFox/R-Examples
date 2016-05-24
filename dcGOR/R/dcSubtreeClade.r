#' Function to extract a subtree under a given clade from a phylo-formatted phylogenetic tree
#'
#' \code{dcSubtreeClade} is supposed to extract a subtree under a given clade from a phylo-formatted phylogenetic tree. In addition to the tree in subject, another input is a built-in integer specifying an internal node/clade of interest. Alternatively, the internal node of interest can be given by its label (if there are internal node labels). As a result, a subtree under a given clade is also represented as an object of class 'phylo'.
#'
#' @param phy an object of class 'phylo'
#' @param choose.node an integer specifying which internal node is chosen. For an object of class 'phylo', the tree has built-in ID for internal nodes, ranging from \eqn{Ntip+1} to \eqn{Ntip+Nnode}, where \eqn{Ntip} and \eqn{Nnode} are the number of tips and internal nodes. Internal nodes are indexed in a pre-ordered manner. The subtree under the given interna node will be extracted
#' @param choose.node.label a character specifying which internal node is chosen. For the tree with internal node labels, the extraction of subtree can be done in this way
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to TRUE for display
#' @return
#' an object of class 'phylo'
#' @note
#' If a valid 'choose.node' is given, then 'choose.node.label' will be ignored.
#' @export
#' @seealso \code{\link{dcTreeConnectivity}}
#' @include dcSubtreeClade.r
#' @examples
#' # 1) a newick tree without internal node labels
#' tree <- "(((t1:5,t2:5):2,(t3:4,t4:4):3):2,(t5:4,t6:4):6);"
#' phy <- ape::read.tree(text=tree)
#' phy
#' Ntip <- ape::Ntip(phy)
#' Nnode <- ape::Nnode(phy)
#' ape::plot.phylo(phy, type="p", use.edge.length=TRUE)
#' ape::nodelabels(node=Ntip+1:Nnode, col="red", bg="white")
#' # a subtree specified via a built-in internal node ID
#' subphy <- dcSubtreeClade(phy, choose.node=Ntip+2)
#' subphy
#' ape::plot.phylo(subphy, type="p", use.edge.length=TRUE)
#'
#' # 2) a newick tree with internal node labels
#' tree <- "(((t1:5,t2:5)i3:2,(t3:4,t4:4)i4:3)i2:2,(t5:4,t6:4)i5:6)i1;"
#' phy <- ape::read.tree(text=tree)
#' phy
#' ape::plot.phylo(phy, type="p", use.edge.length=TRUE, show.node.label=TRUE)
#' # a subtree specified via an internal node label
#' subphy <- dcSubtreeClade(phy, choose.node.label='i2')
#' subphy
#' ape::plot.phylo(subphy, type="p", use.edge.length=TRUE, show.node.label=TRUE)

dcSubtreeClade <- function(phy, choose.node=NULL, choose.node.label=NULL, verbose=T)
{
    
    if (class(phy) != "phylo"){
        stop("The input 'phy' must belong to the class 'phylo'!")
    }
    
    Ntip <- ape::Ntip(phy)
    Nnode <- ape::Nnode(phy)
    Ntot <- Ntip+Nnode
    
    flag_k <- TRUE
    if(!is.null(choose.node)){
        k <- as.integer(choose.node)
        if(k<Ntip+1 | k>Ntot){
            flag_k <- FALSE
        }
    }else{
        flag_k <- FALSE
    }
    
    if(flag_k==FALSE){
        if(!is.null(phy$node.label)){
            if(!is.null(choose.node.label)){
                ind <- match(choose.node.label, phy$node.label)
                if(!is.na(ind)){
                    k <- ind + Ntip
                }else{
                    stop(sprintf("Please specify either 'choose.node' (between %d and %d) or 'choose.node.label' (from 'phy$node.label').", Ntip+1, Ntot))
                }
            }else{
                stop(sprintf("Please specify either 'choose.node' (between %d and %d) or 'choose.node.label' (from 'phy$node.label').", Ntip+1, Ntot))
            }
        }else{
            stop(sprintf("Please specify either 'choose.node' (between %d and %d).", Ntip+1, Ntot))
        }
    }
    
    
    ######################################################################################
    
    if(k == (Ntip+1)){
        subphy <- phy
    }else{
        ## create a temporary 'phy_tmp' and append node.label (if not there)
        phy_tmp <- phy
        if(is.null(phy_tmp$node.label)){
            phy_tmp$node.label <- (Ntip+1):Ntot
        }
        ## extract all children
        connectivity <- suppressMessages(dcTreeConnectivity(phy_tmp, verbose=verbose))
        all_children <- which(connectivity[k-Ntip,]==1)
        ## get all tips that need to remove
        len.tip.remove <- setdiff(1:Ntip, all_children)
        ## get subtree
        subtree <- ape::drop.tip(phy_tmp, len.tip.remove)
        ## get which internal nodes are extracted
        node.extracted <- match(subtree$node.label, phy_tmp$node.label)
        
        ## finalise the extracted subtree
        subphy <- subtree
        subphy$Nnode <- ape::Nnode(subphy)
        if(!is.null(phy$node.label)){
            subphy$node.label <- phy$node.label[node.extracted]
        }else{
            subphy$node.label <- NULL
        }
    }
    
    if(verbose){
        message(sprintf("From the input tree (with %d tips) and under the internal node '%s', a subtree (with %d tips) has been extracted.", Ntip, k, ape::Ntip(subphy)), appendLF=T)
    }
    
    invisible(subphy)
}
