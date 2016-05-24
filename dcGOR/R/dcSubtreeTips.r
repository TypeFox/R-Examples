#' Function to extract a tip-induced subtree from a phylo-formatted phylogenetic tree
#'
#' \code{dcSubtreeTips} is supposed to extract a tip-induced subtree from a phylo-formatted phylogenetic tree. In addition to the tree in subject, another input is a vector containing tip labels of interest. From valid tip lables, there are two types of subtree to extract. One is first induce clade (an internal node) from tip labels, and then the subtree is extracted under the induced clade. Another type is to extract a subtree only containing given tip labels; in this situation, some internal nodes perhaps need to further trimmed. The resulting subtree is also represented as an object of class 'phylo'.
#'
#' @param phy an object of class 'phylo'
#' @param choose.tip.labels a character specifying which tips are chosen
#' @param subtree.type a character specifying how to extract subtree from given tips. It can be 'clade' or 'tips_only'. The former is first induce clade (an internal node) from tip labels, and then to extract the subtree under the induced clade. The latter is to directly extract the subtree only containing given tip labels, (if necessary), after trimming out unnecessary internal nodes
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to TRUE for display
#' @return
#' an object of class 'phylo'
#' @note
#' nonde
#' @export
#' @seealso \code{\link{dcTreeConnectivity}}, \code{\link{dcSubtreeClade}}
#' @include dcSubtreeTips.r
#' @examples
#' # 1) with internal node labels
#' tree <- "(((t1:5,t2:5)i3:2,(t3:4,t4:4)i4:3)i2:2,(t5:4,t6:4)i5:6)i1;"
#' phy <- ape::read.tree(text=tree)
#' ape::plot.phylo(phy, type="p", use.edge.length=TRUE, show.node.label=TRUE)
#'
#' # 2) tip labels of interest
#' choose.tip.labels <- c('t1','t2','t3')
#' # 2a) extract subtree via an induced clade 
#' subphy <- dcSubtreeTips(phy, choose.tip.labels, subtree.type="clade")
#' ape::plot.phylo(subphy, type="p", use.edge.length=TRUE, show.node.label=TRUE)
#' # 2b) extract subtree containing only tips
#' subphy <- dcSubtreeTips(phy, choose.tip.labels, subtree.type="tips_only")
#' ape::plot.phylo(subphy, type="p", use.edge.length=TRUE, show.node.label=TRUE)

dcSubtreeTips <- function(phy, choose.tip.labels=NULL, subtree.type=c("clade","tips_only"), verbose=T)
{
    
    ## match.arg matches arg against a table of candidate values as specified by choices, where NULL means to take the first one
    subtree.type <- match.arg(subtree.type)
    
    if (class(phy) != "phylo"){
        stop("The input 'phy' must belong to the class 'phylo'!")
    }
    
    ind <- match(choose.tip.labels, phy$tip.label)
    ind_query <- ind[!is.na(ind)]
    if(length(ind_query)<=1){
        stop(sprintf("Please choose at least two tip labels (from 'phy$node.label')."))
    }
    
    Ntip <- ape::Ntip(phy)
    Nnode <- ape::Nnode(phy)
    Ntot <- Ntip+Nnode
    
    ######################################################################################
    
    if(subtree.type=="tips_only"){
        ## create a temporary 'phy_tmp' and append node.label (if not there)
        phy_tmp <- phy
        if(is.null(phy_tmp$node.label)){
            phy_tmp$node.label <- (Ntip+1):Ntot
        }
    
        ## get all tips that need to remove
        len.tip.remove <- setdiff(1:Ntip, ind_query)
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
    }else if(subtree.type=="clade"){
        ## extract all parents
        connectivity <- suppressMessages(dcTreeConnectivity(phy, verbose=verbose))
        parents_flag <- apply(connectivity[,ind_query]==1,1,sum)==length(ind_query)
        all_parents <- as.numeric(names(which(parents_flag)))
        ## get most recent common ancestor
        clade <- max(all_parents)
        
        ## extract a subtree under a given clade
        subphy <- suppressMessages(dcSubtreeClade(phy, choose.node=clade, verbose=verbose))
    }
        
    if(verbose){
        message(sprintf("From the input tree (with %d tips) and %d valid tip labels, a subtree (with %d tips) has been extracted.", Ntip, length(ind_query), ape::Ntip(subphy)), appendLF=T)
    }
    
    invisible(subphy)
}
