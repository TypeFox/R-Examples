#' Function to find common ancestors of two terms/nodes from a direct acyclic graph (DAG)
#'
#' \code{dDAGancestor} is supposed to find a list of common ancestors shared by two terms/nodes, given a direct acyclic graph (DAG; an ontology). If two terms are given as NULL, then a sparse matrix of children x ancestors is built for all terms. If one of them is null, then a sparse matrix of children x ancestors is built but only for non-null input terms.
#'
#' @param g an object of class "igraph" or "graphNEL"
#' @param term1 the first term/node as input
#' @param term2 the second term/node as input
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to true for display
#' @return 
#' \itemize{
#'  \item{When two terms are given: a list of terms/nodes that are common ancestors for two input terms/nodes}
#'  \item{When two terms are given as NULL: a sparse matrix of children x ancestors is built for all terms, with '1' for the reachable and otherwise '0'.}
#'  \item{When one of terms is given as NULL: a sparse matrix of children x ancestors is built but only for non-null input terms, with '1' for the reachable and otherwise '0'.}
#' }
#' @note none
#' @export
#' @seealso \code{\link{dDAGinduce}}
#' @include dDAGancestor.r
#' @examples
#' # 1) load HPPA as igraph object
#' ig.HPPA <-dRDataLoader(RData='ig.HPPA')
#' g <- ig.HPPA
#'
#' # 2) randomly give two terms
#' term1 <- sample(V(g)$name,1)
#' term2 <- sample(V(g)$name,1)
#'
#' # 3) find common ancestors
#' dDAGancestor(g, term1, term2)

dDAGancestor <- function (g, term1=NULL, term2=NULL, verbose=T)
{
    
    if(class(g)=="graphNEL"){
        ig <- igraph.from.graphNEL(g)
    }else{
        ig <- g
    }
    if (class(ig) != "igraph"){
        stop("The function must apply to either 'igraph' or 'graphNEL' object.\n")
    }
 
    allterms <- V(ig)$name
    if(is.null(term1) & is.null(term2)){
    
        if(verbose){
            message(sprintf("Build a sparse matrix of children x ancestors (with %d rows and %d columns (%s)...", length(allterms), length(allterms), as.character(Sys.time())), appendLF=T)
        }

        ## find all ancestors for any node
        neighs.in <- igraph::neighborhood(ig, order=vcount(ig), nodes=V(ig), mode="in")
        ## store in a sparse matrix of children X ancestors
        sCP <- Matrix::Matrix(0, nrow=length(allterms), ncol=length(allterms), sparse=T)
        for(i in 1:length(neighs.in)){
            sCP[i,neighs.in[[i]]] <- 1
        }
        rownames(sCP) <- colnames(sCP) <- allterms
          
        return(sCP)
    
    }else if( (!is.null(term1) & is.null(term2)) | (is.null(term1) & !is.null(term2)) ){
        
        if(( !is.null(term1) & is.null(term2)) ){
            terms <- term1
        }else if( is.null(term1) & !is.null(term2) ){
            terms <- term2
        }
        
        ## checking input term1
        terms <- terms[!is.na(terms)]
        flag <- terms %in% V(ig)$name
        if(sum(flag)!=0){
            terms <- terms[flag]
        }else{
            terms <- V(ig)$name
        }
        
        if(verbose){
            message(sprintf("Build a sparse matrix of children x ancestors (with %d rows and %d columns (%s)...", length(terms), length(allterms), as.character(Sys.time())), appendLF=T)
        }
        
        ## find all ancestors for any node
        neighs.in <- igraph::neighborhood(ig, order=vcount(ig), nodes=terms, mode="in")
        ## store in a sparse matrix of children X ancestors
        sCP <- Matrix::Matrix(0, nrow=length(terms), ncol=length(allterms), sparse=T)
        for(i in 1:length(neighs.in)){
        	sCP[i,unlist(neighs.in[i])] <- 1       	
            #sCP[i,neighs.in[[i]]] <- 1
        }
        rownames(sCP) <- terms
        colnames(sCP) <- allterms
        
        return(sCP)
        
    }else{
    
        if(sum(c(term1,term2) %in% allterms)!=2){
            stop("The function requires that both input terms can be found in the input graph.\n")
        }
    
        if(verbose){
            message(sprintf("Find common ancestors shared by %s and %s (%s)...", term1, term2, as.character(Sys.time())), appendLF=T)
        }
    
        subg1 <- dDAGinduce(g=ig, nodes_query=term1, path.mode="all_paths")
        subg2 <- dDAGinduce(g=ig, nodes_query=term2, path.mode="all_paths")
    
        ancestors <- intersect(V(subg1)$name, V(subg2)$name)
        
        return(ancestors)
    }
}


   
