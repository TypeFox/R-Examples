#' Function to convert an object between graph classes
#'
#' \code{xConverter} is supposed to convert an object between classes 'dgCMatrix' and 'igraph'.
#'
#' @param obj an object of class "dgCMatrix" or "igraph"
#' @param from a character specifying the class converted from. It can be one of "dgCMatrix" and "igraph"
#' @param to a character specifying the class converted to. It can be one of "dgCMatrix" and "igraph"
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to true for display
#' @return an object of class "dgCMatrix" or "igraph"
#' @note Conversion is also supported between classes 'dgCMatrix' and 'igraph'
#' @export
#' @seealso \code{\link{xRDataLoader}}
#' @include xConverter.r
#' @examples
#' # generate a ring graph
#' g <- make_ring(10, directed=TRUE)
#' 
#' # convert the object from 'igraph' to 'dgCMatrix' class
#' xConverter(g, from='igraph', to='dgCMatrix')
#' 
#' \dontrun{
#' # Conversion between 'dgCMatrix' and 'igraph'
#' # ig.EF (an object of class "igraph" storing as a directed graph)
#' g <- xRDataLoader('ig.EF')
#' g
#' 
#' # convert the object from 'igraph' to 'dgCMatrix' class
#' s <- xConverter(g, from='igraph', to='dgCMatrix')
#' s[1:10,1:10]
#'
#' # convert the object from 'dgCMatrix' to 'igraph' class
#' ig <- xConverter(s, from="dgCMatrix", to="igraph")
#' ig
#' }

xConverter <- function (obj, from=c("dgCMatrix","igraph"), to=c("igraph","dgCMatrix"), verbose=TRUE)
{
    
    from <- match.arg(from)
    to <- match.arg(to)
    
    if (class(obj) != from){
        #stop(sprintf("The class of your input object '%s' is '%s', mismatched as you intended (from='%s').\n", deparse(substitute(obj)), class(obj), from))
    }
    
    if(from!="igraph" & to!="igraph"){
        stop(sprintf("Conversion between '%s' and '%s' is not supported.\n", from, to))
    }
    
    if(from==to){
        warnings(sprintf("Since the class '%s' converted from is the same as the class '%s' converted to, it will return exactly what you input.\n", from, to))
        return(obj)
    }
    
    if(from=="igraph"){
        
        ## get node data frame
        data <- igraph::get.data.frame(obj, what="vertices")
        
        ## get adjacency matrix
        if ("weight" %in% list.edge.attributes(obj)){
            objConverted <- igraph::get.adjacency(obj, type="both", attr="weight", edges=F, names=T, sparse=getIgraphOpt("sparsematrices"))
        }else{
            objConverted <- igraph::get.adjacency(obj, type="both", attr=NULL, edges=F, names=T, sparse=getIgraphOpt("sparsematrices"))
        }
        
    }else if(from=="dgCMatrix"){
        
        ## node info
        nodes <- data.frame(name=rownames(obj))
        nodenames <- rownames(obj)
        
        ## adjacency matrix
        adjM <- obj
        tmp <- which(as.matrix(adjM!=0), arr.ind=T)
        
        ## un-direct graph
        if(from=="dgCMatrix"){
        	ind <- which(tmp[,1]<tmp[,2])
        	ttmp <- matrix(0, nrow=length(ind), ncol=2)
            ttmp[1:length(ind),] <- tmp[ind,]
            tmp <- ttmp
        }
		
		        
        ## weighted or not
        weight_flag <- T
        if(all(adjM[tmp]==1)){
            weight_flag <- F
        }
        if(weight_flag){
            relations <- data.frame(from=nodenames[tmp[,1]], to=nodenames[tmp[,2]], weight=adjM[tmp])
        }else{
            relations <- data.frame(from=nodenames[tmp[,1]], to=nodenames[tmp[,2]])
        }
        
        ## convert to "igraph"
        if(from=="dgCMatrix"){
            objConverted <- igraph::graph.data.frame(d=relations, directed=F, vertices=nodes)
        }
        
    }
    
    if(verbose){
        message(sprintf("Your input object '%s' of class '%s' has been converted into an object of class '%s'.", deparse(substitute(obj)), from, to), appendLF=T)
    }
    
    return(objConverted)
}
