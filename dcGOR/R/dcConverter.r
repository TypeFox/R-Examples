#' Function to convert an object between graph classes
#'
#' \code{dcConverter} is supposed to convert an object between classes 'Onto' and 'igraph', or between 'Dnetwork' and 'igraph', or between 'Cnetwork' and 'igraph'.
#'
#' @param obj an object of class "Onto", "igraph", "Dnetwork" or "Cnetwork"
#' @param from a character specifying the class converted from. It can be one of "Onto", "igraph", "Dnetwork" and "Dnetwork"
#' @param to a character specifying the class converted to. It can be one of "Onto", "igraph", "Dnetwork" and "Dnetwork"
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to true for display
#' @return an object of class "Onto", "igraph", "Dnetwork" or "Cnetwork"
#' @note Conversion is also supported between classes 'Onto' and 'igraph', or between 'Dnetwork' and 'igraph', or between 'Cnetwork' and 'igraph'
#' @export
#' @seealso \code{\link{dcRDataLoader}}, \code{\link{Onto-class}}, \code{\link{Dnetwork-class}}, \code{\link{Cnetwork-class}}
#' @include dcConverter.r
#' @examples
#' \dontrun{
#' # 1) conversion between 'Onto' and 'igraph'
#' # 1a) load onto.GOMF (as 'Onto' object)
#' on <- dcRDataLoader('onto.GOMF')
#' on
#' # 1b) convert the object from 'Onto' to 'igraph' class
#' ig <- dcConverter(on, from='Onto', to='igraph')
#' ig
#' # 1c) convert the object from 'igraph' to 'Onto' class
#' dcConverter(ig, from='igraph', to='Onto')
#'
#' # 2) conversion between 'Dnetwork' and 'igraph'
#' # 2a) computer a domain semantic network (as 'Dnetwork' object)
#' g <- dcRDataLoader('onto.GOMF')
#' Anno <- dcRDataLoader('SCOP.sf2GOMF')
#' dag <- dcDAGannotate(g, annotations=Anno, path.mode="shortest_paths", verbose=FALSE)
#' alldomains <- unique(unlist(nInfo(dag)$annotations))
#' domains <- sample(alldomains,5) # randomly sample 5 domains
#' dnetwork <- dcDAGdomainSim(g=dag, domains=domains, method.domain="BM.average", method.term="Resnik", parallel=FALSE, verbose=FALSE)
#' dnetwork
#' # 2b) convert the object from 'Dnetwork' to 'igraph' class
#' ig <- dcConverter(dnetwork, from='Dnetwork', to='igraph')
#' ig
#' # 2c) convert the object from 'igraph' to 'Dnetwork' class
#' dcConverter(ig, from='igraph', to='Dnetwork')
#' }

dcConverter <- function (obj, from=c("Onto","igraph","Dnetwork","Cnetwork"), to=c("igraph","Onto","Dnetwork","Cnetwork"), verbose=TRUE)
{
    
    from <- match.arg(from)
    to <- match.arg(to)
    
    if (class(obj) != from){
        stop(sprintf("The class of your input object '%s' is '%s', mismatched as you intended (from='%s').\n", deparse(substitute(obj)), class(obj), from))
    }
    
    if(from!="igraph" & to!="igraph"){
        stop(sprintf("Conversion between '%s' and '%s' is not supported.\n", from, to))
    }
    
    if ((from=="Onto" & to=="Dnetwork") | (from=="Dnetwork" & to=="Onto")){
        #stop(sprintf("Conversion between '%s' and '%s' is not supported.\n", from, to))
    }
    
    if(from==to){
        warnings(sprintf("Since the class '%s' converted from is the same as the class '%s' converted to, it will return exactly what you input.\n", from, to))
        return(obj)
    }
    
    if(from=="igraph"){
        
        ## get node data frame
        data <- igraph::get.data.frame(obj, what="vertices")
        
        ## make sure it is data.frame after removing the first column
        if(ncol(data)==2){
            df <- data.frame(data[,2])
            colnames(df) <- colnames(data)[2]
            rownames(df) <- rownames(data)
            nodeI <- new("InfoDataFrame", data=df)
        }else if(ncol(data)==1){
            nodeI <- new("InfoDataFrame", data=data)
        }else{
            nodeI <- new("InfoDataFrame", data=data.frame(data[,-1]))
        }
        
        ## get adjacency matrix
        if ("weight" %in% list.edge.attributes(obj)){
            adjM <- igraph::get.adjacency(obj, type="both", attr="weight", edges=F, names=T, sparse=getIgraphOpt("sparsematrices"))
        }else{
            adjM <- igraph::get.adjacency(obj, type="both", attr=NULL, edges=F, names=T, sparse=getIgraphOpt("sparsematrices"))
        }
        
        ## convert to either "Onto" or "Dnetwork"
        if(to=="Onto"){
            ## for Onto
            objConverted <- new("Onto", adjMatrix=adjM, nodeInfo=nodeI)
        }else if(to=="Dnetwork"){
            ## for Dnetwork
            objConverted <- new("Dnetwork", adjMatrix=adjM, nodeInfo=nodeI)
        }else if(to=="Cnetwork"){
            ## for Cnetwork
            objConverted <- new("Cnetwork", adjMatrix=adjM, nodeInfo=nodeI)
        }
        
    }else if(from=="Onto" | from=="Dnetwork" | from=="Cnetwork"){
        
        ## node info
        nodes <- nInfo(obj)
        nodes <- data.frame(name=rownames(nodes), nodes)
        nodenames <- nodeNames(obj)
        
        ## adjacency matrix
        adjM <- adjMatrix(obj)
        tmp <- which(adjM!=0, arr.ind=T)
        
        ## un-direct graph for "Dnetwork"
        if(from=="Dnetwork" | from=="Cnetwork"){
            tmp <- tmp[tmp[,1]<tmp[,2],]
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
        if(from=="Onto"){
            ## for Onto
            objConverted <- igraph::graph.data.frame(d=relations, directed=T, vertices=nodes)
        }else if(from=="Dnetwork" | from=="Cnetwork"){
            ## for Dnetwork
            objConverted <- igraph::graph.data.frame(d=relations, directed=F, vertices=nodes)
            
            #objConverted <- igraph::graph.adjacency(as.matrix(adjM), mode="undirected", weighted=TRUE, diag=FALSE)
        }
        
    }
    
    if(verbose){
        message(sprintf("Your input object '%s' of class '%s' has been converted into an object of class '%s'.", deparse(substitute(obj)), from, to), appendLF=T)
    }
    
    return(objConverted)
}