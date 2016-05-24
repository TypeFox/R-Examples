#' Graph a wiki page
#' 
#' Creating a graph objet in witch an edge is a wiki page and a vertice is a
#' link beetween two wiki pages.
#' 
#' @param x Can be either a string defining a page title or a character vector,
#'   and the function will return the graph of all the links beetween al the
#'   pages included in the vector
#' @param domain The domain where is located the wiki
#' @param namespace The namespace pages the function will graph
#'   
#' @return An igraph object
#' @export
#' 
#' @import igraph
#'   
#' @family graph functions
#'   
#' @examples
#' \donttest{
#' # Graph of the 'Action' article in the french wiki.
#' graph_page('Action')
#' }
#' 
#' # Graphing a group of page
#' page <- c('Karl Marx','Classe sociale','Industrie')
#' # Return a graph where the 3 edges represents 'Karl Marx', 'Classe sociale' and 'Industrie', 
#' # and the vertices the link present or not beetween this pages.
#' g <- graph_page(page)

graph_page <- function(x, domain = "fr", namespace = "0") {
    if (length(x) == 1) {
        
        listLinks <- page_links(x, domain, namespace)
        listLinks <- c(listLinks$list, x)
        graph <- graph_page(listLinks, domain = domain, namespace = namespace)
        
    } else {
        
        graph <- make_empty_graph()
        graph <- add_vertices(graph, length(x), attr = list(title = x, type = "articles"))
        
        edgelist <- pbsapply(x, page_islink, x)
        edgelist <- sapply(edgelist, match, x)
        names(edgelist) <- NULL
        
        firstRow <- vector()
        
        for (i in 1:length(edgelist)) {
            firstRow <- c(firstRow, rep(i, length(edgelist[[i]])))
        }
        secondRow <- unlist(edgelist)
        edgelist <- matrix(c(firstRow, secondRow), ncol = length(firstRow), byrow = TRUE)
        
        
        graph <- add.edges(graph, edgelist)
        
    }
    
    return(graph)
} 
