#' Create single-gene FGS
#'
#' Each network node X becomes an FGS itself, so that e.g. each ID 'X' generates a list entry named 'X', the content of which is c('X'), i.e. this same gene/protein. The total length of the FGS list then equals the number of nodes in the network. This is a way to create single-gene "quasi-pathways" for a more specific network enrichment analysis.

#' @param Net.list  the global network object, pre-created using \code{\link{import.net}}
#' @param Lowercase render gene/protein IDs lower-case (Default:1)
#' @seealso \code{\link{import.net}}, \code{\link{import.gs}}
#' @examples
#' data(net.kegg)
#' netpath <- net.kegg
#' net <- import.net(netpath)
#' fgs.genes <- as_genes_fgs(net)
#' print(fgs.genes[1:10])
#' @export

as_genes_fgs <- function(Net.list, Lowercase = 1) {
if (is.null(Net.list)) {stop("No list given...");}
if (Lowercase > 0) {n1 <- tolower(names(Net.list$links));} else {n1 <- names(Net.list$links);}
fgs.list 	<- as.list(n1);
names(fgs.list) <- n1;
return(fgs.list);
}
