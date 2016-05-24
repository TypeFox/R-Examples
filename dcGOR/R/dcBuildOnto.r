#' Function to build an object of the S4 class Onto from input files
#'
#' \code{dcBuildOnto} is supposed to build an object of of the S4 class \code{\link{Onto}}, given input files. These input files include 1) a file containing term relations, and 2) a file containing term/node information. 
#'
#' @param relations.file an input file containing term relations (i.e. edges from parent terms to child terms). For example, a file containing relations between GO Molecular Function (GOMF) terms can be found in \url{http://dcgor.r-forge.r-project.org/data/onto/igraph_GOMF_edges.txt}. As seen in this example, the input file must contain the header (in the first row) and two columns: 1st column for parent term ID, and 2nd column for child term ID. Note: the file should use the tab delimiter as the field separator between columns
#' @param nodes.file an input file containing term/node information. For example, a file containing GO Molecular Function (GOMF) terms can be found in \url{http://dcgor.r-forge.r-project.org/data/onto/igraph_GOMF_nodes.txt}. As seen in this example, the input file must contain the header (in the first row) and five columns: 1st column 'name' for node names (actually term ID; must be unique), 2nd column 'term_id' for term ID, 3rd 'term_name' for term name, 4th column 'term_namespace' for term namespace, and 5th column 'term_distance' for term distance. These five columns must be provided, the content in the first two columns are identical, and the content for the last column can be arbitrary (if it is hard to prepare). Note: the file should use the tab delimiter as the field separator between columns
#' @param output.file an output file used to save the built object as an RData-formatted file. If NULL, this file will be saved into "Onto.RData" in the current working local directory
#' @return 
#' Any use-specified variable that is given on the right side of the assigement sign '<-', which contains the built \code{Onto} object.
#' Also, an RData file specified in "output.file" is saved in the local directory.  
#' @note If there are no use-specified variable that is given on the right side of the assigement sign '<-', then no object will be loaded onto the working environment.
#' @export
#' @seealso \code{\link{Onto}}
#' @include dcBuildOnto.r
#' @examples
#' \dontrun{
#' # build an "Onto" object for GO Molecular Function
#' onto.GOMF <- dcBuildOnto(relations.file="http://dcgor.r-forge.r-project.org/data/onto/igraph_GOMF_edges.txt", nodes.file="http://dcgor.r-forge.r-project.org/data/onto/igraph_GOMF_nodes.txt", output.file="onto.GOMF.RData")
#' onto.GOMF
#' }

dcBuildOnto <- function(relations.file, nodes.file, output.file="Onto.RData")
{
    
    if(is.null(relations.file) | is.na(relations.file)){
        stop("The input file 'relations.file' must be provided!\n")
    }
    
    if(is.null(nodes.file) | is.na(nodes.file)){
        stop("The input file 'nodes.file' must be provided!\n")
    }
        
    if(is.null(output.file)){
        warnings("Since the output file is not provided, the function will use the default output file 'Onto.RData'!\n")
        output.file <- "Onto.RData"
    }
    
    relations <- utils::read.table(relations.file,sep="\t",quote="",header=T)
    nodes <- utils::read.table(nodes.file,sep="\t",quote="",header=T)
    
    ###########################
    # remove all those no-ASCII
    f <- function(x) iconv(x, "latin1", "ASCII", sub="")
    nodes <- apply(nodes,1:2,f)
    ###########################
    
    ig.obo <- igraph::graph.data.frame(d=relations, directed=T, vertices=nodes)
    x <- dcConverter(ig.obo, from='igraph', to='Onto', verbose=F)
    
    # remove the RData extension 
    output.var <- gsub(".RData$", "", output.file, ignore.case=T, perl=T)
    output.var <- gsub(".RDat$", "", output.var, ignore.case=T, perl=T)
    output.var <- gsub(".RDa$", "", output.var, ignore.case=T, perl=T)
    
    do.call(assign, list(output.var, x))
    save(list=output.var, file=output.file)
    
    
    if(file.exists(output.file)){
        message(sprintf("An object of S4 class 'Onto' has been built and saved into '%s'.", file.path(getwd(),output.file)), appendLF=T)
    }
    
    invisible(x)
}