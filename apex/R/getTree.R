#' Build phylogenies from multiple gene data
#'
#' This function builds separate phylogenetic trees for each gene in a \linkS4class{multidna} object, specifying a method for computing pairwise distances between individuals, and a method to build the tree from the distance matrix. By default, procedures from ape are used.
#'
#' @export
#'
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#'
#' @param x a \linkS4class{multidna} object.
#' @param pool a logical indicating if all genes should be pooled (concatenated) to obtain a single tree; defaults to FALSE.
#' @param genes an optional vector indicating the genes to retain for the concatenation; any way to subset the list in x@@dna is acceptable; by default, all genes are used.
#' @param model a character string passed to \code{\link[ape]{dist.dna}} describing the model to be used to compute genetic distances; defaults to 'N', the absolute number of mutations separating sequences.
#' @param pairwise.deletion a logical passed to \code{\link[ape]{dist.dna}} indicating if pairwise deletions should be used; the alternative is to remove all sites for which at least one missing value is present.
#' @param method a function building a tree from a matrix of pairwise genetic distances.
#' @param ladderize a logical indicating if the tree should be ladderized; defaults to TRUE.
#' @param negative.branch.length a logical indicating if negative branch lengths should be allowed (e.g. in the case of Neighbor-Joining reconstruction), or not, in which case they are set to 0 (FALSE, default).
#' @param ... further arguments passed to the tree reconstruction method defined by 'method'.
#'
#'
#' @aliases getTree
#'
#' @seealso \code{\link{dist.multidna}}
#'
#' @return a \code{multiPhylo} object
#'
#' @examples
#'
#' ## simple conversion with nicely ordered output
#' data(woodmouse)
#' genes <- list(gene1=woodmouse[,1:500], gene2=woodmouse[,501:965])
#' x <- new("multidna", genes)
#' x
#' plot(x)
#'
#' ## make trees, default parameters
#' trees <- getTree(x)
#' trees
#' plot(trees, type="unrooted")
#'
#' ## make one single tree based on concatenated genes
#' tree <- getTree(x, pool=TRUE)
#' tree
#' plot(tree, type="unrooted")
#'
getTree <- function(x, pool=FALSE, genes=TRUE, model = "N",
                    pairwise.deletion = TRUE, method=nj,
                    ladderize=TRUE, negative.branch.length=FALSE, ...){
    ## subset data
    x <- x[,genes]

    ## define tree building function
    f1 <- function(e){
        res <- method(dist.dna(e, model=model, pairwise.deletion=pairwise.deletion), ...)
        if(ladderize) res <- ladderize(res)
        if(!negative.branch.length) res$edge.length[res$edge.length<0] <- 0
        return(res)
    }

    ## if genes are  pooled
    if(pool){
        x <- concatenate(x)
        return(f1(x))
    }

    ## otherwise: one tree per gene
    out <- lapply(x@dna, f1)
    names(out) <- names(x@dna)
    class(out) <- "multiPhylo"
    return(out)
} # end getTree
