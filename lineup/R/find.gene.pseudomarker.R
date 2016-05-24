## find.gene.pseudomarker.R
## Karl W Broman

# find.gene.pseudomarker:
#' Find nearest peudomarker to each gene
#'
#' Pull out the pseudomarker that is closest to the position of each of a
#' series of genes.
#'
#' We first convert positions (by interpolation) from those contained within
#' \code{cross} to physical coordinates contained in \code{pmap}.  We then use
#' \code{\link[qtl]{find.pseudomarker}} to identify the closest pseudomarker to
#' each gene location.
#'
#' We also include the positions of the pseudomarkers, and we print a warning
#' message if pseudomarkers are > 2 Mbp from the respective gene.
#'
#' @param cross An object of class \code{"cross"} containing data for a QTL
#' experiment.  See the help file for \code{\link[qtl]{read.cross}} in the
#' R/qtl package (\url{http://www.rqtl.org}).
#' @param pmap A physical map of the markers in \code{cross}, with locations in
#' Mbp.  This is a list whose components are the marker locations on each
#' chromosome.
#' @param geneloc A data frame specifying the physical locations of the genes.
#' There should be two columns, \code{chr} for chromosome and \code{pos} for
#' position in Mbp.  The rownames should indicate the gene names.
#' @param where Indicates whether to pull pseudomarkers from the genotype
#' probabilities (produced by \code{\link[qtl]{calc.genoprob}}) or from the
#' imputed genotypes (produced by \code{\link[qtl]{sim.geno}}).
#' @return A data frame with columns \code{chr} (the chromosome) and
#' \code{pmark} (the name of the pseudomarker).  The third column \code{pos}
#' contains the Mbp position of the pseudomarker.  The final column is the
#' signed distance between the gene and the pseudomarker.  The rownames
#' indicate the gene names.
#' @author Karl W Broman, \email{kbroman@@biostat.wisc.edu}
#' @seealso \code{\link[qtl]{find.pseudomarker}},
#' \code{\link[qtl]{find.pseudomarkerpos}}, \code{\link{plotEGclass}},
#' \code{\link{disteg}}, \code{\link{calc.locallod}}
#' @keywords utilities
#' @examples
#' data(f2cross, expr1, genepos, pmap)
#' library(qtl)
#' \dontshow{
#' n_ind <- 20
#' n_genes <- 5
#' f2cross <- f2cross[,1:n_ind]
#' expr1 <- expr1[1:n_ind,1:n_genes]
#' genepos <- genepos[1:n_genes,]}
#' # calc QTL genotype probabilities
#' f2cross <- calc.genoprob(f2cross, step=1)
#'
#' # find nearest pseudomarkers
#' pmark <- find.gene.pseudomarker(f2cross, pmap, genepos, "prob")
#'
#' @export
find.gene.pseudomarker <-
    function(cross, pmap, geneloc, where=c("prob", "draws"))
{
    where <- match.arg(where)
    if(!(where %in% names(cross$geno[[1]])))
        stop("You first need to run ", ifelse(where=="prob", "calc.genoprob", "sim.geno"), ".")

    cross <- qtl::replacemap(cross, pmap)
    res <- data.frame(chr=geneloc$chr,
                      pmark=qtl::find.pseudomarker(cross, geneloc$chr, geneloc$pos, where, addchr=FALSE),
                      stringsAsFactors=FALSE)

    rownames(res) <- rownames(geneloc)

    pmark <- res$pmark
    gr <- grep("^loc[0-9]+\\.*[0-9]*(\\.[0-9]+)*$", pmark)
    if(length(gr)>0)
        pmark[gr] <- paste("c", res$chr[gr], ".", pmark[gr], sep="")
    upmark <- unique(pmark)
    thepos <- qtl::find.pseudomarkerpos(cross, upmark, where)
    res$pos <- thepos[match(pmark, rownames(thepos)),2]

    res <- cbind(res, dist.from.gene=(d <- geneloc$pos - res$pos))
    d <- abs(d)
    if(any(d > 2)) {
        ngap <- sum(d>2)
        maxd <- max(d)
        warning(ngap, " genes differ from pseudomarker pos by > 2 Mbp, with gaps as big as ", round(maxd, 1), " Mbp")
    }

    res
}
