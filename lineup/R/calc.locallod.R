## calc.locallod.R
## Karl W Broman

# calc.locallod
#
#' Calculate LOD score at physical position of each gene
#'
#' For gene expression data with physical positions of the genes, calculate the
#' LOD score at those positions to assess evidence for local eQTL.
#'
#' \code{cross} and \code{pheno} must contain exactly the same individuals in
#' the same order.  (Use \code{\link{findCommonID}} to line them up.)
#'
#' We consider the expression phenotypes in batches: those whose closest
#' pseudomarker is the same.
#'
#' We use Haley-Knott regression to calculate the LOD scores.
#'
#' Actually, we use a bit of a contortion of the data to force the
#' \code{\link[qtl]{scanone}} function in R/qtl to calculate the LOD score at a
#' single position.
#'
#' We omit any transcripts that map to the X chromosome; we can only handle
#' autosomal loci for now.
#'
#' @param cross An object of class \code{"cross"} containing data for a QTL
#' experiment.  See the help file for \code{\link[qtl]{read.cross}} in the
#' R/qtl package (\url{http://www.rqtl.org}).  There must be a phenotype named
#' \code{"id"} or \code{"ID"} that contains the individual identifiers.
#' @param pheno A data frame of phenotypes (generally gene expression data),
#' stored as individuals x phenotypes.  The row names must contain individual
#' identifiers.
#' @param pmark Pseudomarkers that are closest to the genes in \code{pheno}, as
#' output by \code{\link{find.gene.pseudomarker}}.
#' @param addcovar Additive covariates passed to \code{\link{scanone}}.
#' @param intcovar Interactive covariates passed to \code{\link{scanone}}.
#' @param verbose If TRUE, print tracing information.
#' @param n.cores Number of CPU cores to use in the calculations. With
#' \code{n.cores=0}, \code{\link[parallel]{detectCores}} is used to
#' detect the number of available cores.
#'
#' @return A vector of LOD scores.  The names indicate the gene names (columns in
#' \code{pheno}).
#' @author Karl W Broman, \email{kbroman@@biostat.wisc.edu}
#' @seealso \code{\link{find.gene.pseudomarker}}, \code{\link{plotEGclass}},
#' \code{\link{findCommonID}}, \code{\link{disteg}}
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
#' # line up f2cross and expr1
#' id <- findCommonID(f2cross, expr1)
#'
#' # calculate LOD score for local eQTL
#' locallod <- calc.locallod(f2cross[,id$first], expr1[id$second,], pmark)
#'
#' @export
calc.locallod <-
    function(cross, pheno, pmark, addcovar=NULL, intcovar=NULL, verbose=TRUE,
             n.cores=1)
{
    if(any(pmark$chr == "X")) {
        warning("Dropping X chr loci; we can only handle autosomes for now.")
        pmark <- pmark[pmark$chr != "X",]
    }

    if(qtl::nind(cross) != nrow(pheno))
        stop("cross and pheno have incompatible numbers of individuals.")
    m <- match(colnames(pheno), rownames(pmark))
    if(any(is.na(m))) pheno <- pheno[,!is.na(m),drop=FALSE]
    m <- match(rownames(pmark), colnames(pheno))
    if(any(is.na(m))) pmark <- pmark[!is.na(m),,drop=FALSE]
    pheno <- pheno[,match(rownames(pmark), colnames(pheno)),drop=FALSE]

    cpmark <- paste(pmark$chr, pmark$pmark, sep=":")
    upmark <- unique(cpmark)

    lod <- rep(NA, ncol(pheno))
    temp <- cross
    n.ind <- qtl::nind(cross)

    # function to do calculation at unique pseudomarkers
    tmpf <- function(i) {
        wh <- which(cpmark == upmark[i])

        y <- pheno[,wh,drop=FALSE]
        gp <- cross$geno[[pmark$chr[wh[1]]]]$prob[,pmark$pmark[wh[1]],,drop=FALSE]

        # create dummy cross
        temp$pheno <- cbind(y, cross$pheno)
        temp$geno <- list("1"=list(data=cbind("m1"=rep(1, n.ind)), map=c("m1"=1),prob=gp))
        attr(temp$geno[[1]]$prob, "map") <- c("m1"=1)
        class(temp$geno[[1]]) <- "A"
        lod[wh] <- unlist(qtl::scanone(temp, method="hk", addcovar=addcovar, intcovar=intcovar,
                                       pheno.col=1:ncol(y))[-(1:2)])
    }

    # if n.cores == 0, detect available cores
    if(n.cores == 0) n.cores <- parallel::detectCores()

    if(n.cores > 1) {
        if(Sys.info()[1] == "Windows") { # Windows doesn't support mclapply
            cl <- parallel::makeCluster(n.cores)
            on.exit(parallel::stopCluster(cl))
            result <- parallel::clusterApply(cl, seq(along=upmark), tmpf)
        } else {
            result <- parallel::mclapply(seq(along=upmark), tmpf, mc.cores=n.cores)
        }
    } else {
        result <- lapply(seq(along=upmark), tmpf)
    }

    # fill in results
    for(i in seq(along=upmark)) {
        wh <- which(cpmark == upmark[i])
        lod[wh] <- result[[i]]
    }

    names(lod) <- colnames(pheno)
    lod
}
