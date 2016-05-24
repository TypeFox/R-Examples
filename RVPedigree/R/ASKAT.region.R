##' Runs the ASKAT method on a given genomic region
##'
##' @title Run the ASKAT method on a genomic region defined by a start
##'     and a stop base pair coordinate
##' @param regionname (character) Name of the region/gene on which you
##'     are running the association test. This name is used in the
##'     output of this function and can be used to distinguish
##'     different regions if this function is run multiple times.
##' @param U (optional) Matrix of Eigenvectors of the relationship matrix
##'     obtained from spectral decomposition of the relationship matrix:
##'     \eqn{\Phi = U S U^T}. If this parameter is not given, it will
##'     be computed, so when running this function for many regions
##'     time can be saved by specifying not only \code{Phi}, but also
##'     \code{S} and \code{U}.
##' @param S (optional) Matrix of Eigenvalues of the relationship matrix
##'     obtained from spectral decomposition of the relationship matrix:
##'     \eqn{\Phi = U S U^T}. If this parameter is not given, it will
##'     be computed, so when running this function for many regions,
##'     time can be saved by specifying not only \code{Phi}, but also
##'     \code{S} and \code{U}.
##' @param RH.Null (optional) output of \code{\link{Estim.H0.ASKAT}}
##'     function. In analyses of many regions, it is not necessary to calculate the the null
##'     hypothesis for each region. One estimation per trait is
##'     enough.
##' @inheritParams read.haplo
##' @inheritParams RVPedigree
##' @return A data frame containing the results of the association
##'     test. The data frame contains the following columns:
##'     \itemize{
##'     \item \code{Score.Test}: the score of the given association test
##'     \item \code{P.value}: the p-value of the association test
##'     \item \code{N.Markers}: the number of markers in the region
##'     \item \code{regionname}: Name of the region/gene on which you
##'     are running the association test
##'     }
##' @author Lennart C. Karssen, Sodbo Sharapov
##' @export
ASKAT.region <- function(y=NULL,
                         X=NULL,
                         Phi=NULL,
                         type="bed",
                         filename=NULL,
                         map=NULL,
                         chr=0,
                         startpos=0,
                         endpos=0,
                         regionname=NULL,
                         U=NULL,
                         S=NULL,
                         RH.Null=NULL,
                         weights=NULL){
    y <- check_pheno(y)
    check_covariates(X, y)
    check_relmatrix(Phi)
    check_files(filename, type)
    check_positions(startpos, endpos)
    check_weights(weights)

    if (is.null(U) | is.null(S)) {
        ## The eigen vectors and/or eigen values of the kinship matrix
        ## have not been given so we need to compute them.
        Phi.eig <- eigen(Phi)
        U       <- Phi.eig$vectors
        S       <- Phi.eig$values
    }

    ## Read genotypes in the selected genomic region
    haplotypes <- read.haplo(type=type,
                          filename=filename,
                          map=map,
                          chr=chr,
                          startpos=startpos,
                          endpos=endpos)

    if (ncol(haplotypes) == 0) {
        warning("No genotypes available in the region from ",
                      startpos, " to ", endpos, " on chromosome ", chr)
        result.df <- data.frame(Score.Test=NA,
                                P.value=NA,
                                N.Markers=NA)
        if (!is.null(regionname)) {
            rownames(result.df) <- regionname
        }
        return(result.df)
    }

    G <- Get.G(haplotypes)

    weights <- compute.weights(G, weights)

    if(is.null(RH.Null)) {
        RH.Null <- Estim.H0.ASKAT(y=y, X=X, S=S, U=U)
    }

    pval <- pvalue.ASKAT(RH0=RH.Null,
                         y=y,
                         X=X,
                         Phi=Phi,
                         W=weights,
                         G=G)

    result.df <- data.frame(Score.Test=pval$score,
                            P.value=pval$p.value,
                            N.Markers=ncol(haplotypes))

    if (!is.null(regionname)) {
        rownames(result.df) <- regionname
    }

    return(result.df)
}
