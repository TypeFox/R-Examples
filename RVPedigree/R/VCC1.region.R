##' Runs the VC-C1 method on a given genomic region
##'
##' @title Run the VC-C1 method on a genomic region defined by a start
##'     and a stop base pair coordinate
##' @inheritParams ASKAT.region
##' @param RH.Null (optional) output of \code{\link{Estim.H0.VCC}}
##'     function. Practically, you don't need to calculate the null
##'     hypothesis for every region. One estimation per trait is
##'     enough.
##' @return A data frame containing the results of the association
##'     test. The data frame contains the following columns:
##'     \itemize{
##'     \item \code{Score.Test}: the score of the given association test
##'     \item \code{P.value}: the p-value of the association test
##'     \item \code{N.Markers}: the number of markers in the region
##'     \item \code{regionname}: Name of the region/gene on which you
##'     are running the association test
##'     }
##' @author Sodbo Sharapov, Lennart C. Karssen
##' @export
VCC1.region <- function(y=NULL,
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

    n <- length(y)

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
        RH.Null <- Estim.H0.VCC(y=y, X=X, U=U, S=S)
    }

    Prep = Preparation.VCC(n=n, RH0=RH.Null, Phi=Phi)

    pval <- pvalue.VCC1(Prep=Prep,
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
