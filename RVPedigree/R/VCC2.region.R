##' Runs the VC-C2 method on a given genomic region
##'
##' @title Run the VC-C2 method on a genomic region defined by a start
##' and a stop base pair coordinate
##' @inheritParams read.haplo
##' @inheritParams VCC1.region
##' @param Nperm Integer, number of permutations to use for empirical
##' p-value estimation (default: 100).
##' @param Ncores (integer) Number of processor (CPU) cores to be used
##'     in parallel when doing the permutations to determine the
##'     p-value (default: 1).
##' @return A data frame containing the results of the association
##'     test. The data frame contains the following columns:
##'     \itemize{
##'     \item \code{Score.Test}: the score of the given association test
##'     \item \code{P.value}: the p-value of the association test
##'     \item \code{N.Markers}: the number of markers in the region
##'     \item \code{regionname}: Name of the region/gene on which you
##'     are running the association test
##'     }
##' @author Lennart C. Karssen
##' @importFrom kinship2 kindepth
##' @export
VCC2.region <- function(y=NULL,
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
                        weights=NULL,
                        Nperm=100,
                        Ncores=1){
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

    ## Read the pedigree structure from the Plink file
    pedigree <- read.pedigree(type=type, filename=filename)

    ## Calculate the generation ID/pedigree depth of each individual.
    ## This needs to be done on a per-family basis
    famids <- unique(pedigree[, 1])
    gen.id <- NULL
    for ( fam in 1:length(famids) ) {
        subped <- pedigree[pedigree[, 1] == famids[fam], 1:4]
        gen.id <- append(gen.id,
                         kinship2::kindepth(id=subped[, 2],
                                            dad.id=subped[, 3],
                                            mom.id=subped[, 4],
                                            align=TRUE))
    }

    ## Compute the residuals like in VC-C1

    if(is.null(RH.Null)) {
        RH.Null <- Estim.H0.VCC(y=y, X=X, U=U, S=S)
    }

    n <- length(y)
    prep.VCC  <- Preparation.VCC(n=n, RH0=RH.Null, Phi=Phi)

    pval <- pvalue.VCC2(P=prep.VCC$P,
                        G=G,
                        W=weights,
                        Nperm=Nperm,
                        n=n,
                        pedigree=pedigree,
                        haplotypes=haplotypes,
                        generation.id=gen.id,
                        Ncores=Ncores)

    result.df <- data.frame(Score.Test=pval$score,
                            P.value=pval$p.value,
                            N.Markers=ncol(haplotypes))

    if (!is.null(regionname)) {
        rownames(result.df) <- regionname
    }

    return(result.df)
}
