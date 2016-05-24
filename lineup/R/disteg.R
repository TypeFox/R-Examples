## disteg.R
## Karl W Broman

# disteg
#
#' Calculate distance between two gene expression data sets
#'
#' Calculate a distance between all pairs of individuals for two gene
#' expression data sets
#'
#' We consider the expression phenotypes in batches, by which pseudomarker they
#' are closest to.  For each batch, we pull the genotype probabilities at the
#' corresponding pseudomarker and use the individuals that are in common
#' between \code{cross} and \code{pheno} and whose maximum genotype probability
#' is above \code{min.genoprob}, to form a classifier of eQTL genotype from
#' expression values, using k-nearest neighbor (the function
#' \code{\link[class]{knn}}). The classifier is applied to all individuals with
#' expression data, to give a predicted eQTL genotype. (If the proportion of
#' the k nearest neighbors with a common class is less than
#' \code{min.classprob}, the predicted eQTL genotype is left as \code{NA}.)
#'
#' If \code{repeatKNN} is TRUE, we repeat the construction of the k-nearest
#' neighbor classifier after first omitting individuals whose proportion of
#' mismatches between observed and inferred eQTL genotypes is greater than
#' \code{max.selfd}.
#'
#' Finally, we calculate the distance between the observed eQTL genotypes for
#' each individual in \code{cross} and the inferred eQTL genotypes for each
#' individual in \code{pheno}, as the proportion of mismatches between the
#' observed and inferred eQTL genotypes.
#'
#' If \code{weightByLinkage} is \code{TRUE}, we use weights on the mismatch
#' proportions for the various eQTL, taking into account their linkage. Two
#' tightly linked eQTL will each be given half the weight of a single isolated
#' eQTL.
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
#' @param min.genoprob Threshold on genotype probabilities; if maximum
#' probability is less than this, observed genotype taken as \code{NA}.
#' @param k Number of nearest neighbors to consider in forming a k-nearest
#' neighbor classifier.
#' @param min.classprob Minimum proportion of neighbors with a common class to
#' make a class prediction.
#' @param classprob2drop If an individual is inferred to have a genotype
#' mismatch with classprob > this value, treat as an outlier and drop from the
#' analysis and then repeat the KNN construction without it.
#' @param repeatKNN If TRUE, repeat k-nearest neighbor a second time, after
#' omitting individuals who seem to not be self-self matches
#' @param max.selfd Min distance from self (as proportion of mismatches between
#' observed and predicted eQTL genotypes) to be excluded from the second round
#' of k-nearest neighbor.
#' @param phenolabel Label for expression phenotypes to place in the output
#' distance matrix.
#' @param weightByLinkage If TRUE, weight the eQTL to account for their
#' relative positions (for example, two tightly linked eQTL would each count
#' about 1/2 of an isolated eQTL)
#' @param map.function Used if \code{weightByLinkage} is TRUE
#' @param verbose if TRUE, give verbose output.
#' @return A matrix with \code{nind(cross)} rows and \code{nrow(pheno)}
#' columns, containing the distances.  The individual IDs are in the row and
#' column names.  The matrix is assigned class \code{"lineupdist"}.
#'
#' The names of the genes that were used to construct the classifier are saved
#' in an attribute \code{"retained"}.
#'
#' The observed and inferred eQTL genotypes are saved as attributes
#' \code{"obsg"} and \code{"infg"}.
#'
#' The denominators of the proportions that form the inter-individual distances
#' are in the attribute \code{"denom"}.
#' @author Karl W Broman, \email{kbroman@@biostat.wisc.edu}
#' @seealso \code{\link{distee}}, \code{\link{summary.lineupdist}},
#' \code{\link{pulldiag}}, \code{\link{omitdiag}}, \code{\link{findCommonID}},
#' \code{\link{find.gene.pseudomarker}}, \code{\link{calc.locallod}},
#' \code{\link{plot.lineupdist}}, \code{\link[class]{knn}},
#' \code{\link{plotEGclass}}
#' @keywords utilities
#' @examples
#' library(qtl)
#'
#' # load example data
#' data(f2cross, expr1, pmap, genepos)
#' \dontshow{
#' keep <- c(1:20, 197, 553, 573, 740, 794, 822, 1474, 1522,
#'           1591, 1645, 2080, 2643, 2984, 3089, 3672, 4010, 4039,
#'           4159, 4191, 4198, 4213, 4401, 4544, 4593, 4925)
#' expr1 <- expr1[,keep]
#' genepos <- genepos[keep,]}
#'
#' # calculate QTL genotype probabilities
#' f2cross <- calc.genoprob(f2cross, step=1)
#'
#' # find nearest pseudomarkers
#' pmark <- find.gene.pseudomarker(f2cross, pmap, genepos)
#'
#' # line up individuals
#' id <- findCommonID(f2cross, expr1)
#'
#' # calculate LOD score for local eQTL
#' locallod <- calc.locallod(f2cross[,id$first], expr1[id$second,], pmark)
#'
#' # take those with LOD > 25
#' expr1s <- expr1[,locallod>25,drop=FALSE]
#'
#' # calculate distance between individuals
#' #     (prop'n mismatches between obs and inferred eQTL geno)
#' d <- disteg(f2cross, expr1s, pmark)
#'
#' # plot distances
#' plot(d)
#'
#' # summary of apparent mix-ups
#' summary(d)
#'
#' # plot of classifier for and second eQTL
#' par(mfrow=c(2,1), las=1)
#' plotEGclass(d)
#' plotEGclass(d, 2)
#'
#' @useDynLib lineup
#' @export
disteg <-
    function(cross, pheno, pmark, min.genoprob=0.99,
             k=20, min.classprob=0.8, classprob2drop=1, repeatKNN=TRUE,
             max.selfd=0.3, phenolabel="phenotype",
             weightByLinkage=FALSE,
             map.function=c("haldane", "kosambi", "c-f", "morgan"),
             verbose=TRUE)
{
    # individuals in common between two data sets
    theids <- findCommonID(cross, pheno)
    if(length(theids$first) < k)
        stop("You need at least ", k, " individuals in common between the data sets.")

    # make sure pheno and pmark line up
    m <- findCommonID(colnames(pheno), rownames(pmark))
    pheno <- pheno[,m$first,drop=FALSE]
    pmark <- pmark[m$second,,drop=FALSE]

    # drop X chromosome from cross
    chrtype <- sapply(cross$geno, class)
    if(any(chrtype=="X")) {
        warning("Dropping X chromosome")
        cross <- subset(cross, names(cross$geno)[chrtype=="A"])
    }
    crosschr <- names(cross$geno)

    # find transcript chr in cross
    m <- match(pmark$chr, crosschr)
    if(any(is.na(m))) {
        warning("Dropping ", sum(is.na(m)), " transcripts with unknown chromosome assignment.")
        pheno <- pheno[,!is.na(m), drop=FALSE]
        pmark <- pmark[!is.na(m),, drop=FALSE]
    }
    if(ncol(pheno) < 1)
        stop("Need at least one expression phenotype.")

    # unique eQTL locations
    cpmark <- paste(pmark$chr, pmark$pmark, sep=":")
    upmark <- unique(cpmark)

    m <- match(upmark, cpmark)
    thechr <- pmark$chr[m]
    thepos <- pmark$pos[m]

    # construct weights to account for linkage
    linkwts <- rep(1, length(thechr))
    if(weightByLinkage) {
        uchr <- unique(thechr)

        map.function <- match.arg(map.function)
        mf <- switch(map.function,
                     "haldane" = qtl::mf.h,
                     "kosambi" = qtl::mf.k,
                     "c-f" = qtl::mf.cf,
                     "morgan" = qtl::mf.m)

        for(i in uchr) {
            if(sum(thechr==i)==1) next # just one eQTL on this chromosome

            d <- thepos[thechr==i]
            D <- matrix(1-2*mf(abs(outer(d, d, "-"))), ncol=length(d))
            linkwts[thechr==i] <- (length(d) + 1 - colSums(D))/length(d)
        }
    }

    # to contain observed and inferred eQTL genotypes
    obsg <- matrix(ncol=length(upmark), nrow=qtl::nind(cross))
    infg <- matrix(ncol=length(upmark), nrow=nrow(pheno))
    colnames(obsg) <- colnames(infg) <-
        apply(pmark[match(upmark, cpmark),c(1,3)], 1, paste, collapse="@")
    rownames(obsg) <- qtl::getid(cross)
    rownames(infg) <- rownames(pheno)

    # loop over eQTL; use k-nearest neighbor to classify
    if(verbose) cat("First pass through knn\n")
    ysave <- vector("list", length(upmark))
    names(ysave) <- colnames(obsg)
    for(i in seq(along=upmark)) {
        wh <- which(cpmark == upmark[i])
        pmarkchr <- pmark$chr[wh[1]]
        pmarkpmark <- pmark$pmark[wh[1]]

        ysave[[i]] <- y <- pheno[,wh,drop=FALSE]
        gp <- cross$geno[[pmarkchr]]$prob[,pmarkpmark,,drop=FALSE]
        gi <- apply(gp, 1, function(a) which(a==max(a, na.rm=TRUE)))
        gmx <- apply(gp, 1, max, na.rm=TRUE)
        gi[gmx < min.genoprob] <- NA
        obsg[,i] <- gi

        ysub <- y[theids$second,,drop=FALSE]
        gisub <- gi[theids$first]
        keep <- !is.na(gisub) & (rowSums(is.na(ysub)) == 0) # have genotype and all phenotypes
        keep2 <- rowSums(is.na(y)) == 0 # have all phenotypes

        knnout <- class::knn(ysub[keep,,drop=FALSE], ysub[keep,,drop=FALSE], gisub[keep],
                             k=k, l=ceiling(k*min.classprob), prob=TRUE)
        pr <- attr(knnout, "prob")
        okeep <- keep
        keep[gisub[keep] != knnout & pr >= classprob2drop] <- FALSE
        if(verbose && sum(okeep) > sum(keep))
            cat(" -- Classifier ", i, ": dropping ", sum(okeep) - sum(keep), " outliers.\n", sep="")

        infg[keep2,i] <- class::knn(ysub[keep,,drop=FALSE], y[keep2,,drop=FALSE], gisub[keep],
                                    k=k, l=ceiling(k*min.classprob))

    }


    if(repeatKNN) {
        if(verbose) cat("Calculate self-self distances\n")
        # calculate self-self distances
        pd <- rep(NA, nrow(theids$mat))
        names(pd) <- rownames(theids$mat)
        for(i in rownames(theids$mat)[theids$mat$inBoth])
            pd[i] <- mean(obsg[i,] != infg[i,], na.rm=TRUE)

        # bad individuals
        bad <- names(pd)[!is.na(pd) & pd>= max.selfd]

        # repeat the k-nearest neighbor classification without the bad individuals
        if(verbose) cat("Second pass through knn\n")
        for(i in seq(along=upmark)) {
            wh <- which(cpmark == upmark[i])

            y <- pheno[,wh,drop=FALSE]
            gi <- obsg[,i]

            ysub <- y[theids$second,,drop=FALSE]
            gisub <- gi[theids$first]
            keep <- !is.na(gisub) & (rowSums(is.na(ysub)) == 0) & is.na(match(rownames(ysub), bad))
            keep2 <- rowSums(is.na(y)) == 0

            infg[!keep2,i] <- NA

            knnout <- class::knn(ysub[keep,,drop=FALSE], ysub[keep,,drop=FALSE], gisub[keep],
                                 k=k, l=ceiling(k*min.classprob), prob=TRUE)
            pr <- attr(knnout, "prob")
            okeep <- keep
            keep[gisub[keep] != knnout & pr >= classprob2drop] <- FALSE
            if(verbose && sum(okeep) > sum(keep))
                cat(" -- Classifier ", i, ": dropping ", sum(okeep) - sum(keep), " outliers.\n", sep="")

            infg[keep2,i] <- class::knn(ysub[keep,,drop=FALSE], y[keep2,,drop=FALSE], gisub[keep],
                                        k=k, l=ceiling(k*min.classprob))
        }
    }

    # calculate final distance
    if(verbose) cat("Calculate distance matrix\n")

    z <- .C("R_propmismatch",
            as.integer(ncol(obsg)),
            as.integer(nrow(obsg)),
            as.integer(t(obsg)),
            as.integer(nrow(infg)),
            as.integer(t(infg)),
            as.double(linkwts),
            prop=as.double(rep(NA, nrow(obsg)*nrow(infg))),
            denom=as.double(rep(NA, nrow(obsg)*nrow(infg))),
            PACKAGE="lineup",
            NAOK=TRUE)

    d <- matrix(z$prop, ncol=nrow(infg))
    denom <- matrix(z$denom, ncol=nrow(infg))
    dimnames(denom) <- dimnames(d) <- list(rownames(obsg), rownames(infg))

    attr(d, "d.method") <- "prop.mismatch"
    attr(d, "labels") <- c("genotype", phenolabel)
    if(repeatKNN) {
        attr(d, "orig.selfd") <- pd
        attr(d, "badind") <- bad
    }
    attr(d, "obsg") <- obsg
    attr(d, "infg") <- infg
    attr(d, "y") <- ysave
    attr(d, "denom") <- denom
    attr(d, "linkwts") <- linkwts
    attr(d, "genonames") <- qtl::getgenonames(class(cross)[1], "A", "simple",
                                              qtl::getsex(cross), attributes(cross))
    class(d) <- c("eg.lineupdist", "lineupdist")

    d
}
