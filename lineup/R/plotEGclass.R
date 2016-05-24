## plotEGclass.R
## Karl W Broman

# plotEGclass
#
#' Plot classifier of eQTL genotype from expression data
#'
#' Diagnostic plot of one of the eQTL classifiers from the results of
#' \code{\link{disteg}}: generally expression phenotype against observed eQTL
#' genotype, colored by inferred eQTL genotype.
#'
#' The function produces a diagnostic plot for studying one of the k-nearest
#' neighbor classifiers underlying the output from \code{\link{disteg}}.
#'
#' In the case of one expression phenotype attached to the selected eQTL, the
#' plot is a dot plot of gene expression against observed eQTL genotype.
#'
#' In the case of two expression phenotypes, the plot is a scatterplot of the
#' two expression phenotypes against each other.
#'
#' In the case of more than two expression phenotypes, we use
#' \code{\link[graphics]{pairs}} to produce a matrix of scatterplots.
#'
#' @param d Output of \code{\link{disteg}}.
#' @param eqtl Numeric index or a character vector (of the form "1@@102.35")
#' indicating the eQTL to consider.
#' @param outercol Indicates how to color the outer edge of the points:
#' \code{"observed"} indicates to color based on observed genotypes;
#' \code{"inferred"} indicates to color based on inferred genotypes; otherwise,
#' give a color.
#' @param innercol Like \code{outercol}, but indicating the interior of the
#' points.
#' @param thecolors The colors to use in the plot.  The last element (after the
#' number of genotypes) indicates the color to use for missing values.
#' @param \dots Passed to \code{\link[graphics]{plot}} and
#' \code{\link[graphics]{points}}.
#' @return None.
#' @author Karl W Broman, \email{kbroman@@biostat.wisc.edu}
#' @seealso \code{\link{disteg}}, \code{\link{plot.lineupdist}},
#' \code{\link{plot2dist}}, \code{\link[class]{knn}}
#' @keywords graphics
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
#' # plot of classifier for and second eQTL
#' par(mfrow=c(2,1), las=1)
#' plotEGclass(d)
#' plotEGclass(d, 2)
#'
#' @importFrom graphics plot par pairs points mtext abline axis
#' @importFrom stats runif
#' @export
plotEGclass <-
    function(d, eqtl=1, outercol="inferred", innercol="observed",
             thecolors=c("#7B68ED", "#1B9E78", "#CA3767", "#E59E00"), ...)
{
    if(!("eg.lineupdist" %in% class(d)))
        stop("Input d must be as produced by disteg().")

    # inner and outer colors based on...
    choices <- c("observed", "inferred")
    if(!is.na(pm <- pmatch(outercol, choices)))
        outercol <- choices[pm]
    if(!is.na(pm <- pmatch(innercol, choices)))
        innercol <- choices[pm]

    ginf <- attr(d, "infg")
    gobs <- attr(d, "obsg")
    y <- attr(d, "y")
    gnames <- attr(d, "genonames")

    if(length(eqtl) > 1) {
        warning("eqtl should have length 1; using first element.")
        eqtl <- eqtl[1]
    }
    if(is.character(eqtl)) { # name of an eQTL
        eqtlnam <- eqtl
        eqtl <- match(eqtl, colnames(ginf))
        if(is.na(eqtl))
            stop("Can't find eQTL ", eqtl)
    }
    else {
        if(eqtl < 1 || eqtl > ncol(ginf))
            stop("eqtl must be between 1 and ", ncol(ginf))

        eqtlnam <- colnames(ginf)[eqtl]
    }

    # clean up eQTL name
    if(length(grep("@", eqtlnam)) == 1) {
        spl <- unlist(strsplit(eqtlnam, "@"))
        splnum <- unlist(strsplit(as.character(round(as.numeric(spl[2]),2)), "\\."))
        if(length(splnum)==1 || nchar(splnum[2]) == 0)
            splnum[2] <- "00"
        else if(nchar(splnum[2]) == 1)
            splnum[2] <- paste(splnum[2], "0", sep="")
        eqtlnam <- paste(spl[1], " @ ", splnum[1], ".", splnum[2], sep="")
    }

    ids <- findCommonID(rownames(gobs), rownames(ginf))

    ginf <- ginf[,eqtl]
    gobs <- gobs[,eqtl]
    y <- y[[eqtl]]
    hasy <- rowSums(is.na(y)) == 0

    gobs.sub <- gobs[ids$first]
    ginf.sub <- ginf[ids$second]
    y.sub <- y[ids$second,,drop=FALSE]

    ynog <- y[-ids$second,,drop=FALSE]
    ynog <- ynog[rowSums(is.na(ynog))==0,,drop=FALSE]

    if(ncol(y)==1) { # dot plot
        if(nrow(ynog) > 0) gnames <- c(gnames, "NA")

        xjit <- 0.2
        u <- runif(length(gobs.sub), -xjit, xjit)

        # point colors
        obscol <- rep(NA, length(gobs.sub))
        gobs.sub[is.na(gobs.sub)] <- 4
        for(i in seq(along=gnames))
            obscol[gobs.sub==i] <- thecolors[i]
        infcol <- rep(NA, length(ginf.sub))
        ginf.sub[is.na(ginf.sub)] <- 4
        for(i in seq(along=gnames))
            infcol[ginf.sub==i] <- thecolors[i]

        col <- switch(outercol, "observed"=obscol, "inferred"=infcol, rep(outercol, length(obscol)))
        bg <- switch(innercol, "observed"=obscol, "inferred"=infcol, rep(innercol, length(obscol)))

        plot(gobs.sub+u, y.sub, xlab="Observed genotype", ylab=colnames(y),
             main=eqtlnam, xaxt="n", xlim=c(0.5, length(gnames)+0.5),
             pch=21, col=col, bg=bg, ylim=range(y, na.rm=TRUE), type="n", ...)
        abline(v=seq(along=gnames), lty=2, col="gray60")
        axis(side=1, at=seq(along=gnames), gnames)

        # make sure the mismatches are on top
        wh <- !is.na(gobs.sub) & !is.na(ginf.sub) & gobs.sub==ginf.sub
        points((gobs.sub+u)[wh], y.sub[wh],
               pch=21, col=col[wh], bg=bg[wh], lwd=2, ...)
        points((gobs.sub+u)[!wh], y.sub[!wh],
               pch=21, col=col[!wh], bg=bg[!wh], lwd=2, ...)

        if(nrow(ynog) > 0) { # phenotype no genotype
            theinf <- ginf[rownames(ynog)]

            infcol <- rep("", length(nrow(ynog)))
            theinf[is.na(theinf)] <- 4
            for(i in seq(along=gnames))
                infcol[theinf==i] <- thecolors[i]
            obscol <- rep("black", length(infcol))

            col <- switch(outercol, "observed"=obscol, "inferred"=infcol, rep(outercol, length(obscol)))
            bg <- switch(innercol, "observed"=obscol, "inferred"=infcol, rep(innercol, length(obscol)))

            points(runif(nrow(ynog), 4-xjit, 4+xjit), ynog,
                   pch=21, col=col, bg=bg, lwd=2, ...)
        }

    }
    else if(ncol(y) == 2) { # scatter plot

        if(nrow(ynog) > 0) gnames <- c(gnames, "NA")

        # point colors
        obscol <- rep("", length(gobs.sub))
        gobs.sub[is.na(gobs.sub)] <- 4
        for(i in seq(along=gnames))
            obscol[gobs.sub==i] <- thecolors[i]
        infcol <- rep("", length(ginf.sub))
        ginf.sub[is.na(ginf.sub)] <- 4
        for(i in seq(along=gnames))
            infcol[ginf.sub==i] <- thecolors[i]

        col <- switch(outercol, "observed"=obscol, "inferred"=infcol, rep(outercol, length(obscol)))
        bg <- switch(innercol, "observed"=obscol, "inferred"=infcol, rep(innercol, length(obscol)))

        plot(y.sub[,1], y.sub[,2], xlab=colnames(y)[1], ylab=colnames(y)[2],
             main=eqtlnam, pch=21, col=col, bg=bg,
             xlim=range(y[,1], na.rm=TRUE), ylim=range(y[,2], na.rm=TRUE), type="n", ...)

        # make sure the mismatches are on top
        wh <- !is.na(gobs.sub) & !is.na(ginf.sub) & gobs.sub==ginf.sub
        points(y.sub[wh,1], y.sub[wh,2],
               pch=21, col=col[wh], bg=bg[wh], lwd=2, ...)
        points(y.sub[!wh,1], y.sub[!wh,2],
               pch=21, col=col[!wh], bg=bg[!wh], lwd=2, ...)

        if(nrow(ynog) > 0) { # phenotype no genotype
            theinf <- ginf[rownames(ynog)]

            infcol <- rep("", length(nrow(ynog)))
            theinf[is.na(theinf)] <- 4
            for(i in seq(along=gnames))
                infcol[theinf==i] <- thecolors[i]
            obscol <- rep("black", length(infcol))

            col <- switch(outercol, "observed"=obscol, "inferred"=infcol, rep(outercol, length(obscol)))
            bg <- switch(innercol, "observed"=obscol, "inferred"=infcol, rep(innercol, length(obscol)))

            points(ynog[,1], ynog[,2],
                   pch=21, col=col, bg=bg, lwd=2, ...)
        }
    }
    else { # pairs plot

        if(nrow(ynog) > 0) gnames <- c(gnames, "NA")

        # point colors
        obscol <- rep("", length(gobs.sub))
        gobs.sub[is.na(gobs.sub)] <- 4
        for(i in seq(along=gnames))
            obscol[gobs.sub==i] <- thecolors[i]
        infcol <- rep("", length(ginf.sub))
        ginf.sub[is.na(ginf.sub)] <- 4
        for(i in seq(along=gnames))
            infcol[ginf.sub==i] <- thecolors[i]

        if(nrow(ynog) > 0) { # phenotype no genotype
            theinf <- ginf[rownames(ynog)]

            infadd <- rep("", length(nrow(ynog)))
            theinf[is.na(theinf)] <- 4
            for(i in seq(along=gnames))
                infadd[theinf==i] <- thecolors[i]

            y.sub <- rbind(y.sub, ynog)
            obscol <- c(obscol, rep("black", nrow(ynog)))
            infcol <- c(infcol, infadd)
        }

        col <- switch(outercol, "observed"=obscol, "inferred"=infcol, rep(outercol, length(obscol)))
        bg <- switch(innercol, "observed"=obscol, "inferred"=infcol, rep(innercol, length(obscol)))

        o <- order(as.numeric((col!=bg)))
        y.sub <- y.sub[o,,drop=FALSE]
        col <- col[o]
        bg <- bg[o]

        par(oma=c(0,0,1.5,0))
        pairs(y.sub, pch=21, col=col, bg=bg, ...)
        mtext(side=3, outer=TRUE, eqtlnam)
        wh <- (col == bg)
    }
}
