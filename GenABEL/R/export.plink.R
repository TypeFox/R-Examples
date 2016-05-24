#' Export GenABEL data in PLINK format
#'
#' Export GenABEL data in PLINK format. This function is
#' a simple wrapper to the \code{\link{export.merlin}} function
#' with specific arguments + few lines of code to
#' export phenotypes
#'
#' @param data GenABEL data object of 'gwaa.data'-class to
#' be exported.
#'
#' @param filebasename base file name for exported data,
#' extensions '.ped', '.map' and '.phe' (for phenotype file)
#' are added for specific output files.
#'
#' @param phenotypes NULL (no phenotypes exported), "all" (default) for
#' all phenotypes or a vector of character with names of phenotypes
#' to be exported.
#'
#' @param transpose if TRUE (default), 'tped' files will be produced, else
#' 'ped' files are produced.
#'
#' @param export012na if TRUE, export in numeric (0, 1, 2, NA) format,
#' as opposed to ATGC format (default: FALSE).
#'
#' @param ... arguments passed to \code{\link{export.merlin}}.
#'
#' @author Yurii Aulchenko
#'
#' @keywords IO
#'

"export.plink" <- function(data, filebasename="plink", phenotypes= "all",
                           transpose=TRUE, export012na=FALSE, ...)
{

  if (!is.null(phenotypes)) {
    phef <- paste(filebasename, ".phe", sep="")
    phed <- phdata(data)
    phed <- data.frame(FID=seq(1:dim(phed)[1]),
                       IID=phed[, "id"],
                       phed[, which(names(phed) != "id")])

    if (phenotypes != "all") {
      phed <- phed[, c("FID", "IID", phenotypes)]
    }

    write.table(phed,
                file=phef,
                row.names=FALSE,
                col.names=TRUE,
                quote=FALSE,
                sep=" ")
  }

  if (!transpose) {
    ## Export to .ped and .map
    pedf <- paste(filebasename, ".ped", sep="")
    mapf <- paste(filebasename, ".map", sep="")

    export.merlin(data,
                  pedfile=pedf,
                  datafile=NULL,
                  mapfile=mapf,
                  format="plink", ... )
  } else {
    ## export TFAM
    sx <- male(data)
    sx[sx==0] <- 2
    tfam <- data.frame(FID=c(1:nids(data)),
                       IID=idnames(data),
                       father=0,
                       mother=0,
                       sex=sx,
                       trait=-9,
                       stringsAsFactors=FALSE)
    write.table(tfam,
                file=paste(filebasename, ".tfam", sep=""),
                row.names=FALSE,
                col.names=FALSE,
                quote=FALSE)
    ## export genotypic data (TPED)
    pedfilename <- paste(filebasename, ".tped", sep="")
    tmp <- .Call("export_plink_tped",
                 as.character(snpnames(data)),
                 as.character(chromosome(data)),
                 as.double(map(data)),
                 as.raw(gtdata(data)@gtps),
                 as.integer(nsnps(data)),
                 as.integer(nids(data)),
                 as.character(coding(data)),
                 as.character(pedfilename),
                 as.logical(export012na))
  }
}
