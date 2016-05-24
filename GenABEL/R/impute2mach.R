#' converts IMPUTE to MACH files
#'
#' function to convert IMPUTE files to MACH format
#'
#' @param genofile IMPUTE genotype file name
#' @param infofile IMPUTE info file name
#' @param samplefile IMPUTE sample file name
#' @param machbasename base name for MACH-formatted outputs
#' @param maketextdosefile whether a text dosefile is to
#' be generated (if not, only filevector (*.fvi / *.fvd) files, usable
#' with ProbABEL/DatABEL, will be generated). Default: TRUE
#' @param ... arguments passed to \link{extract.annotation.impute}
#' (DO CHECK the documentation, otherwise your annotation may be
#' skrewed up!)
#'
#' @author Yurii Aulchenko
#'
#' @keywords IO manip
#'
#' @return nothing returned except files generated on the disk
#'
#'

impute2mach <- function(genofile, infofile, samplefile, machbasename,
                        maketextdosefile = TRUE, ... )
{
    if (!require(DatABEL))
        stop("this function requires DatABEL package to be installed")

    if (!is.character(machbasename)) stop("machbasename must be character")

    if (length(machbasename) == 1) {
        machdose <- paste(machbasename, ".machdose", sep="")
        machprob <- paste(machbasename, ".machprob", sep="")
        machinfo <- paste(machbasename, ".machinfo", sep="")
        machlegend <- paste(machbasename, ".machlegend", sep="")
    } else if (length(machbasename) == 4) {
        if (anyDuplicated(machbasename)) stop("names must be unique")
        machdose <- machbasename[1]
        machprob <- machbasename[2]
        machinfo <- machbasename[3]
        machlegend <- machbasename[4]
    } else stop("machbasename must be character of length 1 or 4")

    # create temporary DA file

      #print("before impute2databel");
    dfo <- impute2databel(genofile=genofile, samplefile=samplefile, outfile=genofile)
        #print("after impute2databel");

        #if (maketextdosefile) {
        #tmpname2 <- get_temporary_file_name()
        # transpose file
        # ...
        #dfo <- as(dfo, "matrix")
    #}
    # get annotattion
    #print("calling extract.annotation.impute")
    annot <- extract.annotation.impute(genofile=genofile, infofile=infofile, ... )
    #print(annot[1:5, ])

    # arrange MLINFO and legend file
    #SNP    Al1    Al2    Freq1    MAF    Quality    Rsq
    annot$MAF <- pmin(annot$Freq1, (1. - annot$Freq1))
    #print(annot[1:5, ])
    info_annot <- annot[, c("name", "A1", "A0", "Freq1", "MAF", "Quality", "Rsq")]
    #print(info_annot[1:5, ])
    write.table(info_annot, file=machinfo,
                row.names=FALSE,
                col.names=TRUE,
                quote=FALSE,
                sep="\t")
    legend_annot <- annot[, c("name", "pos", "A1", "A0")]
    #print(legend_annot[1:5, ])
    write.table(legend_annot, file=machlegend,
                row.names=FALSE,
                col.names=TRUE,
                quote=FALSE,
                sep="\t")

    if (maketextdosefile) {
        # arrange MLDOSE file
        ids <- get_dimnames(dfo)[[1]]
        if (file.exists(machdose)) unlink(machdose)

        outfile <- file(machdose, open="wt")

        for (i in 1:dim(dfo)[1])
        {
            # when using transposed DA object, use as.vector(dfo[, i]) (COLUMN!!!)
            outline <- c(ids[i], "MLDOSE", as(dfo[i, ], "vector"))
            #print(outline)
            #print(i)
            #print(class(outline))
            write(x=outline, file=outfile,
                  append=TRUE,
                  ncolumns=length(outline),
                  sep=" ")

            if ((i %% 100)==0 || i==dim(dfo)[1]) print(i)
        }
        close(outfile)
        #unlink(paste(tmpname2, "*", sep=""))
    }
    #print("AAA")
    disconnect(dfo)
    rm(dfo)
    gc()
}
