##' Meta-analysis of a series of results from GWA scans
##'
##' Performes meta-analysis of results of multiple GWA studies
##' stored in files.
##'
##' The function looks for files named "POPSextens" in the directory
##' "dir", reads them and pool consecutively using
##' \code{\link{metagwa.tables}} function.
##'
##' The source files should contain a number of variables, such as allelic
##' coding, code for the effective allele, etc. Please refer to the
##' formetascore function of GenABEL package for details.
##'
##' When the file is not available for the first population, or there are
##' non-unique population names, the function stops with an error message.
##'
##' If there are no files corresponding to next populations, these are
##' skipped in analysis with a warning message (see details in the
##' description of returned values).
##'
##' @param dir Path to directory containing files with GWA results
##' @param pops A vector specifying study names, which make the first
##' part of the file names
##' @param extens Fixed extension to the file name
##' @param maf Filter threshold for the absolute number/frequency of
##' minor allele.  If the parameter is < 1, this MAF threshold is used
##' for filtering before meta-analysis; if this parameter is >1, n*MAF
##' -- absolute number of the minor allele copies -- is used
##' @param call Filter threshold for SNP call rate
##' @param phwe Filter threshold for SNP P-value for HWE
##' @param precorrect Should GC be applied to the original data before pooling
##' @param correct.pooled Whether to apply Genomic Control correction to
##'    the study named "POOLED"
##' @return   Results are dumped to file named "POOLEDextens" in the directory
##'  "dir".  Additionally, a list is returned with elements
##'  \item{analysed.pops}{list of populations successfully analysed}
##' @author Yurii Aulchenko
##' @seealso \code{\link{metagwa.tables}}
##' @keywords htest
##' @export
"metagwa.files" <-
    function(dir=".", pops, extens, maf=5, call=0.95, phwe=1.e-8,
             precorrect=TRUE, correct.pooled=FALSE) {

        analysed.pops <- rep(TRUE, length(pops))
        if (length(unique(pops)) != length(pops)) {
            stop("Some population names are not unique!")
        }

        cpop <- 1
        tryerr <- "try-error"
        while (tryerr == "try-error") {
            fname <- paste(dir, "/", pops[cpop], extens, sep="")
            cat("Population", pops[cpop], ", reading", fname)
            df <- try(read.csv(file=fname, header=TRUE, stringsAsFactors=FALSE))
            if (class(df) == "try-error") {
                warning(paste("File", fname,
                              "can not be open. Skipping population",
                              pops[cpop]), immediate.=TRUE)
                tryerr <- "try-error"
                analysed.pops[cpop] <- FALSE
                cpop <- cpop + 1
                if (cpop > length(pops)) stop("Ran out of populations...")
                next;
            }

            if (maf>1) {
                if (is.na(match("n", names(df)))) {
                    stop("Number of people (n) required to filter on number of rare allele copies")
                }

                norac <- pmin(2 * df$n * df$effallelefreq, 2 * df$n *
                              (1. - df$effallelefreq))
                goods <- which(df$call >= call &
                               norac >= maf &
                               df$pexhwe >= phwe)
                rm(norac)
            } else {
                goods <- which(df$call >= call &
                               pmin(df$effallelefreq,
                                    1. - df$effallelefreq) >= maf &
                               df$pexhwe >= phwe)
            }

            df <- df[goods, ]
            df <- df[!(df$allele1=="1" | df$allele2=="1"), ]
            cat(" done\n")
            cat("Dimesions after filters are", dim(df), "\n")

            if (dim(df)[1] <= 0) {
                tryerr <- "try-error"
                analysed.pops[cpop] <- FALSE
                cpop <- cpop + 1
                if (cpop > length(pops))
                    stop("Ran out of populations...")

                warning(paste("Data in file", fname,
                              ", population", pops[cpop],
                              ": no SNPs left after filtering. Skipping"),
                        immediate.=TRUE)
                next;
            }
            tryerr <- "no"
        }                               # END while (try-err)


        df$npops <- rep(0, dim(df)[1])
        df$npops[!is.na(df$beta)] <- 1
        names1 <- df$name
        for (pop in pops[(which(analysed.pops==TRUE)[2]):length(pops)]) {
            cat("population ", pop, ", ", sep="")
            fname <- paste(dir, "/", pop, extens, sep="")
            cat("reading", fname)
            gc()
            dft <- try(read.csv(file=fname, header=TRUE, stringsAsFactors=FALSE))
            if (class(dft) == "try-error") {
                warning(paste("File", fname,
                              "can not be open. Skipping population", pop),
                        immediate.=TRUE)
                analysed.pops[which(pops == pop)] <- FALSE
                next;
            }

            cat(" done\n")
            if (maf > 1) {
                if (is.na(match("n", names(dft))))
                    stop("Number of people (n) required to filter on number of rare allele copies")
                norac <- pmin(2 * dft$n*dft$effallelefreq, 2 * dft$n *
                              (1. - dft$effallelefreq))
                goods <- which(dft$call>=call & norac >= maf & dft$pexhwe >= phwe)
                rm(norac)
            } else {
                goods <- which(dft$call>=call &
                               pmin(dft$effallelefreq,
                                    1. - dft$effallelefreq) >= maf &
                               dft$pexhwe >= phwe)
            }

            dft <- dft[goods, ]
            dft <- dft[which(!(dft$allele1=="1" | dft$allele2=="1")), ]
            cat("Dimesions after filters are", dim(dft), "\n")

            if (dim(dft)[1] <= 0) {
                warning(paste("Data in file", fname,
                              ", population", pop,
                              ": no SNPs left after filtering. Skipping"),
                        immediate.=TRUE)
                analysed.pops[which(pops == pop)] <- FALSE
                next;
            }
            gc()

            if (pop==(pops[analysed.pops])[2]) {
                basepop <- (pops[analysed.pops])[1]
            } else {
                basepop <- "POOLED"
            }

            df <- metagwa.tables(df, dft, name.x=basepop, name.y=pop,
                                 precorrect=precorrect,
                                 correct.pooled=correct.pooled)
            cat("Dimesions after pooling are", dim(df), "\n")
        }

        write.csv(df, file=paste(dir, "/POOLED", extens, sep=""),
                  row.names=FALSE, quote=FALSE)
        out <- list()
        out$analysed.pops <- pops[analysed.pops]
        out
    }
