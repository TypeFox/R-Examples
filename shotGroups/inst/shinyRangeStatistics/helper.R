library(shotGroups)
library(shiny)

#####---------------------------------------------------------------------------
## option sets and their respective inverse
#####---------------------------------------------------------------------------

rangeStat    <- c("Extreme spread"="1", "Figure of Merit"="2", "Bounding Box Diagonal"="3")
rangeStatInv <- c("1"="ES", "2"="FOM", "3"="D")

rangeStatSig    <- c("Rayleigh \\(\\sigma\\)"="1", "Extreme spread"="2",
                     "Figure of Merit"="3", "Bounding Box Diagonal"="4")
rangeStatSigInv <- c("1"="Rayleigh", "2"="ES", "3"="FOM", "4"="D")

unitsDst    <- c("m"="1",  "yard"="2", "feet"="3")
unitsDstInv <- c("1"="m",  "2"="yd", "3"="ft")

unitsXY     <- c("cm"="1", "mm"="2", "inch"="3")
unitsXYInv  <- c("1"="cm", "2"="mm", "3"="in")

unitsAbs    <- c("m"="1", "cm"="2", "mm"="3", "yard"="4", "feet"="5", "inch"="6")
unitsAbsInv <- c("1"="m", "2"="cm", "3"="mm", "4"="yd", "5"="ft", "6"="in")

unitsAng    <- c("MOA"="1", "SMOA"="2", "mrad"="3", "NATO mil"="4")
unitsAngInv <- c("1"="MOA", "2"="SMOA", "3"="mrad", "4"="mil")

#####---------------------------------------------------------------------------
## function to write list to tab-delimited text file
#####---------------------------------------------------------------------------

textOut <- function(x, name, f) {
    if(is.list(x) && !is.data.frame(x) && !inherits(x, "htest")) {
        Map(textOut, x, names(x), f)
    } else if(!inherits(x, "htest") && !isS4(x)) {
        name <- gsub("^([[:digit:]]+)(.*)", "group\\1\\2", name)
        suppressWarnings(write.table(data.frame(x=name, stringsAsFactors=FALSE),
                                     file=f, row.names=FALSE,
                                     col.names=FALSE, quote=FALSE, append=TRUE))
        if(is.matrix(x) || is.data.frame(x)) {
            colnames(x) <- gsub("^([[:digit:]]+)(.*)", "group\\1\\2", colnames(x))
            colnames(x) <- gsub(" \\(", "_CIlo", colnames(x))
            colnames(x) <- gsub(" \\)", "_CIup", colnames(x))
            if(is.matrix(x)) {
                x <- if(all(c("x", "y") %in% rownames(x))) {
                    data.frame(coord=rownames(x), x, stringsAsFactors=FALSE)
                } else {
                    data.frame(unit=rownames(x), x, stringsAsFactors=FALSE)
                }
            }
            x <- suppressWarnings(rbind(x, ""))
        } else if(is.vector(x)) {
            x <- rbind(t(x), "")
        }

        suppressWarnings(write.table(x, file=f, sep="\t", row.names=FALSE,
                                     quote=FALSE, append=TRUE))
    }

    return(invisible(NULL))
}
