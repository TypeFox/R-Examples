library(shiny)
library(shotGroups)

#####---------------------------------------------------------------------------
## option sets and their respective inverse
#####---------------------------------------------------------------------------

dataBuiltIn <- c("1"="DF300BLK", "2"="DF300BLKhl", "3"="DFcciHV", "4"="DFcm",
                 "5"="DFsavage", "6"="DFscar17", "7"="DFtalon")

dataBuiltInInv <- as.list(names(dataBuiltIn))
dataBuiltInInv <- setNames(dataBuiltInInv, dataBuiltIn)

unitsDst    <- c("m"="1",  "yard"="2", "feet"="3")
unitsDstInv <- c("1"="m",  "2"="yd", "3"="ft")

unitsXY     <- c("cm"="1", "mm"="2", "inch"="3")
unitsXYInv  <- c("1"="cm", "2"="mm", "3"="in")

unitsAbs    <- c("m"="1", "cm"="2", "mm"="3", "yard"="4", "feet"="5", "inch"="6")
unitsAbsInv <- c("1"="m", "2"="cm", "3"="mm", "4"="yd", "5"="ft", "6"="in")

unitsAng    <- c("degree"="1", "MOA"="2", "SMOA"="3", "radian"="4", "milliradian"="5", "NATO mil"="6")
unitsAngInv <- c("1"="deg",    "2"="MOA", "3"="SMOA", "4"="rad",    "5"="mrad",        "6"="mil")

unitsPlot    <- c("cm"="1", "mm"="2", "inch"="3", "degree"="4", "MOA"="5", "SMOA"="6", "radian"="7", "milliradian"="8", "NATO mil"="9")
unitsPlotInv <- c("1"="cm", "2"="mm", "3"="in",   "4"="deg",    "5"="MOA", "6"="SMOA", "7"="rad",    "8"="mrad",        "9"="mil")

hitpRUnits   <- c('same as xy-coords', 'm', 'cm', 'mm', 'yd', 'ft', 'in', 'deg', 'MOA', 'SMOA', 'rad', 'mrad', 'mil')
hitpRUnit    <- as.list(seq_along(hitpRUnits))
hitpRUnit    <- setNames(hitpRUnit, hitpRUnits)
hitpRUnitInv <- as.list(c("unit", hitpRUnits[-1]))
hitpRUnitInv <- setNames(hitpRUnitInv, seq_along(hitpRUnit))

CEPtypes <- c("Correlated Normal"="1", "Grubbs-Pearson"="2", "Grubbs-Patnaik"="3",
              "Grubbs-Liu"="4", "Rayleigh"="5", "Krempasky"="6", "Ignani"="7",
              "RMSE"="8", "Ethridge"="9", "RAND-234"="10", "Valstar"="11")

CEPtypesInv <- c("1"="CorrNormal", "2"="GrubbsPearson", "3"="GrubbsPatnaik",
                 "4"="GrubbsLiu", "5"="Rayleigh", "6"="Krempasky", "7"="Ignani",
                 "8"="RMSE", "9"="Ethridge", "10"="RAND", "11"="Valstar")

CItypes    <- c("basic"="1", "percentile"="2", "BCa"="3", "normal"="4")
CItypesInv <- c("1"="basic", "2"="perc", "3"="bca", "4"="norm")

## function to extract group names and make them an option set
getGroups <- function(x, choices=FALSE) {
    x       <- droplevels(x)
    groups  <- levels(x$series)
    groups  <- setNames(groups, seq_along(groups))
    groupsC <- seq_along(groups)
    if(choices) {
        return(setNames(groupsC, groups))
    } else {
        return(groups)
    }
}

#####---------------------------------------------------------------------------
## output components that need a CI / CEP / confidence level
#####---------------------------------------------------------------------------

CEPOutComps <- c("CEP", "confEll", "confEllRob")
CIOutComps  <- c("sigmaCI", "sigmaMRci", "sdXci", "sdYci", "RSDci", "MRci", "ctrXci", "ctrYci")

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
