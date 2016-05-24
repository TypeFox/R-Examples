library(shiny)
library(shotGroups)
#library(mvoutlier)
#library(energy)

#####---------------------------------------------------------------------------
## number of diagrams for each panel
#####---------------------------------------------------------------------------

nShapePlots   <- 7
nSpreadPlots  <- 3
nComparePlots <- 4

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

unitsPlot    <- c("cm"="1", "mm"="2", "inch"="3", "degree"="4", "MOA"="5", "SMOA"="6", "radian"=7, "milliradian"="8", "NATO mil"="9")
unitsPlotInv <- c("1"="cm", "2"="mm", "3"="in",   "4"="deg",    "5"="MOA", "6"="SMOA", "7"="rad",  "8"="mrad",        "9"="mil")

CEPtypes <- c("Correlated Normal"="1", "Grubbs-Pearson"="2", "Grubbs-Patnaik"="3",
              "Grubbs-Liu"="4", "Rayleigh"="5", "Krempasky"="6", "Ignani"="7",
              "RMSE"="8", "Ethridge"="9", "RAND-234"="10", "Valstar"="11")

CEPtypesInv <- c("1"="CorrNormal", "2"="GrubbsPearson", "3"="GrubbsPatnaik",
                 "4"="GrubbsLiu", "5"="Rayleigh", "6"="Krempasky", "7"="Ignani",
                 "8"="RMSE", "9"="Ethridge", "10"="RAND", "11"="Valstar")

CItypes    <- c("basic"="1", "percentile"="2", "BCa"="3", "normal"="4")
CItypesInv <- c("1"="basic", "2"="perc", "3"="bca", "4"="norm")

shapeOut <- c("cor XY"="1", "cor XY robust"="2", "Outliers"="3",
              "Shapiro X"="4", "Shapiro Y"="5", "mult normality"="6")

shapeOutInv <- c("1"="corXY", "2"="corXYrob", "3"="Outliers",
                 "4"="ShapiroX", "5"="ShapiroY", "6"="multNorm")

spreadOut <- c("sd X, Y"="1", "sd X CI"="2", "sd Y CI"="3", "sd X, Y robust"="4",
               "cov XY"="5",
               "cov XY robust"="6", "dist to center"="7", "sigma CI"="8", "RSD CI"="9",
               "MR CI"="10", "max spread"="11", "group rect"="12",
               "min group rect"="13", "min circle"="14", "conf ellipse"="15",
               "conf ellipse robust"="16", "conf ellipse shape"="17",
               "conf ellipse shape robust"="18",
               "CEP"="19")

spreadOutInv <- c("1"="sdXY", "2"="sdXci", "3"="sdYci", "4"="sdXYrob",
                  "5"="covXY", "6"="covXYrob", "7"="distToCtr", "8"="sigmaCI",
                  "9"="RSDci", "10"="MRci", "11"="maxPairDist", "12"="groupRect",
                  "13"="groupRectMin", "14"="minCircleRad", "15"="confEll",
                  "16"="confEllRob", "17"="confEllShape", "18"="confEllShapeRob",
                  "19"="CEP")

locationOut <- c("center"="1", "center robust"="2", "distance POA"="3" ,
                 "distance POA robust"="4",
                 "Hotelling"="5", "center X CI"="6", "center Y CI"="7")

locationOutInv <- c("1"="ctr", "2"="ctrRob", "3"="distPOA", "4"="distPOArob",
                    "5"="Hotelling", "6"="ctrXci", "7"="ctrYci" )

compOut2 <- c("center"="1", "dist POA"="2", "MANOVA"="3", "cor XY"="4", "sd XY"="5",
              "sd XY CI"="6", "mean dist to center"="7", "max spread"="8",
              "bound box fig of merit"="9",
              "bound box diagonal"="10", "min circle radius"="11", "sigma"="12",
              "MR"="13", "sigma + MR CI"="14", "CEP"="15", "AnsariX"="16",
              "AnsariY"="17", "Wilcoxon"="18")

compOut2Inv <- c("1"="ctr", "2"="distPOA", "3"="MANOVA", "4"="corXY", "5"="sdXY",
                 "6"="sdXYci", "7"="meanDistToCtr", "8"="maxPairDist", "9"="bbFoM",
                 "10"="bbDiag", "11"="minCircleRad", "12"="sigma", "13"="MR",
                 "14"="sigmaMRci", "15"="CEP", "16"="AnsariX",
                 "17"="AnsariY", "18"="Wilcoxon")

compOut3plus <- c("ctr"="1", "distPOA"="2", "MANOVA"="3", "corXY"="4", "sdXY"="5",
                  "sdXYci"="6", "meanDistToCtr"="7", "maxPairDist"="8", "bbFoM"="9",
                  "bbDiag"="10", "minCircleRad"="11", "sigma"="12", "MR"="13",
                  "sigmaCI"="14", "CEP"="15", "FlignerX"="16",
                  "FlignerY"="17", "Kruskal"="18")

compOut3plusInv <- c("1"="ctr", "2"="distPOA", "3"="MANOVA", "4"="corXY", "5"="sdXY",
                     "6"="sdXYci", "7"="meanDistToCtr", "8"="maxPairDist", "9"="bbFoM",
                     "10"="bbDiag", "11"="minCircleRad", "12"="sigma", "13"="MR",
                     "14"="sigmaMRci", "15"="CEP", "16"="FlignerX",
                     "17"="FlignerY", "18"="Kruskal")

targetL    <- as.list(c(1, seq_along(targets)+1))
targetL    <- setNames(targetL, c("none", names(targets)))
targetLinv <- as.list(c(NA, names(targets)))
targetLinv <- setNames(targetLinv, c(1, seq_along(targets)+1))

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

# ## print list of objects to tab-delimited file
# l <- list(a=1:5,
#           b=data.frame(x=rnorm(3), y=LETTERS[1:3]),
#           c=list(p=matrix(1:6, nrow=3),
#                  q=LETTERS[1:4]))
# names(l$a) <- LETTERS[1:5]
# colnames(l$c$p) <- c("V1", "V2")
# names(l$c$q) <- paste0("X", 1:4)
# f <- "textOut.txt"
# Map(textOut, l, names(l), f)
