#' Estimate relationship matrix based on pedigree or genomic data
#'
#' Estimate relationship matrix based on pedigree or genomic data.
#' This function can use either Plink files as input (for both
#' pedigree-based and genomic relationship matrix calculation), or
#' genetic data in GenABEL format (for genomic relationships).
#'
#' Note that, with respect to the input parameters it is important to
#'     distinguish the two options of \code{datatype}. If
#'     \code{datatype} is equal to \code{"pedigree"} the user should
#'     specify the options \code{plinkbasefile}, \code{is.binary},
#'     \code{transpose}, \code{header} and \code{path2Plink}. In that
#'     case, the kinship will be estimated based on the pedigree data
#'     in the Plink files. If, however, \code{datatype} is equal to
#'     \code{"genomic"}, there are two options: either one uses the
#'     parameter \code{gwaa.data} to tell this function to compute the
#'     genomics relationship matrix from a previously stored GenABEL
#'     data object. Or, if \code{gwaa.data} is empty, the
#'     Plink-related parameters have to be specified and the genomics
#'     relationship matrix will be computed based on the genetic data
#'     in the Plink file.
#' @param datatype character, \code{"pedigree"} or \code{"genomic"}.
#'     Estimate the relationship matrix using either pedigree
#'     information or genotype data
#' @param plinkbasefile character, path to files in plink format. E.g.
#'     if you have files \code{test.ped} and \code{test.map},
#'     plinkbasefile should be \code{test}. More details are also given in the vignette.
#' @param is.binary logical, indicate whether the plink files are in
#'     binary format (\code{.bed}/\code{.bim}/\code{.fam})
#' @param transpose logical, indicate whether the plink text files are
#'     transposed or not (\code{.tped} and \code{.tfam} files)
#' @param header logical, indicate whether the input text files have
#'     header or not
#' @param path2Plink character, path to the binary (executable file)
#'     of plink_1.90 or later (the Plink_1.90 binary is used for
#'     efficient computation of the relationship matrix). More details are also given in the 
#'     vignette.
#' @param weight character, either \code{"no"} or \code{"freq"}. We
#'     suggest to use \code{weight="freq"}, which weighs by allelic
#'     frequency assuming HWE. See help for the \code{ibs()} function
#'     of the
#'     \href{https://cran.r-project.org/packages=GenABEL}{GenABEL
#'     package}
#' @param gwaa.data object, name of object gwaa.data-class, which is
#'     GenABEL genotype/phenotype data format
#' @param pedigreefile reserved for future development
#' @author Sodbo Sharapov
#' @return matrix object which contains the estimated relationship
#'     matrix
#' @examples
#' system.file("extdata", "2012.csv", package = "testdat")
#' \dontrun{
#' pedRel <- GetRelMatrix(datatype="pedigree",
#'                        plinkbasefile=system.file("extdata",
#'                                                  "data",
#'                                          package="RVPedigree"),
#'                        transpose=FALSE)
#' pedRel <- GetRelMatrix(datatype="pedigree",
#'                        plinkbasefile=system.file("extdata",
#'                                                  "dataT",
#'                                          package="RVPedigree"),
#'                        transpose=TRUE)
#' pedRel <- GetRelMatrix(datatype="pedigree",
#'                        plinkbasefile=system.file("extdata",
#'                                                  "dataB",
#'                                          package="RVPedigree"),
#'                        is.binary=TRUE)
#' pedRel <- GetRelMatrix(datatype="pedigree",
#'                        plinkbasefile=system.file("extdata",
#'                                                  "OneFamilyExample",
#'                                          package="RVPedigree"),
#'                        transpose=FALSE, header=TRUE)
#' pedRel <- GetRelMatrix(datatype="pedigree",
#'                        plinkbasefile=system.file("extdata",
#'                                                  "TwoFamilyExample",
#'                                          package="RVPedigree"),
#'                        transpose=FALSE, header=TRUE)
#'
#' load(system.file("extdata", "gwaa.data.RData")
#' genRel <- GetRelMatrix(datatype="genomic", gwaa.data=data1, weight="no")
#' genRel <- GetRelMatrix(datatype="genomic", gwaa.data=data1, weight="freq")
#' genRelError <- GetRelMatrix(datatype="genomic", gwaa.data=pedRel, weight="freq")
#'
#' genPlinkRel <- GetRelMatrix(data="genomic",
#'                             path2Plink="plink_1.90",
#'                             system.file("extdata",
#'                                         "data",
#'                                          package="RVPedigree"),
#' genPlinkRel[1:10, 1:10]
#'
#' genPlinkRel <- GetRelMatrix(data="genomic",
#'                             path2Plink="plink_1.90",
#'                             system.file("extdata",
#'                                         "dataT",
#'                                          package="RVPedigree"),
#'                             transpose=TRUE)
#' genPlinkRel[1:10, 1:10]
#'
#' genPlinkRel <- GetRelMatrix(data="genomic",
#'                             path2Plink="plink_1.90",
#'                             system.file("extdata",
#'                                         "dataB",
#'                                          package="RVPedigree"),
#'                             is.binary=TRUE)
#' genPlinkRel[1:10, 1:10]
#'
#' # GetRelMatrix(file="OneFamilyExample.ped", datatype="genomic")
#' # GetRelMatrix(file="OneFamilyExample.ped", datatype="pedasdfa")
#' # GetRelMatrix(file="OneFamilyExaample", datatype="pedigree")
#' }
#'
#' @export
GetRelMatrix <- function(datatype = NULL, plinkbasefile = NULL,
                         is.binary = FALSE, transpose = FALSE,
                         header = FALSE, path2Plink = "plink",
                         weight = "freq", gwaa.data = NULL,
                         pedigreefile = NULL){
    if (is.null(datatype)){
        stop("Data type was not specified. Please, choose 'pedigree' or 'genomic'")
    }

    if (datatype != "pedigree" & datatype != "genomic"){
        stop("Incorrect data type. Data type should be either 'pedigree' or 'genomic'")
    }

    if (datatype == "pedigree"){
        if (!is.binary){
            if (!transpose){
                message("Trying to read file in regular (untransposed) PLINK text format")

                if (!file.exists(paste(plinkbasefile, "ped", sep="."))){
                    stop(paste(plinkbasefile,
                               ".ped doesn't exist. Please check the path",
                               " to input file in PLINK format",
                               sep=""))
                }
                PedTmp <- read.table(
                    file=paste(plinkbasefile, "ped", sep="."),
                    colClasses=c(rep("character", 4),
                        "integer",
                        rep("NULL", count.fields(paste(plinkbasefile,
                                                       "ped",
                                                       sep="."))[1]-5)
                                 ),
                    header=header)
            }else{
                message("Trying to read file in PLINK text transposed format")

                if (!file.exists(paste(plinkbasefile, "tfam", sep="."))){
                    stop(paste(plinkbasefile,
                               ".tfam doesn't exist. Please check the path",
                               " to input file in PLINK format",
                               sep=""))
                }
                PedTmp <- read.table(
                    file=paste(plinkbasefile, "tfam", sep="."),
                    colClasses=c(rep("character", 4),
                        "integer",
                        rep("NULL", count.fields(paste(plinkbasefile,
                                                       "tfam",
                                                       sep="."))[1]-5)
                                 ),
                    header=header)
            }
        }else{
            message("Trying to read file in PLINK binary format")
            if (!file.exists(paste(plinkbasefile, "fam", sep="."))){
                stop(paste(plinkbasefile,
                           ".fam doesn't exist. Please check the path",
                           " to input file in PLINK format",
                           sep=""))
            }
            PedTmp <- read.table(
                file=paste(plinkbasefile, "fam", sep="."),
                colClasses=c(rep("character", 4),
                    "integer",
                    rep("NULL",
                        count.fields(paste(plinkbasefile,
                                           "fam",
                                           sep="."))[1]-5)
                             ),
                header=header)
        }
        message("Pedigree data was read successfully")
        message("Using kinship2 R package to estimate pedigree relationship")

        ManyFamilies <- length(unique(PedTmp[, 1])) != 1
        if (!ManyFamilies){
            tmpTped <- with(PedTmp, pedigree(id = PedTmp[, 2],
                                             mom = PedTmp[, 3],
                                             dad = PedTmp[, 4],
                                             sex = PedTmp[, 5]))
            return(kinship2::kinship(tmpTped) * 2)
        }
        if (ManyFamilies){
            return(as.matrix(kinship2::makekinship(famid=PedTmp[, 1],
                                                   id=PedTmp[, 2],
                                                   mother.id=PedTmp[, 3],
                                                   father.id=PedTmp[, 4])) * 2)
        }
    }
    if (datatype == "genomic"){
        if (!is.null(gwaa.data)){
            if (class(gwaa.data)[1] == "gwaa.data"){
                message("GenABEL gwas.data was found. Using GenABEL::ibs() to estimate genomic relationship matrix")

                if (weight == "no"){
                    message("Your weight option is 'no'. Estimated genomic relationship will be biased.")
                    message("Try to use weight='freq'")
                }
                gRelTmp = GenABEL::ibs(data=gwaa.data, weight = weight)
                gRelTmp[upper.tri(gRelTmp)] <- gRelTmp[lower.tri(gRelTmp)]
                return(gRelTmp * 2)
            }else{
                stop(paste(deparse(substitute(gwaa.data)),
                           "object is not a gwaa.data-class object. Abort",
                           sep=" "))
            }
        }
        if (!is.binary){
            if (!transpose){
                message("Trying to read file in regular (untransposed) PLINK text format")
                if (!file.exists(paste(plinkbasefile, "ped", sep="."))){
                    stop(paste(plinkbasefile,
                               ".ped doesn't exist. Please check the path",
                               " to input file in PLINK format",
                               sep=""))
                }
                system(paste(path2Plink,
                             "--noweb --file",
                             plinkbasefile,
                             "--make-rel square --out RelMat",
                             sep=" "))
                gRel <- read.table(file="RelMat.rel",
                                   stringsAsFactors=FALSE,
                                   header=FALSE)
                system("rm RelMat.log RelMat.rel RelMat.rel.id")
                return(gRel)

            }else{
                message("Trying to read file in PLINK text transposed format")
                if (!file.exists(paste(plinkbasefile, "tped", sep="."))){
                    stop(paste(plinkbasefile,
                               ".tped doesn't exist. Please check the path",
                               " to input file in PLINK format",
                               sep=""))
                }
                system(paste(path2Plink, "--noweb --tfile",
                             plinkbasefile,
                             "--make-rel square --out RelMat",
                             sep=" "))
                gRel <- read.table(file="RelMat.rel",
                                   stringsAsFactors=FALSE,
                                   header=FALSE)
                system("rm RelMat.log RelMat.rel RelMat.rel.id")
                return(gRel)
            }
        }else{
            message("Genotype data in PLINK binary format")
            message("Trying to read file in PLINK binary format")
            if (!file.exists(paste(plinkbasefile, "bed", sep="."))){
                stop(paste(plinkbasefile,
                           ".bed doesn't exist. Please check the path",
                           " to input file in PLINK format",
                           sep=""))
            }
            system(paste(path2Plink,
                         "--noweb --bfile",
                         plinkbasefile,
                         "--make-rel square --out RelMat",
                         sep=" "))
            gRel <- read.table(file="RelMat.rel",
                               stringsAsFactors=FALSE,
                               header=FALSE)
            system("rm RelMat.log RelMat.rel RelMat.rel.id")
            return(gRel)
        }
    }
}
