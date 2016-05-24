#' Calculate figure of merit
#' 
#' Caclulate the figure of merit (an objective measure of observer performance) for each treatment-reader combination. 
#' 
#' @param dataset The dataset to be analyzed, see \link{RJafroc-package}.
#' @param fom The figure of merit to be used in the calculation. The default is \code{"wJAFROC"}. See "Details".
#' 
#' @return An \code{c(I, J)} array, where the row names are the IDs of the treatments and column names are the IDs of the readers.
#' 
#' @details Allowed figures of merit are: (1) \code{"Wilcoxon"} for ROC data; (2) \code{"JAFROC1"},
#' \code{"JAFROC"}, \code{"wJAFROC1"}, \code{"wJAFROC"} (the default), \code{"HrAuc"}, \code{"SongA1"}, \code{"SongA2"}\bold{**} , 
#' \code{"HrSe"}, \code{"HrSp"}, \code{"MaxLLF"}, \code{"MaxNLF"}, \code{"MaxNLFAllCases"}, \code{"ExpTrnsfmSp"}, 
#' for free-response data and (3) \code{"ROI"} for ROI data. 
#' The JAFROC FOMs are described in the paper by Chakraborty and Berbaum. The Song FOMs are described in the paper by Song et al.
#' The \code{"MaxLLF"}, \code{"MaxNLF"} and \code{"MaxNLFAllCases"} FOMs correspond to ordinate, abscissa and abscissa, respectively, of the highest 
#' point on the FROC operating characteristic obtained by counting all the LL marks on diseased, all NL marks on non-diseased cases, 
#' and all NL marks on all cases, respectively). The \code{"ExpTrnsfmSp"} FOM is described in the paper by Popescu. 
#' The \code{"ROI"} FOM is described in the paper by Obuchowski et al.
#' 
#' ** \bold{The Song A2 figure of merit is computationally very intensive.}   
#'
#' @examples
#' FigureOfMerit(dataset = rocData, fom = "Wilcoxon")
#' 
#' FigureOfMerit(dataset = frocData)
#' 
#' @references
#' Chakraborty, D. P., & Berbaum, K. S. (2004). Observer studies involving detection and localization: modeling, analysis, and validation. 
#' Medical Physics, 31(8), 1-18.
#' 
#' Song T, Bandos AI, Rockette HE, Gur D (2008) On comparing methods for discriminating between actually negative and actually positive subjects 
#' with FROC type data. Medical Physics 35: 1547-1558.
#' 
#' Popescu, L. M. (2011). Nonparametric signal detectability evaluation using an exponential transformation of the FROC curve. 
#' Medical Physics, 38(10), 5690. 
#' 
#' Obuchowski, N. A., Lieber, M. L., & Powell, K. A. (2000). Data Analysis for Detection and Localization of Multiple Abnormalities 
#' with Application to Mammography. Academic Radiology, 553-554.
#'  
#' @export
#'  
FigureOfMerit <- function(dataset, fom = "wJAFROC") {
  NL <- dataset$NL
  LL <- dataset$LL
  lesionNum <- dataset$lesionNum
  lesionID <- dataset$lesionID
  lesionWeight <- dataset$lesionWeight
  maxNL <- dim(NL)[4]
  dataType <- dataset$dataType
  modalityID <- dataset$modalityID
  readerID <- dataset$readerID
  I <- length(modalityID)
  J <- length(readerID)
  K <- dim(NL)[3]
  K2 <- dim(LL)[3]
  K1 <- K - K2  
  
  if (dataType != "ROC" && fom == "Wilcoxon") 
    stop("Cannot use Wilcoxon with FROC or ROI data")
  
  if (dataType == "ROI" && fom != "ROI") {
    errMsg <- paste0("ROI dataset cannot be analyzed using ", fom, " figure of merit.")
    stop(errMsg)
  }
  
  if (K1 == 0 && (fom != "JAFROC1" && fom != "wJAFROC1")) {
    errMsg <- paste0("Only JAFROC1 or wJAFROC1 figures of merit are allowed for datasets with no non-diseased cases.")
    stop(errMsg)
  }
  
  fomArray <- array(dim = c(I, J))
  for (i in 1:I) {
    for (j in 1:J) {
      nl <- NL[i, j, , ]
      ll <- LL[i, j, , ]
      dim(nl) <- c(K, maxNL)
      dim(ll) <- c(K2, max(lesionNum))
      fomArray[i, j] <- MyFOM(nl, ll, lesionNum, lesionID, lesionWeight, maxNL, fom)
    }
  }
  rownames(fomArray) <- paste("Trt", modalityID, sep = " - ")
  colnames(fomArray) <- paste("Rdr", readerID, sep = " - ")
  return(fomArray)
} 
