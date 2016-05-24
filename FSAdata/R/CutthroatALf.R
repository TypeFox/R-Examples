#' @title Capture histories (9 samples) of Cutthroat Trout from Auke Lake.
#'
#' @description Summarized (\dQuote{RMark} format) capture histories of Cutthroat Trout (\emph{Oncorhynchus clarki}) in Auke Lake, Alaska, from samples taken in 1998-2006.
#'
#' @name CutthroatALf
#' 
#' @docType data
#' 
#' @format A data frame with 47 observations on the following 2 variables.
#'  \describe{
#'    \item{ch}{Unique capture history (as a character string)}
#'    \item{freq}{Frequency of fish with that capture history}
#'  }
#'
#' @section Topic(s):
#'  \itemize{
#'    \item Population Size
#'    \item Abundance
#'    \item Mark-Recapture
#'    \item Capture-Recapture
#'    \item Jolly-Seber
#'    \item Capture History 
#'  }
#' 
#' @concept Abundance 'Population Size' 'Mark-Recapture' 'Capture-Recapture' 'Jolly-Seber' 'Capture History'
#' 
#' @source Entered from Appendix A.3 of Harding, R.D., C.L. Hoover, and R.P. Marshall. 2010.  Abundance of Cutthroat Trout in Auke Lake, Southeast Alaska, in 2005 and 2006.  Alaska Department of Fish and Game Fisheries Data Series No. 10-82.  Accessed from http://www.sf.adfg.state.ak.us/FedAidPDFs/FDS10-82.pdf.
#' 
#' @seealso See \code{CutthroatAL} for the same data in \dQuote{individual} fish format (i.e., the data in this file were converted using \code{capHistConvert} from \pkg{FSA}).  See \code{mrOpen} from \pkg{FSA} for an example analysis.
#' 
#' @keywords datasets
#' 
#' @examples
#' data(CutthroatALf)
#' str(CutthroatALf)
#' head(CutthroatALf)
#'
NULL
