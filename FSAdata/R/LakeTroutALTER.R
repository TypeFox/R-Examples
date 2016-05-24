#' @title Biological data for Lake Trout from the Arctic LTER (AK).
#' 
#' @description Biological data (lengths, weight, age, and sex) of Lake Trout (\emph{Salvelinus namaycush}) sampled from Lake NE12 of the Arctic Long Term Ecological Research location.
#' 
#' @details Lake trout were removed from Lake NE12 in the summers of 1986, 1988, and 1989 using five-panel experimental gill nets (mesh size of 0.75, 1, 1.5, 2, and 2.5 inches).  Lengths, weights, and sex were recorded from the fish while otoliths, and if possible, stomachs and gonads were removed for future analysis.  A check was performed on several otoliths by an independent colleague and prevents introduction of bias due to familiarity with the samples.  The original file was \dQuote{cleaned} in the following ways:
#' \enumerate{
#'   \item Only Lake Trout were kept in the data file.
#'   \item All unknown sex fish were removed.
#'   \item Fish with missing data (length, weight,age, or sex) were removed.
#'   \item Decimals were removed from the ages.
#'   \item The unique IDs for fish from 1989 were changed to start at 500.
#'   \item The weight of fish number 509 was changed from 100 to 1100.
#' }
#' 
#' @name LakeTroutALTER
#' 
#' @docType data
#' 
#' @format A data frame of 86 observations on the following 6 variables:
#'  \describe{
#'    \item{id}{A unique identification number.} 
#'    \item{tl}{Total Length (nearest mm) at capture.} 
#'    \item{fl}{Fork Length (nearest mm) at capture.}
#'    \item{sl}{Standard Length (nearest mm) at capture.} 
#'    \item{w}{Weight (nearest g) at capture.} 
#'    \item{otorad}{Total otolith radius (mm) at capture.}
#'    \item{age}{Age (completed growing seasons) at capture.} 
#'    \item{sex}{Sex of the fish (\code{F}=female and \code{M}=male).} 
#'  }
#' 
#' @section Topic(s):
#'  \itemize{
#'    \item Length Frequency 
#'    \item Weight-Length 
#'    \item Length Conversion 
#'    \item Growth
#'    \item von Bertalanffy 
#'    \item Size Structure
#'  }
#'  
#' @concept 'Length Frequency' 'Weight-Length' 'Growth' 'von Bertalanffy' 'Size Structure' 'Length Conversion'
#' 
#' @source Was (does not appear to be available there now) from http://ecosystems.mbl.edu/ARC/lakes/fish/89mcne12.html.  It seems like it should still be available from the Arctic LTER site at http://ecosystems.mbl.edu/ARC/lakes/fish/index.shtml.
#' 
#' @keywords datasets
#' 
#' @examples
#' data(LakeTroutALTER)
#' str(LakeTroutALTER)
#' head(LakeTroutALTER)
#' op <- par(mfrow=c(2,2),pch=19)
#' ## Four (of many possible) examples
#' hist(LakeTroutALTER$tl,main="")
#' plot(w~tl,data=LakeTroutALTER)
#' plot(tl~fl,data=LakeTroutALTER)
#' plot(tl~age,data=LakeTroutALTER)
#' par(op)
#' 
NULL