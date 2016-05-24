#' @title Biological data for Slimy Sculpin from the Arctic LTER (AK).
#' 
#' @description Biological data (length, weight, and age) for Slimy Sculpin (\emph{Cottus cognatus}) sampled from Lake S-6 of the Arctic Long Term Ecological Research location.
#' 
#' @details Slimy Sculpin were captured from Lake S-6 during the summers of 1988 and 1989.  Sculpins were captured using canning jars with acetate funnels inserted in them (similar in design to some minnow traps). The sculpin traps were set by raft and placed on the bottom of the lake in a specific habitat (inshore, rocky, rock-mud interface, and mud bottom). (Roughly) every three days the traps were retrieved and the sculpins, if any, were removed.  The ages of the fish were determined by counting the rings in their otoliths.
#' 
#' @name SculpinALTER
#' 
#' @docType data
#' 
#' @format A data frame of 117 observations on the following 3 variables:
#'  \describe{
#'    \item{w}{Weight (nearest 0.1 g) at capture} 
#'    \item{tl}{Total length (nearest mm) at capture}
#'    \item{age}{Age (completed growing seasons from otoliths) at capture}
#'  }
#'  
#' @section Topic(s):
#'  \itemize{
#'    \item Length Frequency 
#'    \item Weight-Length
#'    \item Growth
#'    \item von Bertalanffy
#'    \item Size Structure
#'  }
#'  
#' @concept Growth 'von Bertalanffy' 'Length Frequency' 'Weight-Length' 'Size Structure'
#' 
#' @source Was (does not appear to be available there now) from http://ecosystems.mbl.edu/ARC/lakes/fish/89hes6ag.html.  It seems like it should still be available from the Arctic LTER site at http://ecosystems.mbl.edu/ARC/lakes/fish/index.shtml.
#' 
#' @keywords datasets
#' 
#' @examples
#' data(SculpinALTER)
#' str(SculpinALTER)
#' head(SculpinALTER)
#' op <- par(mfrow=c(2,2),pch=19)
#' hist(SculpinALTER$tl,main="")
#' hist(SculpinALTER$age,main="")
#' plot(w~tl,data=SculpinALTER)
#' plot(tl~age,data=SculpinALTER)
#' par(op)
#' 
NULL