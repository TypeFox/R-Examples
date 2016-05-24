#' @title Ages, lengths, and maturity for Yelloweye Rockfish.
#' 
#' @description Ages, lengths, and maturity for female Yelloweye Rockfish (\emph{Sebastes rubberimus}) from Oregon.
#' 
#' @name YERockfish
#' 
#' @docType data
#' 
#' @format A data frame with 159 observations on the following 5 variables.
#'  \describe{
#'    \item{date}{Date fish was collected} 
#'    \item{length}{Total length (cm)} 
#'    \item{age}{Otolith age} 
#'    \item{maturity}{Maturity state (\code{Immature} or \code{Mature})} 
#'    \item{stage}{Stage of maturity (\code{1}:Immature, \code{2}:Maturing, \code{3}:Mature, \code{4}:Fertilized, \code{5}:Ripe, \code{6}:Spent, \code{7}:Resting} 
#'  }
#'  
#' @section Topic(s):
#'  \itemize{
#'    \item Growth
#'    \item Maturity
#'    \item von Bertalanffy 
#'  }
#'  
#' @concept Growth 'von Bertalanffy' 'Maturity'
#' 
#' @source Obtained directly (from Bob Hannah).  Date were used in Hannah, R.W, M.T.O. Blume, and J.E. Thompson.  2009.  Length and age at maturity of female yelloweye rockfish (\emph{Sebastes rubberimus}) and cabezon (\emph{Scorpaenichthys marmoratus}) from Oregon waters based on histological evaluation of maturity.  Oregon Department of Fish and Wildlife, Information Reports 2009-04.  [Was (is?) from http://www.dfw.state.or.us/mrp/publications/docs/Info200904_YlwEyeRF_Maturity.pdf]
#' 
#' @keywords datasets
#' 
#' @examples
#' data(YERockfish)
#' str(YERockfish)
#' head(YERockfish)
#' op <- par(mfrow=c(2,2),pch=19)
#' plot(length~age,data=YERockfish,ylab="Total Length (cm)",xlab="Age")
#' hist(YERockfish$length,xlab="Total Length (cm)",main="")
#' tbl1 <- with(YERockfish,table(age,maturity))
#' (ptbl1 <- prop.table(tbl1,margin=1))
#' plot(ptbl1[,2]~as.numeric(row.names(ptbl1)),type="l",xlab="Age",ylab="Proportion Mature")
#' tbl2 <- with(YERockfish,table(length,maturity))
#' (ptbl2 <- prop.table(tbl2,margin=1))
#' plot(ptbl2[,2]~as.numeric(row.names(ptbl2)),type="l",xlab="Age",ylab="Proportion Mature")
#' par(op)
#' 
NULL