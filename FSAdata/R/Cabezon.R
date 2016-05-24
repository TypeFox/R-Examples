#' @title Ages, lengths, and maturity for female Cabezon from Oregon.
#' 
#' @description Ages, lengths, and maturity for female Cabezon (\emph{Scorpaenichthys marmoratus}) from Newport and Depoe Bay, Oregon.
#' 
#' @name Cabezon
#' 
#' @docType data
#' 
#' @format A data frame with 525 observations on the following 5 variables.
#'  \describe{
#'    \item{date}{Date fish was collected} 
#'    \item{length}{Total length (cm)} 
#'    \item{age}{Otolith age} 
#'    \item{maturity}{Maturity state (Immature or Mature)} 
#'    \item{stage}{Stage of maturity (1:Immature, 2:Maturing, 3:Mature, 4:Fertilized, 5:Ripe, 6:Spent, 7:Resting}
#'  }
#' 
#' @section Topic(s):
#'  \itemize{
#'    \item Maturity 
#'    \item Growth 
#'    \item von Bertalanffy 
#'  }
#' 
#' @concept Growth 'von Bertalanffy' Maturity
#' 
#' @source Actual data obtained directly (from Bob Hanna) from Hannah, R.W, M.T.O. Blume, and J.E. Thompson.  2009.  Length and age at maturity of female yelloweye rockfish (\emph{Sebastes rubberimus}) and cabezon (\emph{Scorpaenichthys marmoratus}) from Oregon waters based on histological evaluation of maturity.  Oregon Department of Fish and Wildlife, Information Reports 2009-04.  [Was (is?) from http://www.dfw.state.or.us/mrp/publications/docs/Info200904_YlwEyeRF_Maturity.pdf.]
#' 
#' @keywords datasets
#' 
#' @examples
#' data(Cabezon)
#' str(Cabezon)
#' head(Cabezon)
#' op <- par(mfrow=c(2,2),pch=19)
#' plot(length~age,data=Cabezon,ylab="Total Length (cm)",xlab="Age")
#' hist(Cabezon$length,xlab="Total Length (cm)",main="")
#' tbl1 <- xtabs(~age+maturity,data=Cabezon)
#' (ptbl1 <- prop.table(tbl1,margin=1))
#' plot(ptbl1[,2]~as.numeric(row.names(ptbl1)),type="l",xlab="Age",ylab="Proportion Mature")
#' tbl2 <- xtabs(~length+maturity,data=Cabezon)
#' (ptbl2 <- prop.table(tbl2,margin=1))
#' plot(ptbl2[,2]~as.numeric(row.names(ptbl2)),type="l",xlab="Length",ylab="Proportion Mature")
#' par(op)
#' 
NULL
