#' @title Stock and recruitment data for Alaskan Pink Salmon, 1960-1990.
#' 
#' @description Harvest, escapement, and return of northern Southeast Alaska Pink Salmon (\emph{Oncorhynchus gorbuscha}), 1960-1991, and average sea surface temperature between June and November off Sitka, Alaska, 1960-1990.
#' 
#' @name PSalmonAK
#' 
#' @docType data
#' 
#' @format A data frame of 34 rows on the following 5 variables:
#'  \describe{
#'    \item{year}{Year of data}
#'    \item{harvest}{Harvest (thousands of fish)}
#'    \item{escapement}{Escapement (thousands of fish)}
#'    \item{return}{Returns (thousands of fish) as sum of harvest and escapement from two years later (lagging is for proper brood year correspondence)}
#'    \item{SST}{Average sea surface temperature (C) between June and November off Sitka, AK from one year latter (lagging is for matching when the salmon are actually in the ocean)}
#'  }
#'  
#' @section Topic(s):
#'  \itemize{
#'    \item Stock-Recruit
#'    \item Recruitment
#'  }
#' 
#' @concept 'Stock-Recruit' Recruitment
#' 
#' @source From Table 3.1 in Quinn, T.J. and R.B. Deriso.  1999.  Quantitative Fish Dynamics.  Oxford University Press, New York, New York.
#' 
#' @keywords datasets
#' 
#' @examples
#' data(PSalmonAK)
#' str(PSalmonAK)
#' head(PSalmonAK)
#' op <- par(mfrow=c(1,2))
#' plot(return~year,data=PSalmonAK)
#' plot(return~escapement,data=PSalmonAK)
#' par(op)
#' 
NULL