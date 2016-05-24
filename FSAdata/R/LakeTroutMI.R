#' @title Stock and recruitment data for Lake Trout in Lake Superior, 1971-1991.
#' 
#' @description Stock and recruitment data for the 1971-1991 year-classes of Lake Trout (\emph{Salvelinus namaycush}) in Michigan waters of Lake Superior.
#' 
#' @name LakeTroutMI
#' 
#' @docType data
#' 
#' @format A data frame of 105 observations on the following 5 variables:
#'  \describe{
#'    \item{year}{Year of data}
#'    \item{recruits}{Recuit index -- geometric mean number of age-7 fish/km/net-night}
#'    \item{wild}{Wild fish spawning stock index -- geometric mean number of wild age-8 and older fish/km/net-night}
#'    \item{stocked}{Stocked fish spawning stock index -- geometric mean number of stocked age-8 and older fish/km/net-night}
#'    \item{area}{Lake Superior management unit} 
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
#' @source From (approximately) figures in Richards, J.M., M.J. Hansen, C.R. Bronte, and S.P. Sitar.  2004.  Recruitment dynamics of the 1971-1991 year-classes of Lake Trout in Michigan waters of Lake Superior.  North American Journal of Fisheries Management.  24:475-489.
#' 
#' @keywords datasets
#' 
#' @examples
#' data(LakeTroutMI)
#' LakeTroutMI$stock <- LakeTroutMI$wild+LakeTroutMI$stocked
#' str(LakeTroutMI)
#' head(LakeTroutMI)
#' op <- par(mfrow=c(3,2),pch=19)
#' plot(recruits~year,data=LakeTroutMI,subset=area=="MI3",type="l",ylim=c(0,max(recruits,na.rm=TRUE)))
#' lines(recruits~year,data=LakeTroutMI,subset=area=="MI4",col="blue")
#' lines(recruits~year,data=LakeTroutMI,subset=area=="MI5",col="green")
#' lines(recruits~year,data=LakeTroutMI,subset=area=="MI6",col="red")
#' lines(recruits~year,data=LakeTroutMI,subset=area=="MI7",col="yellow")
#' plot(recruits~stock,data=LakeTroutMI,subset=area=="MI3",main="MI3")
#' plot(recruits~stock,data=LakeTroutMI,subset=area=="MI4",col="blue",main="MI4")
#' plot(recruits~stock,data=LakeTroutMI,subset=area=="MI5",col="green",main="MI5")
#' plot(recruits~stock,data=LakeTroutMI,subset=area=="MI6",col="red",main="MI6")
#' plot(recruits~stock,data=LakeTroutMI,subset=area=="MI7",col="yellow",main="MI7")
#' par(op)
#' 
NULL