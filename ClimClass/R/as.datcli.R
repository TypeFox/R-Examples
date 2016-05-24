NULL
#' as.datcli
#' 
#' @description Tranforms a data frame (see example dataset) into a data frame format like  '\code{datcli}' in '\code{climatol}' package
#' 
#' @param df data frame or list 
#' @param station name
#' @param MonthField character string for month field in \code{df}. Default is \code{"month"}.
#' @param PrecField character string for Mean Precipitation field in \code{df}. Default is \code{"P"}.
#' @param MinTempField character string  for Mean Daily Minimum Temperature field in \code{df}. Default is \code{"Tn"}.
#' @param MaxTempField character string  for Mean Daily Maximum Temperature field in \code{df}. Default is \code{"Tx"}.
#' @param MeanTempField character string for Mean Daily Maximum Temperature field in \code{df}. Default is \code{"Tm"}.
#' @param AbsMinTempField character string for Absolute Monthly Minimum Temperature field in \code{df}. Default is \code{"AbsTn"}. 
#' @param AbsMinTempOffset estimated offset between Average Min Temperature and  Absolute Min Temperature.
#' @param StationField  character string for Station  field in \code{df}. Default is \code{"station"}.
#'    
#' @author Emanuele Cordano
#' 
#' @export
#' @seealso \url{www.climatol.eu}, \url{http://www.zoolex.org/walter.html}
#' 
#' @examples 
#' 
#' ### Not Run!! 
#' # Install 'climatol' from 'http://www.climatol.eu/' first
#' ### Then load the package, uncomment and run the following line
#' # library(climatol)
#' 	library(stringr)
#'  data(Trent_climate)
#' 
#'  TrentinoClimateDf <- do.call(rbind,clima_81_10)
#'  names <- rownames(TrentinoClimateDf)
#'  TrentinoClimateDf$station <- 
#'  unlist(lapply(X=str_split(names,pattern="[.]"),FUN=function(x) {x[1]}))
#'  
#' 
#'  station <- "T0129"
#' datcli <- as.datcli(TrentinoClimateDf,station=station)
#' 
#' ### Not Run!! 
#' # Install 'climatol' from 'http://www.climatol.eu/' first
#' ### Then load the package, uncomment and run the following line
#' # diagwl(datcli,est=station,alt=100,per="Period",mlab="en") ## plots a Walter-Lieth's climograph
#' 


as.datcli <- function(df,station,MonthField="month",
		     PrecField="P",
			 MinTempField="Tn",
			 MaxTempField="Tx",
			 MeanTempField="Tm",
			 AbsMinTempField="AbsTn",
			 AbsMinTempOffset=4,
			 StationField="station") {
	
		 	 
		 if (length(station)>1)	 {
			 
			 station <- station[1]
			 warning("Only first station is considered!")
		 }
		 
		 if (is.data.frame(df)) {
		 	out <- df[df[,StationField]==station,]
		 } else if (is.list(df)) {
			 
			 out <- df[[station]]
		 }
		 if  (!(AbsMinTempField %in% names(out))) {
			 
			 warning("Missing absolute Monthly Minimum Temperature: calculated by AbsMinTempOffset!")
			 out[,AbsMinTempField] <- out[,MinTempField]-AbsMinTempOffset
		 }
		 
		 row.names(out) <- sprintf("2000-%02d-01",out[,MonthField])
		 
		 
		 out <- out[,c(PrecField,MaxTempField,MinTempField,AbsMinTempField)]
		 
		 sort <- sort(row.names(out))
		 
		 
		 out <- out[sort,]
		 
		 row.names(out) <- months(as.Date(row.names(out)),abbreviate=TRUE)
	
		 out <- t(as.matrix(out))
		 
	     return(out)	
	
	
	
}
