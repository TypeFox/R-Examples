#' @title Catch-at-age for Yellow Perch from Chequamegon Bay, Lake Superior.
#' 
#' @description Catch-at-age for Yellow Perch (\emph{Perca flavescens}) from Chequamegon Bay, Lake Superior, 1973-1988.
#' 
#' @name YPerchCB1
#' 
#' @docType data
#' 
#' @format A data frame with 39 observations on the following 5 variables.
#'  \describe{
#'    \item{age}{Age in that capture year (1973-1988).}
#'    \item{year73}{Number of fish at each age in capture year 1973.}
#'    \item{year74}{Number of fish at each age in capture year 1974.}
#'    \item{year75}{Number of fish at each age in capture year 1975.}
#'    \item{year76}{Number of fish at each age in capture year 1976.}
#'    \item{year77}{Number of fish at each age in capture year 1977.}
#'    \item{year78}{Number of fish at each age in capture year 1978.}
#'    \item{year79}{Number of fish at each age in capture year 1979.}
#'    \item{year80}{Number of fish at each age in capture year 1980.}
#'    \item{year81}{Number of fish at each age in capture year 1981.}
#'    \item{year82}{Number of fish at each age in capture year 1982.}
#'    \item{year83}{Number of fish at each age in capture year 1983.}
#'    \item{year84}{Number of fish at each age in capture year 1984.}
#'    \item{year85}{Number of fish at each age in capture year 1985.}
#'    \item{year86}{Number of fish at each age in capture year 1986.}
#'    \item{year87}{Number of fish at each age in capture year 1987.}
#'    \item{year88}{Number of fish at each age in capture year 1988.} 
#'  }
#'  
#' @section Topic(s):
#'  \itemize{
#'    \item Mortality
#'    \item Catch curve 
#'  }
#'  
#' @concept Mortality 'Catch Curve'
#' 
#' @source From Table 1 in Bronte, C.R., J.H. Selgeby, D.V.Swedberg.  1993.  Dynamics of a yellow perch population in Western Lake Superior.  North American Journal of Fisheries Management 13:511-523.
#' 
#' @keywords datasets
#' 
#' @examples
#' data(YPerchCB1)
#' str(YPerchCB1)
#' head(YPerchCB1)
#' max.n <- max(YPerchCB1[,-1])  # maximum catch
#' op <- par(mfrow=c(4,4),mar=c(3.5,3.5,1,1),mgp=c(2,0.75,0))
#' plot(log(year73)~age,data=YPerchCB1,ylim=c(0,log(max.n)),main="1973")
#' plot(log(year74)~age,data=YPerchCB1,ylim=c(0,log(max.n)),main="1974")
#' plot(log(year75)~age,data=YPerchCB1,ylim=c(0,log(max.n)),main="1975")
#' plot(log(year76)~age,data=YPerchCB1,ylim=c(0,log(max.n)),main="1976")
#' plot(log(year77)~age,data=YPerchCB1,ylim=c(0,log(max.n)),main="1977")
#' plot(log(year78)~age,data=YPerchCB1,ylim=c(0,log(max.n)),main="1978")
#' plot(log(year79)~age,data=YPerchCB1,ylim=c(0,log(max.n)),main="1979")
#' plot(log(year80)~age,data=YPerchCB1,ylim=c(0,log(max.n)),main="1980")
#' plot(log(year81)~age,data=YPerchCB1,ylim=c(0,log(max.n)),main="1981")
#' plot(log(year82)~age,data=YPerchCB1,ylim=c(0,log(max.n)),main="1982")
#' plot(log(year83)~age,data=YPerchCB1,ylim=c(0,log(max.n)),main="1983")
#' plot(log(year84)~age,data=YPerchCB1,ylim=c(0,log(max.n)),main="1984")
#' plot(log(year85)~age,data=YPerchCB1,ylim=c(0,log(max.n)),main="1985")
#' plot(log(year86)~age,data=YPerchCB1,ylim=c(0,log(max.n)),main="1986")
#' plot(log(year87)~age,data=YPerchCB1,ylim=c(0,log(max.n)),main="1987")
#' plot(log(year88)~age,data=YPerchCB1,ylim=c(0,log(max.n)),main="1988")
#' par(op)
#' op <- par(mfrow=c(1,2),mar=c(3.5,3.5,1,1),mgp=c(2,0.75,0))
#' 
#' ## plot for 1973 and 1982 year-classes (be very careful)
#' # get 1973 year-class as diagonal for ages 0-9 and years 1973-1982
#' yc73 <- diag(as.matrix(YPerchCB1[,2:11]))
#' plot(log(yc73)~YPerchCB1$age,main="1973 Year-Class")
#' yc82 <- diag(as.matrix(YPerchCB1[1:7,11:17]))
#' plot(log(yc82)~YPerchCB1$age[1:7],main="1982 Year-Class")
#' 
NULL