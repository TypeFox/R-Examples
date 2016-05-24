#' superposition of discharge, unit hydrograph
#' 
#' superposition of precipitation along unit hydrograph (to simulate Q from P)
#' 
#' @return list with optimized n and k, Nash-Sutcliffe Index, and simulated discharge
#' @author Berry Boessenkool, \email{berry-b@@gmx.de}, July 2013
#' @seealso \code{\link{lsc}} where superPos is used, \code{\link{unitHydrograph}}
#' @keywords hplot ts
#' @export
#' @examples 
#' 
#' N <- c(9,5,2,14,1,3) # [mm/hour]
#' UH <- c(0.1, 0.4, 0.3, 0.1, 0.1) # [1/h]
#' sum(UH) # sum must be 1
#' 
#' superPos(N, UH)
#' # If catchment area = 34 km^2 and precipitation is homogenous:
#' superPos(N/10^3, UH) * 34*10^6 / 3600  # m^3/s  # Add baseflow and you're done...
#' 
#' SP <- data.frame(Prec=c(N, 0,0,0,0),
#'           P1=c( UH*N[1], 0,0,0,0,0),
#'           P2=c(0, UH*N[2], 0,0,0,0),
#'           P3=c(0,0, UH*N[3], 0,0,0),
#'           P4=c(0,0,0, UH*N[4], 0,0),
#'           P5=c(0,0,0,0, UH*N[5], 0),
#'           P6=c(0,0,0,0,0, UH*N[6] ),
#'           runoff=superPos(N, UH))
#' SP # SuperPosition
#' 
#' SPcum <- t( apply(SP[2:7], 1, cumsum) )
#' 
#' plot(N, type="h", col=2:7, lwd=3, xlim=c(1, 10), ylim=c(30,0), lend=1)
#' par(new=TRUE)
#' plot(1, type="n", ylim=c(0, 15), xlim=c(1, 10), axes=FALSE, ann=FALSE)
#' axis(4, las=1)
#' polygon(x=c(1:10, 10:1), y=c(SPcum[,1], rep(0, 10)), col=2)
#' for(i in 2:6) polygon(x=c(1:10, 10:1), y=c(SPcum[,i], rev(SPcum[,i-1])), col=i+1)
#' text(2.5, 1, "Shape of UH")
#' 
#' lines( superPos(N, UH), lwd=3)
#' 
#' plot(UH, type="o", ylim=c(0, 0.4), las=1)
#' lines(UH, type="h" )
#' 
#' 
#' # Effect of distribution of Prec:
#' P_a <- c(1,2,3,4,5,6,7,8)
#' P_b <- c(4,4,4,4,4,4,4,4,4)
#' P_c <- c(8,7,6,5,4,3,2,1)
#' sum(P_a) ; sum(P_b) ; sum(P_c)
#' 
#' UH_1 <- unitHydrograph(n=2, k=2.3, t=1:25)
#' UH_2 <- unitHydrograph(n=5.5, k=1.8, t=1:25)
#' 
#' par(mfrow=c(2,3), mar=c(2,3,2,1), las=1)
#' plot(P_a, type="h", col=3, lwd=3, ylim=c(0,8), main="Precipitation a")
#' plot(P_b, type="h", col=4, lwd=3, ylim=c(0,8), main="Precipitation b")
#' plot(P_c, type="h", col=5, lwd=3, ylim=c(0,8), main="Precipitation c")
#' #
#' plot(UH_1, type="l", main="unit hydrograph", ylab="",xlab="Zeit")
#' lines(UH_2, col=2)
#' text(c(7,14), c(0.12, 0.07), c("UH_1","UH_2"), col=1:2)
#' abline(h=0)
#' #
#' plot( superPos(P=P_a, UH=UH_1), col=3, ylim=c(0,5), type="l",
#'       main="Discharge", ylab="Q [m^3/s]")
#' lines(superPos(P=P_b, UH=UH_1), col=4)
#' lines(superPos(P=P_c, UH=UH_1), col=5)
#' legend("topright", c("P a","P b", "P c"), title="with UH_1", col=3:5, lty=1)
#' #
#' plot( superPos(P=P_a, UH=UH_2), col=3, ylim=c(0,5), type="l", 
#'       main="Discharge", ylab="Q [m^3/s]")
#' lines(superPos(P=P_b, UH=UH_2), col=4)
#' lines(superPos(P=P_c, UH=UH_2), col=5)
#' legend("topright", c("P a","P b", "P c"), title="with UH_2", col=3:5, lty=1)
#' 
#' @param P Vector with precipitation values
#' @param UH Vector with discrete values of the Unit Hydrograph. This can be any UH summing to one, not just the storage cascade model.
#' 
superPos <- function(
P,
UH)
{
added <- length(UH)-1
qsim <- rep(0, length(P)+added )
for(i in 1:length(P) ) qsim[i:(i+added)] <- P[i]*UH + qsim[i:(i+added)]
qsim
}

