#' Obtain VT table from Track
#' 
#' The VT table containes speeds, steplengths, orientations and other summaries derived from a track. The output of this function is (typically) meant to feed the \code{\link{WindowSweep}} function. 
#' 
#' @param Data a track to analyze.  Must contain columns: X, Y and Time (as a POSIX object).  The \code{track} class is a robust entry.  
#' @param units the time units for the analysis; one of \code{sec, min, hour, day}.
#' @param skiplast filters away last step.
#' 
#' @return a data frame with the following columns:
#' \item{Z.start, Z.end}{the start and end locations (as complex coordinates)}
#' \item{S}{the step length}
#' \item{Phi, Theta}{absolute and turning angle, respectively}
#' \item{T.start, T.end}{start and time of steps (numeric - in given units)}
#' \item{T.mid}{temporal midpoint of the step }
#' \item{dT}{duration of the step}
#' \item{V}{approximate speed (S/dT)}
#' \item{T.POSIX}{the temporal midpoint of the step as a POSIX objects.}
#' 
#' @author Eliezer Gurarie
#' @seealso \code{\link{WindowSweep}}
#' @examples
#' data(Simp)
#' plot(Simp)
#' Simp.VT <- GetVT(Simp)
#' head(Simp.VT)
#' # Distribution of estimated speeds
#' hist(Simp.VT$V, col="grey", breaks=20)
#' # Distribution of turning angles
#' require(circular)
#' rose.diag(Simp.VT$Theta, bins=24)

GetVT <- function(Data, units = "hour", skiplast=TRUE)
{
  if (!"Z" %in% names(Data)) 
    Data$Z <- Data$X + (0+1i) * Data$Y
  Z.start <- Data$Z[-nrow(Data)]
  Z.end <- Data$Z[-1]
  S <- Mod(diff(Data$Z))
  Phi <- Arg(diff(Data$Z))
  Theta <- c(NA, diff(Phi))
  
  
  T.POSIX <- Data$Time[-nrow(Data)] + diff(Data$Time)/2
  
  if(inherits(Data$Time, "POSIXt"))
  {  
    Data$Time <- as.numeric(Data$Time-Data$Time[1])
    Data$Time <- Data$Time/ifelse(units == "sec", 1, 
                                  ifelse(units == "min", 60, 
                                         ifelse(units == "hour", 60*60, 
                                                ifelse(units == "day", 60*60*24, 
                                                       stop("Invalid time unit.")))))
  }
  
  T.start <- Data$Time[-nrow(Data)]
  T.end <- Data$Time[-1]
  dT <- T.end-T.start
  V <- S/as.vector(dT)
  
  T.mid <- (T.start + T.end)/2
  VT.table <- data.frame(Z.start, Z.end, S, Phi, Theta, T.start, T.end, T.mid, dT, V, T.POSIX)
  if(skiplast) VT.table <- VT.table[-1,]
  return(VT.table)
}
