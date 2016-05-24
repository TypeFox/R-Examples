#' @title Compute METE species area relationship (SAR)
#'
#' @description Uses raw data or state variables to calculate METE SAR 
#' and EAR (endemics area relatiohsip) as well as compute the observed 
#' SAR or EAR from data, if provided
#'
#' @details Currently only doublings of area are supported. Predictions 
#' and comparison to data can be made via several options. If \code{spp} 
#' and \code{abund} are not provided then only theoretical predictions 
#' are returned without emperical SAR or EAR results. In this case areas 
#' can either be specified by providing \code{Amin} and \code{A0} from 
#' which a vector of doubling areas is computed, or my providing \code{row}, 
#' \code{col} and \code{A0} in which case \code{row} and \code{col} are 
#' taken to be the number of desired rows and columns used to construct 
#' a grid across the landscape. If data are provided in the form of 
#' \code{spp} and \code{abund} then either \code{row} and \code{col} or 
#' \code{x} and \code{y} must be provided for each data entry (i.e. the 
#' length of \code{row} and \code{col} or \code{x} and \code{y} must equal
#' the length of \code{spp} and \code{abund}). If \code{x} and \code{y}
#' are provided then the landscape is gridded either by specifying 
#' \code{Amin} (the size of the smallest grid cell) or by providing the
#' number or desired rows and columns via the \code{row} and \code{col}
#' arguments.
#' 
#' 
#' @param spp vector of species identities
#' @param abund numberic vector abundances associated with each record
#' @param row identity of row in a gridded landscape associated with each record, or desired number of rows to divide the landcape into
#' @param col identity of column in a gridded landscape associated with each recod, or desired number of columns to divide the landcape into
#' @param x the x-coordinate of an individual if recorded
#' @param y the y-coordinate of an individual if recorded
#' @param S0 total number of species
#' @param N0 total abundance
#' @param Amin the smallest area, either the anchor area for upscaling or the desired area to downscale to
#' @param A0 the largest area, either the area to upscale to or the total area from which to downscale
#' @param upscale logical, should upscaling or downscaling be carried out
#' @param EAR logical, should the EAR or SAR be computed
#' 
#' @export
#' 
#' @examples
#' \dontrun{
#' data(anbo)
#' 
#' ## using row and col from anbo dataset
#' anbo.sar1 <- meteSAR(anbo$spp, anbo$count, anbo$row, anbo$col, Amin=1, A0=16)
#' plot(anbo.sar1)
#' 
#' ## using simulated x, y data
#' anbo.sar2 <- meteSAR(anbo$spp, anbo$count, x=anbo$x, y=anbo$y, row=4, col=4)
#' plot(anbo.sar2)
#' 
#' ## using just state variable
#' thr.sar <- meteSAR(Amin=1, A0=16, S0=50, N0=500) 
#' }
#' @return an object of class \code{meteRelat} with elements
#' \describe{
#'   \item{\code{pred}}{predicted relationship; an object of class \code{sar}}
#'   \item{\code{obs}}{observed relationship; an object of class\code{sar}}
#' }
#'
#' @author Andy Rominger <ajrominger@@gmail.com>, Cory Merow
#' @seealso sad, meteESF, metePi
#' @references Harte, J. 2011. Maximum entropy and ecology: a theory of abundance, distribution, and energetics. Oxford University Press.
# @aliases - a list of additional topic names that will be mapped to
# this documentation when the user looks them up from the command
# line.
# @family sar

meteSAR <- function(spp, abund, row, col, x, y, S0 = NULL, N0 = NULL,
                    Amin, A0, upscale=FALSE, EAR=FALSE) {    
  ## figure out vector of sizes in units of cells; right now only doublings supported
  ## not needed if upscale is TRUE
  if(!upscale) {
    areaInfo <- .findAreas(
      spp=if(missing(spp)) NULL else spp,
      abund=if(missing(abund)) NULL else abund, 
      row=if(missing(row)) NULL else row, 
      col=if(missing(col)) NULL else col, 
      x=if(missing(x)) NULL else x, 
      y=if(missing(y)) NULL else y, 
      Amin=if(missing(Amin)) NULL else Amin, 
      A0=if(missing(A0)) NULL else A0)
    areas <- areaInfo$areas
    row <- areaInfo$row
    col <- areaInfo$col
    nrow <- areaInfo$nrow
    ncol <- areaInfo$ncol
    Amin <- areaInfo$Amin
    A0 <- areaInfo$A0
  }
  
  if(upscale & EAR) stop('upscaling EAR not currently supported')
  
  ## the ESF
  if(!missing(spp) & !missing(abund)) {
    S0 <- length(unique(spp))
    N0 <- sum(abund)
  }
  if(is.null(S0) | is.null(N0)) stop('must provide spp and abund data or state variables S0 and N0')
  thisESF <- meteESF(S0=S0, N0=N0)
  
  ## calculate empirical SAR
  if(!missing(spp) & !missing(abund)) {
    eSAR <- empiricalSAR(spp, abund, row=row, col=col, Amin=Amin, A0=A0, EAR=EAR)
  } else {
    eSAR <- NULL
  }
  
  ## calculate theoretical SAR
  if(upscale) {
    thrSAR <- upscaleSAR(thisESF, Amin, A0, EAR)
  } else {
    thrSAR <- downscaleSAR(thisESF, areas*Amin, A0, EAR)
  }
  
  out <- list(obs=eSAR, pred=thrSAR)
  class(out) <- 'meteRelat'
  
  return(out)
}



#================================================================
#' @title Empirical SAR or EAR
#'
#' @description computes observed SAR or EAR from raw data
#'
#' @details Currently only doublings of area are supported. There are 
#' several options for specifying areas. Either \code{row} and \code{col} or 
#' \code{x} and \code{y} must be provided for each data entry (i.e. the 
#' length of \code{row} and \code{col} or \code{x} and \code{y} must equal
#' the length of \code{spp} and \code{abund}). If \code{x} and \code{y}
#' are provided then the landscape is gridded either by specifying 
#' \code{Amin} (the size of the smallest grid cell) or by providing the
#' number or desired rows and columns via the \code{row} and \code{col}
#' arguments. If only \code{row} and \code{col} are provided these are taken
#' to be the row and column identities of each data entry
#' 
#' 
#' 
#' @param spp vector of species identities
#' @param abund numberic vector abundances associated with each record
#' @param row identity of row in a gridded landscape associated with each record, or desired number of rows to divide the landcape into
#' @param col identity of column in a gridded landscape associated with each recod, or desired number of columns to divide the landcape into
#' @param x the x-coordinate of an individual if recorded
#' @param y the y-coordinate of an individual if recorded
#' @param Amin the smallest area, either the anchor area for upscaling or the desired area to downscale to
#' @param A0 the largest area, either the area to upscale to or the total area from which to downscale
#' @param EAR logical, should the EAR or SAR be computed
#' 
#' @export
#' 
#' @examples
#' data(anbo)
#' anbo.obs.sar <- empiricalSAR(anbo$spp, anbo$count, anbo$row, anbo$col, Amin=1, A0=16)
#' plot(anbo.obs.sar)
#' anbo.obs.ear <- empiricalSAR(anbo$spp, anbo$count, anbo$row, anbo$col, Amin=1, A0=16, EAR=TRUE)
#' plot(anbo.obs.ear)
#' 
#' ## empirical SAR from simulated x, y data
#' anbo$x <- runif(nrow(anbo), 0, 1) + anbo$column
#' anbo$y <- runif(nrow(anbo), 0, 1) + anbo$row
#' meteSAR(anbo$spp, anbo$count, x=anbo$x, y=anbo$y, row=4, col=4)
#' 
#' @return an object of class \code{sar} inheriting from \code{data.frame} with 
#' columns \code{A} and \code{S} giving area and species richness, respectively
#'
#' @author Andy Rominger <ajrominger@@gmail.com>, Cory Merow
#' @seealso meteESF, meteSAR, downscaleSAR, upscaleSAR
#' @references Harte, J. 2011. Maximum entropy and ecology: a theory of abundance, distribution, and energetics. Oxford University Press.
# @aliases - a list of additional topic names that will be mapped to
# this documentation when the user looks them up from the command
# line.
# @family sar


empiricalSAR <- function(spp, abund, row, col, x, y, Amin, A0, EAR=FALSE) {
  ## figure out vector of sizes in units of cells; right now only doublings supported
  areaInfo <- .findAreas(
    spp=if(missing(spp)) NULL else spp,
    abund=if(missing(abund)) NULL else abund, 
    row=if(missing(row)) NULL else row, 
    col=if(missing(col)) NULL else col, 
    x=if(missing(x)) NULL else x, 
    y=if(missing(y)) NULL else y, 
    Amin=if(missing(Amin)) NULL else Amin, 
    A0=if(missing(A0)) NULL else A0)
  areas <- areaInfo$areas
  row <- areaInfo$row
  col <- areaInfo$col
  nrow <- areaInfo$nrow
  ncol <- areaInfo$ncol
  Amin <- areaInfo$Amin
  A0 <- areaInfo$A0
  
  ## loop over areas
  out <- lapply(areas, function(a) {
    nspp <- .getSppInGroups(spp, abund, row, col, .getNeighbors(a, nrow, ncol), EAR)
    data.frame(A=a*Amin, S=nspp)
  })
  out <- do.call(rbind, out)
  
  ## make output of class `sar' and tell it about empirical v. theoretical and ear v. sar
  attr(out, 'source') <- 'empirical'
  attr(out, 'type') <- ifelse(EAR, 'ear', 'sar')
  class(out) <- 'sar'
  
  return(out)
}



#================================================================
#' @title Downscale the species area relationship (SAR) or endemics area relationship (EAR)
#'
#' @description Compute METE SAR by downscaling from some larger area \code{A0} to a smaller areas.
#'
#' @details Unlike the other SAR functions, downscaling can be computed for any arbitrary scale
#' \eqn{\leq A_0}. 
#' 
#' @param x an object of class meteESF
#' @param A numerical vector of areas (<= \code{A0}) for which the METE prediction is desired
#' @param A0 total study area
#' @param EAR logical. TRUE computes the endemics area relatinship
#'  
#' @export
#' 
#' @examples
#' data(anbo)
#' anbo.esf <- meteESF(spp=anbo$spp, abund=anbo$count)
#' anbo.thr.downscale <- downscaleSAR(anbo.esf, 2^(seq(-3, 4, length=7)), 16)
#' plot(anbo.thr.downscale)
#' 
#' ## theoretical SARs from state variables only
#' thr.downscale <- downscaleSAR(meteESF(S0=40, N0=400), 2^seq(-1,4,by=1), 16)
#' thr.downscaleEAR <- downscaleSAR(meteESF(S0=40, N0=400), 2^seq(-1, 4, by=1), 16, EAR=TRUE)
#' plot(thr.downscale, ylim=c(0, 40), col='red')
#' plot(thr.downscaleEAR, add=TRUE, col='blue')
#'              
#' @return an object of class \code{sar} inheriting from \code{data.frame} with 
#' columns \code{A} and \code{S} giving area and species richness, respectively
#'
#' @author Andy Rominger <ajrominger@@gmail.com>, Cory Merow
#' @seealso meteESF, meteSAR, empiricalSAR, upscaleSAR
#' @references Harte, J. 2011. Maximum entropy and ecology: a theory of abundance, distribution, and energetics. Oxford University Press.
# @aliases - a list of additional topic names that will be mapped to
# this documentation when the user looks them up from the command
# line.
# @family sar

downscaleSAR <- function(x, A, A0, EAR=FALSE) {
  n0 <- 1:x$state.var['N0']
  
  ## difference between EAR and SAR is for EAR we get Pi(n0) [fun .getPin0]
  ## and for SAR we get 1 - Pi(0) [1 - .getPi0]
  if(EAR) {
    piFun <- function(a) .getPin0(n0, a, A0)
  } else {
    piFun <- function(a) 1 - .getPi0(n0, a, A0)
  }
  
  ## function to get species number at scale `a'
  getspp <- function(a) {
    probs <- piFun(a) * 
      with(x, 
           metePhi(n0, La[1], La[2], Z, 
                   state.var['S0'], state.var['N0'], 
                   ifelse(is.na(state.var['E0']), 1e+06, state.var['E0'])))
    
    return(x$state.var['S0'] * sum(probs))
  }
  
  ## loop over A
  nspp <- sapply(A, getspp)
  
  ## should return matrix with column for area and column for spp
  out <- data.frame(A=A, S=nspp)
  attr(out, 'source') <- 'theoretical'
  attr(out, 'type') <- ifelse(EAR, 'ear', 'sar')
  class(out) <- 'sar'
  
  return(out)
}




#================================================================
#' @title upscale SAR
#'
#' @description Based on information at an anchor scale (\code{A0}) 
#' calcuate predicted species area relationship at larger scales
#'
#' @details Currently only doublings of area are supported and only 
#' the SAR (not EAR) is supported. Upscaling works by iteratively 
#' solving for the constraints (\eqn{S} and \eqn{N} at larger scales)
#' that would lead to the observed data at the anchor scale. See 
#' references for more details on this approach.
#' 
#' 
#' @param x an object of class meteESF
#' @param A0 the anchor scale at which community data are availible.
#' @param Aup the larges area to which to upscale
#' @param EAR logical. TRUE computes the endemics area relatinship; currently not supported
#' 
#' @export
#' 
#' @examples
## combine SAR for scales at which we have data with upscaled SAR
#' data(anbo)
#' anbo.sar <- meteSAR(anbo$spp, anbo$count, anbo$row, anbo$col, Amin=1, A0=16)
#' anbo.sar
#' plot(anbo.sar, xlim=c(1, 2^10), ylim=c(3, 50), log='xy')
#' 
#' ## get upscaled SAR and add to plot
#' anbo.esf <- meteESF(spp=anbo$spp, abund=anbo$count) # need ESF for upscaling
#' anbo.sarUP <- upscaleSAR(anbo.esf, 16, 2^10)
#' plot(anbo.sarUP, add=TRUE, col='blue')
#' 
#'              
#' @return an object of class \code{sar} inheriting from \code{data.frame} with 
#' columns \code{A} and \code{S} giving area and species richness, respectively
#'
#' @author Andy Rominger <ajrominger@@gmail.com>, Cory Merow
#' @seealso meteESF, meteSAR, empiricalSAR, downscaleSAR
#' @references Harte, J. 2011. Maximum entropy and ecology: a theory of abundance, distribution, and energetics. Oxford University Press.
# @aliases - a list of additional topic names that will be mapped to
# this documentation when the user looks them up from the command
# line.
# @family sar

upscaleSAR <- function(x, A0, Aup, EAR=FALSE) {
  ## vector of areas starting with anchor area A0
  Aups <- A0 * 2^(0:ceiling(log(Aup/A0)/log(2)))
  
  ## vector of abundances at each area
  N0s <- x$state.var['N0'] * 2^(0:ceiling(log(Aup/A0)/log(2)))
  
  ## vector of number of species at each area
  S0s <- numeric(length(Aups))
  S0s[1] <- x$state.var['S0']
  
  ## vector to hold termination codes from nleqslv about whether optimization succeeded
  termcodes <- numeric(length(Aups))
  
  ## need to recursively solve constraint fun (solution in `.solveUpscale') up to Aup
  for(i in 2:length(Aups)) {
    S0s[i] <- .solveUpscale(S0s[i-1], N0s[i-1])
  }
  
  ## should return matrix with column for area and column for spp
  out <- data.frame(A=Aups, S=S0s)
  attr(out, 'source') <- 'theoretical'
  attr(out, 'type') <- ifelse(EAR, 'ear', 'sar')
  class(out) <- 'sar'
  
  return(out)
}
