#' @title meteSSF
#'  
#' @description \code{meteSSF} calculates the ``spatial structure
#' function'' \eqn{\Pi(n)} (analogous to the ecosystem structure function). From 
#' the SSF the spatial abundance distribution can be calculated.
#'
#' @details
#' Uses either data or state variables to calculate the Spatial Structure 
#' Function (SSF). Uses internal code to determine when computation-saving approximations 
#' can be safely made
#' 
#' @param spp A vector of species names
#' @param sppID  A character giving the name of the desired species (as it appears in `spp') 
#' @param abund A vector of abundances 
#' @param row A vector of row IDs for each observation
#' @param col A vector of column IDs for each observation
#' @param x A vector of x coordinates for each observation
#' @param y A vector of y coordinates for each observation
#' @param n0 Total abundance in area A0
#' @param A The area at which abundances were recorded
#' @param A0 Total study area
#' 
#' 
#' @keywords lagrange multiplier, METE, MaxEnt
#'
#'  @export
#' 
#' @examples
#' data(anbo)
#' ## calculate SSF Pi
#' pi1 <- meteSSF(anbo$spp, 'crcr', anbo$count, row=anbo$row, col=anbo$column, A=1, A0=16)
#' pi1
#' 
#' @return An object of class \code{meteSSF} with elements
#' \describe{
#'   \item{\code{data}}{The data used to construct the SSF}
#'   \item{\code{La}}{Vector of Lagrange multipliers}
#'   \item{\code{La.info}}{Termination information from optimization procedure}
#'   \item{\code{state.var}}{State variables used to constrain entropy maximization}
#' }
#'
#' @author Andy Rominger <ajrominger@@gmail.com>, Cory Merow
#' @seealso metePi
#' @references Harte, J. 2011. Maximum entropy and ecology: a theory of abundance, distribution, and energetics. Oxford University Press.

meteSSF <- function(spp, sppID, abund, row, col, x, y, n0=sum(abund), A, A0) {
  ## get grid regardless of starting input but only if spp and abund are given
  if(!missing(spp) & !missing(abund)) {
    areaInfo <- .findAreas(
      spp=if(missing(abund)) NULL else spp,
      abund=if(missing(abund)) NULL else abund, 
      row=if(missing(row)) NULL else row, 
      col=if(missing(col)) NULL else col, 
      x=if(missing(x)) NULL else x, 
      y=if(missing(y)) NULL else y, 
      Amin=if(missing(A)) NULL else A, 
      A0=if(missing(A0)) NULL else A0)
    areas <- areaInfo$areas
    row <- areaInfo$row
    col <- areaInfo$col
    nrow <- areaInfo$nrow
    ncol <- areaInfo$ncol
    Amin <- areaInfo$Amin
    A0 <- areaInfo$A0
    
    ## focal area might have changed slightly
    A <- areas[which.min(abs(A - areas))]
    
    ## get abundance in grid
    abund[spp != sppID] <- 0
    abund <- .getAbundInGroups(abund, row, col, .getNeighbors(A, nrow, ncol))
  } else {
    abund <- NULL
  }
  
  ## make SSF for state variables
	out <- .makeSSF(n0, A, A0)
	
	out$data$n <- abund
	class(out) <- c('meteSSF', 'meteESF')
	return(out)
}

## helper function to get abundances in all cells
.getAbundInGroups <- function(abund, row, col, groups) {
  cellID <- paste(row, col, sep=',')
  cellGroup <- as.factor(groups$group[match(cellID, groups$cells)])
  levels(cellGroup) <- unique(groups$group)
  out <- tapply(abund, cellGroup, sum)
  out[is.na(out)] <- 0
  
  return(out)
}

#==============================================================================
#' @title Equation of the PMF of the METE spatial species abundance distribution
#'
#' @description
#' \code{metePi} is a low level function that returns the spatial species abundance 
#' distribution \eqn{Pi(n)} predicted by METE; vectorized in n
#'
#' @details
#' See Examples
#' 
#' @param n A vector giving abundances of each entry
#' @param la The spatial Lagrange multiplier returned by \code{meteSSF}
#' @param n0 Total abundance in area A0
#' @export
#' 
#' @examples
#' metePi(0:10, 0.01, 100)
#' 
#' @return a numeric vector giving the probability of each entry in \code{n}
#'
#' @author Andy Rominger <ajrominger@@gmail.com>, Cory Merow
#' @seealso metePi
#' @references Harte, J. 2011. Maximum entropy and ecology: a theory of abundance, distribution, and energetics. Oxford University Press.

metePi <- function(n,la,n0) {
	1/.mete.Pi.Z(la, n0) * exp(-la*n)
}



##==========================================
## helper functions for meteSSF and metePi

##	normalization constant for Pi
.mete.Pi.Z <- function(la,n0) {
	if(la != 0) {
		return((1-exp(-la*(n0+1)))/(1-exp(-la)))
	} else {
		return(n0+1)
	}
}

##	constraint function for Pi lagrange multiplier
##	function: x = e^(-la)
##  vectorized over `x'
.pi.cons <- function(x, n0, A, A0) {
	lhs <- rep(NA,length(x))
	
	case1 <- x > 1 & (n0+1)*log(x) <= log(2e+64)
	case2 <- x > 1 & (n0+1)*log(x) > log(2e+64)
	case3 <- x < 1
	case4 <- x == 1
	case0 <- case1 | case3
	
	lhs[case0] <- x[case0]/(1-x[case0]) - ((n0 + 1)*x[case0]^(n0+1))/(1-x[case0]^(n0+1))
	lhs[case2] <- x[case2]/(1-x[case2]) + n0 + 1
	lhs[case4] <- n0/2
	
	return(lhs - n0*A/A0)
}


##	make *S*patial *S*tructure *F*unction, like the ESF slot in mete class
#' @importFrom stats uniroot
.makeSSF <- function(n0, A, A0, eq52=.useEq52(n0,A,A0)) {
	if(A/A0 == 0.5) {
		return(list(La=0, La.info='analytic solution',
                state.var=c(n0=n0,A=A,A0=A0)))
	} else if(eq52) {
		La <- -log((n0*A/A0)/(1+n0*A/A0))
		return(list(La=La, La.info='analytic solution',
                state.var=c(n0=n0,A=A,A0=A0)))
	} else {
		if(.pi.cons(1-.Machine$double.eps^0.45,n0,A,A0) > 0) {
			upper <- 1-.Machine$double.eps^0.45
		} else {
			upper <- 1.01*(n0*(A/A0-1)-1)/(n0*(A/A0-1))
		}
		sol <- uniroot(.pi.cons,c(0,upper),
                   n0, A, A0,
                   tol=.Machine$double.eps^0.75)
		
		La <- -log(sol$root)
		
		return(list(La=La,La.info=sol[-1],
                state.var=c(n0=n0,A=A,A0=A0)))
	}
}


## anything with n0 > 8000 and A/A0 < 0.5 use approx,
## otherwise follow this eq (validation in `check_eq52.R')
## vectorized over `n0'
.useEq52 <- function(n0,A,A0) {
	test <- exp(4.9 -1.36*log(n0) + 
                0.239*log(n0)^2 -0.0154*log(n0)^3)
	res <- A0/A >= test
	
	res[n0 > 2^16 & A/A0 < 0.5] <- TRUE
	res[A/A0 == 0.5] <- TRUE
	
	return(res)
}

