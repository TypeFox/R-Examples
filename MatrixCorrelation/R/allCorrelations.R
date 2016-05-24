#' @title All correlations
#'
#' @description Compare all correlation measures in the package (or a subset)
#'
#' @param X1 first \code{matrix} to be compared (\code{data.frames} are also accepted).
#' @param X2 second \code{matrix} to be compared (\code{data.frames} are also accepted).
#' @param ncomp1 maximum number of subspace components from the first \code{matrix}.
#' @param ncomp2 maximum number of subspace components from the second \code{matrix}.
#' @param methods \code{character} vector containing a subset of the supported methods: "SMI", "RV", "RV2", "RVadj", "r1", "r2", "r3", "r4", "GCD".
#' @param digits number of digits for numerical output.
#' @param plot logical indicating if plotting should be performed (default = TRUE).
#' @param xlab optional x axis label.
#' @param ylab optional y axis label.
#' @param ... additional arguments for \code{SMI} or \code{plot}.
#'
#' @details For each of the three coefficients a single scalar is computed to describe
#' the similarity between the two input matrices.
#'
#' @return A single value measuring the similarity of two matrices.
#'
#' @author Kristian Hovde Liland
#'
#' @references
#' \itemize{
#'  \item{RV:}{ Robert, P.; Escoufier, Y. (1976). "A Unifying Tool for Linear Multivariate
#'   Statistical Methods: The RV-Coefficient". Applied Statistics 25 (3): 257-265.}
#'  \item{RV2:}{ Smilde, AK; Kiers, HA; Bijlsma, S; Rubingh, CM; van Erk, MJ (2009). "Matrix correlations
#'  for high-dimensional data: the modified RV-coefficient". Bioinformatics 25(3): 401-5.}
#'  \item{Adjusted RV:}{ Maye, CD; Lorent, J; Horgan, GW. (2011). "Exploratory analysis of multiple omics
#'  datasets using the adjusted RV coefficient". Stat Appl Genet Mol Biol. 10(14).}
#' }
#'
#' @seealso \code{\link{SMI}}, \code{\link{RV}} (RV2/RVadj), \code{\link{r1}} (r2/r3/r4/GCD).
#'
#' @examples
#' X1  <- scale( matrix( rnorm(100*300), 100,300), scale = FALSE)
#' usv <- svd(X1)
#' X2  <- usv$u[,-3] %*% diag(usv$d[-3]) %*% t(usv$v[,-3])
#'
#' allCorrelations(X1,X2, ncomp1 = 5,ncomp2 = 5)
#'
#' @export
allCorrelations <- function(X1,X2,
                            ncomp1, ncomp2,
                            methods = c("SMI","RV","RV2","RVadj","r1","r2","r3","r4","GCD"),
                            digits = 3, plot = TRUE, xlab = '', ylab = '', ...){
  # Handle components
  if(missing(ncomp1) && ("SMI" %in% methods || "GCD" %in% methods)){
    ncomp1 <- min(5, Rank(X1)-1)
    warning(paste('ncomp1 not specified, defaulting to ', ncomp1, sep=""))
  }
  if(missing(ncomp2) && ("SMI" %in% methods || "GCD" %in% methods)){
    ncomp2 <- min(5, Rank(X2)-1)
    warning(paste('ncomp2 not specified, defaulting to ', ncomp2, sep=""))
  }

  # Compute all correlations
  n.met   <- length(methods)
  results <- numeric(n.met)
  for(i in 1:n.met){
    results[i] <- switch(methods[i],
                         SMI   = SMI(X1,X2, ncomp1,ncomp2, ...)[ncomp1,ncomp2],
                         RV    = RV(X1,X2),
                         RV2   = RV2(X1,X2),
                         RVadj = RVadj(X1,X2),
                         r1    = r1(X1,X2),
                         r2    = r2(X1,X2),
                         r3    = r3(X1,X2),
                         r4    = r4(X1,X2),
                         GCD   = GCD(X1,X2, ncomp1,ncomp2))
  }

  # Optionally plot
  if(plot == TRUE){
	  plot(1:n.met,results, type='h', axes=FALSE, # ylim = c(ifelse(any(results<0),-1,1),1.1),
		   xlab = xlab, ylab = ylab, panel.first = grid(nx=NA, ny=NULL), ...)
	  points(which(results>=0),results[results>=0],col='black')
	  if(any(results<0))
		points(which(results<0),results[results<0],col='red')
	  # grid(nx=NA, ny=NULL)
	  axis(3,1:n.met,methods)
	  axis(1,1:n.met,methods)
	  axis(2)
	  box()
  }

  # Process return values
  names(results) <- methods
  if(!is.null(digits) && !is.na(digits)){
    results <- round(results,digits)
  }
  results
}
