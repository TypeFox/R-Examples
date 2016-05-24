#' \code{q} Nearest Neighbors Test
#' 
#' \code{qnn.test} calculates statistics related to the q nearest neighbors method of 
#' comparing case and control point patterns under the random labeling hypothesis.
#' 
#' @param x A \code{ppp} object from the \code{spatstat} package with marks for the case and control groups.
#' @param q A vector of positive integers indicating the values of \code{q} for which to do the q nearest neighbors test.
#' @param nsim The number of simulations from which to compute p-value.
#' @param case The position of the name of the "case" group in levels(x$marks).  The default is 2.
#' 
#' @return Returns a list with the following components: 
#' \item{qsum}{A dataframe with the number of neighbors (q), test statistic (Tq), and p-value for each test.}
#' \item{consum}{A dataframe with the contrasts (contrast), test statistic (Tcon), and p-value (pvaluecon) for the test of contrasts.}
#' @author Joshua French
#' @import spatstat
#' @importFrom stats dist
#' @importFrom utils combn
#' @export
#' @references Waller, L.A., and Gotway, C.A. (2005).  Applied Spatial Statistics for Public Health Data.  Hoboken, NJ: Wiley.  Cuzick, J., and Edwards, R. (1990). Spatial clustering for inhomogeneous populations. Journal of the Royal Statistical Society. Series B (Methodological), 73-104.  Alt, K.W., and Vach, W. (1991). The reconstruction of "genetic kinship" in prehistoric burial complexes-problems and statistics. Classification, Data Analysis, and Knowledge Organization, 299-310.
#' @examples 
#' data(grave)
#' qnn.test(grave, q = c(3, 5, 7, 9, 11, 13, 15))
#' 
qnn.test = function(x, q = 5, nsim = 499, case = 2)
{
  
  if(!is.element("ppp", class(x))) stop("x must be a ppp object")
  if(is.null(x$marks)) stop("x must be marked as cases or controls")
  if(!is.factor(x$marks)) stop("The marks(x) must be a factor")
  nlev = length(levels(x$marks))
  if(case < 1 || case > nlev) stop("case must be an integer between 1 and length(levels(x$marks))")
  if(nsim < 0 | !is.finite(nsim)) stop("nsim must be a non-negative integer")
  
  # the case and index of cases
  casename = levels(x$marks)[case]
  idxcase = which(x$marks == casename)
  
  # number of observations
  N = x$n
  N1 = length(idxcase) # number of cases
  
  q = sort(q) # sorted q values
  maxq = max(q) # largest q
  
  # distances, then ordered
  D = as.matrix(stats::dist(cbind(x$x, x$y))) #calculate pairwise distances
  diag(D) = Inf #set diagonal to infinity
  odistq = t(apply(D, 1, order)[1:maxq, ]) #Find maxq nearest neighbors
  
  rcases = matrix(0, ncol = N1, nrow = nsim) #matrix to store index of random cases
  for(i in 1:nsim) rcases[i,] = sample(1:N, N1) #sample "cases"
  Delta = matrix(0, nrow = N, ncol = nsim + 1) # storage of indicator for cases for each nsim and observed
  for(i in 1:nsim) Delta[rcases[i,],i] <- 1 #indicate whether observation i in sample i is a case
  Delta[idxcase,nsim + 1] = 1 # cases for observed data
  Tq = matrix(0, nrow = length(q), ncol = nsim + 1) # store test statistics
  con = utils::combn(1:length(q), 2) #pairwise contrasts
  Tcon = matrix(0, nrow = ncol(con), ncol = nsim + 1) # store constrast test statistics
  
  for(i in 1:length(q))
  {
    Wq = matrix(0, nrow = N, ncol = N)
    for(j in 1:N)
    {
      nn = odistq[j, 1:q[i]]
      Wq[j, nn] = 1
    }
    Tq[i,] = rowSums(crossprod(Delta, Wq) * t(Delta)) # calculate test statistics
  }
  
  p = rowMeans(Tq >= Tq[, nsim + 1]) #calculate p-values
  
  # calculate test statistics for pairwise contrasts
  for(i in 1:ncol(con)) Tcon[i, ] = Tq[con[2, i], ] - Tq[con[1,i], ]
  
  connames = vector("character", ncol(con))
  for(i in 1:ncol(con)) connames[i] = paste("T",q[con[2,i]]," - T", q[con[1,i]], sep="")
  pcon = rowMeans(Tcon >= Tcon[, nsim + 1]) #calculate p-values
  
  return(list(qsum = data.frame(q = q, Tq = Tq[,nsim + 1], pvalue = p), 
              consum = data.frame(contrast = connames, 
                                  Tcontrast = Tcon[,nsim + 1], pvaluecon = pcon)))
  
}



