
#'Generate weight matrix
#'
#'Generates a weight matrix given a group of ages and years
#'and a set of cohorts which are to be given zero weight. This
#'is useful for excluding some data points when fitting a 
#'Stochastic Mortality Model (see \code{\link{fit.StMoMo}}).
#'
#'@param ages vector of ages.
#'@param years vector of years.
#'@param clip number of cohorts in the boundary to assign a zero 
#'weight. This can be be used to zero weigh some of the first and 
#'last cohorts in the data.
#'@param zeroCohorts other cohort for which a zero weight is to be assigned.
#'   
#'@return A 0-1 matrix with 0 for the zero-weighed cohorts.
#'
#'@seealso \code{\link{fit.StMoMo}}
#'
#'@examples
#' #Zero-weight the first three and last three cohorts
#' wxt1 <- genWeightMat(55:89,  EWMaleData$years, clip = 3)
#' APCfit1 <- fit(apc(), Dxt = EWMaleData$Dxt, Ext = EWMaleData$Ext, 
#'                ages = EWMaleData$ages, years = EWMaleData$years, 
#'                ages.fit = 55:89, wxt = wxt1)
#' plot(APCfit1, parametricbx = FALSE, nCol = 3)
#' 
#' #Also Zero-weight the 1886 cohort
#' wxt2 <- genWeightMat(55:89,  EWMaleData$years, clip = 3, 
#'                      zeroCohorts = 1886)
#' APCfit2 <- fit(apc(), Dxt = EWMaleData$Dxt, Ext = EWMaleData$Ext, 
#'                ages = EWMaleData$ages, years = EWMaleData$years, 
#'                ages.fit = 55:89, wxt = wxt2)
#' plot(APCfit2, parametricbx = FALSE, nCol = 3)
#'
#'@export
genWeightMat <- function(ages, years, clip = 0, zeroCohorts = NULL) {
  nAges <- length(ages)
  nYears <- length(years)  
  cohorts <- (years[1] - ages[nAges]):(years[nYears] - ages[1])
  nCohorts <- length(cohorts)  
  wxt <- matrix(1, nrow = nAges, ncol = nYears)
  rownames(wxt) <- ages
  colnames(wxt) <- years
  if (clip > 0) {
    zeroCohorts <- c(zeroCohorts, head(cohorts, clip), tail(cohorts, clip))
  }
  for (c in zeroCohorts){
    h <- c - cohorts[1] + 1 - nAges
    if (h <= 0){
      col <- 1
      row <- -h + 1
    } else {
      row <- 1
      col <- h + 1
    }
    while (col <= nYears && row <= nAges) {
      wxt[row, col] <- 0
      row <- row + 1
      col <- col + 1          
    }  
  }
  wxt
}
