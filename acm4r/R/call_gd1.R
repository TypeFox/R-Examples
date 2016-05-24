#' Genetic Distance 1
#'
#' Genetic Distance 1 which includes a lambda parameter
#'
#' @param intable is the name of the acm table.
#' @param gcrit is the threshold of genetic distance such that two lanes are matching if and only if the genetic distance between them is less than gcrit. Is used to get a list of ids for subjects who have gd distance less than gcrit
#' @param lambda is the minimum number of bands that need to match
#' @return A matrix with two colums of fingerprint IDs. Each line represents a match according to \code{gcrit}.
#' @references 
#' Salamon et. al (1998) Accommodating Error Analysis in Comparison and Clustering of Molecular Fingerprints. Emerging Infectious Diseases Vol. 4, No. 2, April-June 1998
#' @references Abasci LLC. JAMES v1.0 User Documentation. 2002. 
#' @author Andrea Benedetti \email{andrea.benedetti@@mcgill.ca}
#' @author Sahir Rai Bhatnagar
#' @author XiaoFei Zhao
#' @note Pass the result of this function to and pass this result to \code{clusters} to synthesize the results
#' @examples 
#' \dontshow{
#' (WD <- getwd())
#' if (!is.null(WD)) setwd(WD)
#' data(replicates.in)
#' write.table(replicates.in,  "replicates.in", quote=FALSE, row.names=FALSE)
#' data(experiments.in)
#' write.table(experiments.in, "experiments.in", quote=FALSE, row.names=FALSE)
#' call_erra("replicates.in",  dnum=1, sd=1, delete=TRUE)
#' res1<-call_acm("experiments.in")
#' }
#' #matching by GD1
#' res_gd1<-call_gd1(res1)
#' @export 
call_gd1<-function (intable, gcrit = 5, lambda = 0.05) 
{
  if (ncol(intable) < 5) {
    fprintf(stderr(), "Error in call_gd1: the input intable should have at least 5 columns")
    stop(1)
  }
  numout = 0
  outtable = matrix(data = "", nrow = nrow(intable), ncol = 2)
  for (i in 1:nrow(intable)) {
    scan_uli_1 = as.character(intable[i, 1])
    scan_uli_2 = as.character(intable[i, 2])
    m = as.integer(intable[i, 3])
    n = as.integer(intable[i, 4])
    k = as.integer(intable[i, 5])
    if (!is.integer(m) | !is.integer(n) | !is.integer(k)) {
      fprintf(stderr(), "Error in call_gd1: The %dth row of the intable table is not of the format \n\t\t\t\t\t[string, string, integer, integer, integer]", 
              i)
      stop(2)
    }
    d <- rawgd1(n, m, k)/lambda
    if (d <= gcrit) {
          numout = numout + 1
      outtable[numout, 1] = scan_uli_1
      outtable[numout, 2] = scan_uli_2
    }
   
  }
  outtable<-outtable[outtable[,1]!="",]
  return(outtable)
}


rawgd1<-
function (n, m, k) 
{
  if (n < 1 | m < 1) {
    return - 1
  }
  if (k > n | k < 0) {
    return - 2
  }
  if (m < n) {
    i <- n
    n <- m
    m <- i
  }
  sum <- 0
  if (!(n == m & k == m)) {
    if (n <= m - 1) {
      for (i in n:(m - 1)) {
        sum <- sum + 1/i
      }
    }
    sum2 <- 0
    if (k > 0) {
      if (k <= n - 1) {
        for (i in k:(n - 1)) {
          sum2 <- sum2 + 1/i
        }
      }
    }
    else {
      sum2 <- sum2 + 1
      if (k + 1 <= n - 1) {
        for (i in (k + 1):(n - 1)) {
          sum2 <- sum2 + 1/i
        }
      }
    }
    sum <- sum + 2 * sum2
  }
  else {
    sum <- 0.5/m
  }
  return(sum)
}

fprintf<-
function (filename = "", format = "", ...) 
{
  write(sprintf(format, ...), file = filename, append = TRUE)
}