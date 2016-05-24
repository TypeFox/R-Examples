#' BABLBS (Band Added Band Lost Band Substituted)
#'
#' BABLBS clustering means that two fingerprints match if they 
#' differ only through a lost band, a gained band, or a substituted band
#'
#' @param intable is the name of the acm table.
#' @param minmatch is the threshold of band number such that two lanes are matching if and only if the number of matching bands between them is more than minmatch.
#' @return a matrix of 2 columns where each row represents a pair of matching lanes.
#' @references 
#' Salamon et. al (1998) Accommodating Error Analysis in Comparison and Clustering of Molecular Fingerprints. Emerging Infectious Diseases Vol. 4, No. 2, April-June 1998
#' @references Abasci LLC. JAMES v1.0 User Documentation. 2002. 
#' @author Andrea Benedetti \email{andrea.benedetti@@mcgill.ca}
#' @author Sahir Rai Bhatnagar
#' @author XiaoFei Zhao
#' @note \code{call_bablbs} doesn't work for minmatch<5
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
#' #bablbs, GD1, and GD2 all work with the results from call_acm
#' #now check matching by bablbs 
#' res_bab<-call_bablbs(res1)
#' @export 
call_bablbs<-function (intable, minmatch = 2) 
{
  if (ncol(intable) < 5) {
    fprintf(stderr(), "Error in call_bablbs: the input intable should have at least 5 columns")
    stop(1)
  }
  numout = 0
  outtable = matrix(data = NA, nrow = nrow(intable), ncol = 2)
  for (i in 1:nrow(intable)) {
    scan_uli_1 = as.character(intable[i, 1])
    scan_uli_2 = as.character(intable[i, 2])
    m = as.integer(intable[i, 3])
    n = as.integer(intable[i, 4])
    k = as.integer(intable[i, 5])
    if (!is.integer(m) | !is.integer(n) | !is.integer(k)) {
      fprintf(stderr(), "Error in bablbs: The %dth row of the intable table is not of the format \n\t\t\t\t\t[string, string, integer, integer, integer]",               i)
      stop(2)
    }
    if (k >= minmatch & ((n - m == 1 & k == m) | (m - n ==1 & k == n) | (m == n & k == m - 1) | (m == n & k == 
                                                                                             m))) {
      numout = numout + 1
      outtable[numout, 1] = scan_uli_1
      outtable[numout, 2] = scan_uli_2
    }
  }
  outtable = outtable[!is.na(outtable[, 1]) & !is.na(outtable[, 2]),]
  #fprintf(stdout(), "%d pairs read of which %d were ba-bl-bs matches with at least %d bands", nrow(intable), numout, minmatch)
  outtable}
