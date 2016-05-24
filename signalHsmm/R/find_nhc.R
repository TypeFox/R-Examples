#' Localize n-, h- and c-region in signal peptide
#'
#' Finds borders between distinct regions constituting signal peptides using a heuristic
#' algorithm.
#'
#' @param protein a vector of amino acids or object of class 
#' \code{\link[seqinr]{SeqFastaAA}}.
#' @param signal range of signal peptide. If \code{NULL}, the attribute \code{sig}
#' of \code{protein} will be used.
#' @return a vector of length 4 containing positions of:
#' \enumerate{
#'   \item start of n-region,
#'   \item start of h-region,
#'   \item start of c-region,
#'   \item cleavage site.
#' }
#' @references Henrik Nielsen, Anders Krogh (1998). Prediction of signal peptides
#' and signal anchors by a hidden Markov model. \emph{Proc. Sixth Int. Conf. on 
#' Intelligent Systems for Molecular Biology}.
#' @export

find_nhc <- function(protein, signal = NULL) {
  protein <- toupper(protein)
  if (is.null(signal)) 
    signal <- attr(protein, "signal")
  
  sig <- protein[signal[1]:signal[2]]
  start_c <- length(sig) - 2
  
  #noh number of hydrophobic residues
  noh <- 0
  while(noh < 2 && start_c > 1) {
    start_c <- start_c - 1
    noh <- noh + ifelse(sig[start_c] %in% c("A", "I", "L", "M", "F", "W", "V"), 1, -1)
    noh <- ifelse(noh < 0, 0, noh)
  }
  start_c <- start_c + 2
  
  start_h <- start_c - 6
  #if statement to prevent negative values
  if (start_h > 1) {
    #nonh number of nonhydrophobic residues
    nonh <- 0
    #noc number of charged
    noc <- 0
    while(nonh < 3 && noc == 0 && start_h > 1) {
      start_h <- start_h - 1
      nonh <- nonh + ifelse(sig[start_h] %in% c("A", "I", "L", "M", "F", "W", "V"), -1, 1)
      nonh <- ifelse(nonh < 0, 0, nonh)
      noc <- ifelse(sig[start_h] %in% c("R", "H", "K", "D", "E"), 1, 0)
    }
  } else {
    start_h <- 1
  }
  
  prestart_c <- start_c - 1
  noh <- 0
  while(noh == 0 && start_h < prestart_c) {
    start_h <- start_h + 1
    noh <- noh + ifelse(sig[start_h] %in% c("A", "I", "L", "M", "F", "W", "V"), 1, 0)
  }
  #c(start_n = signal[1], start_h = start_h, start_c = start_c, cs = signal[2])
  c(start_n = 1, start_h = start_h, start_c = start_c, cs = signal[2])
}
