dif <- function(seq1, seq2, gap = FALSE, aa.strict = FALSE) {

  if (length(seq1) !=  length(seq2))
    stop("seq1 and seq2 are not of the same length")

  #turn !aa and !gap in NA
  seq1[!as.logical(is.aa(seq1, aa.strict) + is.gap(seq1))] <- NA
  seq2[!as.logical(is.aa(seq2, aa.strict) + is.gap(seq2))] <- NA

  #turn gap into NA
  if (!gap) {
    seq1[is.gap(seq1)] <- NA
    seq2[is.gap(seq2)] <- NA
  }

  #complete.cases takes NA into account
  res <- sum(seq1 != seq2, na.rm = TRUE)/sum(complete.cases(cbind(seq1, seq2)))
  return (res)
}