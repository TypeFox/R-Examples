#' Degenerate protein sequence
#'
#' 'Degenerates' amino acid or nucleic sequence by aggregating 
#' elements to bigger groups.
#' 
#' @param seq \code{character} vector or matrix representing single sequence.
#' @param element_groups list of groups to which elements of sequence should be aggregated.
#' @keywords manip
#' @return a \code{character} vector or matrix (if input is a matrix) 
#' containing aggregated elements.
#' @note Both sequence and \code{element_groups} should contain lower-case letters.
#' Upper-case will be automatically converted without a message.
#' 
#' Characters not present in the \code{element_groups} will be converted to NA with a 
#' warning.
#' @export
#' @seealso \code{\link{l2n}} to easily convert information stored in biological sequences from 
#' letters to numbers.
#' @keywords manip
#' @examples
#' sample_seq <- c(1, 3, 1, 3, 4, 4, 3, 1, 2)
#' table(sample_seq)
#' 
#' #aggregate sequence to purins and pyrimidines
#' deg_seq <- degenerate(sample_seq, list(w = c(1, 4), s = c(2, 3)))
#' table(deg_seq)

degenerate <- function(seq, element_groups) {
  tmp_seq <- tolower(seq)
  if (!all(unique(tmp_seq) %in% unlist(element_groups))) {
    warning("'seq' contains elements not present in 'element_groups'. Such elements will be replaced by NA.")
    tmp_seq[!(tmp_seq %in% unlist(element_groups))] <- NA
  }
  
  
  for (i in 1L:length(element_groups)) {
    tmp_seq[tmp_seq %in% element_groups[[i]]] <- names(element_groups)[i]
  }
  
  if(class(seq) == "matrix")
    dim(tmp_seq) <- dim(seq)
  
  tmp_seq
}

#' Convert letters to numbers
#'
#' Converts biological sequence from letter to number notation.
#' @param seq \code{character} vector representing single sequence.
#' @param seq_type the type of sequence. Can be \code{rna}, \code{dna} or \code{prot}.
#' @keywords manip
#' @return a \code{numeric} vector containing converted elements.
#' @export
#' @keywords manip
#' @seealso 
#' \code{l2n} is a wrapper around \code{\link{degenerate}}.
#' 
#' Inverse function: \code{\link{n2l}}.
#' @examples
#' sample_seq <- c("a", "d", "d", "g", "a", "g", "n", "a", "l")
#' l2n(sample_seq, "prot")

l2n <- function(seq, seq_type) {
  if (!(seq_type %in% c("prot", "dna", "rna")))
    stop("The value of 'what' must be: 'dna', 'rna' or 'prot'.")
  elements_list <- switch(seq_type,
                          rna = c("a", "c", "g", "u"),
                          dna = c("a", "c", "g", "t"),
                          prot = c("a", "c", "d", "e", "f", 
                                   "g", "h",  "i", "k", "l", 
                                   "m", "n", "p", "q", "r", 
                                   "s", "t", "v", "w", "y"))
  names(elements_list) <- 1L:length(elements_list)
  as.numeric(degenerate(seq, elements_list))
}


#' Convert numbers to letters
#'
#' Converts biological sequence from number to letter notation.
#' @param seq \code{numeric} vector representing single sequence.
#' @param seq_type the type of sequence. Can be \code{rna}, \code{dna} or \code{prot}.
#' @keywords manip
#' @return a \code{numeric} vector containing converted elements.
#' @export
#' @keywords manip
#' @seealso 
#' \code{n2l} is a wrapper around \code{\link{degenerate}}.
#' 
#' Inverse function: \code{\link{l2n}}.
#' @examples
#' sample_seq <- c(1, 3, 3, 6, 1, 6, 12, 1, 10)
#' n2l(sample_seq, "prot")

n2l <- function(seq, seq_type) {
  if (!(seq_type %in% c("prot", "dna", "rna")))
    stop("The value of 'what' must be: 'dna', 'rna' or 'prot'.")
  names_list <- switch(seq_type,
                          rna = c("a", "c", "g", "u"),
                          dna = c("a", "c", "g", "t"),
                          prot = c("a", "c", "d", "e", "f", 
                                   "g", "h",  "i", "k", "l", 
                                   "m", "n", "p", "q", "r", 
                                   "s", "t", "v", "w", "y"))
  elements_list <- 1L:length(names_list)
  names(elements_list) <- names_list
  degenerate(seq, elements_list)
}