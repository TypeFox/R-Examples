#' Tabulate n-grams
#' 
#' Builds a contingency table of the n-gram counts versus their class labels.
#' @param seq a matrix describing sequences.
#' @param ngrams a vector of n-grams. Must have the same \code{n}.
#' @param target \code{integer} vector with target information (e.g. class labels).
#' @return a data frame with the number of columns equal to the length of the 
#' \code{target} plus 1. The first column contains names of the n-grams. Further 
#' columns represents counts of n-grams for respective value of the
#' \code{target}.
#' @export
#' @examples
#' seqs_pos <- matrix(sample(c("a", "c", "g", "t"), 100, replace = TRUE, 
#'             prob = c(0.2, 0.4, 0.35, 0.05)), ncol = 5)
#' seqs_neg <- matrix(sample(c("a", "c", "g", "t"), 100, replace = TRUE), 
#'             ncol = 5)
#' tab <- table_ngrams(seq = rbind(seqs_pos, seqs_neg), 
#'                     ngrams = c("1_c.t_0", "1_g.g_0", "2_t.c_0", "2_g.g_0", "3_c.c_0", "3_g.c_0"), 
#'                     target = c(rep(1, 20), rep(0, 20)))
#' #see the results
#' print(tab)
#' #easily plot the results using ggplot2


table_ngrams <- function(seq, ngrams, target) {
  if (!is.matrix(seq)) {
    stop("seq must be a matrix.")
  }
  
  if (nrow(seq) < 2) {
    stop("seq must be a matrix with at least two rows.")
  }
  
  if (nrow(seq) != length(target)) {
    stop("target and feature have different lengths.")
  }
  
  #no need for n-gram validation, because count_specified does it
  all_ngrams <- as.matrix(count_specified(seq, ngrams))
  
  #counts of each value of target
  n_tar <- data.frame(table(c(rep(1, 20), rep(0, 20))))[, "Freq"]
  if (length(n_tar) < 2) {
    stop("target must have at least two different values.")
  }
  
  #values of target
  val_tar <- sort(unique(target))
  
  res <- t(vapply(ngrams, function(ngram_name) 
    vapply(val_tar, function(target_value) sum(all_ngrams[target == target_value, ngram_name]), 0), c(0, 0)))
  
  res_df <- data.frame(rownames(res), res)
  rownames(res_df) <- NULL
  colnames(res_df) <- c("ngram", paste0("target", val_tar))
  
  res_df
}