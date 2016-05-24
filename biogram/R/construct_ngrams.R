#' Construct and filter n-grams
#'
#' Step-by-step builds and selects important n-grams. \code{construct_ngrams} extracts 
#' unigrams from the sequences and filters significant features (with p-value below 
#' \code{conf_level}).
#'
#' @param target \code{integer} vector with target information (e.g. class labels).
#' @inheritParams count_ngrams
#' @param n_max size of constructed n-grams.
#' @param conf_level confidence level.
#' @param gap \code{logical}, if \code{TRUE} gap are used. See Details.
#' @return a vector of n-grams.
#' @details Gap parameter determines if \code{construct_ngrams} performs feature 
#' selection on exact n-grams (\code{gap} equal to FALSE) or on all features in the 
#' Hamming distance 1 from the n-gram (\code{gap} equal to TRUE).
#' 
#' @export
#' @seealso
#' Feature filtering method: \code{\link{test_features}}.
#' @examples
#' \donttest{
#' deg_seqs <- degenerate(human_cleave[, 1L:9],
#' list(`1` = c(1, 6, 8, 10, 11, 18),
#'      `2` = c(2, 13, 14, 16, 17),
#'      `3` = c(5, 19, 20),
#'      `4` = c(7, 9, 12, 15),
#'      '5' = c(3, 4)))
#' bigrams <- construct_ngrams(human_cleave[, "tar"], deg_seqs, 1L:5, 2)
#' }

construct_ngrams <- function(target, seq, u, n_max, conf_level = 0.95, gap = TRUE) {
  #build unigrams
  unigrams <- count_ngrams(seq, 1, u, pos = TRUE)
  #filter unigrams
  test_unigrams <- test_features(target, unigrams)
  signif_ngrams <- cut(test_unigrams, breaks = c(0, 1 - conf_level, 1))[[1]]
  
  res <- list()
  
  res[[1]] <- signif_ngrams
  
  for(i in 2L:n_max) {
    if(length(signif_ngrams) != 0) {
      signif_ngrams <- if(i > 2 && gap) {
        build_and_test_gaps(signif_ngrams, seq, i, u, target)
      } else {
        build_and_test(signif_ngrams, seq, i, u, target, conf_level)
      }
          
      res[[i]] <- signif_ngrams
    } else {
      NULL
    }
  }
  res
}

build_and_test <- function(signif_ngrams, seq, n, u, target, conf_level) {
  nplusgrams <- unique(add_1grams(signif_ngrams))
  new_counts <- count_specified(seq, nplusgrams)
  new_test <- test_features(target, new_counts)
  cut(new_test, breaks = c(0, conf_level, 1))[[1]]
}

build_and_test_gaps <- function(signif_ngrams, seq, n, u, target) {
  nplusgrams <- unique(add_1grams(signif_ngrams))
  gap_nplusgrams <- gap_ngrams(nplusgrams)
  new_counts <- count_specified(seq, gap_nplusgrams)
  
  #artificial counts created by sum of counts of all gapped n-grams belonging to n-gram
  art_counts <- do.call(cbind, lapply(1L:length(nplusgrams), function(id_ngram) {
    as.numeric(row_sums(new_counts[, n*(id_ngram - 1) + 1:n]) > 0)
  }))
  
  colnames(art_counts) <- nplusgrams
  
  new_test <- test_features(target, art_counts)
  cut(new_test, breaks = c(0, 0.05, 1))[[1]]
}