#' Count specified n-grams
#'
#' Counts specified n-grams in the input sequence(s).
#'
#' @param seq a vector or matrix describing sequence(s). 
#' @param ngrams a vector of n-grams. Must have the same \code{n}.
#' @return a \code{\link[slam]{simple_triplet_matrix}} where columns represent
#' n-grams and rows sequences.
#' @export
#' @details \code{\link{count_specified}} counts only selected n-grams declared by
#' user in the \code{ngrams} parameter. Declared n-grams must be written using the
#' \code{biogram} notation.
#' @seealso Count all possible n-grams: \code{\link{count_ngrams}}.
#' @examples
#' seqs <- matrix(c(1, 2, 2, 1, 2, 1, 2, 1, 1, 2, 3, 4, 1, 2, 2, 4), nrow = 2)
#' count_specified(seqs, ngrams = c("1.1.1_0.0", "2.2.2_0.0", "1.1.2_0.0"))
#' 
#' seqs <- matrix(sample(1L:5, 200, replace = TRUE), nrow = 20)
#' count_specified(seqs, ngrams = c("2_4.2_0", "2_1.4_0", "3_1.3_0",
#'                                  "2_4.2_1", "2_1.4_1", "3_1.3_1",
#'                                  "2_4.2_2", "2_1.4_2", "3_1.3_2"))

count_specified <- function(seq, ngrams) {
  #validate n-grams
  validated_ngram <- sapply(ngrams, is_ngram)
  if(!all(validated_ngram))
    stop("Improper n-grams: ", paste(names(which(!validated_ngram)), collapse = ", "))
  
  #if sequence is not a matrix (single sequence), convert it to matrix with 1 row
  if (class(seq) != "matrix")
    seq <- matrix(seq, nrow = 1)
  
  #length of sequence
  len_seq <- ncol(seq)
  #number of sequences
  n_seqs <- nrow(seq)
  
  df <- ngrams2df(ngrams)
  
  #splitted ngrams
  sn_grams <- strsplit(df[, "ngram"], ".", fixed = TRUE)
  
  #n in ngram
  n <- unique(vapply(sn_grams, length, 1))
  if(length(n) > 1)
    stop("n-grams must have the same n.")
  
  #n-grams are grouped by their unique distance to speed up function
  #when ncol(df) == 3 n-grams are positioned
  res <- if(ncol(df) == 3) {
    do.call(cbind, lapply(unique(df[, "distance"]), function(unique_dist) {
      dist_df <- df[df[, "distance"] == unique_dist, ]
      #all possible n-gram positions
      all_ngram_pos <- get_ngrams_ind(len_seq, n, 
                                      as.numeric(strsplit(df[, "distance"], ".", fixed = TRUE)[[1]]))
      vapply(1L:nrow(dist_df), function(ngram_id)
        vapply(1L:n_seqs, function(single_seq) {
          #positions of the n-gram of interest
          single_ngram_pos <- sapply(all_ngram_pos, function(single_pos) 
            single_pos[dist_df[ngram_id, "position"]])
          as.numeric(all(as.character(seq[single_seq, single_ngram_pos]) == sn_grams[[ngram_id]]))
        }, 0), rep(0, n_seqs))
    }))
  } else {
    do.call(cbind, lapply(unique(df[, "distance"]), function(unique_dist) {
      #unpositioned n-grams
      dist_df <- df[df[, "distance"] == unique_dist, ]
      #all possible n-gram positions
      all_ngram_pos <- get_ngrams_ind(len_seq, n, 
                                      as.numeric(strsplit(df[, "distance"], ".", fixed = TRUE)[[1]]))
      vapply(1L:nrow(dist_df), function(ngram_id)
        vapply(1L:n_seqs, function(single_seq) {
          #positions of the n-gram of interest
          all_ngram_pos <- do.call(rbind, all_ngram_pos)
          sum(apply(all_ngram_pos, 2, function(single_ngram_pos)
            as.numeric(all(as.character(seq[single_seq, single_ngram_pos]) == sn_grams[[ngram_id]]))), na.rm = TRUE)
        }, 0), rep(0, n_seqs))
    }))
  }
  
  #reoder results - were shuffled because of indexing on unqiue distance
  res <- res[, order(unlist(lapply(unique(df[, "distance"]), function(unique_dist) 
    which(df[, "distance"] == unique_dist, ))))]
  
  if(class(res) == "numeric") {
    res <- matrix(res, ncol = 1)
  }
  #name columns
  colnames(res) <- ngrams
  
  as.simple_triplet_matrix(res)
}
