#' Get all possible n-Grams
#'
#' Creates vector of all posible n_grams (for given \code{n}).
#'
#' @inheritParams count_ngrams
#' @param possible_grams number of possible n-grams. If not \code{NULL} n-grams do not
#' contain information about position
#' @return a character vector. Elements of n-gram are separated by dot.
#' @note Input data must be a matrix or data frame of numeric elements.
#' @details See Details section of \code{\link{count_ngrams}} for more 
#' information about n-grams naming convention. The possible information about distance 
#' must be added by hand (see examples).
#' @export
#' @examples 
#' #bigrams for standard aminoacids
#' create_ngrams(2, 1L:20)
#' #bigrams for standard aminoacids with positions, 10 amino acid long sequence, so 
#' #only 9 bigrams can be located in sequence
#' create_ngrams(2, 1L:20, 9)
#' #bigrams for DNA with positions, 10 nucleotide long sequence, distance 1, so only 
#' #8 bigrams in sequence
#' #paste0 adds information about distance at the end of n-gram
#' paste0(create_ngrams(2, 1L:4, 8), "_0")


create_ngrams <- function(n, u, possible_grams = NULL) {
  grid_list <- lapply(1L:n, function(i) u)
  res <- apply(expand.grid(grid_list), 1, function(x)
    paste(x, collapse = "."))
  if (!is.null(possible_grams))
    res <- as.vector(sapply(res, function(i) paste(1L:possible_grams, i, sep = "_")))
  res
}


#' Extract n-grams from sequence
#'
#' Extracts vector of n-grams present in sequence(s).
#'
#' @inheritParams count_ngrams
#' @details A format of \code{d} vector is discussed in Details of 
#' \code{\link{count_ngrams}}.
#' @return A \code{character} matrix of n-grams, where every row corresponds to a
#' different sequence.
#' @export
#' @examples 
#' #trigrams from multiple sequences
#' seqs <- matrix(sample(1L:4, 600, replace = TRUE), ncol = 50)
#' seq2ngrams(seqs, 3, 1L:4)

seq2ngrams <- function(seq, n, u, d = 0, pos = FALSE) {
  if (!(is.matrix(seq) || is.vector(seq)))
    stop("'seq' must be vector or matrix.")
  
  #if sequence is not a matrix (single sequence), convert it to matrix with 1 row
  if (class(seq) != "matrix")
    seq <- matrix(seq, nrow = 1)
  
  #length of sequence
  len_seq <- ncol(seq)
  #number of sequences
  n_seqs <- nrow(seq)
  
  #look for n-gram indices for d
  ngram_ind <- get_ngrams_ind(len_seq, n, d)
  
  max_grams <- calc_max_grams(len_seq, n, ngram_ind)
  
  #extract n-grams from sequene
  res <- t(vapply(1L:n_seqs, function(i) {
    grams <- seq2ngrams_helper(seq[i, ], ind = ngram_ind, max_grams)
    paste(grams, paste0(attr(ngram_ind, "d"), collapse = "."), 
          sep = "_")
  }, rep("a", max_grams)))
  if (max_grams == 1)
    res <- t(res)
  
  #add position information if requested
  if(pos)
    res <- do.call(cbind, lapply(1L:ncol(res), function(pos_id)
      paste0(pos_id, "_", res[, pos_id])))
  
  res
}


#' Gap n-grams
#'
#' Introduces gaps in the n-grams.
#'
#' @inheritParams position_ngrams
#' @return A \code{character} vector of (n-1)-grams with introduced gaps.
#' @details A single element of the input n-gram at a time will be replaced 
#' by a gap. For example, introducing gaps in n-gram \code{2_1.1.2_0.1} 
#' will results in three n-grams: \code{3_1.2_1} (where the \code{2_1_0} unigram 
#' was replaced by a gap), \code{2_1.2_2} and \code{2_1.1_0}.
#' @export
#' @examples 
#' gap_ngrams(c("2_1.1.2_0.1", "3_1.1.2_0.0", "3_2.2.2_0.0"))

gap_ngrams <- function(ngrams) {
  #check if unigrams are there 
  
  #no need to validate n-grams, decode does it for us
  decoded <- decode_ngrams(ngrams)
  
  df <- ngrams2df(ngrams)
  
  #splitted ngrams
  sn_grams <- strsplit(df[, "ngram"], ".", fixed = TRUE)
  distances <- strsplit(df[, "distance"], ".", fixed = TRUE)
  
  unlist(lapply(1L:length(ngrams), gap_single_ngram, sn_grams, distances, df, decoded))
}


gap_single_ngram <- function(ngram_id, sn_grams, distances, df, decoded) {
  sn_gram <- sn_grams[[ngram_id]]
  distance <- as.numeric(distances[[ngram_id]])
  pos_start <- df[ngram_id, "position"]
  #single decoded
  s_decoded <- strsplit(decoded[ngram_id], "")[[1]]
  
  ids <- which(s_decoded != "_")
  res <- unlist(lapply(ids, function(id) {
    s_decoded[id] <- "_"
    code_ngrams(paste0(s_decoded, collapse = ""))
  }))
  
  #positions of gapped n-grams
  #first position is increased, because first element of the n-gram becomes a gap
  res_positions <- c(pos_start + ids[2] - ids[1], rep(pos_start, length(res) - 1))
  
  unlist(lapply(1L:length(ids), function(i)
    paste0(res_positions[i], "_", res[i])))
}