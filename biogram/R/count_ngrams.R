#' Count n-grams in sequences
#'
#' Counts all n-grams or position-specific n-grams present in the input sequence(s).
#'
#' @param seq a vector or matrix describing sequence(s). 
#' @param n \code{integer} size of n-gram.
#' @param u \code{integer}, \code{numeric} or \code{character} vector of all
#' possible unigrams.
#' @param d \code{integer} vector of distances between elements of n-gram (0 means 
#' consecutive elements). See Details.
#' @param pos \code{logical}, if \code{TRUE} position-specific n_grams are counted.
#' @param scale \code{logical}, if \code{TRUE} output data is normalized. Should be
#' used only for n-grams without position information. See \code{Details}.
#' @param threshold \code{integer}, if not equal to 0, data is binarized into
#' two groups (larger or equal to threshold vs. smaller than threshold).
#' 
#' @details A \code{distance} vector should be always \code{n} - 1 in length.
#' For example when \code{n} = 3, \code{d} = c(1,2) means A_A__A. For \code{n} = 4, 
#' \code{d} = c(2,0,1) means A__AA_A. If vector \code{d} has length 1, it is recycled to
#' length \code{n} - 1.
#' 
#' n-gram names follow a specific convention and have three parts for position-specific
#' n-grams and two parts otherwise. The parts are separated by \code{_}. The \code{.} symbol
#' is used to separate elements within a part. The general naming scheme is 
#' \code{POSITION_NGRAM_DISTANCE}. The optional \code{POSITION} part of the name indicates
#' the actual position of the n-gram in the sequence(s) and will be present 
#' only if \code{pos} = \code{TRUE}. This part is always a single integer. The \code{NGRAM}
#' part of the name is a sequence of elements in the n-gram. For example, \code{4.2.2}
#' indicates the n-gram 422 (e.g. TCC). The \code{DISTANCE} part of the name is a vector of
#' distance(s). For example, \code{0.0} indicates zero distances (continuous n-grams), while
#' \code{1.2} represents distances for the n-gram A_A__A.
#' 
#' Examples of n-gram names:
#' \itemize{
#' \item{46_4.4.4_0.1 : trigram 44_4 on position 46}
#' \item{12_2.1_2     : bigram 2__1 on position 12}
#' \item{8_1.1.1_0.0  : continuous trigram 111 on position 8}
#' \item{1.1.1_0.0    : continuous trigram 111 without position information}
#' }
#' 
#' @return a \code{\link[slam]{simple_triplet_matrix}} where columns represent
#' n-grams and rows sequences. See \code{Details} for specifics of the naming convention.
#' 
#' @note By default, the counted n-gram data is stored in a memory-saving format.
#' To convert an object to a 'classical' matrix use the \code{\link[base]{as.matrix}}
#' function. See examples for further information.
#' @export
#' @seealso 
#' Create vector of possible n-grams: \code{\link{create_ngrams}}.
#' 
#' Extract n-grams from sequence(s): \code{\link{seq2ngrams}}.
#' 
#' Get indices of n-grams: \code{\link{get_ngrams_ind}}.
#' 
#' Count n-grams for multiple values of n: \code{\link{count_multigrams}}.
#' 
#' Count only specified n-grams: \code{\link{count_specified}}.
#' @examples 
#' #count trigrams without position information for nucleotides
#' count_ngrams(sample(1L:4, 50, replace = TRUE), 3, 1L:4, pos = FALSE)
#' #count position-specific trigrams from multiple nucleotide sequences
#' seqs <- matrix(sample(1L:4, 600, replace = TRUE), ncol = 50)
#' ngrams <- count_ngrams(seqs, 3, 1L:4, pos = TRUE)
#' #output results of the n-gram counting to screen
#' as.matrix(ngrams)

count_ngrams <- function(seq, n, u, d = 0, pos = FALSE, 
                         scale = FALSE, threshold = 0) {
  
  #if sequence is not a matrix (single sequence), convert it to matrix with 1 row
  if (class(seq) != "matrix")
    seq <- matrix(seq, nrow = 1)
  
  if (scale && pos)
    stop("Cannot scale positioned n-grams (scaling a sparse matrix).")
  
  u_seq <- unique(as.vector(seq))
  if (!(all(u_seq %in% u)))
    warning(paste0("'seq' contains following unigrams not present in 'u' parameter:\n",
                   paste(u_seq[!(u_seq %in% u)], collapse = ", ")))
  
  #length of sequence
  len_seq <- ncol(seq)
  #number of sequences
  n_seqs <- nrow(seq)
  
  #create list of n-grams for n
  possib_ngrams <- create_ngrams(n, u)  
  
  #look for n-gram indices for d
  ngram_ind <- get_ngrams_ind(len_seq, n, d)
  max_grams <- calc_max_grams(len_seq, n, ngram_ind)
  
  
  #extract n-grams from sequence
  grams <- vapply(1L:n_seqs, function(i)
    seq2ngrams_helper(seq[i, ], ind = ngram_ind, max_grams), rep("a", max_grams))
  
  if (pos) {    
    #get positioned possible n-grams
    pos_possib_ngrams <- create_ngrams(n, u, max_grams)
    
    res <- do.call(cbind, lapply(possib_ngrams, function(current_ngram)
      as.simple_triplet_matrix(t(vapply(1L:n_seqs, function(current_sequence)
        grams[, current_sequence] == current_ngram, rep(0, max_grams))))))
    
    colnames(res) <- pos_possib_ngrams
    
  } else {
    
    res <- do.call(cbind, lapply(possib_ngrams, function(current_ngram)
      as.simple_triplet_matrix(vapply(1L:n_seqs, function(current_sequence)
        sum(grams[, current_sequence] == current_ngram, na.rm = TRUE), 0))))
    colnames(res) <- possib_ngrams
  }
  
  if (threshold > 0) {
    ind_sums <- rowSums(res)
    res <- res[ind_sums >= threshold, ]
  }
  
  if (scale)
    res <- res/max_grams
  
  colnames(res) <- paste(colnames(res), paste0(attr(ngram_ind, "d"), collapse = "."), 
                         sep = "_")
  res
}

count_ngrams_helper <- function(seq, feature_list, n, ind, pos) {
  #feature list(list of possible n-grams) is outside, because count_ngrams is meant to
  #be inside the loop
  #same for indices
  if (n > 1) {
    element_matrix <- do.call(cbind, lapply(ind, function(i) seq[i]))
    grams <- apply(element_matrix, 1, function(x) 
      paste(x, collapse="."))
  } else {
    grams <- seq
  }
  
  if (pos) {
    grams <- paste(1L:length(grams), grams, sep = "_")
  } 
  
  sapply(feature_list, function(i)
    sum(i == grams))
}

#ultrafast function for n-gram extraction
seq2ngrams_helper <- function(seq, ind, max_grams) {
  #get all consecutive n-grams from sequence
  #length(ind) > 1 - equivalent of n > 1 
  if (length(ind) > 1) {
    #element_matrix contains elements of n-gram in matrix structure
    element_matrix <- vapply(ind, function(i) seq[i], rep(seq[1], max_grams))
    
    #rare situation, only one n-gram in sequence
    if(max_grams == 1) {
      grams <- paste(element_matrix, collapse=".")
    } else {
      grams <- vapply(1L:nrow(element_matrix), function(ith_row)
        paste(element_matrix[ith_row, ], collapse="."), "a")
    }
    
    
  } else {
    #in case of unigrams take sequence - output must be character
    grams <- as.character(seq)
  }
  
  #all n-grams from seuqence - character vector
  grams
}

