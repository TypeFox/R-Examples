#' Get indices of n-grams
#'
#' Computes list of n-gram elements positions in sequence.
#'
#' @param len_seq \code{integer} value describing sequence's length.
#' @inheritParams count_ngrams
#' @return A list with number of elements equal to \code{n}. Every element is a 
#' vector containing locations of given n-gram letter. For example, first element of
#' list contain indices of first letter of all n-grams. The attribute \code{d}
#' of output contains distances between letter used to compute locations 
#' (see Details).
#' @details A format of \code{d} vector is discussed in Details of 
#' \code{\link{count_ngrams}}.
#' @export
#' @examples 
#' #positions trigrams in sequence of length 10
#' get_ngrams_ind(10, 9, 0)

get_ngrams_ind <- function(len_seq, n, d) {
  #n - size of gram (i.e. 2-grams: "AA", "AC", ...)
  #d - distance between two consecutive letter (a vector of distances)
  
  #calculate indices of n-grams elements
  ind <- lapply(1L:n, function(i) 
    (1 + i - 1):(len_seq - n + i))
  
  if(length(d) != 1 && length(d) != n - 1)
    stop("Length of d must be 1 or n - 1")
  
  if(n > 1) {
    #if distance vector is too short, recycle it
    if(length(d) == 1 && n > 2)
      d <- rep(d, n - 1)
    
    if(sum(d) > 0) {
      ind[-1] <- lapply(1L:length(d), function(i)
        ind[[i + 1]] + sum(d[1L:i]))
      not_taken <- ind[[1]][(length(ind[[1]]) - sum(d) + 1):length(ind[[1]])]
      ind <- lapply(ind, function(i) i[-not_taken])
    }
    
    attr(ind, "d") <- d
  } else {
    #distance is a nonsense for unigrams
    attr(ind, "d") <- 0
  }
  
  ind
}


#' Count total number of n-grams
#'
#' Computes total number of n-grams that can be extracted from sequences taking 
#' into account their length (even or uneven).
#'
#' @inheritParams count_ngrams
#' @return A number of n-grams. 
#' @details A format of \code{d} vector is discussed in Details of 
#' \code{\link{count_ngrams}}.
#' @export
#' @examples 
#' seqs <- matrix(sample(1L:4, 600, replace = TRUE), ncol = 50)
#' #make several sequences shorter by replacing them partially with NA
#' seqs[8L:11, 46L:50] <- NA
#' seqs[1L, 31L:50] <- NA
#' count_total(seqs, 3, c(1, 0))

count_total <- function(seq, n, d) {
  seq_lengths <- ncol(seq) - apply(seq, 1, function(i) sum(is.na(i)))
  
  #unique lengths of sequences
  tablel <- data.frame(table(seq_lengths))
  
  #calculate number of posible n-grams for each unique length
  tablel <- cbind(tablel, totals = vapply(as.numeric(as.character(tablel[["seq_lengths"]])), function(single_length) {
    ind <- get_ngrams_ind(single_length, n, d)
    calc_max_grams(single_length, n, ind)
  }, 0))
  
  #multiply number of possible n-grams by time the sequence of this length occur in data
  sum(apply(tablel[, -1], 1, function(i)
    i["Freq"] * i["totals"]))
}


#helper function calculating maximum number of n-grams possible. Throws an
#error if there is no possibility of extracting n-gram from a sequence 
#(when result is negative)

calc_max_grams <- function(len_seq, n, ngram_ind){
  #use attr(ngram_ind, "d") instead of d because of distance recycling
  max_grams <- len_seq - n - sum(attr(ngram_ind, "d")) + 1
  if (max_grams < 1)
    stop("n-gram too long.")
  max_grams
}