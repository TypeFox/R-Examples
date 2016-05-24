#' Decode n-grams
#'
#' Transforms a vector of n-grams into a human-friendly form.
#'
#' @param ngrams a \code{character} vector of n-grams.
#' @return a \code{character} vector of length equal to the number of n-grams. 
#' @note Decoded n-grams lose the position information.
#' @export
#' @seealso
#' Validate n-gram structure: \code{\link{is_ngram}}.
#' 
#' Inverse function: \code{\link{code_ngrams}}.
#' @examples
#' decode_ngrams(c("2_1.1.2_0.1", "3_1.1.2_2.0", "3_2.2.2_0.0"))
decode_ngrams <- function(ngrams) {
  validated_ngram <- sapply(ngrams, is_ngram)
  if(!all(validated_ngram))
    stop("Improper n-grams: paste(names(which(!validated_ngram)), collapse = ", ").")
  
  sngrams <- strsplit(ngrams, "_")
  vapply(sngrams, decode_single_ngrams, "a")
}

decode_single_ngrams <- function(splitted_ngram) {
  pos_inf <- ifelse(length(splitted_ngram) == 3, TRUE, FALSE)
  seq <- strsplit(splitted_ngram[1 + pos_inf], ".", fixed = TRUE)[[1]]
  if(length(seq) > 1) {
    dists <- strsplit(splitted_ngram[2 + pos_inf], ".", fixed = TRUE)[[1]]
    #distances in bar form
    bar_dists <- vapply(dists, function(i) 
      paste(rep("_", i), collapse = ""), "a")
    paste(c(vapply(1L:(length(seq) - 1), function(i)
      c(seq[i], bar_dists[i]), c("a", "a")), seq[length(seq)]), collapse = "")
  } else {
    seq
  }
}


#' Code n-grams
#'
#' Code human-friendly representation of n-grams into a biogram format.
#'
#' @param decoded_ngrams a \code{character} vector of decoded n-grams.
#' @return a \code{character} vector of n-grams. 
#' @export
#' @seealso Inverse function: \code{\link{decode_ngrams}}.
#' @examples
#' code_ngrams(c("11_2", "1__12", "222"))
#' code_ngrams(c("aaa_b", "d__aa", "abd"))

code_ngrams <- function(decoded_ngrams) {
  #ad some checks for decoded n-grams (allow only letters, numbers and underscores)
  as.vector(sapply(decoded_ngrams, function(decoded_ngram) {
    sn <- strsplit(decoded_ngram, "")[[1]]
    
    #get indices of elements
    id_elements <- which(sn != "_")
    
    #calculate distances between elements
    dists <- sapply(2L:length(id_elements), function(id) 
      id_elements[id] - id_elements[id - 1] - 1)
    paste0(paste(sn[sn != "_"], collapse = "."), "_", 
           paste(dists, collapse = "."))
  }))
}


#' n-grams to data frame
#'
#' Tranforms a vector of n-grams into a data frame.
#' 
#' @inheritParams decode_ngrams
#' @return a \code{data.frame} with 2 (in case of n-grams without known position) or
#' three columns (n-grams with position information). 
#' @export
#' @seealso
#' Decode n-grams: \code{\link{decode_ngrams}}.
#' @examples
#' ngrams2df(c("2_1.1.2_0.0", "3_1.1.2_0.0", "3_2.2.2_0.0", "2_1.1_0"))

ngrams2df <- function(ngrams) {
  sngrams <- strsplit(ngrams, "_")
  df <- data.frame(do.call(rbind, strsplit(ngrams, "_")), stringsAsFactors = FALSE)
  if(ncol(df) == 2) {
    colnames(df) <- c("ngram", "distance")
  } else {
    colnames(df) <- c("position", "ngram", "distance")
    df[, "position"] <- as.numeric(df[, "position"])
  }
  df
}

         
