#' Position n-grams
#'
#' Tranforms a vector of positioned n-grams into a list of positions filled with n-grams 
#' that start on them.
#'
#' @param ngrams a vector of positioned n-grams (as created by \code{\link{count_ngrams}}).
#' @param df logical, if \code{TRUE} returns a data frame, if \code{FALSE} returns a list.
#' @param unigrams_output logical, if \code{TRUE} extracts unigrams from the data and returns
#' information about their position. 
#' @return if \code{df} is \code{FALSE}, returns a list of length equal to the number of unique 
#' n-gram starts present in n-grams. Each element of the list contains n-grams that start on 
#' this position. If \code{df} is \code{FALSE}, returns a data frame where first column contains 
#' n-grams and the second column represent their start positions. 
#' @export
#' @seealso
#' Transform n-gram name to human-friendly form: \code{\link{decode_ngrams}}.
#' 
#' Validate n-gram structure: \code{\link{is_ngram}}.
#' @examples
#' #position data in the list format
#' position_ngrams(c("2_1.1.2_0.1", "3_1.1.2_0.0", "3_2.2.2_0.0"))
#' #position data in the data frame format
#' position_ngrams(c("2_1.1.2_0.1", "3_1.1.2_0.0", "3_2.2.2_0.0"), df = TRUE)

position_ngrams <- function(ngrams, df = FALSE, unigrams_output = TRUE) {
  validated_ngram <- sapply(ngrams, is_ngram)
  if(!all(validated_ngram))
    stop("Improper n-grams: ", paste(names(which(!validated_ngram)), collapse = ", "))
  
  sngrams <- strsplit(ngrams, "_")
  #check if there is information about position
  if(length(sngrams[[1]]) != 3)
    stop("n-grams do not have position information.")
  
  #table of positions
  pos_table <- if(unigrams_output) {
    do.call(rbind, lapply(sngrams, function(single_ngram) {
      unigrams <- strsplit(single_ngram[[2]], ".", fixed = TRUE)[[1]]
      dists <- strsplit(single_ngram[[3]], ".", fixed = TRUE)[[1]]
      #positions of unigrams
      uni_positions <- as.numeric(single_ngram[[1]])
      #for loop suitable only for n-grams bigger than unigrams
      if(nchar(single_ngram[[2]]) > 1)
        for (next_unigram in as.numeric(dists))
          uni_positions[length(uni_positions) + 1] <- next_unigram + uni_positions[length(uni_positions)] + 1
      data.frame(ngrams = paste0(unigrams, "_0"), pos = uni_positions)
    }))
  } else {
    pos_table <- do.call(rbind, lapply(sngrams, function(single_ngram) {
      ngram_name <- paste0(single_ngram[2], "_", single_ngram[3])
      #positions of unigrams
      ngram_position <- as.numeric(single_ngram[[1]])
      data.frame(ngrams = ngram_name, pos = ngram_position)
    }))
  }

  if(df) {
    res <- pos_table
    colnames(res) <- c("ngram", "position")
    #duplicates shouldn't be removed
    #res <- res[!duplicated(res), ]
    res <- res[order(res[["ngram"]]), ]
    res <- res[order(res[["position"]]), ]
    rownames(res) <- NULL
  } else {
    
    res <- lapply(sort(unique(pos_table[["pos"]])), function(unique_pos) 
      sort(pos_table[pos_table[["pos"]] == unique_pos, "ngrams"]))
    names(res) <- sort(unique(pos_table[["pos"]]))
  }
  
  res
}


