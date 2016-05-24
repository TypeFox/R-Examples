#' Add 1-grams
#'
#' Builds (n+1)-grams from n-grams.
#'
#' @inheritParams position_ngrams
#' @return a vector of n-grams (where \code{n} is equal to the \code{n} of the input plus one) 
#' with position information.
#' @details n-grams are built by pasting existing n-grams with unigrams extracted 
#' from them.
#' @export
#' @seealso
#' Function used by \code{add_1grams} to extract unigrams: \code{\link{position_ngrams}}.
#' @note All n-grams must have the same length (\code{n}).
#' @examples
#' add_1grams(c("1_1_0", "2_1_0", "5_1_0", "7_1_0", "4_2_0", 
#' "5_2_0", "7_2_0", "8_5_0"))
#' add_1grams(c("1_2.3.4_0.0", "4_1.1.1_0.0"))

add_1grams <- function(ngrams) {
  validated_ngram <- sapply(ngrams, is_ngram)
  if(!all(validated_ngram))
    stop("Improper n-grams: ", paste(names(which(!validated_ngram)), 
                                     collapse = ", "))
  
  splitted_ngrams <- strsplit(ngrams, "_")
  if(unique(sapply(splitted_ngrams, length)) != 3)
    stop("Use only n-grams with position information.")
  n <- sapply(strsplit(sapply(splitted_ngrams, function(ngram) 
    ngram[2]), ".", fixed = TRUE), length)
  if(length(unique(n)) != 1)
    stop("Unequal n-gram size. Use n-grams with the same size (n).")
  
  positioned_ngrams <- position_ngrams(ngrams, df = TRUE, 
                                       unigrams_output = FALSE)
  
  #calculate end of the n-gram
  ngrams_ends <- if(n[1] == 1) {
    #if unigram, end is equal to the start of n-gram
    positioned_ngrams[["position"]]
  } else {
    #if n > 1, end is equal to the start of n-gram plus n plus sum of distances
    positioned_ngrams[["position"]] + 
      apply(positioned_ngrams, 1, function(ngram) {
        sn <- strsplit(as.character(ngram[1]), "_")[[1]]
        ngram_end <- sum(as.numeric(strsplit(sn[2], ".", fixed = TRUE)[[1]])) +
          length(strsplit(sn[2], ".", fixed = TRUE)[[1]])
        ngram_end
      })
  }
  
  position_data <- cbind(positioned_ngrams, ngrams_ends)
  #ngram, start position of n-gram, end position of n-gram
  colnames(position_data) <- c("ngram", "pstart", "pend")
  
  #create table for unigrams that will be added to existing n-grams
  positioned_ugrams <- position_ngrams(ngrams, df = TRUE, unigrams_output = TRUE)
  positioned_ugrams <- positioned_ugrams[!duplicated(positioned_ugrams), ]
  positioned_ugrams[["ngram"]] <- as.character(positioned_ugrams[["ngram"]])
  
  u_positions <- unique(positioned_ugrams[["position"]])
  # n-grams to which we cannot add anything on the right side
  # position_data[["pend"]] < max(u_positions)
  # n-grams to which we cannot add anything on the left side
  # position_data[["pstart"]] > min(u_positions)
  #add unigrams on the right side
  
  res_right <- add_unigrams_right(position_data[position_data[["pend"]] < max(u_positions), ], 
                                  positioned_ugrams, n = n[1])
  res_left <- add_unigrams_left(position_data[position_data[["pend"]] > min(u_positions), ], 
                                positioned_ugrams, n = n[1])
  #work with long n-grams (positions bigger than 9)
  res <- c(res_left, res_right)
  
  names(res) <- NULL
  res
}


#TODO - add_unigrams_right and add_unigrams_left should be collapsed into a single function 
add_unigrams_right <- function(position_data, positioned_ugrams, n) 
  unlist(apply(position_data, 1, function(single_row) {
    chosen_ngram <- strsplit(single_row["ngram"], "_")[[1]]
    single_position <- as.numeric(single_row["pend"])
    #u_grams that may be pasted
    other_ugrams <- positioned_ugrams[positioned_ugrams[["position"]] > single_position, ]
    #position in other_ugrams is now distance between single_position and their position
    other_ugrams[["position"]] <- other_ugrams[["position"]] - single_position - 1
    #remaining distance - cut redundant 0 when n = 1
    remain_distance <- ifelse(n == 1, "", paste0(chosen_ngram[[2]], "."))
    apply(other_ugrams, 1, function(other_ugram)
      paste0(as.numeric(single_row["pstart"]), "_", #position 
             chosen_ngram[[1]], ".", substr(other_ugram[1], 1, 1), #ngram 
             "_", remain_distance, as.numeric(other_ugram[2])))})) #distance



add_unigrams_left <- function(position_data, positioned_ugrams, n) 
  unlist(apply(position_data, 1, function(single_row) {
    chosen_ngram <- strsplit(single_row["ngram"], "_")[[1]]
    single_position <- as.numeric(single_row["pstart"])
    #u_grams that may be pasted
    other_ugrams <- positioned_ugrams[positioned_ugrams[["position"]] < single_position, ]
    #position in other_ugrams is now distance between single_position and their position
    other_ugrams[["position"]] <- single_position - other_ugrams[["position"]] - 1
    #remaining distance - cut redundant 0 when n = 1
    remain_distance <- ifelse(n == 1, "", paste0(".", chosen_ngram[[2]]))
    apply(other_ugrams, 1, function(other_ugram)
      paste0(as.numeric(single_row["pstart"]) - as.numeric(other_ugram[2]) - 1, "_", #position 
             substr(other_ugram[1], 1, 1), ".", chosen_ngram[[1]], #ngram 
             "_", as.numeric(other_ugram[2]), remain_distance)) #distance
  }))


