is.gap <- function(seq) {

  #characters allowed to be a gap
  gap <- c("-", ".", "~")

  #NA not considered as gap
  return(seq %in% gap)
}