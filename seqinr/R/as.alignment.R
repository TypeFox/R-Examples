#
# Constructor for class alignment
#
as.alignment <- function(nb = NULL, nam = NULL, seq = NULL, com = NULL){
  ali <- list(nb = as.numeric(nb), nam = nam, seq = seq, com = com)
  class(ali) <- "alignment"
  return(ali)
}

