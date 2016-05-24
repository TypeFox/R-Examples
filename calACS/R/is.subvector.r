# Internal function:
# Check if the input vector subvec is a subvector of the other input vector vec
#
# param subvec the input subvector
# param vec the input vector that may contain the subvector

is.subvector <- function(subvec, vec){
  #subvec <- c('a')
  #vec <- c('a','b')

  substr<-paste(subvec, collapse = ', ')
  fullstr<-paste(vec, collapse = ', ')

  grepl(pattern = substr,x = fullstr)
}
