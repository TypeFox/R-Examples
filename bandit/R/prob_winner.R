prob_winner <- function(post){
  k = ncol(post)
  w = table(factor(max.col(post), levels = 1:k))
  return(w/sum(w))
}
