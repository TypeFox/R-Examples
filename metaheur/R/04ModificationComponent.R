
getnewcandidate <- function(grid, taboo, taboolistlength, uniquepreprocessors, copyofcurrentbest){

set.seed(as.numeric(Sys.time()))
repeat{
  candidate_phase <- sample(1:ncol(grid@grid),1)
  candidate_preprocessor <- sample(1:length(uniquepreprocessors[[candidate_phase]]), 1)

  candidate_preprocessor <- unlist(grid@grid[candidate_preprocessor, candidate_phase])
  candidate_new <- copyofcurrentbest
  candidate_new[,candidate_phase] <- candidate_preprocessor

  condition1 <- lapply(utils::tail(taboo,taboolistlength), function(x) identical(unname(unlist(candidate_new)),unname(unlist(x))))
  condition2 <- all(unlist(condition1)==FALSE)

  if(condition2==TRUE) {
    break}
}

  return(candidate_new)
}





