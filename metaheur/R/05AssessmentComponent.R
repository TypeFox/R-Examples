
getfirstassessment <- function(numberofcurrentbest, grid, predictors, nholdout){
  dat <- grid@data[[numberofcurrentbest]]
  dat1 <- data.frame(y=dat@y, x=dat@x)
  s2 <- preprocomb::getprogrammaticprediction(dat1, predictors, nholdout)[1]
  s2 <- apply(s2, 2, mean)
}

getcandidatedata <- function(grid, candidate_new, returntype){

  res <- logical()
  for (i in 1:nrow(grid@grid))
  {
    a1 <- unname(unlist(grid@grid[i,]))
    b1 <- as.character(unlist(candidate_new[1,]))
    res[i] <- (identical(a1,b1))
  }

  temp <- which(res==TRUE)
  dat <- grid@data[[temp]]
  dat1 <- data.frame(y=dat@y, x=dat@x)

  list(dat1, temp)

}

getconsequentassessment <- function(dat1, predictors, nholdout) {
  r <- preprocomb::getprogrammaticprediction(dat1, predictors, nholdout)[1]
  r <- apply(r, 2, mean)
}
