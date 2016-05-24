#'The function takes in multiple vectors of any length, and returns the one with the longest length.
#'The tieBreaker variable controls if the first or the last of the longest vectors gets returned in case
#'there are multiple
#'
#'@param ... vectors of any length
#'@param tieBreaker decides if the first or the last longest vector gets returned if there are multiple longest vectors.
#'Can be either 'first' or 'last'. Default to 'last'.
#'@examples
#'longestVec(1:5, c('a','b'))
#'
#'@export

longestVec <- function(... , tieBreaker='last'){
  vectors<-list(...)
  output<-NULL

  if (tieBreaker=='first'){
    for(i in 1:length(vectors)){
      if(length(vectors[[i]])>(length(output))){
        output <- vectors[[i]]
      }
    }
  } else if (tieBreaker=='last'){
    for(i in 1:length(vectors)){
      if(length(vectors[[i]])>=length(output)){
        output <- vectors[[i]]
      }
    }
  }

  output
}

# longestVec(c(1,1,2,1),c(1,1,1))
