un <- unname





rtpois <- function(n, lambda, atLeast = -1){
  x <- rpois(n, lambda)
  if(all(x >= atLeast)) return(x)
  while(any(x < atLeast)){
    x[x < atLeast] <- rpois(sum(x < atLeast), lambda)
  }
  x
}
#rtpois(10, 1, 1)




rtbinom <- function(n, size, prob, atLeast = 0, atMost = size){
  x <- rbinom(n, size, prob)
  if(all(x >= atLeast & x <= atMost)) return(x)
  while(any(x < atLeast | x > atMost)){
    x[x < atLeast | x > atMost] <- 
      rbinom(sum(x < atLeast | x > atMost), size, prob)
  }
  x  
}
#rtbinom(10, 5, .5)
#rtbinom(10, 5, .5, 1, 4)




