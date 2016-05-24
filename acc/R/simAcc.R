#' @export
#' @importFrom utils head tail 
#' @importFrom graphics par axis title plot rect legend
#' @importFrom mhsmm simulate.hmmspec hmmspec dnorm.hsmm rnorm.hsmm
#' @importFrom zoo rollmean rollsum rollmedian
#' @importFrom PhysicalActivity dataCollapser
simAcc <- function(minutes,mvpaLevel,seedset=1234){
  
  randomTime <- seq(ISOdate(2015,1,1),ISOdate(2020,1,1),"min")
  J <- 3; initial <- rep(1/J, J)
  P <- matrix(rep(NA,9),byrow='TRUE',nrow=J)
  
  if(mvpaLevel=='low'){
    P <- matrix(c(0.95, 0.04, 0.01, 
                  0.09, 0.9, 0.01, 
                  0.1, 0.2, 0.7), byrow='TRUE',nrow = J)
  }
  
  if(mvpaLevel=='moderate'){
    P <- matrix(c(0.95, 0.04, 0.01, 
                  0.09, 0.8, 0.11, 
                  0.1, 0.1, 0.8), byrow='TRUE',nrow = J)
  }
  
  if(mvpaLevel=='high'){
    P <- matrix(c(0.95, 0.04, 0.01, 
                  0.09, 0.7, 0.21, 
                  0.1, 0.1, 0.8), byrow='TRUE',nrow = J)
  }
  
  b <- list(mu = c(0, 30, 2500), sigma = c(0, 30, 1000))
  model <- hmmspec(init = initial, trans = P, parms.emission = b,dens.emission = dnorm.hsmm)
  train <- simulate.hmmspec(model, nsim = (minutes), seed = seedset, rand.emission = rnorm.hsmm)
  counts <- data.frame(TimeStamp = randomTime[1:minutes], counts = round(train$x,0))
  counts
}