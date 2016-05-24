updateLambdas <- function(N.k, n.k, b.k, I.t){
  lambda.k <- rep(0, length(N.k))
  for(i in seq_along(N.k)){
    # i <- 1
    if(N.k[i]!=0){
      potentials <- N.k[i]-n.k[i]
      if(any(potentials<0)) stop('Impossible n.k value')
      lambda <- b.k[i] * potentials * I.t
      
      lambda.k[i] <- ifelse(is.na(lambda), 0, lambda)
    }
  }
  return(lambda.k)
}
## Testing:
# example(makeRdsSample)
# lambda.k <- chords:::updateLambdas(
#   N.k =true.Nks ,
#   n.k = rep(0, length.out = length(true.Nks)), 
#   b.k =true.log.bks+log(100) ,
#   I.t=1)
# plot(lambda.k, type='h')



### Generate sample:
## Needs to return an rds.object including:
# I.t
# degree.in
# degree.out
# rds.sample$NS1
# rds.sample$interviewDt

makeRdsSample <- function (N.k, b.k, sample.length) {
  n.ks <- rep(0, length(N.k))
  event.time <- 0
  I.t <- 1
  rds.object.simulated <- rdsObjectConstructor(
    rds.sample = data.frame(NS1=rep(NA, sample.length), interviewDt=rep(NA, sample.length)),
    I.t = rep(NA, sample.length),
    degree.in = rep(0, sample.length),
    degree.out = rep(0, sample.length))
  for(period in seq_len(sample.length)){
    lambda.k <- updateLambdas(N.k, n.ks, b.k, I.t)
    ## FIXME: floating point errors when summing rate?
    time.to.event <- rexp(1, sum(lambda.k))
    event.time <- event.time + time.to.event
    rds.object.simulated$rds.sample$interviewDt[period] <- event.time
    
    event.type <- rmultinom(1, 1, prob = lambda.k)
    in.degree <-  which(as.logical(event.type))
    rds.object.simulated$degree.in[period] <- in.degree
    n.ks <- n.ks + event.type  
    
    rds.object.simulated$I.t[period] <- I.t
    I.t <- I.t+1
    
    rds.object.simulated$rds.sample$NS1[period] <- in.degree
  }
  return(rds.object.simulated)
}
## Testing:
# rds.simulated.object <- makeRdsSample(N.k =N.k , sample.length = 900L)
# plot(rds.simulated.object$rds.object$rds.sample$interviewDt)


compareNkEstimate <- function(Nk1, Nk2){
#     object1 <- nk.estimates.2
#     object2 <- nk.estimates
  
  if(length(Nk1) > length(Nk2)) {
    .temp <- Nk2
    Nk2 <- Nk1
    Nk1 <- .temp
  }
  len1 <- length(Nk1)
  len2 <- length(Nk2)
  
  # Padding if lengths differ:
  if(len1!= len2){
    Nk1[(length(Nk1)+1):len2] <- 0
  }

  y.lim <- max(c(Nk2,Nk1))
  plot(Nk2, type='h', lwd=2, main='N_k',ylim=c(0,y.lim))
  points(Nk1, col='red', type='h')
  
  return(list(max1=max(Nk2/Nk1, na.rm=TRUE),
              min=min(Nk2/Nk1, na.rm=TRUE)))
}