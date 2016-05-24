rsplitnorm <-
function(n, mode = 0, sd = 1, skew = 0, sd1 = NULL, sd2 = NULL) {
  n <- max(n,length(mode),length(sd),length(skew),length(sd1),length(sd2))
  if(length(mode)<n)
    mode[1:n]<-mode
  if(length(sd)<n)
    sd[1:n]<-sd
  if(length(skew)<n)
    skew[1:n]<-skew
  if(length(sd1)<n)
    sd1[1:n]<-sd1
  if(length(sd2)<n)
    sd2[1:n]<-sd2
  var0 <- sd^2
  if (!is.null(sd1)) 
    var1 <- sd1^2
  if (!is.null(sd2)) 
    var2 <- sd2^2
  if (sum(is.null(sd1), is.null(sd2)) == 1) 
    stop("give either sd of both sides (sd1 and sd2) or neither (sd only)")
  if (is.null(sd1) & is.null(sd2)) {
    var1 <- var0/(1 + skew)
    var2 <- var0/(1 - skew)
    sd1 <- sqrt(var1)
    sd2 <- sqrt(var2)
  }
  if (any(findInterval(skew, c(-1,1), rightmost.closed=TRUE)!=1) )
    stop("skew must be between -1 and 1")
  u <- runif(n)
  f <- qsplitnorm(u, mode = mode, sd = sd, skew = skew, sd1 = sd1, sd2 = sd2)
  f
}
