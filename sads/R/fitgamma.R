fitgamma <- function(x, trunc, start.value, ...){
  dots <- list(...)
  if (!missing(trunc)){
    if (min(x)<=trunc) stop("truncation point should, be lower than the lowest data value")
}
  if(missing(start.value)){
      if(missing(trunc)){
          ka <- (mean(x)/sd(x))^2
          theta <- var(x)/mean(x)
          kahat <- function(k, dados){
              eq <- length(dados)*(log(k) - log(mean(dados)) - digamma(k)) + sum(log(dados))
          }
          ka <- uniroot(kahat, interval = c(min(theta, ka), max(theta, ka)), dados = x)$root
          theta <- mean(x)/ka
      }
      else{
          xh <- hist(x, plot=FALSE)
          xbr <- xh$breaks
          Eh <- matrix(ncol=2, nrow=length(xbr)-1)
          for(i in 1:(length(xbr)-1)){
              m1 <- matrix(c(1,1,-1,1), ncol=2)
              Eh[i,] <- solve(m1, xbr[i:(i+1)])
          }
          P <- xh$counts/sum(xh$counts)
          P[P==0] <- min(P[P>0])
          Y <- log(P[-length(P)])-log(P[-1]) -(log(Eh[-nrow(Eh),2])-log(Eh[-1,2]))
          X1 <- Eh[-1,1]-Eh[-nrow(Eh),1]
          X2 <- log(Eh[-nrow(Eh),1])-log(Eh[-1,1])
          st1 <- unname(coef(lm(Y~X1+X2-1)))
          ka <- st1[2]+1
          theta <- 1/st1[1]
      }
  }
  else{
      ka <- start.value[1]
      theta <-start.value[2]
  }
  if (missing(trunc)){
      LL <- function(shape, rate) -sum(dgamma(x, shape, rate, log = TRUE))
  }
  else {
      LL <- function(shape, rate) -sum(dtrunc("gamma", x = x, coef = list(shape = shape, rate = rate), trunc = trunc, log = TRUE))
  }  
  result <- do.call("mle2", c(list(LL, start = list(shape = ka, rate = 1/theta), data = list(x = x)), dots))
  new("fitsad", result, sad="gamma", distr = distr.depr, trunc = ifelse(missing(trunc), NaN, trunc)) 
}
