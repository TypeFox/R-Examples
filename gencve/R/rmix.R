rmix <- function(n=100) {
  nHalf <- ceiling(n/2)
  n <- 2*nHalf
  y <- as.factor(c(rep("green", nHalf),rep("red", nHalf)))
  gm <- matrix(c(-0.253433,1.74148,0.266693,0.371234,2.09647,1.23336,-0.0612727,
                 -0.208679,2.70354,0.596828,2.37721,-1.18641,1.05691,-0.683894,
                 0.578884,-0.0683458,0.624252,0.598738,1.67335,-0.289316),
               byrow=TRUE, ncol=2)
  rm <- matrix(c(1.19937,0.248409,-0.302561,0.945419,0.0572723,2.41973,1.32932,
                 0.819226,-0.0793842,1.6138,3.50793,1.05299,1.61392,0.671738,
                 1.00754,1.36831,-0.454621,1.08607,-1.79802,1.92978),
               byrow=TRUE, ncol=2)
  k <- ceiling(10*runif(n))
  greenMean <- gm[k[1:nHalf],]
  redMean <- rm[k[1:nHalf],]
  mns <- rbind(greenMean, redMean)
  m<-as.data.frame(t(apply(mns, 1, function(x)
    c(rnorm(1, mean=x[1], sd=1/sqrt(5)), rnorm(1, mean=x[2], sd=1/sqrt(5))))))
  names(m) <- c("x1", "x2")
  cbind(m, y=y)
}
