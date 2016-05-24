library(testthat)
library(mcmcse)

context("validResults")

fun_test <- function(x, method = "bm", size = "sqroot", g = NULL, 
              level = 0.95, large = FALSE)
{ 
  chain <- as.matrix(x)
  if(!is.matrix(x) && !is.data.frame(x))
    stop("'x' must be a matrix or data frame.")

  if (is.function(g)) 
  {
    chain <- apply(x, 1, g)

    if(is.vector(chain))
    {
      chain <- as.matrix(chain)
    }else
    {
      chain <- t(chain)
    }
  }

  ## Setting dimensions on the mcmc output. 
  n = dim(chain)[1]
  p = dim(chain)[2]


  ## Initializing b_n 
  if(size == "sqroot")
  {
    b = floor(sqrt(n))
  } 
  else if(size == "cuberoot") {
    b = floor(n^(1/3))
  }
  else {
    if (!is.numeric(size) || size <= 1 || size > n) 
        stop("'size' must be a numeric quantity not larger than n.")
    b = floor(size)
  }
  a = floor(n/b)
  m <- 0
  
  ## Overall means of the mcmc output
  mu.hat <- colMeans(chain)

  ## Batch means
  batch_means <- matrix(0, nrow = a, ncol = p)
  
  ## Setting matrix sizes to avoid dynamic memory 
  sig.mat = matrix(0, nrow = p, ncol = p)
  
  ## only does bm and obm
  if(method != "bm" && method != "bartlett" && method != "tukey")
  {
    stop("No such method available")
  }
  
  ## Batch Means
  if(method == "bm")
  {
    
  	for(i in 1:a)
  	{
  		batch_means[i,] = colMeans(chain[((i-1)*b+1): (i*b),])
  	}

  	for(j in 1:a)
  	{
  		sig.mat <- sig.mat + tcrossprod(batch_means[j,] - mu.hat)
  	}

  	sig.mat <- b*sig.mat/(a-1)
    m <- a - 1
  }
  
  ## Modified Bartlett Window

  if(method == "bartlett" || method == "tukey")
  {
    m <- n - b
  	lag <- function(x)
  	{ 
  		ifelse(method == "bartlett",  (1 - x/b), (1 + cos(pi * x/b))/2 ) 
  	}

  	## centering
  	for(i in 1:n)
  	{
  		chain[i,] <- chain[i,] - mu.hat
  	}

  	for(s in 1:(b-1))
  	{
  		for(j in 1: (n-s))
  		{
  			foo <- chain[j,]%*%t(chain[j+s,])
  			sig.mat <- sig.mat + lag(s)*(foo + t(foo))
  		}
  	}
    sig.mat <- (sig.mat + t(chain)%*%chain)/n
  }

  crit <- ifelse(large, qchisq(level, df = p)/n,
              exp(log(p) + log(m) - log(n) - log(m-p+1) + log(qf(level, p, m-p+1))) )
  dummy <- log(2) + (p/2)*log(pi*crit) - log(p) - lgamma(p/2)
  # log.det.sig <- sum(log(eigen(sig.mat, only.values = TRUE)$values))
  det.sig <- det(sig.mat)
  volume.sig <- exp(dummy + (1/2)*log(det.sig))
  pth.vol <- (volume.sig)^(1/p)
  return(list("cov" = sig.mat, "vol" = pth.vol, "est" = mu.hat, "nsim" = n, 
      "method" = method, "large" = large, "size" = b))
}

n <- 1e3
p <- 3
b <- floor(sqrt(n))
out <- matrix(rnorm(n*p), nrow = n, ncol = p)

mbatch <- fun_test(out)
mbart <- fun_test(out, method = "bartlett")
mtukey <- fun_test(out, method = "tukey")

test_that("test if functions return correct values", {

expect_equal(mcse.multi(out), mbatch)
expect_equal(mcse.multi(out, method = "bartlett"), mbart)
expect_equal(mcse.multi(out, method = "tukey"), mtukey)

})

