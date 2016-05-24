require(mvtnorm)
set.seed(0)
N <- 500


## Uncorrelated Normal ##

mu <- c(0, 0)
Sigma <- matrix(c(1/4, 0, 0, 1), 2, 2)

Uncorrelated.Normal <- rmvnorm(N, mean = mu, sigma = Sigma)



## Correlated Normal ##

mu <- c(0, 0)
Sigma <- matrix(c(1, 0.7, 0.7, 1), 2, 2)

Correlated.Normal <- rmvnorm(N, mean = mu, sigma = Sigma)



## Skewed ##

mu <- list()
Sigma <- list()

mu[[1]] <- c(0, 0)
Sigma[[1]] <- matrix(c(1, 0, 0, 1), 2, 2)

mu[[2]] <- c(0.5, 0.5)
Sigma[[2]] <- matrix(c(4/9, 0, 0, 4/9), 2, 2)

mu[[3]] <- c(13/12, 13/12)
Sigma[[3]] <- matrix(c(25/81, 0, 0, 25/81), 2, 2)

Skewed <- matrix(0.0, N, 2)
for(i in 1:N) {
	component <- sample(1:3, 1, prob = c(1/5, 1/5, 3/5))
	Skewed[i, ] <- rmvnorm(1, mean = mu[[component]], sigma = Sigma[[component]])
}



## Kurtotic ##

mu <- list()
Sigma <- list()

mu[[1]] <- c(0, 0)
Sigma[[1]] <- matrix(c(1, 1, 1, 4), 2, 2)

mu[[2]] <- c(0, 0)
Sigma[[2]] <- matrix(c(4/9, -1/9, -1/9, 1/9), 2, 2)

Kurtotic <- matrix(0.0, N, 2)
for(i in 1:N) {
	component <- sample(1:2, 1, prob = c(2/3, 1/3))
	Kurtotic[i, ] <- rmvnorm(1, mean = mu[[component]], sigma = Sigma[[component]])
}



## Bimodel I ##

mu <- list()
Sigma <- list()

mu[[1]] <- c(-1, 0)
Sigma[[1]] <- matrix(c(4/9, 0, 0, 4/9), 2, 2)

mu[[2]] <- c(1, 0)
Sigma[[2]] <- matrix(c(4/9, 0, 0, 4/9), 2, 2)

Bimodal1 <- matrix(0.0, N, 2)
for(i in 1:N) {
	component <- sample(1:2, 1, prob = c(1/2, 1/2))
	Bimodal1[i, ] <- rmvnorm(1, mean = mu[[component]], sigma = Sigma[[component]])
}



## Bimodel II ##

mu <- list()
Sigma <- list()

mu[[1]] <- c(-3/2, 0)
Sigma[[1]] <- matrix(c(1/16, 0, 0, 1), 2, 2)

mu[[2]] <- c(3/2, 0)
Sigma[[2]] <- matrix(c(1/16, 0, 0, 1), 2, 2)

Bimodal2 <- matrix(0.0, N, 2)
for(i in 1:N) {
	component <- sample(1:2, 1, prob = c(1/2, 1/2))
	Bimodal2[i, ] <- rmvnorm(1, mean = mu[[component]], sigma = Sigma[[component]])
}



## Bimodel III ##

mu <- list()
Sigma <- list()

mu[[1]] <- c(-1, 1)
Sigma[[1]] <- matrix(c(4/9, 4/15, 4/15, 4/9), 2, 2)

mu[[2]] <- c(1, -1)
Sigma[[2]] <- matrix(c(4/9, 4/15, 4/15, 4/9), 2, 2)

Bimodal3 <- matrix(0.0, N, 2)
for(i in 1:N) {
	component <- sample(1:2, 1, prob = c(1/2, 1/2))
	Bimodal3[i, ] <- rmvnorm(1, mean = mu[[component]], sigma = Sigma[[component]])
}



## Bimodel IV ##

mu <- list()
Sigma <- list()

mu[[1]] <- c(1, -1)
Sigma[[1]] <- matrix(c(4/9, 14/45, 14/45, 4/9), 2, 2)

mu[[2]] <- c(-1, 1)
Sigma[[2]] <- matrix(c(4/9, 0, 0, 4/9), 2, 2)

Bimodal4 <- matrix(0.0, N, 2)
for(i in 1:N) {
	component <- sample(1:2, 1, prob = c(1/2, 1/2))
	Bimodal4[i, ] <- rmvnorm(1, mean = mu[[component]], sigma = Sigma[[component]])
}



## Trimodal I ##

mu <- list()
Sigma <- list()

mu[[1]] <- c(-6/5, 6/5)
Sigma[[1]] <- matrix(c(9/25, 27/250, 27/250, 9/25), 2, 2)

mu[[2]] <- c(6/5, -6/5)
Sigma[[2]] <- matrix(c(9/25, -27/125, -27/125, 9/25), 2, 2)

mu[[3]] <- c(0, 0)
Sigma[[3]] <- matrix(c(1/16, 1/80, 1/80, 1/16), 2, 2)

Trimodal1 <- matrix(0.0, N, 2)
for(i in 1:N) {
	component <- sample(1:3, 1, prob = c(9/20, 9/20, 1/10))
	Trimodal1[i, ] <- rmvnorm(1, mean = mu[[component]], sigma = Sigma[[component]])
}



## Trimodal II ##

mu <- list()
Sigma <- list()

mu[[1]] <- c(-6/5, 0)
Sigma[[1]] <- matrix(c(9/25, 63/250, 63/250, 9/25), 2, 2)

mu[[2]] <- c(6/5, 0)
Sigma[[2]] <- matrix(c(9/25, 63/250, 63/250, 9/25), 2, 2)

mu[[3]] <- c(0, 0)
Sigma[[3]] <- matrix(c(9/25, -63/250, -63/250, 9/25), 2, 2)

Trimodal2 <- matrix(0.0, N, 2)
for(i in 1:N) {
	component <- sample(1:3, 1, prob = c(1/3, 1/3, 1/3))
	Trimodal2[i, ] <- rmvnorm(1, mean = mu[[component]], sigma = Sigma[[component]])
}



## Trimodal III ##

mu <- list()
Sigma <- list()

mu[[1]] <- c(-1, 0)
Sigma[[1]] <- matrix(c(9/25, 63/250, 63/250, 49/100), 2, 2)

mu[[2]] <- c(1, 2*sqrt(3)/3)
Sigma[[2]] <- matrix(c(9/25, 0, 0, 49/100), 2, 2)

mu[[3]] <- c(1, -2*sqrt(3)/3)
Sigma[[3]] <- matrix(c(9/25, 0, 0, 49/100), 2, 2)

Trimodal3 <- matrix(0.0, N, 2)
for(i in 1:N) {
	component <- sample(1:3, 1, prob = c(3/7, 3/7, 1/7))
	Trimodal3[i, ] <- rmvnorm(1, mean = mu[[component]], sigma = Sigma[[component]])
}



## Quadrimodal ##

mu <- list()
Sigma <- list()

mu[[1]] <- c(-1, 1)
Sigma[[1]] <- matrix(c(4/9, 8/45, 8/45, 4/9), 2, 2)

mu[[2]] <- c(-1, -1)
Sigma[[2]] <- matrix(c(4/9, 4/15, 4/15, 4/9), 2, 2)

mu[[3]] <- c(1, -1)
Sigma[[3]] <- matrix(c(4/9, -14/45, -14/45, 4/9), 2, 2)

mu[[4]] <- c(1, 1)
Sigma[[4]] <- matrix(c(4/9, -2/9, -2/9, 4/9), 2, 2)

Quadrimodal <- matrix(0.0, N, 2)
for(i in 1:N) {
	component <- sample(1:4, 1, prob = c(1/8, 3/8, 1/8, 3/8))
	Quadrimodal[i, ] <- rmvnorm(1, mean = mu[[component]], sigma = Sigma[[component]])
}


## clean up ##

rm(component, i, mu, N, Sigma)




