# tennis

Pmat <- function(p) {
  # transition matrix where p is prob A wins a point

  # states
  #01 0,0   #02 15,0   #03 30,0   #04 40,0
  #05 0,15  #06 15,15  #07 30,15  #08 40,15
  #09 0,30  #10 15,30  #11 30,30  #12 40,30
  #13 0,40  #14 15,40  #15 30,40  #16 Deuce 
  #17 advA
  #18 advB
  #19 Awin
  #20 Bwin
  
  P <- matrix(0, 20, 20)
  
  P[1,2] <- P[2,3] <- P[3,4] <- P[4,19]<- p
  P[1,5] <- P[2,6] <- P[3,7] <- P[4,8] <- 1-p
  P[5,6] <- P[6,7] <- P[7,8] <- P[8,19]<- p
  P[5,9] <- P[6,10]<- P[7,11] <-P[8,12]<- 1-p

  P[9,10] <- P[10,11] <- P[11,12] <-P[12,19] <- p
  P[9,13] <- P[10,14] <- P[11,15] <-P[12,16] <- 1-p

  P[13,14] <-P[14,15] <- P[15,16] <-P[16,17] <- p
  P[13,20] <-P[14,20] <- P[15,20] <-P[16,18] <- 1-p
 
  P[17,19] <- p ; P[17,16] <- 1-p
  P[18,16] <- p ; P[18,20] <- 1-p
  
  P[19,19] <- 1
  P[20,20] <- 1
  
  return(P)
}


# probability A wins
pvec <- seq(0,1,.05)
pwin <- rep(0, length(pvec))
for (i in 1:length(pvec)) {
  P <- Pmat(pvec[i])
  T <- P[1:18,1:18]
  S <- P[1:18,19:20]
  A <- solve(diag(rep(1,18)) - T, S)
  pwin[i] <- A[1,1]
}
plot(pvec, pwin, type="l")


# length of a game when p = 0.5
P <- Pmat(0.5)
n <- 0                      # time
pi_n <- c(1, rep(0, 19))    # pi_n[i] = Pr(X(n) = i)
cdf <- pi_n[19] + pi_n[20]  # cdf[i+1] = Pr(absorbtion time <= i)
max_n <- 30
for (n in 1:max_n) {
  pi_n <- pi_n %*% P
  cdf[n+1] <- pi_n[19] + pi_n[20]
}
pdf <- cdf - c(0, cdf[1:max_n])  # pdf[i+1] = Pr(absorbtion time = i)
plot(0:max_n, pdf,type="h")

# check expected absorbtion time
sum(0:max_n*pdf)
T <- P[1:18,1:18]
U1 <- solve(diag(rep(1,18)) - T, rep(1,18))
U1[1]


# length of a game by simulation
n_fin <- function(P) {
  # time until absorbtion
  x <- 1  # current state
  n <- 0  # current time
  while (x <= 18) {
    x <- sample(1:20, 1, prob=P[x,])  # update state
    n <- n + 1                        # update time
  }
  return(n)
}

N <- 10000  # num of simulations
n_vec <- rep(0, N)
for (i in 1:N) n_vec[i] <- n_fin(P)
for (n in 0:max_n) points(n, mean(n_vec==n), pch=2, col="red")
mean(n_vec)
