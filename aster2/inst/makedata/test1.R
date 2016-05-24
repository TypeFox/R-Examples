
 ##### this is the script used to make the data set test1.rda #####

 RNGversion("2.11.1")
 set.seed(42)

 n <- 100

 # multinomial
 m1 <- as.numeric(runif(n) < 1 / 3)
 m2 <- as.numeric(runif(n) < 1 / 2)
 m2 <- m2 * (1 - m1)
 m3 <- 1 - m1 - m2

 # normal
 x <- rnorm(n, 5, 1)
 n1 <- x * m1
 n2 <- x^2 * m1

 # bernoulli
 b1 <- as.numeric(runif(n) < 1 / 2)
 b1 <- b1 * m2

 # poisson
 p1 <- rpois(n, 3)
 p1 <- p1 * m3

 # zero-truncated poisson
 z1 <- rpois(10 * n, 4)
 z1 <- z1[z1 > 0]
 z1 <- z1[1:n]
 z1 <- z1 * b1

 test1 <- data.frame(m1, m2, m3, n1, n2, b1, p1, z1)

 names(test1)
 dim(test1)
 sapply(test1, class)

 save(test1, file = "test1.rda")

