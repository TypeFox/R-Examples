
 library(trust)

 options(digits = 3)

 ##### four-way contingency table with all two-way interactions

 d <- c(3, 4, 5, 6)
 n <- 1000

 ##### model matrix
 m <- NULL
 for (i in 1:d[1]) {
     for (j in 1:d[2]) {
         mfoo <- array(0, dim = d)
         mfoo[i, j, , ] <- 1
         mfoo <- as.vector(mfoo)
         m <- cbind(m, mfoo)
     }
 }
 for (i in 1:d[1]) {
     for (j in 1:d[3]) {
         mfoo <- array(0, dim = d)
         mfoo[i, , j, ] <- 1
         mfoo <- as.vector(mfoo)
         m <- cbind(m, mfoo)
     }
 }
 for (i in 1:d[1]) {
     for (j in 1:d[4]) {
         mfoo <- array(0, dim = d)
         mfoo[i, , , j] <- 1
         mfoo <- as.vector(mfoo)
         m <- cbind(m, mfoo)
     }
 }
 for (i in 1:d[2]) {
     for (j in 1:d[3]) {
         mfoo <- array(0, dim = d)
         mfoo[ , i, j, ] <- 1
         mfoo <- as.vector(mfoo)
         m <- cbind(m, mfoo)
     }
 }
 for (i in 1:d[2]) {
     for (j in 1:d[4]) {
         mfoo <- array(0, dim = d)
         mfoo[ , i, , j] <- 1
         mfoo <- as.vector(mfoo)
         m <- cbind(m, mfoo)
     }
 }
 for (i in 1:d[3]) {
     for (j in 1:d[4]) {
         mfoo <- array(0, dim = d)
         mfoo[ , , i, j] <- 1
         mfoo <- as.vector(mfoo)
         m <- cbind(m, mfoo)
     }
 }
 dimnames(m) <- NULL
 foo <- qr(m)
 m <- m[ , foo$pivot[seq(1, foo$rank)]]

 ##### true parameter value
 set.seed(42)
 theta.true <- 0.25 * rnorm(ncol(m))
 theta.true <- round(theta.true, 5)

 ##### simulate data
 eta <- as.numeric(m %*% theta.true)
 p <- exp(eta)
 p <- p / sum(p)
 x <- sample(nrow(m), n, replace = TRUE, prob = p)
 x <- tabulate(x, nbins = nrow(m))

 ##### save data
 iffy <- try(read.table("fred.txt"), silent = TRUE)
 if (inherits(iffy, "try-error")) {
     data <- data.frame(x = x, m = m)
     write.table(data, file = "fred.txt", row.names = FALSE)
 }
 data <- read.table(file = "fred.txt", header = TRUE)
 x <- data$x
 data$x <- NULL
 m <- as.matrix(data)
 dimnames(m) <- NULL

 ##### log likelihood
 objfun <- function(theta) {
     eta <- as.numeric(m %*% theta)
     p <- exp(eta)
     f <- sum(x * eta - p)
     g <- as.numeric(t(x - p) %*% m)
     B <- sweep(- m, 1, p, "*")
     B <- t(m) %*% B
     list(value = f, gradient = g, hessian = B)
 }

 ##### check it
 sally <- objfun(theta.true)
 epsilon <- 1e-8
 mygrad <- double(length(theta.true))
 for (i in 1:length(mygrad)) {
     theta.eps <- theta.true
     theta.eps[i] <- theta.true[i] + epsilon
     sally.eps <- objfun(theta.eps)
     mygrad[i] <- (sally.eps$value - sally$value) / epsilon
 }
 all.equal(sally$gradient, mygrad, tolerance = length(mygrad) * epsilon)
 myhess <- matrix(NA, length(theta.true), length(theta.true))
 for (i in 1:length(mygrad)) {
     theta.eps <- theta.true
     theta.eps[i] <- theta.true[i] + epsilon
     sally.eps <- objfun(theta.eps)
     myhess[i, ] <- (sally.eps$gradient - sally$gradient) / epsilon
 }
 all.equal(sally$hessian, myhess, tolerance = length(mygrad) * epsilon)

 fred <- trust(objfun, theta.true, 1, sqrt(ncol(m)), blather = TRUE)
 fred$converged
 max(abs(fred$gradient))
 length(fred$r)
 data.frame(type = fred$steptype, rho = fred$rho, change = fred$preddiff,
     accept = fred$accept, r = fred$r)
 (fred$stepnorm / fred$r)[fred$accept & fred$steptype != "Newton"]

 ##### note: FAILS to converge because function is unbounded below -- minimum
 #####     value does not exist
 #####
 ##### this is what happens when a luser forgets minimize = FALSE in a
 ##### maximization problem.
 #####
 ##### also it revealed a bug in trust (now fixed), where it used to
 ##### test whether fred(beta.up) == 0 before calling uniroot

