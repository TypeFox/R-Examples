
 library(rcdd)

 ### optimal solution exists -- file samplelp1.ine in cddlib ###

 hrep <- rbind(
 c(0, 1,  1,  0,  0),
 c(0, 1,  0,  1,  0),
 c(0, 1,  0,  0,  1),
 c(0, 1, -1,  0,  0),
 c(0, 1,  0, -1,  0),
 c(0, 1,  0,  0, -1))

 a <- c(1, 1, 1)

 lpcdd(hrep, a, minimize = FALSE)

 ### optimal solution exists -- file samplelp2.ine in cddlib ###

 hrep <- rbind(
 c("0",   "0", "1",  "1",  "0",  "0"),
 c("0",   "0", "0",  "2",  "0",  "0"),
 c("1",   "3", "0", "-1",  "0",  "0"),
 c("1", "9/2", "0",  "0", "-1", "-1"))

 a <- c("2", "3/5", "0", "0")

 lpcdd(hrep, a)

 ### primal inconsistent problem ###

 hrep <- rbind(
 c("0",  "0",  "1",  "0"),
 c("0",  "0",  "0",  "1"),
 c("0", "-2", "-1", "-1"))

 a <- c("1", "1")

 lpcdd(hrep, a)

 lpcdd(q2d(hrep), q2d(a))

 ### dual inconsistent problem ###

 hrep <- rbind(
 c("0", "0", "1", "0"),
 c("0", "0", "0", "1"))

 a <- c("1", "1")

 lpcdd(hrep, a, minimize = FALSE)

 lpcdd(q2d(hrep), q2d(a), minimize = FALSE)

 ### negative Lagrange multipliers

 set.seed(42)
 d <- 20
 k <- 6
 foo <- matrix(sample(seq(-1000, 1000), k * d, replace = TRUE), ncol = d)
 foo <- rbind(foo, diag(d))
 foo <- rbind(foo, - diag(d))
 foo <- cbind(c(rep(0, k), rep(1, 2 * d)), foo)
 foo <- cbind(c(rep(1, k), rep(0, 2 * d)), foo)

 w <- sample(seq(-1000, 1000), d, replace = TRUE)

 out <- lpcdd(d2q(foo), d2q(w))

 out$solution.type

 q2d(out$primal.solution)

 x <- out$primal.solution
 lambda <- qneg(out$dual.solution)   ### see tutorial
 l <- foo[ , 1]
 b <- foo[ , 2]
 v <- foo[ , - c(1, 2)]

 ##### check gradient of Lagrangian function zero
 all(qsign(qmq(w, qmatmult(rbind(lambda), v))) == 0)

 ##### check primal feasibility
 slack <- qpq(b, qmatmult(v, cbind(x)))
 all(qsign(slack) >= 0)

 ##### check dual feasibility
 all(qsign(lambda) * (1 - l) >= 0)

 ##### check complementary slackness
 all(qsign(slack) * qsign(lambda) == 0)

 ##### number of negative lagrange multipliers (shows exercised relevant code)
 sum(qsign(lambda) < 0)
