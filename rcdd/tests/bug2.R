
 library(rcdd)

 A <- rbind(
    c(0,  1.000,   -1,    0,    0,    0),
    c(0,  1.000,   -1,    0,    0,    0),
    c(0, -0.001,    1,    0,    0,    0),
    c(0, -0.001,    0,    1,    0,    0),
    c(0, -0.001,    0,    0,    1,    0),
    c(0, -0.001,    0,    0,    0,    1),
    c(1,  1.000,   -1,   -1,   -1,   -1),
    c(1,  0.300,   -1,   -1,   -1,    0),
    c(1,  7.990,   -1,   -3,   -5,   -8))

 b <- c(0, 0, 0, 0)

 out <- lpcdd(A, b)
 out

 all(A[out$dual.direction < 0, 1] == 1)

 fred <- rbind(out$dual.direction) %*% A
 fred <- as.numeric(fred)
 sally <- fred[2]
 fred <- fred[- c(1, 2)]
 all.equal(max(abs(fred)), 0)
 
 sally >= 0
 ### if FALSE proves the constraints cannot be satisfied


 Arat <- 1000 * A
 Arat <- round(Arat)
 Arat <- z2q(Arat, 0 * Arat + 1000)
 Arat

 brat <- d2q(b)
 brat

 out.rat <- lpcdd(Arat, brat)
 out.rat

 out$solution.type == out.rat$solution.type
 all.equal(out$dual.direction, q2d(out.rat$dual.direction))

 fred <- qmatmult(rbind(out.rat$dual.direction), Arat)
 sally <- fred[2]
 fred <- fred[- c(1, 2)]
 all(fred == "0")
 
 qsign(sally) >= 0
 ### if FALSE proves the constraints cannot be satisfied

