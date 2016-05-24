# ################################################
# ------------------------------------------------
# QR
# ------------------------------------------------
# ################################################

library(pbdDMAT , quiet=T)

init.grid()


comm.set.seed(seed=1234, diff=F)

mean <- 10
sd <- 100

n <- 250
p <- 100

n <- 10
p <- 6

x <- matrix(rnorm(n*p), n, p)
y <- matrix(rnorm(n), nrow=n)

dx <- as.ddmatrix(x)
dy <- as.ddmatrix(y)


# In general, the "QR" returns will be different because the algorithms
# differ.  However, Q and R must agree, so they are tested below.
qr1 <- qr(dx)
qr2 <- qr(x)

# Q
out1 <- qr.Q(qr1)
out1 <- as.matrix(out1)
out2 <- qr.Q(qr2)

test <- all.equal(out1, out2)
comm.print(test, quiet=T)

# R
out1 <- qr.R(qr1)
out1 <- as.matrix(out1)
out2 <- qr.R(qr2)

test <- all.equal(out1, out2)
comm.print(test, quiet=T)

# Qy / Q^Ty
out1 <- qr.qy(qr1, dy)
out1 <- as.matrix(out1)
out2 <- qr.qy(qr2, y)

test <- all.equal(out1, out2)
comm.print(test, quiet=T)


out1 <- qr.qty(qr1, dy)
out1 <- as.matrix(out1)
out2 <- qr.qty(qr2, y)

test <- all.equal(out1, out2)
comm.print(test, quiet=T)



finalize()
