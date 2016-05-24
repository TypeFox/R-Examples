# 1. Data input
library(erer); data(daBed); colnames(daBed)
y <- daBed 
lab8 <- c("CN", "VN", "ID", "MY", "CA", "BR", "IT")
share <- paste("s", c(lab8, "RW"), sep = "")
price <- paste("lnp", c(lab8, "RW"), sep = "")
shift <- c("dum1", "dum2", "dum3")
expen <- "rte" 

# 2. Determining dimensions
nShare <- length(share)                      # nShare = 8
nExoge <- 1 + length(expen) + length(shift)  # nExoge = 5
nParam <- nExoge + nShare                    # nParam = 13
nTotal <- (nShare - 1) * nParam              # nTotal = 91

# 3. Restriction for homogeneity
# 3.1. Approach A: one matrix
m.h <- matrix(data = 0, nrow = nShare - 1, ncol = nTotal)
for (i in 1:(nShare - 1)) {
  for (j in 1:nShare) {
    m.h[i, (i - 1) * nParam + nExoge + j] <- 1 
  }
}

# 3.2. Approach B: two matrices
da <- matrix(data = 0, nrow = nShare - 1, ncol = nExoge)
db <- matrix(data = 0, nrow = nShare - 1, ncol = nShare)
m.h2 <- NULL
for (i in 1:(nShare - 1)) {
  db2 <- db
  db2[i, ] <- 1
  dc <- cbind(da, db2)
  m.h2 <- rbind(m.h2, c(t(dc)))
}

# 4. Restriction for symmetry
# 4.1. Approach A: one matrix
m.s <- matrix(0, nrow = (nShare - 2) * (nShare - 1) / 2, ncol = nTotal)
k <- 0
for (i in 1:(nShare - 2)) {
  for (j in (i + 1):(nShare - 1)) {
    k <- k + 1
    m.s[k, (i - 1) * nParam + nExoge + j] <- 1
    m.s[k, (j - 1) * nParam + nExoge + i] <- -1
  }
}

# 5. Restrictions: right side
r.h  <- rep(x = 0, times = nrow(m.h))
r.s  <- rep(x = 0, times = nrow(m.s))

# 6. Output
identical(m.h, m.h2)
str(m.h); str(m.s); str(r.h); str(r.s) 
m.h[, 1:24]
m.s[1:3, 1:24] 