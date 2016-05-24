##
##  r a n k . r  Test suite
##


Rank <- pracma::Rank
nullspace <- pracma::nullspace

magic <- pracma::magic
all.equal(Rank(c()), 0)
r <- numeric(8)
for (i in 3:10){
    r[i-2] <- Rank(magic(i))
}
identical(r, c(3, 3, 5, 5, 7, 3, 9, 7))

hilb <- pracma::hilb
identical(Rank(hilb(6)), 6L)

# Vorzeichenwechsel moeglich
# N <- nullspace(magic(4))
# all.equal(as.numeric(N), c(0.2236, 0.6708, -0.6708, -0.2236),
#           tolerance = 1e-5)
