x <- as.matrix(c(3,1))            # vector as column matrix
A <- rbind(c(5,2),c(3,1))
Ainv <- solve(A); Ainv            # solve() computes inverse
A %*% Ainv
Ainv %*% A
Ainv %*% x                        # solution to system
