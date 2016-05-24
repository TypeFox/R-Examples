`mNNinfo2` <-
function (n, R, Q) 
{
    N <- sum(n)
    k <- length(n)
    l <- names(n)
    EN <- matrix(0, nrow = k, ncol = k)
    VN <- VarN <- matrix(0, nrow = k * k, ncol = k * k)
    for (i in 1:k) {
        for (j in 1:k) {
            EN[i, j] <- n[i] * (n[j] - (i == j))/(N - 1)
        }
    }
    for (l1 in 1:(k * k)) {
        i <- 1 + (l1 - 1)%/%k
        j <- 1 + (l1 - 1)%%k
        for (l2 in l1:(k * k)) {
            i2 <- 1 + (l2 - 1)%/%k
            j2 <- 1 + (l2 - 1)%%k
            if ((i == i2) & (j == j2)) {
                if (i == j) {
                  p2 <- n[i] * (n[i] - 1)/(N * (N - 1))
                  p3 <- p2 * (n[i] - 2)/(N - 2)
                  p4 <- p3 * (n[i] - 3)/(N - 3)
                  VN <- check(1, VN, l1, l2)
                  VarN[l1, l2] <- (N + R) * p2 + (2 * N - 2 * 
                    R + Q) * p3 + (N * (N - 3) - Q + R) * p4 - 
                    EN[i, j] * EN[i, j]
                }
                else {
                  p2 <- n[i] * n[j]/(N * (N - 1))
                  p3 <- p2 * (n[i] - 1)/(N - 2)
                  p4 <- p3 * (n[j] - 1)/(N - 3)
                  VN <- check(2, VN, l1, l2)
                  VarN[l1, l2] <- N * p2 + Q * p3 + (N * (N - 
                    3) - Q + R) * p4 - EN[i, j] * EN[i, j]
                }
            }
            else if ((i == j) & (i == i2) & (j != j2)) {
                p3 <- n[i] * (n[i] - 1) * n[j2]/(N * (N - 1) * 
                  (N - 2))
                p4 <- p3 * (n[i] - 2)/(N - 3)
                VN <- check(3, VN, l1, l2)
                VarN[l1, l2] <- VarN[l2, l1] <- (N - R) * p3 + 
                  (N * (N - 3) - Q + R) * p4 - EN[i, j] * EN[i2, 
                  j2]
            }
            else if ((i2 == j2) & (i == i2) & (j != j2)) {
                p3 <- n[i2] * (n[i2] - 1) * n[j]/(N * (N - 1) * 
                  (N - 2))
                p4 <- p3 * (n[i] - 2)/(N - 3)
                VN <- check(3, VN, l1, l2)
                VarN[l1, l2] <- VarN[l2, l1] <- (N - R) * p3 + 
                  (N * (N - 3) - Q + R) * p4 - EN[i, j] * EN[i2, 
                  j2]
            }
            else if ((i2 == j2) & (j == j2) & (i != i2)) {
                p3 <- n[j] * (n[j] - 1) * n[i]/(N * (N - 1) * 
                  (N - 2))
                p4 <- p3 * (n[j] - 2)/(N - 3)
                VN <- check(4, VN, l1, l2)
                VarN[l1, l2] <- VarN[l2, l1] <- (N - R + Q) * 
                  p3 + (N * (N - 3) - Q + R) * p4 - EN[i, j] * 
                  EN[i2, j2]
            }
            else if ((i == j) & (i == j2) & (i != i2)) {
                p3 <- n[i] * (n[i] - 1) * n[i2]/(N * (N - 1) * 
                  (N - 2))
                p4 <- p3 * (n[i] - 2)/(N - 3)
                VN <- check(14, VN, l1, l2)
                VarN[l1, l2] <- VarN[l2, l1] <- (N - R + Q) * 
                  p3 + (N * (N - 3) - Q + R) * p4 - EN[i, j] * 
                  EN[i2, j2]
            }
            else if ((i == j) & (i2 == j2) & (i != i2)) {
                p4 <- n[i] * (n[i] - 1) * n[i2] * (n[i2] - 1)/(N * 
                  (N - 1) * (N - 2) * (N - 3))
                VN <- check(5, VN, l1, l2)
                VarN[l1, l2] <- VarN[l2, l1] <- (N * (N - 3) - 
                  Q + R) * p4 - EN[i, j] * EN[i2, j2]
            }
            else if ((i == j) & (i2 != i) & (j2 != j) & (i2 != 
                j2)) {
                p4 <- n[i] * (n[i] - 1) * n[i2] * n[j2]/(N * 
                  (N - 1) * (N - 2) * (N - 3))
                VN <- check(6, VN, l1, l2)
                VarN[l1, l2] <- VarN[l2, l1] <- (N * (N - 3) - 
                  Q + R) * p4 - EN[i, j] * EN[i2, j2]
            }
            else if ((i2 == j2) & (i2 != i) & (j2 != j) & (i != 
                j)) {
                p4 <- n[i2] * (n[i2] - 1) * n[i] * n[j]/(N * 
                  (N - 1) * (N - 2) * (N - 3))
                VN <- check(6, VN, l1, l2)
                VarN[l1, l2] <- VarN[l2, l1] <- (N * (N - 3) - 
                  Q + R) * p4 - EN[i, j] * EN[i2, j2]
            }
            else if ((i == i2) & (i != j) & (i2 != j2) & (j != 
                j2)) {
                p4 <- n[i] * (n[i] - 1) * n[j] * n[j2]/(N * (N - 
                  1) * (N - 2) * (N - 3))
                VN <- check(7, VN, l1, l2)
                VarN[l1, l2] <- VarN[l2, l1] <- (N * (N - 3) - 
                  Q + R) * p4 - EN[i, j] * EN[i2, j2]
            }
            else if ((i == j2) & (i2 == j) & (i != j)) {
                p2 <- n[i] * n[j]/(N * (N - 1))
                p3 <- p2 * (n[i] - 1 + n[j] - 1)/(N - 2)
                p4 <- p2 * (n[i] - 1) * (n[j] - 1)/((N - 2) * 
                  (N - 3))
                VN <- check(8, VN, l1, l2)
                VarN[l1, l2] <- VarN[l2, l1] <- R * p2 + (N - 
                  R) * p3 + (N * (N - 3) - Q + R) * p4 - EN[i, 
                  j] * EN[i2, j2]
            }
            else if ((i != j) & (j == i2) & (i2 != j2) & (i != 
                j2)) {
                p3 <- n[i] * n[j] * n[j2]/(N * (N - 1) * (N - 
                  2))
                p4 <- p3 * (n[j] - 1)/(N - 3)
                VN <- check(9, VN, l1, l2)
                VarN[l1, l2] <- VarN[l2, l1] <- (N - R) * p3 + 
                  (N * (N - 3) - Q + R) * p4 - EN[i, j] * EN[i2, 
                  j2]
            }
            else if ((i != j) & (j == j2) & (i2 != j2) & (i != 
                i2)) {
                p3 <- n[i] * n[j] * n[i2]/(N * (N - 1) * (N - 
                  2))
                p4 <- p3 * (n[j] - 1)/(N - 3)
                VN <- check(10, VN, l1, l2)
                VarN[l1, l2] <- VarN[l2, l1] <- Q * p3 + (N * 
                  (N - 3) - Q + R) * p4 - EN[i, j] * EN[i2, j2]
            }
            else if ((i != j) & (i == j2) & (i2 != j2) & (j != 
                i2)) {
                p3 <- n[i] * n[j] * n[i2]/(N * (N - 1) * (N - 
                  2))
                p4 <- p3 * (n[i] - 1)/(N - 3)
                VN <- check(11, VN, l1, l2)
                VarN[l1, l2] <- VarN[l2, l1] <- (N - R) * p3 + 
                  (N * (N - 3) - Q + R) * p4 - EN[i, j] * EN[i2, 
                  j2]
            }
            else if ((i != j) & (i != i2) & (i != j2) & (j != 
                i2) & (j != j2) & (i2 != j2)) {
                p4 <- n[i] * n[j] * n[i2] * n[j2]/(N * (N - 1) * 
                  (N - 2) * (N - 3))
                VN <- check(12, VN, l1, l2)
                VarN[l1, l2] <- VarN[l2, l1] <- (N * (N - 3) - 
                  Q + R) * p4 - EN[i, j] * EN[i2, j2]
            }
        }
    }
    v <- as.vector(t(outer(l, l, paste, sep = "")))
    dimnames(VN) <- dimnames(VarN) <- list(v, v)
    list(EN = EN, VN = VN, VarN = VarN)
}

