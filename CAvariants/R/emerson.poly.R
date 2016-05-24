emerson.poly <-
function (mj, pj) 
{
#################Emerson polynomials, recurrence formulas
#mj: natural scores for example c(1,2,3,4)
#pj
#####################
    nc <- length(mj)
    Dj <- diag(pj)
    B <- matrix(1, (nc + 1), nc)
    B[1, ] <- 0
    Sh <- Th <- Vh <- NULL
    for (i in 3:(nc + 1)) {
        for (j in 1:nc) {
            Th[i] <- mj %*% Dj %*% B[i - 1, ]^2
            Vh[i] <- mj %*% Dj %*% (B[i - 1, ] * B[i - 2, ])
            Sh[i] <- sqrt(mj^2 %*% Dj %*% B[i - 1, ]^2 - Th[i]^2 - 
                Vh[i]^2)^(-1)
            B[i, j] <- Sh[i] * ((mj[j] - Th[i]) * B[i - 1, j] - 
                Vh[i] * B[i - 2, j])
        }
    }
B<-t(B)
B1<-B[,-c(1,2)]
BT<-B[,-c(1)]
    list(B=B1,BT=BT)
}
