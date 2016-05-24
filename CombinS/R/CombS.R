CombS <-
function (n, l, s) 
{
    V <- n * l
    if (s > l || s <= 1) {
        "Choose : 1 < s <= l"
    }
    else {
reso<-(n*l)%%(2*s)
bbo<-ifelse(reso==0,"Yes","No")
        lamda <- NULL
        Op <- function(n, l, s) {
            V <- n * l
            mt <- matrix(1:V, nrow = n, byrow = TRUE)
            n <- dim(mt)[1]
            l <- dim(mt)[2]
            V <- n * l
            a <- combn(1:l, s)
            b <- combn(1:n, 2)
            v <- dim(b)[2]
            vv <- dim(a)[2]
            MAT <- NULL
            A <- 1
            y <- 1
            while (A <= vv) {
                for (k in 1:v) {
                  zz <- a[, A]
                  ss <- b[, k]
                  MAT[[y]] <- as.vector(t(mt[ss, zz]))
                  y <- y + 1
                }
                A <- A + 1
            }
            return(Reduce("rbind", MAT))
        }
        BIB <- Op(n, l, s)

        T <- BIB[1, 1]
        R <- length(which(T == BIB))
        if (l ==s) {
            lamda[[1]] <- (n - 1) * choose(l - 2, s - 2)
            lamda[[2]] <- choose(l - 1, s - 1)
            lamda[[3]] <- choose(l - 2, s - 2)
            return(list(PBIB = BIB, Type = "Singular group divisible design", 
                V = V, B = dim(BIB)[1], R = R, K = 2 * s, lamda = lamda, Resolvable=bbo))
        }
        else {
            lamda[[1]] <- R
            lamda[[2]] <- 1

            return(list(PBIB = BIB, Type = "Rectangular PBIB design", 
                V = V, B = dim(BIB)[1], R = R, K = 2 * s, lamda = lamda, Resolvable=bbo))
        }
    }
}
