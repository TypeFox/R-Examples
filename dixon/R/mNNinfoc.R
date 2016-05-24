`mNNinfoc` <-
function (xy, label, nnid = NULL, splancs = TRUE, fortran) 
{
    if (is.null(nnid)) {
        nnid <- NNid(xy, splancs)
    }
    n <- table(label)
    N <- sum(n)
    l <- names(n)
    k <- length(l)
    R <- sum((1:N) == nnid[nnid])
    Q1 <- table(c(1:6, table(nnid)))
    Q <- 2 * Q1[2] + 6 * Q1[3] + 12 * Q1[4] + 20 * Q1[5] + 30 * 
        Q1[6] - 70
    ON <- matrix(0, nrow = k, ncol = k)
    for (i in 1:k) {
        for (j in 1:k) {
            ON[i, j] <- sum((label == l[i]) & (label[nnid] == 
                l[j]))
        }
    }
    if (fortran==TRUE) temp <- mNNinfo2b(n, R, Q) else  temp <- mNNinfo2(n, R, Q)
    rownames(ON) <- colnames(ON) <- rownames(temp$EN) <- colnames(temp$EN) <- l
    list(ON = ON, EN = temp$EN, VarN = temp$VarN, R = R, Q = Q)
}


