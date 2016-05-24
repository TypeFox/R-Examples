UType <-
function (lst) {
bob<-lst$Resolvable

if (bob=="Yes") {
mat<-lst$PBIB
    w <- mat
    v <- sort(unique(as.vector(w)))
    W <- w
    malist <- NULL
    TRI <- function(a, b) all(a %in% b)
    k <- 1
    while (TRUE) {
        vv <- v
        malist[[k]] <- c(0, 0)
        it <- 0
        while (TRUE) {
            it <- it + 1
            bool <- apply(W, 1, TRI, vv)
            if (!any(bool)) {
                break
            }
            else {
                ind <- which(bool)[1]
                u <- W[ind, ]
                vv <- vv[!(vv %in% u)]
                W <- W[-ind, , drop = FALSE]
                malist[[k]] <- rbind(malist[[k]], u)
            }
        }
        malist[[k]] <- malist[[k]][-1, ]
        if (!all(as.vector(W) %in% v) | length(W) == 0) 
            break
        k <- k + 1
    }
    x <- Reduce("rbind", malist)
    a <- max(x)
    b <- length(malist)
    c <- dim(x)[1]
    lev <- c/b
    ex <- length(v)
    UD <- matrix(nrow = a, ncol = b)
    v <- sort(unique(as.vector(x)))
    for (i in 1:b) {
        q <- malist[[i]]
        e <- c()
        for (j in v) {
            e[j] <- which(q == j)
            if (e[j] > lev) {
                e[j] <- (e[j]%%lev)
            }
            if (e[j] == 0) {
                e[j] <- lev
            }
        }
        UD[, i] <- e
    }
    ud <- na.omit(UD)
    N <- dim(ud)[1]
    F <- dim(ud)[2]
    return(list(v = N, r = F, UtypeDesign = ud))
}
else {"The PBIB is not resolvable"}}
