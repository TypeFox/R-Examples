Aintmap <-
function (L, R, Lin = NULL, Rin = NULL) 

{

    n <- length(L)

    if (is.null(Lin) & is.null(Rin))
    {

        Lin <- rep(FALSE, n)

        Rin <- rep(TRUE, n)

        Lin[L == R] <- TRUE

        Rin[R == Inf] <- FALSE

    }

    else if (length(Lin) == 1 & length(Rin) == 1 & is.logical(Lin) & is.logical(Rin))
    {

        Lin <- rep(Lin, n)

        Rin <- rep(Rin, n)

    }

    else if (length(Lin) != n | length(Rin) != n | !all(is.logical(Lin)) | !all(is.logical(Rin)))
    {

        stop("Lin and Rin should be either NULL, logical length 1 or length same as L,R")

    }

    if (n != length(R))
 
        stop("length of L and R must be the same")

    LRvalues <- sort(unique(c(0, L, R, Inf)))

    eps <- min(diff(LRvalues))/2

    Le <- L

    Re <- R

    Le[!Lin] <- L[!Lin] + eps

    Re[!Rin] <- R[!Rin] - eps

    oLR <- order(c(Le, Re + eps/2))

    Leq1.Req2 <- c(rep(1, n), rep(2, n))

    flag <- c(0, diff(Leq1.Req2[oLR]))

    R.right.of.L <- (1:(2 * n))[flag == 1]

    intmapR <- c(L, R)[oLR][R.right.of.L]

    intmapL <- c(L, R)[oLR][R.right.of.L - 1]

    intmapRin <- c(Lin, Rin)[oLR][R.right.of.L]

    intmapLin <- c(Lin, Rin)[oLR][R.right.of.L - 1]

    intmap <- matrix(c(intmapL, intmapR), byrow = TRUE, nrow = 2)

    attr(intmap, "LRin") <- matrix(c(intmapLin, intmapRin), byrow = TRUE, 
nrow = 2)

    k <- dim(intmap)[[2]]

    Lbracket <- rep("(", k)

    Lbracket[intmapLin] <- "["

    Rbracket <- rep(")", k)

    Rbracket[intmapRin] <- "]"

    intname <- paste(Lbracket, intmapL, ",", intmapR, Rbracket, 
sep = "")

    A <- matrix(0, n, k, dimnames = list(1:n, intname))

    intmapLe <- intmapL

    intmapLe[!intmapLin] <- intmapL[!intmapLin] + eps

    intmapRe <- intmapR

    intmapRe[!intmapRin] <- intmapR[!intmapRin] - eps

    for (i in 1:n) {

        tempint <- Le[i] <= intmapRe & Re[i] >= intmapLe

        A[i, tempint] <- 1

    }

    out <- list(A = A, intmap = intmap)

    out

}

