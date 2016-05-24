sub.pipath <- 
function (lambda, x, y, kernel = "linear", kparam = NULL, eps = 1e-07, Nmoves = 100 * n, ridge = 0) 
{
    K <- Kmat(x, x, kernel, kparam)
    if (is.null(K)) stop("kernel function is not properly defined") 

    n <- length(y)
    index.plus <- which(y == 1)
    index.minus <- which(y == -1)
    nplus <- length(index.plus)
    nminus <- length(index.minus)
    Right <- Elbow <- NULL
    Left <- seq(n)

    Kscript <- K * outer(y, y)
    Kstar <- matrix(0, 1, 1)
    alpha <- matrix(0, n, Nmoves)
    alpha0 <- double(Nmoves)
    Elbow.list <- as.list(seq(Nmoves))
    Left.list <- as.list(seq(Nmoves))
    Pi <- double(Nmoves)
    init <- pi.initialization(x, y, K, kernel, kparam, lambda) ###

    Elbow.list[[1]] <- Elbow <- init$Elbow
    Left.list[[1]] <- Left <- init$Left
    Right <- init$Right
    alpha0[1] <- init$alpha0
    alpha[, 1] <- init$alpha
    Pi[1] <- 0.5
    first <- c(0, y[Elbow])
    rest <- matrix(cbind(y[Elbow], Kscript[Elbow, Elbow]), length(Elbow), 
        length(Elbow) + 1)
    Kstar <- rbind(first, rest)
    fl <- (K %*% (alpha[, 1] * y) + alpha0[1])/lambda
    k <- 1
    for (k in 1:Nmoves) {
        if (length(Elbow) == 0) {
            nlplus <- length(intersect(Left, index.plus))
            nl <- length(Left)
            balanced <- abs(nlplus/nl - Pi[k]) < eps
            if (balanced) {
                empty <- empty.solution(pi = Pi[k], a = alpha[, 
                  k], a0 = alpha0[k], lambda, x, y, K, Left, Right)
                Pi[k + 1] <- empty$pi
                alpha[, k + 1] <- empty$alpha
                alpha0[k + 1] <- empty$alpha0
                Elbow <- empty$Elbow
                Left <- empty$Left
                Right <- empty$Right
                first <- c(0, y[Elbow])
                rest <- matrix(cbind(y[Elbow], Kscript[Elbow, 
                  Elbow]), 1, 2)
                Kstar <- rbind(first, rest)
                fl <- (K %*% (alpha[, k + 1] * y) + alpha0[k + 
                  1])/lambda
            }
            else stop("Empty elbow situation with unbalanced left set")
        }
        else {
            ne <- length(Elbow)
            nl <- length(Left)
            Kle <- matrix(K[Left, Elbow], nl, ne)
            kl <- y[Elbow] * apply(Kle, 2, sum)
            bstar <- PiSolveKstar(Kstar, vec = kl, con = nl, ridge = ridge)
            b0 <- bstar[1]
            b <- bstar[-1]
            gl <- (matrix(K[, Elbow], n, ne) %*% (y[Elbow] * 
                b) - matrix(K[, Left], n, nl) %*% rep(1, nl) + 
                b0)/lambda
            temp <- -alpha[Elbow, k] + Pi[k] * b
            pi.left <- (y[Elbow] == -1) * temp/(b - 1) + (y[Elbow] == 
                1) * (1 + temp)/(1 + b)
            pi.right <- temp/b
            pi01 <- c(pi.right, pi.left)
            pi.exit <- min(pi01[pi01 > (Pi[k] + eps)], 2)
            pi.exit[is.nan(pi.exit) + is.na(pi.exit) != 0] <- 2
            pii <- (y - fl)/gl + Pi[k]
            pii[is.nan(pii)] <- -2
            pi.entry <- min(pii[pii > (Pi[k] + eps)], 2)
            pi.entry[is.nan(pi.entry) + is.na(pi.entry) != 0] <- 2
            new.pi <- min(pi.entry, pi.exit)
            Pi[k + 1] <- new.pi
            w <- wvec(Pi[k + 1], y = y)
            alpha[, k + 1] <- alpha[, k]
            alpha[Elbow, k + 1] <- alpha[Elbow, k] - (Pi[k] - 
                Pi[k + 1]) * b
            alpha[Left, k + 1] <- w[Left]
            alpha0[k + 1] <- alpha0[k] - (Pi[k] - Pi[k + 1]) * 
                b0
            fl <- fl + (Pi[k + 1] - Pi[k]) * gl
            if (abs(sum(fl) + n) < eps) {
                immobile = prod(abs(fl + 1) < eps)
                if (immobile == 1) 
                  break
            }
            if (pi.exit < pi.entry) {
                i1 <- Elbow[pi.left == new.pi]
                i2 <- Elbow[pi.right == new.pi]
                i <- union(i1, i2)
                Elbow <- setdiff(Elbow, i)
                Left <- union(Left, i1)
                Right <- union(Right, i2)
            }
            else {
                i1 <- Left[pii[Left] == new.pi]
                i2 <- Right[pii[Right] == new.pi]
                Elbow <- c(Elbow, i1, i2)
                Left <- setdiff(Left, i1)
                Right <- setdiff(Right, i2)
            }
            first <- c(0, y[Elbow])
            rest <- matrix(cbind(y[Elbow], Kscript[Elbow, Elbow]), 
                length(Elbow), length(Elbow) + 1)
            Kstar <- rbind(first, rest)
        }
        if (Pi[k + 1] > 1 - eps) 
            break
        k <- k + 1
        Elbow.list[[k]] <- Elbow
        Left.list[[k]] <- Left
    }
    obj <- list(Pi = Pi[seq(k)], alpha = alpha[, seq(k)], alpha0 = alpha0[seq(k)], 
        pi = Pi[seq(k)], Elbow = Elbow.list[seq(k)], Left = Left.list[seq(k)], 
        kernel = kernel, K = K,
        param.kernel = kparam, x = x, y = y, lambda = lambda)
    obj
}
