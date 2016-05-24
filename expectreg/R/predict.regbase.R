predict.regbase <-
function (object, newdata = NULL, ...) 
{
    type = object$type
    bnd = object$bnd
    Zspathelp = object$Zspathelp
    phi = object$phi
    center = object$center
    x = object$x
    P = object$P
    if (is.null(newdata)) 
        newdata = x
    else if (all(object$xname %in% names(newdata))) 
        newdata = newdata[, names(newdata) %in% object$xname]
    else if (!is.vector(newdata)) 
        stop("Names of newdata not consistent with original.")
    if (type == "pspline") {
        B.deg = 2
        B.size = 20
        diff.size = 2
        x0 <- min(x, na.rm = TRUE) - 0.001
        x1 <- max(x, na.rm = TRUE) + 0.001
        dx = (x1 - x0)/(B.size - 1)
        B = splineDesign(knots = seq(x0 - dx * B.deg, x1 + dx * 
            B.deg, by = dx), x = newdata, ord = B.deg + 1, outer.ok = TRUE)
        P <- diag(dim(B)[2])
        P <- diff(P, diff = diff.size)
        if (center) {
            tildeU <- matrix(0, dim(B)[2], diff.size)
            for (i in 1:diff.size) tildeU[, i] <- (1:(dim(B)[2]))^(i - 
                1)
            tildeZ <- t(P) %*% solve(P %*% t(P))
            P = t(P) %*% P
            e = eigen(P)
            tildeU = e$vectors[, ncol(e$vectors):(ncol(e$vectors) - 
                diff.size)]
            tildeZ = e$vectors[, (ncol(e$vectors) - diff.size - 
                1):1] %*% diag(1/sqrt(e$values[(ncol(e$vectors) - 
                diff.size - 1):1]))
            U <- B %*% tildeU
            Z <- B %*% tildeZ
            B = cbind(U[, -1], Z)
            P <- diag(c(rep(0, ncol(U) - 1), rep(1, ncol(Z))))
        }
    }
    else if (type == "2dspline") {
        B.deg = 2
        B.size = 20
        diff.size = 2
        x0 <- min(x[, 1], na.rm = TRUE) - 0.001
        x1 <- max(x[, 1], na.rm = TRUE) + 0.001
        dx = (x1 - x0)/(B.size - B.deg)
        Bx = splineDesign(knots = seq(x0 - dx * B.deg, x1 + dx * 
            B.deg, by = dx), x = newdata[, 1], ord = B.deg + 
            1, outer.ok = TRUE)
        y0 <- min(x[, 2], na.rm = TRUE) - 0.001
        y1 <- max(x[, 2], na.rm = TRUE) + 0.001
        dy = (y1 - y0)/(B.size - B.deg)
        By = splineDesign(knots = seq(y0 - dy * B.deg, y1 + dy * 
            B.deg, by = dy), x = newdata[, 2], ord = B.deg + 
            1, outer.ok = TRUE)
        B = matrix(NA, nrow = dim(newdata)[1], ncol = dim(Bx)[2] * 
            dim(By)[2])
        for (i in 1:dim(Bx)[1]) B[i, ] = as.vector(Bx[i, ] %o% 
            By[i, ])
        D = diff(diag(dim(Bx)[2]), diff = diff.size)
        P = t(D) %*% D
        P = diag(dim(Bx)[2]) %x% P + P %x% diag(dim(Bx)[2])
        if (center) {
            tildeX = matrix(1, nrow = dim(Bx)[2] * dim(By)[2], 
                ncol = 4)
            tildeX[, 2] = rep(1:dim(Bx)[2], times = dim(By)[2])
            tildeX[, 3] = rep(1:dim(By)[2], each = dim(Bx)[2])
            tildeX[, 4] = tildeX[, 2] * tildeX[, 3]
            ek = eigen(P)
            ek$values = ek$values[1:(length(ek$values) - 2 * 
                diff.size)]
            ek$vectors = ek$vectors[, 1:(dim(ek$vectors)[2] - 
                2 * diff.size)]
            L = ek$vectors %*% sqrt(diag(ek$values))
            tildeU = L %*% solve(t(L) %*% L)
            X = B %*% tildeX
            U = B %*% tildeU
            B = cbind(X[, -1], U)
            P <- diag(c(rep(0, ncol(X) - 1), rep(1, ncol(U))))
        }
    }
    else if (type == "radial") {
        x = x[order(x[, 1]), ]
        knots = Zspathelp
        B = matrix(NA, nrow = dim(newdata)[1], ncol = dim(knots)[1])
        for (j in 1:dim(knots)[1]) {
            r = sqrt(rowSums((newdata - matrix(unlist(knots[j, 
                ]), nrow = nrow(newdata), ncol = ncol(knots), 
                byrow = T))^2))
            r[r == 0] = 1
            B[, j] = r^2 * log(r)
        }
    }
    else if (type == "krig") {
        c = 9.233
        x = x[order(x[, 1]), ]
        knots = Zspathelp
        B = matrix(NA, nrow = dim(newdata)[1], ncol = dim(knots)[1])
        P = matrix(0, nrow = dim(B)[2], ncol = dim(B)[2])
        for (i in 1:dim(B)[2]) for (j in 1:dim(B)[2]) {
            P[i, j] = sqrt(sum((knots[i, ] - knots[j, ])^2))
        }
        phi = max(P)/c
        for (j in 1:dim(knots)[1]) {
            r = sqrt(rowSums((newdata - matrix(unlist(knots[j, 
                ]), nrow = nrow(newdata), ncol = ncol(knots), 
                byrow = T))^2))
            B[, j] = exp(-r/phi) * (1 + r/phi)
        }
    }
    else if (type == "markov") {
        P2 = bnd2gra(bnd)
        districts = dimnames(P2)[[1]]
        B = matrix(0, nrow = length(newdata), ncol = dim(P2)[2])
        for (i in 1:length(newdata)) B[i, which(districts == 
            newdata[i])] = 1
        if (center) {
            B <- B %*% Zspathelp
        }
    }
    else if (type == "random") {
        districts = sort(unique(x))
        B = matrix(0, nrow = length(newdata), ncol = length(districts))
        for (i in 1:length(newdata)) B[i, which(districts == 
            newdata[i])] = 1
        P = diag(nrow = dim(B)[2])
    }
    else if (type == "ridge") {
        B = matrix(NA, nrow = dim(newdata)[1], ncol = dim(x)[2])
        for (i in 1:dim(x)[2]) B[, i] = x[, i]
        P = diag(nrow = dim(B)[2])
    }
    else if (type == "special") {
        if (is.na(B) || is.na(P)) 
            stop("In 'special' case: Base and Penalty matrix have to be specified.")
        if (is.vector(newdata)) 
            l = length(newdata)
        else if (is.matrix(newdata) || is.data.frame(newdata)) 
            l = nrow(newdata)
        B = as.matrix(object$B)[1:l, , drop = FALSE]
        P = as.matrix(P)
    }
    else if (type == "parametric") {
        x = data.frame(1, x)
        newdata = data.frame(1, newdata)
        B = model.matrix(formula(newdata), newdata)[, -1, drop = FALSE]
        P = matrix(0, nrow = ncol(B), ncol = ncol(B))
    }
    B
}
