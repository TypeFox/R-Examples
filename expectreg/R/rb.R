rb <-
function (x, type = c("pspline", "2dspline", "markov", "radial", 
    "krig", "random", "ridge", "special", "parametric"), B = NA, 
    P = NA, bnd = NA, center = TRUE, by = NA) 
{
    type = match.arg(type)
    Zspathelp = NA
    phi = NA
    xname = deparse(as.list(match.call())$x)
    if (type == "pspline") {
        B.deg = 2
        B.size = 20
        diff.size = 2
        x0 <- min(x, na.rm = TRUE) - 0.001
        x1 <- max(x, na.rm = TRUE) + 0.001
        dx = (x1 - x0)/(B.size - 1)
        B = matrix(0, nrow = length(x), ncol = B.size + B.deg - 
            1)
        notnas = which(!is.na(x))
        B[notnas, ] = splineDesign(knots = seq(x0 - dx * B.deg, 
            x1 + dx * B.deg, by = dx), x = x[notnas], ord = B.deg + 
            1)
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
        Bx = matrix(0, nrow = nrow(x), ncol = B.size)
        notnas = which(!is.na(x[, 1]))
        Bx[notnas, ] = splineDesign(knots = seq(x0 - dx * B.deg, 
            x1 + dx * B.deg, by = dx), x = x[notnas, 1], ord = B.deg + 
            1)
        y0 <- min(x[, 2], na.rm = TRUE) - 0.001
        y1 <- max(x[, 2], na.rm = TRUE) + 0.001
        dy = (y1 - y0)/(B.size - B.deg)
        By = matrix(0, nrow = nrow(x), ncol = B.size)
        notnas = which(!is.na(x[, 2]))
        By[notnas, ] = splineDesign(knots = seq(y0 - dy * B.deg, 
            y1 + dy * B.deg, by = dy), x = x[notnas, 2], ord = B.deg + 
            1)
        B = matrix(NA, nrow = dim(x)[1], ncol = dim(Bx)[2] * 
            dim(By)[2])
        for (i in 1:dim(Bx)[1]) B[i, ] = as.vector(Bx[i, ] %o% 
            By[i, ])
        D = diff(diag(dim(Bx)[2]), diff = diff.size)
        P = t(D) %*% D
        P = diag(dim(Bx)[2]) %x% P + P %x% diag(dim(Bx)[2])
        ek = eigen(P)
        ek$values = ek$values[1:(length(ek$values) - 2 * diff.size)]
        ek$vectors = ek$vectors[, 1:(dim(ek$vectors)[2] - 2 * 
            diff.size)]
        P = t(ek$vectors %*% sqrt(diag(ek$values)))
        if (center) {
            tildeX = matrix(1, nrow = dim(Bx)[2] * dim(By)[2], 
                ncol = 4)
            tildeX[, 2] = rep(1:dim(Bx)[2], times = dim(By)[2])
            tildeX[, 3] = rep(1:dim(By)[2], each = dim(Bx)[2])
            tildeX[, 4] = tildeX[, 2] * tildeX[, 3]
            tildeU = t(P) %*% solve(P %*% t(P))
            X = B %*% tildeX
            U = B %*% tildeU
            B = cbind(X[, -1], U)
            P <- diag(c(rep(0, ncol(X) - 1), rep(1, ncol(U))))
        }
    }
    else if (type == "radial") {
        x = x[order(x[, 1]), ]
        require(fields)
        knots = unique(x)
        knots <- cover.design(R = knots, nd = min(50, nrow(knots)))$design
        Zspathelp = knots
        B = matrix(NA, nrow = dim(x)[1], ncol = dim(knots)[1])
        for (j in 1:dim(knots)[1]) {
            r = sqrt(rowSums((x - matrix(unlist(knots[j, ]), 
                nrow = nrow(x), ncol = ncol(knots), byrow = T))^2))
            r[r == 0] = 1
            B[, j] = r^2 * log(r)
        }
        P = matrix(0, nrow = dim(B)[2], ncol = dim(B)[2])
        for (j in 1:dim(B)[2]) {
            r = sqrt(rowSums((matrix(unlist(knots), nrow = nrow(knots), 
                ncol = ncol(knots)) - matrix(unlist(knots[j, 
                ]), nrow = nrow(knots), ncol = ncol(knots), byrow = T))^2))
            r[r == 0] = 1
            P[, j] = r^2 * log(r)
        }
    }
    else if (type == "krig") {
        cons = 9.233
        x = x[order(x[, 1]), ]
        require(fields)
        knots = unique(x)
        knots <- cover.design(R = knots, nd = min(50, nrow(knots)))$design
        Zspathelp = knots
        B = matrix(NA, nrow = dim(x)[1], ncol = dim(knots)[1])
        P = matrix(0, nrow = dim(B)[2], ncol = dim(B)[2])
        for (j in 1:dim(B)[2]) {
            P[, j] = sqrt(rowSums((matrix(unlist(knots), nrow = nrow(knots), 
                ncol = ncol(knots)) - matrix(unlist(knots[j, 
                ]), nrow = nrow(knots), ncol = ncol(knots), byrow = T))^2))
        }
        phi = max(P)/cons
        P = P/phi
        P = exp(-P) * (1 + P)
        for (j in 1:dim(knots)[1]) {
            r = sqrt(rowSums((x - matrix(unlist(knots[j, ]), 
                nrow = nrow(x), ncol = ncol(knots), byrow = T))^2))
            B[, j] = exp(-r/phi) * (1 + r/phi)
        }
    }
    else if (type == "markov") {
        if (any(!is.na(bnd)) && any(is.na(P))) 
            P = bnd2gra(bnd)
        if (all(is.na(P))) 
            stop("No Neighbourhood defined.")
        if (any(diag(P) < 0) || any(P - diag(P) > 0) || is.null(dimnames(P)[[1]])) 
            stop("Maldefined Neighbourhood.")
        districts = dimnames(P)[[1]]
        B = matrix(0, nrow = length(x), ncol = dim(P)[2])
        for (i in 1:length(x)) B[i, which(districts == x[i])] = 1
        Zspathelp = diag(ncol(P))
        e <- eigen(P)
        if (center) {
            tildeU = e$vectors[, e$values < 1e-05]
            tildeZ = e$vectors[, e$values >= 1e-05] %*% diag(1/sqrt(e$values[e$values >= 
                1e-05]))
            U <- B %*% tildeU
            Z <- B %*% tildeZ
            B = cbind(U[, -1], Z)
            P <- diag(c(rep(0, ncol(U) - 1), rep(1, ncol(Z))))
            Zspathelp = cbind(tildeU, tildeZ)[, -1]
        }
        else P <- t(e$vectors[, -dim(e$vectors)[2]] * sqrt(abs(e$values[-length(e$values)])))
    }
    else if (type == "random") {
        districts = sort(unique(x))
        B = matrix(0, nrow = length(x), ncol = length(districts))
        for (i in 1:length(x)) B[i, which(districts == x[i])] = 1
        P = diag(nrow = dim(B)[2])
    }
    else if (type == "ridge") {
        B = matrix(NA, nrow = dim(x)[1], ncol = dim(x)[2])
        for (i in 1:dim(x)[2]) B[, i] = x[, i]
        P = diag(nrow = dim(B)[2])
    }
    else if (type == "special") {
        if (is.na(B) || is.na(P)) 
            stop("In 'special' case: Base and Penalty matrix have to be specified.")
        B = as.matrix(B)
        P = as.matrix(P)
    }
    else if (type == "parametric") {
        if (class(x) == "matrix") 
            xname = colnames(x)
        xdata = data.frame(1, x)
        names(xdata) = c("X1", xname)
        B = model.matrix(formula(xdata), xdata)[, -1, drop = FALSE]
        P = matrix(0, nrow = ncol(B), ncol = ncol(B))
    }
    constraint = matrix(0, nrow = 2, ncol = ncol(P))
    rb = list(B = B, P = P, x = x, type = type, bnd = bnd, Zspathelp = Zspathelp, 
        phi = phi, center = center, by = by, xname = xname, constraint = constraint)
    class(rb) = c("regbase")
    rb
}
