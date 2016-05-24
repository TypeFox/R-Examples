complement <-
function (x, y = NULL, alpha) 
{
    if (alpha < 0) {
        stop("Parameter alpha must be greater or equal to zero")
    }
    if (!inherits(x, "delvor")) {
        dd.obj <- delvor(x, y)
    }
    else {
        dd.obj <- x
    }
    xy.data <- dd.obj$x
    mesh <- dd.obj$mesh
    dm1 <- sqrt((mesh[, "x1"] - mesh[, "mx1"])^2 + (mesh[, "y1"] - 
        mesh[, "my1"])^2)
    dm2 <- sqrt((mesh[, "x1"] - mesh[, "mx2"])^2 + (mesh[, "y1"] - 
        mesh[, "my2"])^2)
    pm.x <- (mesh[, "x1"] + mesh[, "x2"]) * 0.5
    pm.y <- (mesh[, "y1"] + mesh[, "y2"]) * 0.5
    dm <- sqrt((mesh[, "x1"] - mesh[, "x2"])^2 + (mesh[, "y1"] - 
        mesh[, "y2"])^2) * 0.5
    d.pm.m1 <- sqrt((mesh[, "mx1"] - pm.x)^2 + (mesh[, "my1"] - 
        pm.y)^2)
    d.pm.m2 <- sqrt((mesh[, "mx2"] - pm.x)^2 + (mesh[, "my2"] - 
        pm.y)^2)
    theta.m1 <- atan(dm/d.pm.m1)
    theta.m2 <- atan(dm/d.pm.m2)
    v.x.m1 <- ifelse(d.pm.m1 != 0, (pm.x - mesh[, "mx1"])/d.pm.m1, 
        0)
    v.y.m1 <- ifelse(d.pm.m1 != 0, (pm.y - mesh[, "my1"])/d.pm.m1, 
        0)
    v.x.m2 <- ifelse(d.pm.m2 != 0, (pm.x - mesh[, "mx2"])/d.pm.m2, 
        0)
    v.y.m2 <- ifelse(d.pm.m2 != 0, (pm.y - mesh[, "my2"])/d.pm.m2, 
        0)
    thetam.m1 = acos(v.x.m1)
    thetam.m2 = acos(v.x.m2)
    n.edges <- dim(mesh)[1]
    betw = numeric(length = n.edges)
    vert <- mesh[, "mx1"] == mesh[, "mx2"]
    if (sum(vert) > 0) {
        betw[vert] <- apply(cbind(mesh[vert, "my1"], mesh[vert, 
            "my2"], pm.y[vert]), 1, rank)[3, ] == 2
    }
    if (sum(!vert) > 0) {
        betw[!vert] <- apply(cbind(mesh[!vert, "mx1"], mesh[!vert, 
            "mx2"], pm.x[!vert]), 1, rank)[3, ] == 2
    }
    aux <- alpha^2 - dm^2
    comp1.1 <- NULL
    case1.1 <- mesh[, "bp1"] == 0 & dm1 > alpha
    n.case1.1 <- sum(case1.1)
    if (n.case1.1 > 0) {
        comp1.1 <- cbind(mesh[case1.1, "mx1"], mesh[case1.1, 
            "my1"], dm1[case1.1], matrix(mesh[case1.1, ], nrow = n.case1.1), 
            1, v.x.m1[case1.1], v.y.m1[case1.1], theta.m1[case1.1])
    }
    comp1.2 <- NULL
    case1.2 <- (mesh[, "bp1"] == 0 & dm1 > alpha & aux > 0 & 
        betw == 1) | (mesh[, "bp1"] == 1 & aux > 0 & betw == 
        1)
    n.case1.2 <- sum(case1.2)
    if (n.case1.2 > 0) {
        a1 <- sqrt(aux[case1.2])
        e.x <- pm.x[case1.2] - a1 * v.x.m1[case1.2]
        e.y <- pm.y[case1.2] - a1 * v.y.m1[case1.2]
        theta <- atan(dm[case1.2]/a1)
        comp1.2 <- cbind(e.x, e.y, alpha, matrix(mesh[case1.2, 
            ], nrow = n.case1.2), 1, v.x.m1[case1.2], v.y.m1[case1.2], 
            theta)
    }
    comp1.3.1 <- NULL
    case1.3 <- (mesh[, "bp1"] == 0 & dm1 > alpha & aux > 0 & 
        betw == 0) | (mesh[, "bp1"] == 1 & aux > 0 & betw == 
        0)
    n.case1.3 <- sum(case1.3)
    if (n.case1.3 > 0) {
        a1 <- sqrt(aux[case1.3])
        e.x <- pm.x[case1.3] - a1 * v.x.m1[case1.3]
        e.y <- pm.y[case1.3] - a1 * v.y.m1[case1.3]
        theta <- atan(dm[case1.3]/a1)
        case1.3.1 <- alpha - a1 < dm2[case1.3] - d.pm.m2[case1.3]
        n.case1.3.1 <- sum(case1.3.1)
        if (n.case1.3.1 > 0) {
            comp1.3.1 <- cbind(e.x[case1.3.1], e.y[case1.3.1], 
                alpha, matrix(matrix(mesh[case1.3, ], nrow = n.case1.3)[case1.3.1, 
                  ], nrow = n.case1.3.1), 1, v.x.m1[case1.3][case1.3.1], 
                v.y.m1[case1.3][case1.3.1], theta[case1.3.1])
        }
    }
    comp2.1 <- NULL
    case2.1 <- mesh[, "bp2"] == 0 & dm2 > alpha
    n.case2.1 <- sum(case2.1)
    if (n.case2.1 > 0) {
        comp2.1 <- cbind(mesh[case2.1, "mx2"], mesh[case2.1, 
            "my2"], dm2[case2.1], matrix(mesh[case2.1, ], nrow = n.case2.1), 
            2, v.x.m2[case2.1], v.y.m2[case2.1], theta.m2[case2.1])
    }
    comp2.2 <- NULL
    case2.2 <- (mesh[, "bp2"] == 0 & dm2 > alpha & aux > 0 & 
        betw == 1) | (mesh[, "bp2"] == 1 & aux > 0 & betw == 
        1)
    n.case2.2 <- sum(case2.2)
    if (n.case2.2 > 0) {
        a1 <- sqrt(aux[case2.2])
        e.x <- pm.x[case2.2] - a1 * v.x.m2[case2.2]
        e.y <- pm.y[case2.2] - a1 * v.y.m2[case2.2]
        theta <- atan(dm[case2.2]/a1)
        comp2.2 <- cbind(e.x, e.y, alpha, matrix(mesh[case2.2, 
            ], nrow = n.case2.2), 2, v.x.m2[case2.2], v.y.m2[case2.2], 
            theta)
    }
    comp2.3.1 <- NULL
    case2.3 <- (mesh[, "bp2"] == 0 & dm2 > alpha & aux > 0 & 
        betw == 0) | (mesh[, "bp2"] == 1 & aux > 0 & betw == 
        0)
    n.case2.3 <- sum(case2.3)
    if (n.case2.3 > 0) {
        a1 <- sqrt(aux[case2.3])
        e.x <- pm.x[case2.3] - a1 * v.x.m2[case2.3]
        e.y <- pm.y[case2.3] - a1 * v.y.m2[case2.3]
        theta <- atan(dm[case2.3]/a1)
        case2.3.1 <- alpha - a1 < dm1[case2.3] - d.pm.m1[case2.3]
        n.case2.3.1 <- sum(case2.3.1)
        if (n.case2.3.1 > 0) {
            comp2.3.1 <- cbind(e.x[case2.3.1], e.y[case2.3.1], 
                alpha, matrix(matrix(mesh[case2.3, ], nrow = n.case2.3)[case2.3.1, 
                  ], nrow = n.case2.3.1), 2, v.x.m2[case2.3][case2.3.1], 
                v.y.m2[case2.3][case2.3.1], theta[case2.3.1])
        }
    }
    comp.h <- NULL
    case.h <- mesh[, "bp2"] == 1
    n.case.h <- sum(case.h)
    if (n.case.h > 0) {
        comp.h <- matrix(nrow = n.case.h, ncol = 19)
        index <- which(case.h)
        for (i in 1:n.case.h) {
            if (mesh[index[i], "x1"] == mesh[index[i], "x2"]) {
                a <- mesh[index[i], "x1"]
                if (mesh[index[i], "mx2"] <= a) {
                  sig = -4
                }
                else {
                  sig = -3
                }
                comp.h[i, ] = c(a, 0, sig, mesh[index[i], ], 
                  2, 0, 0, 0)
            }
            else {
                b <- (mesh[index[i], "y2"] - mesh[index[i], "y1"])/(mesh[index[i], 
                  "x2"] - mesh[index[i], "x1"])
                a <- mesh[index[i], "y1"] - mesh[index[i], "x1"] * 
                  b
                if (mesh[index[i], "my2"] <= a + b * mesh[index[i], 
                  "mx2"]) {
                  sig = -2
                }
                else {
                  sig = -1
                }
                comp.h[i, ] = c(a, b, sig, mesh[index[i], ], 
                  2, 0, 0, 0)
            }
        }
    }
    comp <- rbind(comp1.1, comp1.2, comp1.3.1, comp2.1, comp2.2, 
        comp2.3.1, comp.h)
    colnames(comp) <- c("c1", "c2", "r", colnames(mesh), "ind", 
        "v.x", "v.y", "theta")
    invisible(comp)
}
