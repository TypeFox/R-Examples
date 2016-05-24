mono <-
function (x, constraint = c("increase", "decrease", "convex", 
    "concave", "flatend"), by = NA) 
{
    constraint = match.arg(constraint)
    xname = deparse(as.list(match.call())$x)
    B.deg = 2
    B.size = 20
    diff.size = 2
    x0 <- min(x, na.rm = TRUE) - 0.001
    x1 <- max(x, na.rm = TRUE) + 0.001
    dx = (x1 - x0)/(B.size - 1)
    B = matrix(0, nrow = length(x), ncol = B.size + 1)
    notnas = which(!is.na(x))
    B[notnas, ] = splineDesign(knots = seq(x0 - dx * B.deg, x1 + 
        dx * B.deg, by = dx), x = x[notnas], ord = B.deg + 1)
    P <- diag(dim(B)[2])
    P <- diff(P, diff = diff.size)
    if (constraint == "increase") {
        Pc = diff(diag(dim(B)[2]))
    }
    else if (constraint == "decrease") {
        Pc = -diff(diag(dim(B)[2]))
    }
    else if (constraint == "convex") {
        Pc = diff(diag(dim(B)[2]), diff = 2)
    }
    else if (constraint == "concave") {
        Pc = -diff(diag(dim(B)[2]), diff = 2)
    }
    else if (constraint == "flatend") {
        nflat = 2
        D1 = diff(diag(dim(B)[2]))
        v = rep(0, ncol(B) - 1)
        v[1:nflat] = 1
        v[ncol(B) - (1:nflat) + 1] = 1
        V = diag(v)
        Pc = V %*% D1
        Pc = rbind(Pc, -Pc)
    }
    rb = list(B = B, P = P, x = x, type = "pspline", bnd = NA, 
        Zspathelp = NA, phi = NA, center = FALSE, by = by, xname = xname, 
        constraint = Pc)
    class(rb) = c("regbase", "rbmono")
    rb
}
