lcurve <-
function (yy, B, quantile, DD, nb, constmat, types) 
{
    las1 = seq(-2, 4, by = 0.5)
    glatterms = which(types != "parametric")
    lambdas = matrix(las1, nrow = length(las1) * length(glatterms), 
        ncol = 1)
    if (length(glatterms) > 1) 
        for (i in 2:length(glatterms)) lambdas = cbind(lambdas, 
            rep(las1, each = i, times = length(glatterms) - i + 
                1))
    omega = rep(0, nrow(lambdas))
    theta = matrix(NA, nrow(lambdas), ncol(lambdas))
    lambdas = 10^lambdas
    penalty = rep(0, length(types))
    for (i in 1:nrow(lambdas)) {
        penalty[glatterms] = lambdas[i, ]
        aa <- asyregpen.lsfit(yy, B, quantile, abs(penalty), 
            DD, nb, constmat)
        omega[i] = sum(aa$weight * (yy - aa$fitted)^2)
        for (j in 1:ncol(lambdas)) {
            partbasis = (sum(nb[0:(glatterms[j] - 1)]) + 1):(sum(nb[0:glatterms[j]]))
            if (all(B[, 1] == 1)) {
                partDD = DD[, -1, drop = FALSE][-1, , drop = FALSE][, 
                  partbasis, drop = FALSE]
                partaa = aa$a[-1][partbasis]
            }
            else {
                partDD = DD[, partbasis, drop = FALSE]
                partaa = aa$a[partbasis]
            }
            theta[i, j] = t(partaa) %*% t(partDD) %*% partDD %*% 
                partaa
        }
    }
    omega = log(omega)
    theta = log(theta)
    matplot(omega, theta, type = "l")
    deltas = 2 * (sqrt(sum((theta[1, ] - theta[2, ])^2) + (omega[1] - 
        omega[2])^2))
    for (i in 2:(nrow(lambdas) - 1)) {
        deltas[i] = sqrt(sum((theta[i, ] - theta[i - 1, ])^2) + 
            (omega[i] - omega[i - 1])^2)
    }
    deltas[nrow(lambdas)] = 2 * (sqrt(sum((theta[i, ] - theta[i - 
        1, ])^2) + (omega[i] - omega[i - 1])^2))
    print(deltas)
    deltas = sqrt(diff(omega)^2 + rowSums(diff(theta))^2)
    penalty[glatterms] = lambdas[which.min(deltas) + 1, ]
    penalty
}
