TangleInferenceGauss <-
function (data, e, g, ...) 
{
    nc = ncol(data)
    thetaOpt = function(fre, e, g, theta0 = c(0, 10, 10, 1), 
        type = c(1, 1, 1, 1)) {
        cor = cor(fre)
        D_E = as.matrix(dist(e)/max(dist(e)))
        D_G = as.matrix(dist(g)/max(dist(g)))
        n = ncol(cor)
        I = diag(n)
        type = as.logical(type)
        infind = theta0 == Inf
        theta0[infind] = 10^8
        type[infind] = FALSE
        ssc = function(theta) {
            th = theta0
            th[type] = theta
            d = th[1]
            be = th[2]
            bg = th[3]
            g = th[4]
            P = D_E/be + D_G/bg
            tcor = (1 - d) * exp(-P^g) + d * I
            ss = sum((cor - tcor)^2)
            return(ss)
        }
        dssc = function(theta) {
            th = theta0
            th[type] = theta
            d = th[1]
            be = th[2]
            bg = th[3]
            g = th[4]
            P = D_E/be + D_G/bg
            ds1 = (sum(2 * (cor - (1 - d) * exp(-P^g) - d * I) * 
                (exp(-P^g) - I)))
            ds2 = -2 * (cor - (1 - d) * exp(-P^g) - d * I) * 
                (1 - d) * P^g * g * D_E * exp(-P^g)/(be^2 * P)
            ds2 = (sum(ds2[!I]))
            ds3 = -2 * (cor - (1 - d) * exp(-P^g) - d * I) * 
                (1 - d) * P^g * g * D_G * exp(-P^g)/(bg^2 * P)
            ds3 = (sum(ds3[!I]))
            ds4 = 2 * (cor - (1 - d) * exp(-P^g) - d * I) * (1 - 
                d) * P^g * log(P) * exp(-P^g)
            ds4 = (sum(ds4[!I]))
            dss = c(ds1, ds2, ds3, ds4)[type]
            return(dss)
        }
        lb = c(0, 1e-07, 1e-07, 0)[type]
        ub = c(1, Inf, Inf, ifelse(sum(infind) > 0, 2, 1))[type]
        op = optim(theta0[type], ssc, dssc, lower = lb, upper = ub, 
            method = "L-BFGS-B", control = list(factr = 1e-10))
        theta = theta0
        theta[type] = op$par
        theta[infind] = Inf
        names(theta) = c("delta", "beta_E", "beta_G", "gamma")
        return(theta)
    }
    predError = function(data, e, g, theta, pred = 1) {
        order = 1:nc
        order[c(1, pred)] = c(pred, 1)
        data = data[, order]
        e = e[order]
        g = g[order, ]
        D_E = as.matrix(dist(e)/max(dist(e)))
        D_G = as.matrix(dist(g)/max(dist(g)))
        delta = theta[1]
        beta_E = theta[2]
        beta_G = theta[3]
        gamma = theta[4]
        S = (1 - delta) * exp(-(D_E/beta_E + D_G/beta_G)^gamma) + 
            delta * diag(nc)
        x2 = t(as.matrix(data[, -1]))
        S11 = S[1, 1]
        S12 = S[1, -1]
        S21 = S[-1, 1]
        S22 = S[-1, -1]
        mu1 = mu2 = mean(x2)
        x1 = mu1 + S12 %*% solve(S22) %*% (x2 - mu2)
        error = sum((x1 - data[, 1])^2)/nrow(data)
        return(error)
    }
    MSE = rep(0, 3)
    names(MSE) = c("e+g", "e", "g")
    for (model in 1:3) {
        g0 = ifelse(model == 2, Inf, 50)
        e0 = ifelse(model == 3, Inf, 50)
        SE = 0
        for (k in 1:nc) {
            theta = thetaOpt(data[, -k], e[-k], g[-k, ], theta0 = c(0, 
                e0, g0, 1))
            error = predError(data, e, g, theta, pred = k)
            SE = SE + error
        }
        MSE[model] = SE/nc
    }
    return(MSE)
}
