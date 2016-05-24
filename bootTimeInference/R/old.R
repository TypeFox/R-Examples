#' @useDynLib bootTimeInference
#' @importFrom Rcpp sourceCpp
#' @importFrom stats ar lm pnorm rgeom var
NULL


block.size.calibrate <-
function (ret, b.vec = c(1, 3, 6, 10), alpha = 0.05, M = 199,
    K = 1000, b.av = 5, T.start = 50)
{
    b.len = length(b.vec)
    emp.reject.probs = rep(0, b.len)
    Delta.hat = sharpe.ratio.diff(ret)
    ret1 = ret[, 1]
    ret2 = ret[, 2]
    T = length(ret1)
    Var.data = matrix(0, T.start + T, 2)
    Var.data[1, ] = ret[1, ]
    Delta.hat = sharpe.ratio.diff(ret)
    fit1 = lm(ret1[2:T] ~ ret1[1:(T - 1)] + ret2[1:(T - 1)])
    fit2 = lm(ret2[2:T] ~ ret1[1:(T - 1)] + ret2[1:(T - 1)])
    coef1 = as.numeric(fit1$coef)
    coef2 = as.numeric(fit2$coef)
    resid.mat = cbind(as.numeric(fit1$resid), as.numeric(fit2$resid))
    for (k in (1:K)) {
        resid.mat.star = rbind(c(0, 0), resid.mat[sb.sequence(T -
            1, b.av, T.start + T - 1), ])
        # print(resid.mat.star)
        for (t in (2:(T.start + T))) {
            Var.data[t, 1] = coef1[1] + coef1[2] * Var.data[t -
                1, 1] + coef1[3] * Var.data[t - 1, 2] + resid.mat.star[t,
                1]
            Var.data[t, 2] = coef2[1] + coef2[2] * Var.data[t -
                1, 1] + coef2[3] * Var.data[t - 1, 2] + resid.mat.star[t,
                2]
        }
        # print(Var.data)
        Var.data.trunc = Var.data[(T.start + 1):(T.start + T),
            ]
        for (j in (1:b.len)) {
            p.Value = boot.time.inference(Var.data.trunc, b.vec[j],
                M, Delta.hat)$p.Value
            print(p.Value)
            if (p.Value <= alpha) {
                emp.reject.probs[j] = emp.reject.probs[j] + 1
            }
        }
    }
    print(emp.reject.probs)
    emp.reject.probs = emp.reject.probs/K
    print(emp.reject.probs)
    b.order = order(abs(emp.reject.probs - alpha))
    print(b.order)
    b.opt = b.vec[b.order[1]]
    print(b.opt)
    b.vec.with.probs = rbind(b.vec, emp.reject.probs)
    colnames(b.vec.with.probs) = rep("", length(b.vec))
    list(Empirical.Rejection.Probs = b.vec.with.probs, b.optimal = b.opt)
}
boot.time.inference <-
function (ret, b, M, Delta.null = 0, digits = 4)
{
    T = length(ret[, 1])
    l = floor(T/b)
    Delta.hat = sharpe.ratio.diff(ret)
    d = abs(Delta.hat - Delta.null)/compute.se.Parzen.pw(ret)
    p.value = 1
    for (m in (1:M)) {
        ret.star = ret[cbb.sequence(T, b), ]
        Delta.hat.star = sharpe.ratio.diff(ret.star)
        ret1.star = ret.star[, 1]
        ret2.star = ret.star[, 2]
        mu1.hat.star = mean(ret1.star)
        mu2.hat.star = mean(ret2.star)
        gamma1.hat.star = mean(ret1.star^2)
        gamma2.hat.star = mean(ret2.star^2)
        gradient = rep(0, 4)
        gradient[1] = gamma1.hat.star/(gamma1.hat.star - mu1.hat.star^2)^1.5
        gradient[2] = -gamma2.hat.star/(gamma2.hat.star - mu2.hat.star^2)^1.5
        gradient[3] = -0.5 * mu1.hat.star/(gamma1.hat.star -
                                           mu1.hat.star^2)^1.5
        gradient[4] = 0.5 * mu2.hat.star/(gamma2.hat.star - mu2.hat.star^2)^1.5
        y.star = data.frame(ret1.star - mu1.hat.star, ret2.star -
            mu2.hat.star, ret1.star^2 - gamma1.hat.star, ret2.star^2 -
            gamma2.hat.star)
        Psi.hat.star = matrix(0, 4, 4)
        for (j in (1:l)) {
            zeta.star = b^0.5 * colMeans(y.star[((j - 1) * b + 1):(j *
                b), ])
            Psi.hat.star = Psi.hat.star + zeta.star %*% t(zeta.star)
        }
        Psi.hat.star = Psi.hat.star/l
        se.star = as.numeric(sqrt(t(gradient) %*% Psi.hat.star %*%
            gradient/T))
        d.star = abs(Delta.hat.star - Delta.hat)/se.star
        if (d.star >= d) {
            p.value = p.value + 1
        }
    }
    p.value = p.value/(M + 1)
    list(Difference = round(Delta.hat, digits), p.Value = round(p.value,
                                                    digits))
}
cbb.sequence <-
structure(function (T, b)
{
    l = floor(T/b)
    index.sequence = c(1:T, 1:b)
    sequence = rep(0, T)
    start.points = sample(1:T, l, replace = T)
    for (j in (1:l)) {
        start = start.points[j]
        sequence[((j - 1) * b + 1):(j * b)] = index.sequence[start:(start +
            b - 1)]
    }
    sequence
}, source = c("function(T, b){", "    l = floor(T / b)", "    index.sequence = c(1:T, 1:b)",
"    sequence = rep(0, T)", "    start.points = sample(1:T, l, replace = T)",
"    for (j in (1:l)) {", "      start = start.points[j]", "      sequence[((j-1)*b+1):(j*b)] = index.sequence[start:(start+b-1)]",
"    }", "    sequence", "  }"))
compute.alpha.hat <-
structure(function (V.hat)
{
    dimensions = dim(V.hat)
    T = dimensions[1]
    p = dimensions[2]
    numerator = 0
    denominator = 0
    for (i in (1:p)) {
        fit = ar(V.hat[, i], 0, 1, method = "ols")
        rho.hat = as.numeric(fit[2])
        sig.hat = sqrt(as.numeric(fit[3]))
        numerator = numerator + 4 * rho.hat^2 * sig.hat^4/(1 -
            rho.hat)^8
        denominator = denominator + sig.hat^4/(1 - rho.hat)^4
    }
    numerator/denominator
}, source = c("function(V.hat){", "    # see (6.4) of Andrews (1991)",
"    dimensions = dim(V.hat)", "    T = dimensions[1]", "    p = dimensions[2]",
"    numerator = 0", "    denominator = 0", "    for (i in (1:p)) {",
"      fit = ar(V.hat[, i], 0, 1, method = \"ols\")", "      rho.hat = as.numeric(fit[2])",
"      sig.hat = sqrt(as.numeric(fit[3]))", "      numerator = numerator + 4 * rho.hat^2 * sig.hat^4 / (1 - rho.hat)^8",
"      denominator = denominator + sig.hat^4 / (1 - rho.hat)^4",
"    }", "    numerator / denominator", "  }"))
compute.Gamma.hat <-
structure(function (V.hat, j)
{
    dimensions = dim(V.hat)
    T = dimensions[1]
    p = dimensions[2]
    Gamma.hat = matrix(0, p, p)
    if (j >= T)
        stop("j must be smaller than the row dimension!")
    for (i in ((j + 1):T)) Gamma.hat = Gamma.hat + V.hat[i, ] %*%
        t(V.hat[i - j, ])
    Gamma.hat = Gamma.hat/T
    Gamma.hat
}, source = c("function(V.hat, j){", "    # see second part of (2.5) of Andrews (1991)",
"    dimensions = dim(V.hat)", "    T = dimensions[1]", "    p = dimensions[2]",
"    Gamma.hat = matrix(0, p, p)", "    if (j >= T)", "      stop(\"j must be smaller than the row dimension!\")",
"    for (i in ((j+1):T)) ", "      Gamma.hat = Gamma.hat + V.hat[i,] %*% t(V.hat[i - j,])",
"    Gamma.hat = Gamma.hat / T", "    Gamma.hat", "  }"))
compute.Psi.hat <-
structure(function (V.hat)
{
    T = length(V.hat[, 1])
    alpha.hat = compute.alpha.hat(V.hat)
    S.star = 2.6614 * (alpha.hat * T)^0.2
    Psi.hat = compute.Gamma.hat(V.hat, 0)
    j = 1
    while (j < S.star) {
        Gamma.hat = compute.Gamma.hat(V.hat, j)
        Psi.hat = Psi.hat + kernel.Parzen(j/S.star) * (Gamma.hat +
            t(Gamma.hat))
        j = j + 1
    }
    Psi.hat = (T/(T - 4)) * Psi.hat
    Psi.hat
}, source = c("function(V.hat) {", "    # see first part of (2.5) of Andrews (1991)",
"    # except we call Psi.hat what he calls J.hat there", "    T = length(V.hat[,1])",
"    alpha.hat = compute.alpha.hat(V.hat)", "    S.star = 2.6614 * (alpha.hat * T)^.2",
"    Psi.hat = compute.Gamma.hat(V.hat, 0)", "    j = 1", "    while (j < S.star) {",
"      Gamma.hat = compute.Gamma.hat(V.hat, j)", "      Psi.hat = Psi.hat + kernel.Parzen(j / S.star) * (Gamma.hat + t(Gamma.hat))",
"      j = j + 1", "    }", "    Psi.hat = (T / (T - 4)) * Psi.hat",
"    Psi.hat", "  }"))
compute.se.Parzen <-
structure(function (ret)
{
    ret1 = ret[, 1]
    ret2 = ret[, 2]
    T = length(ret1)
    mu1.hat = mean(ret1)
    mu2.hat = mean(ret2)
    gamma1.hat = mean(ret1^2)
    gamma2.hat = mean(ret2^2)
    gradient = rep(0, 4)
    gradient[1] = gamma1.hat/(gamma1.hat - mu1.hat^2)^1.5
    gradient[2] = -gamma2.hat/(gamma2.hat - mu2.hat^2)^1.5
    gradient[3] = -0.5 * mu1.hat/(gamma1.hat - mu1.hat^2)^1.5
    gradient[4] = 0.5 * mu2.hat/(gamma2.hat - mu2.hat^2)^1.5
    V.hat = compute.V.hat(ret)
    Psi.hat = compute.Psi.hat(V.hat)
    se = as.numeric(sqrt(t(gradient) %*% Psi.hat %*% gradient/T))
    se
}, source = c("function(ret){", "    # implements the Parzen kernel estimator with automatic choice of bandwith",
"    ret1 = ret[,1]", "    ret2 = ret[,2]", "    T = length(ret1)",
"    mu1.hat = mean(ret1)", "    mu2.hat = mean(ret2)", "    gamma1.hat = mean(ret1^2)",
"    gamma2.hat = mean(ret2^2)", "    gradient = rep(0, 4)",
"    gradient[1] = gamma1.hat / (gamma1.hat - mu1.hat^2)^1.5",
"    gradient[2] = - gamma2.hat / (gamma2.hat - mu2.hat^2)^1.5",
"    gradient[3] = -0.5 * mu1.hat / (gamma1.hat - mu1.hat^2)^1.5",
"    gradient[4] = 0.5 * mu2.hat / (gamma2.hat - mu2.hat^2)^1.5",
"    V.hat = compute.V.hat(ret)", "    Psi.hat = compute.Psi.hat(V.hat)",
"    se = as.numeric(sqrt(t(gradient) %*% Psi.hat %*% gradient / T))",
"    se", "  }"))
compute.se.Parzen.pw <-
structure(function (ret)
{
    ret1 = ret[, 1]
    ret2 = ret[, 2]
    mu1.hat = mean(ret1)
    mu2.hat = mean(ret2)
    gamma1.hat = mean(ret1^2)
    gamma2.hat = mean(ret2^2)
    gradient = rep(0, 4)
    gradient[1] = gamma1.hat/(gamma1.hat - mu1.hat^2)^1.5
    gradient[2] = -gamma2.hat/(gamma2.hat - mu2.hat^2)^1.5
    gradient[3] = -0.5 * mu1.hat/(gamma1.hat - mu1.hat^2)^1.5
    gradient[4] = 0.5 * mu2.hat/(gamma2.hat - mu2.hat^2)^1.5
    T = length(ret1)
    V.hat = compute.V.hat(ret)
    A.ls = matrix(0, 4, 4)
    V.star = matrix(0, T - 1, 4)
    reg1 = V.hat[1:T - 1, 1]
    reg2 = V.hat[1:T - 1, 2]
    reg3 = V.hat[1:T - 1, 3]
    reg4 = V.hat[1:T - 1, 4]
    for (j in (1:4)) {
        fit = lm(V.hat[2:T, j] ~ -1 + reg1 + reg2 + reg3 + reg4)
        A.ls[j, ] = as.numeric(fit$coef)
        V.star[, j] = as.numeric(fit$resid)
    }
    svd.A = svd(A.ls)
    d = svd.A$d
    d.adj = d
    for (i in (1:4)) {
        if (d[i] > 0.97)
            d.adj[i] = 0.97
        else if (d[i] < -0.97)
            d.adj[i] = -0.97
    }
    A.hat = svd.A$u %*% diag(d.adj) %*% t(svd.A$v)
    D = solve(diag(4) - A.hat)
    reg.mat = rbind(reg1, reg2, reg3, reg4)
    for (j in (1:4)) {
        V.star[, j] = V.hat[2:T, j] - A.hat[j, ] %*% reg.mat
    }
    Psi.hat = compute.Psi.hat(V.star)
    Psi.hat = D %*% Psi.hat %*% t(D)
    se = as.numeric(sqrt(t(gradient) %*% Psi.hat %*% gradient/T))
    se
}, source = c("function(ret){", "    # implements the prewhitened Parzen kernel estimator of A&M (1992)",
"    ret1 = ret[,1]", "    ret2 = ret[,2]", "    mu1.hat = mean(ret1)",
"    mu2.hat = mean(ret2)", "    gamma1.hat = mean(ret1^2)",
"    gamma2.hat = mean(ret2^2)", "    gradient = rep(0, 4)",
"    gradient[1] = gamma1.hat / (gamma1.hat - mu1.hat^2)^1.5",
"    gradient[2] = - gamma2.hat / (gamma2.hat - mu2.hat^2)^1.5",
"    gradient[3] = -0.5 * mu1.hat / (gamma1.hat - mu1.hat^2)^1.5",
"    gradient[4] = 0.5 * mu2.hat / (gamma2.hat - mu2.hat^2)^1.5    ",
"    T = length(ret1)", "    V.hat = compute.V.hat(ret)", "    A.ls = matrix(0, 4, 4)",
"    V.star = matrix(0, T - 1, 4)", "    reg1 = V.hat[1:T-1,1]",
"    reg2 = V.hat[1:T-1,2]", "    reg3 = V.hat[1:T-1,3]", "    reg4 = V.hat[1:T-1,4]",
"    for (j in (1:4)) {", "      fit = lm(V.hat[2:T,j] ~ -1 + reg1 + reg2 + reg3 + reg4)",
"      A.ls[j,] = as.numeric(fit$coef)", "      V.star[,j] = as.numeric(fit$resid)",
"    }", "    # SVD adjustment of A&M (1992, page 957)", "    svd.A = svd(A.ls)",
"    d = svd.A$d", "    d.adj = d", "    for (i in (1:4)) {",
"      if (d[i] > 0.97)", "        d.adj[i] = 0.97", "      else if (d[i] < -0.97)",
"        d.adj[i] = -0.97", "    }", "    A.hat = svd.A$u %*% diag(d.adj) %*% t(svd.A$v)",
"    D = solve(diag(4) - A.hat)", "    reg.mat = rbind(reg1, reg2, reg3, reg4)",
"    for (j in (1:4)) {", "      V.star[,j] = V.hat[2:T,j] - A.hat[j,] %*% reg.mat",
"    }", "    Psi.hat = compute.Psi.hat(V.star)", "    Psi.hat = D %*% Psi.hat %*% t(D)",
"    se = as.numeric(sqrt(t(gradient) %*% Psi.hat %*% gradient / T))",
"    se", "  }"))
compute.V.hat <-
structure(function (ret)
{
    ret1 = ret[, 1]
    ret2 = ret[, 2]
    V.hat = cbind(ret1 - mean(ret1), ret2 - mean(ret2), ret1^2 -
        mean(ret1^2), ret2^2 - mean(ret2^2))
    V.hat
}, source = c("function(ret){", "    # what Andrews (1991) calls V.hat = V(theta.hat) in our context",
"    ret1 = ret[,1]", "    ret2 = ret[,2]", "    V.hat = cbind(ret1 - mean(ret1), ret2 - mean(ret2),",
"      ret1^2 - mean(ret1^2), ret2^2 - mean(ret2^2))", "    V.hat",
"  }"))
hac.inference <-
structure(function (ret, digits = 3)
{
    ret1 = ret[, 1]
    ret2 = ret[, 2]
    mu1.hat = mean(ret1)
    mu2.hat = mean(ret2)
    sig1.hat = var(ret1)^0.5
    sig2.hat = var(ret2)^0.5
    SR1.hat = mu1.hat/sig1.hat
    SR2.hat = mu2.hat/sig2.hat
    SRs = round(c(SR1.hat, SR2.hat), digits)
    diff = SR1.hat - SR2.hat
    names(SRs) = c("SR1.hat", "SR2.hat")
    Delta.hat = SR1.hat - SR2.hat
    se = compute.se.Parzen(ret)
    se.pw = compute.se.Parzen.pw(ret)
    SEs = round(c(se, se.pw), digits)
    names(SEs) = c("HAC", "HAC.pw")
    PV = 2 * pnorm(-abs(diff)/se)
    PV.pw = 2 * pnorm(-abs(diff)/se.pw)
    PVs = round(c(PV, PV.pw), digits)
    names(PVs) = c("HAC", "HAC.pw")
    list(Sharpe.Ratios = SRs, Difference = round(diff, digits),
        Standard.Errors = SEs, p.Values = PVs)
}, source = c("function(ret, digits = 3) {", "    ret1 = ret[,1]",
"    ret2 = ret[,2]", "    mu1.hat = mean(ret1)", "    mu2.hat = mean(ret2)",
"    sig1.hat = var(ret1)^.5", "    sig2.hat = var(ret2)^.5",
"    SR1.hat = mu1.hat / sig1.hat", "    SR2.hat = mu2.hat / sig2.hat",
"    SRs = round(c(SR1.hat, SR2.hat), digits)", "    diff = SR1.hat - SR2.hat",
"    names(SRs) = c(\"SR1.hat\", \"SR2.hat\")", "    Delta.hat = SR1.hat - SR2.hat",
"    se = compute.se.Parzen(ret)", "    se.pw = compute.se.Parzen.pw(ret)",
"    SEs = round(c(se, se.pw), digits)", "    names(SEs) = c(\"HAC\", \"HAC.pw\")",
"    PV = 2 * pnorm(-abs(diff) / se)", "    PV.pw = 2 * pnorm(-abs(diff)/ se.pw)",
"    PVs = round(c(PV, PV.pw), digits)", "    names(PVs) = c(\"HAC\", \"HAC.pw\")",
"    list(Sharpe.Ratios = SRs, Difference = round(diff, digits),",
"         Standard.Errors = SEs, p.Values = PVs)", "    ", "    ",
"  }"))
kernel.Parzen <-
structure(function (x)
{
    if (abs(x) <= 0.5)
        result = 1 - 6 * x^2 + 6 * abs(x)^3
    else if (abs(x) <= 1)
        result = 2 * (1 - abs(x))^3
    else result = 0
    result
}, source = c("function(x)", "{", "  if (abs(x) <= 0.5)", "    result = 1 - 6 * x^2 + 6 * abs(x)^3",
"  else if (abs(x) <= 1)", "    result = 2 * (1 - abs(x))^3",
"  else", "    result = 0", "  result", "}"))
ret.agg <-
structure(list(V1 = c(3.92670407, -2.33048071, -4.85832276, 1.95998592,
-0.43972737, -2.70563325, 2.92574256, 3.87814751, -2.52489141,
1.74590464, -4.06991179, 0.84012606, -0.56170517, 2.84344303,
3.11084836, 1.91652781, 1.2406743, 2.7005608, 4.01653946, 1.19308713,
2.44223324, -1.61503088, 3.6948135, 2.07018962, 2.0208975, 0.83481382,
1.44586084, 1.20524857, 1.51496434, 0.06133578, -4.84827095,
2.44977553, 4.55040666, 1.21709433, 5.4069905, -2.54873859, 3.72536208,
0.54105796, -5.43839809, 5.18967322, 4.95603241, 4.63932427,
7.65874884, -5.72423861, 4.80637372, -3.22547442, 3.90518504,
1.86219433, 0.10582716, 6.23278903, 4.7072399, 0.34272956, -1.02731916,
4.33753278, 0.03111618, -16.62503573, 3.97885471, 6.78552349,
5.69553948, 7.71206682, 3.2070641, -2.46220351, 3.89688785, 1.76306486,
-3.53632421, 5.06137686, -3.53519744, -1.90283675, -2.54655421,
3.85615188, 3.63409054, 9.75479521, -4.83117625, 2.90380609,
3.82397828, -6.45137781, -5.08712293, 5.19073696, -2.53089011,
6.15036353, -5.84708876, -0.82048211, -11.71065812, 2.1564888,
2.76751644, -10.24214404, -8.84161009, 12.06687524, 0.23469865,
-0.65065573, -5.39932906, -5.22087084, -10.37836022, 2.6288369,
5.75676617, 2.07200752, -1.60537242, -2.50126884, 2.24747676,
-6.00599452, -0.49649931, -6.5602266, -8.63037497, -0.04259396,
-10.26085977, 8.24590004, 4.83684118, -5.7832214, -3.33069179,
-0.82456276, 1.73950487, 7.23176297, 3.48935011, 0.94714694,
2.3977011, 1.60398836, -1.8021481, 5.55897751, 0.58601643, 5.52014288
), V2 = c(2.64759888, -0.81627938, -5.72403137, 0.04747615, -4.24156123,
-7.59122784, 3.33909569, 6.7505071, -0.84311657, 4.2812789, -4.4141597,
2.29524585, -3.21205831, 4.18452874, 3.66221478, 2.78664632,
3.5772847, 10.06257052, 10.75298768, 0.90615474, 1.90052051,
-1.6643143, -1.6861009, -5.91624393, 0.99665747, 5.4183424, -0.40340949,
5.81894847, 3.37124409, -4.26584128, -10.75271612, 3.69745056,
7.75711383, -3.00709376, 4.93609027, -3.787201, 6.3226772, -5.64273356,
-7.06975665, 1.7490358, 8.635159, 2.33156526, 9.29364687, -1.61818302,
5.63142548, -8.20496716, 0.50002987, 1.09178582, -1.07181348,
8.75431668, 5.11078556, 1.2007289, -4.2240417, 9.33901336, 0.4632305,
-21.92512391, 11.49516863, 2.9996214, 6.78450801, 12.22539193,
10.21974187, -5.9869284, 12.78140099, 6.38115222, -4.7807636,
8.73816721, -1.53231983, 5.41216931, -2.88542841, 10.47915002,
10.82574708, 16.93645248, -0.69995667, 15.91009419, -2.98171268,
-17.06131019, -12.18635617, 16.94839843, -7.93024496, 14.78038896,
-10.25417492, -12.19078838, -26.51730222, 6.52516368, 1.67180416,
-22.61619191, -23.94040367, 19.0282885, -3.57649912, -11.19311909,
-13.36661023, -15.10994731, -23.94610591, 15.68607237, 10.00690578,
0.01861309, -7.06579239, -10.80461355, 4.75810268, -10.2733745,
-7.13687313, -16.37565909, -7.77480883, -3.71392416, -14.69675717,
13.96299673, 13.48573121, -9.00858074, -0.6360189, -1.54886517,
1.80773149, 4.97737788, 7.06602213, 1.26386946, 2.70543732, 3.73190342,
-1.55509919, 5.91609349, 2.27753758, 1.82369758)), .Names = c("V1",
"V2"), class = "data.frame", row.names = c("1", "2", "3", "4",
"5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15",
"16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26",
"27", "28", "29", "30", "31", "32", "33", "34", "35", "36", "37",
"38", "39", "40", "41", "42", "43", "44", "45", "46", "47", "48",
"49", "50", "51", "52", "53", "54", "55", "56", "57", "58", "59",
"60", "61", "62", "63", "64", "65", "66", "67", "68", "69", "70",
"71", "72", "73", "74", "75", "76", "77", "78", "79", "80", "81",
"82", "83", "84", "85", "86", "87", "88", "89", "90", "91", "92",
"93", "94", "95", "96", "97", "98", "99", "100", "101", "102",
"103", "104", "105", "106", "107", "108", "109", "110", "111",
"112", "113", "114", "115", "116", "117", "118", "119", "120"
))
ret.hedge <-
structure(c(0.6812972, 0.6590867, 1.113184, 0.5253526, 0.8250787,
1.209827, 0.9580597, 0.7780394, 1.692381, 0.6913268, 1.639145,
0.2832465, 0.5421948, 1.11027, 1.30982, 0.2096562, 0.6543112,
1.133846, 1.609258, 2.566037, 1.247054, 1.541927, 0.994488, 1.817061,
2.101866, 1.331784, 1.880532, 1.607264, 2.703011, 1.566412, 1.070305,
0.8612195, 0.8680123, 1.201234, 0.9419864, 0.8329147, 1.590724,
1.660986, 0.9035943, 1.344896, 1.804685, 1.451666, 2.603538,
2.236834, 0.9381944, 1.947973, -0.08834136, 1.354882, 1.932465,
0.877896, 1.375012, 0.7309824, 1.013313, -1.987464, 0.8585362,
-2.286329, -2.110502, -0.04409754, 1.316254, 1.002578, 1.626264,
1.126939, 1.008629, 2.712589, 2.62778, 1.549161, 0.8417694, 0.6403194,
1.156887, 0.7974892, -1.04404, 1.773266, 2.5434, 1.491014, 3.027254,
3.299774, 2.097768, 1.509243, 1.634951, 0.8581319, 1.524086,
0.2827904, -1.000094, 1.512823, 4.973459, 3.659623, 1.236722,
0.7181669, 0.4607545, -0.01620412, 1.663385, 1.685307, 1.841467,
1.487927, 1.383737, 3.16509, 2.20604, 0.6141716, 0.9461478, 2.842833,
-0.1229929, -2.099501, -1.782211, 1.404185, 2.021849, 0.8841528,
3.086808, 2.964086, 3.159487, 1.73617, 1.099206, 2.287825, 1.350948,
-2.127188, -1.677494, 0.4995302, 3.352566, 2.090021, 0.888179,
2.686836, 0.1943571, 0.3118182, 0.79741, -0.06141619, 0.5678554,
0.2104327, 0.2946753, 0.5506035, 0.32118, 0.1259476, 0.3275226,
0.4815175, 0.1553337, 0.07108152, 0.3126823, 0.1501736, 1.127773,
0.7885659, 0.6145729, 0.1698496, -0.05948804, -0.01964648, -0.02591778,
0.240129, 0.2412821, 0.3441884, 0.08540721, -0.06317093, 0.4778987,
0.3831296, 0.0802441, 0.6240868, -0.05434315, 0.2818077, 0.298526,
0.2285909, 0.7542084, 0.300116, 0.2301754, 0.3675995, 0.177732,
0.1374963, 0.09813306, 0.200423, 0.2452992, 0.2735195, 0.2401688,
0.3676544, 0.2083172, 0.2242618, 0.189248, 0.394233, 0.1021266,
0.1932088, 0.2345265, 0.2592063, 0.1928255, 0.3043636, 0.1590346,
0.2699723, 0.3435174, 0.345459, 0.2461398, 0.417002, 0.3009782,
0.4046134, 0.3069134, 0.2635905, 0.236958, 0.2524431, 0.2953638,
0.2651712, 0.2060053, 0.2573521, 0.1668359, 0.1990478, 0.2384473,
0.1966746, 0.1753915, 0.2347065, 0.1921152, 0.1238467, 0.1489801,
0.2201506, 0.2194547, 0.2110882, 0.1880539, 0.2022153, 0.212393,
0.1929784, 0.1913548, 0.2233333, 0.1720814, 0.2011724, 0.1950076,
0.1997783, 0.2030418, 0.1765777, 0.1716257, 0.2371794, 0.1665875,
0.1690597, 0.1706616, 0.1656167, 0.1648148, 0.1589215, 0.1578078,
0.1710534, 0.1627099, 0.1427931, 0.1660376, 0.156027, 0.1709827,
0.1833439, 0.1650404, 0.1808825, 0.241505, 0.2132809, 0.1625625,
0.1750414), .Dim = c(120L, 2L))
sb.sequence <-
structure(function (T, b.av, length = T)
{
    index.sequence = c(1:T, 1:T)
    sequence = rep(0, length + T)
    current = 0
    while (current < length) {
        start = sample(1:T, 1)
        b = rgeom(1, 1/b.av) + 1
        sequence[(current + 1):(current + b)] = index.sequence[start:(start +
            b - 1)]
        current = current + b
    }
    sequence[1:length]
}, source = c("function(T, b.av, length = T){", "    index.sequence = c(1:T, 1:T)",
"    sequence = rep(0, length + T)", "    current = 0", "    while (current < length) {",
"      start = sample(1:T, 1)", "      b = rgeom(1, 1/b.av) + 1",
"      sequence[(current+1):(current+b)] = index.sequence[start:(start+b-1)]",
"      current = current + b", "    }", "    sequence[1:length]",
"  }"))
sharpe.ratio.diff <-
structure(function (ret)
{
    ret1 = ret[, 1]
    ret2 = ret[, 2]
    mu1.hat = mean(ret1)
    mu2.hat = mean(ret2)
    sig1.hat = var(ret1)^0.5
    sig2.hat = var(ret2)^0.5
    SR1.hat = mu1.hat/sig1.hat
    SR2.hat = mu2.hat/sig2.hat
    diff = SR1.hat - SR2.hat
    diff
}, source = c("function(ret){", "    ret1 = ret[,1]", "    ret2 = ret[,2]",
"    mu1.hat = mean(ret1)", "    mu2.hat = mean(ret2)", "    sig1.hat = var(ret1)^.5",
"    sig2.hat = var(ret2)^.5", "    SR1.hat = mu1.hat / sig1.hat",
"    SR2.hat = mu2.hat / sig2.hat", "    diff = SR1.hat - SR2.hat",
"    diff", "  }"))
