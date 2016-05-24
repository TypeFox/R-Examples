compute.sigma.star <-
function (no.bin, no.nor, prop.vec.bin = NULL, corr.vec = NULL, 
    corr.mat = NULL) 
{
    d = no.bin + no.nor
    validation.corr(no.bin, no.nor, prop.vec.bin, corr.vec = corr.vec, 
        corr.mat = corr.mat)
    if ((no.nor < 0) | (floor(no.nor) != no.nor)) {
        warning("Number of normal variables \nmust be a non-negative integer!\n")
    }
    if (is.null(corr.mat)) {
        corr.mat = lower.tri.to.corr.mat(corr.vec, d)
    }
    sigma = corr.mat
    p = prop.vec.bin
    q = 1 - p
    if (no.bin != 0) {
        sigmaBB = diag(no.bin)
        for (i in 1:no.bin) {
            for (j in 1:no.bin) {
                if (i != j) 
                  sigmaBB[i, j] = phi2poly(sigma[i, j], p[i], 
                    p[j])
            }
        }
    }
    if (no.bin > 0 & no.nor > 0) {
        sigmaBN = sigma
        for (i in (no.bin + 1):d) {
            for (j in 1:no.bin) {
                sigmaBN[i, j] = sigmaBN[i, j]/(dnorm(qnorm(p[j]))/sqrt(p[j] * 
                  q[j]))
            }
        }
        sigmaBN = sigmaBN[(no.bin + 1):d, 1:no.bin]
        sigma_star = sigma
        sigma_star[1:no.bin, 1:no.bin] = sigmaBB
        sigma_star[(no.bin + 1):d, 1:no.bin] = sigmaBN
        sigma_star[1:no.bin, (no.bin + 1):d] = t(sigmaBN)
    }
    if (no.bin > 0 & no.nor == 0) {
        sigma_star = sigmaBB
    }
    if (no.bin == 0 & no.nor > 0) {
        sigma_star = sigma
    }
    PD = TRUE
    temp = NULL
    eigenv = eigen(sigma_star)$value
    if (is.positive.definite(sigma_star) == FALSE) {
        temp = sigma_star
        cat("sigma_star before using nearPD is\n")
        print(temp)
        sigma_star = as.matrix(nearPD(sigma_star, corr = TRUE, 
            keepDiag = TRUE)$mat)
        sigma_star = (sigma_star + t(sigma_star))/2
    warning("sigma_star is not positive definite.\nAlgorithm will be using the nearest positive definite matrix which is\n", 
            immediate. = TRUE) 
        print(sigma_star)
        PD = FALSE
    }
    return(list(sigma_star = sigma_star, nonPD = temp, PD = PD, 
        eigenv = eigenv))
}

