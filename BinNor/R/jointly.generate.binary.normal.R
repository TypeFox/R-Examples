jointly.generate.binary.normal <-
function (no.rows, no.bin, no.nor, prop.vec.bin = NULL, mean.vec.nor = NULL, 
    var.nor = NULL, sigma_star = NULL, corr.vec = NULL, corr.mat = NULL, 
    continue.with.warning = TRUE) 
{

if ((no.rows<1)|(floor(no.rows)!=no.rows)){stop("Number of rows must be an integer whose value is at least 1!\n")}
    d = no.bin + no.nor
    validation.bin(no.bin, prop.vec.bin)
    validation.nor(no.nor, mean.vec.nor, var.nor)
    if (is.null(sigma_star)) {
        validation.corr(no.bin, no.nor, prop.vec.bin, corr.vec = corr.vec, 
            corr.mat = corr.mat)
        sig_star = compute.sigma.star(no.bin, no.nor, prop.vec.bin, 
            corr.vec, corr.mat)
        sigma_star = sig_star$sigma_star
        if (sig_star$PD == FALSE & continue.with.warning == FALSE) {
            stop("User has chosen to stop as the final correlation matrix is not positive definite")
        }
    }
    else {
        if (is.positive.definite(sigma_star) == FALSE) {
            print(sigma_star)
            if (continue.with.warning == TRUE) {
                sigma_star = as.matrix(nearPD(sigma_star, corr = TRUE, 
                  keepDiag = TRUE)$mat)
                sigma_star = (sigma_star + t(sigma_star))/2
                warning("sigma_star is not positive definite.\nAlgorithm will be using the nearest positive definite matrix!\nThe nearest positive definite matrix computed is:", 
                  immediate. = TRUE)
            }
            else {
                stop("The final correlation matrix is not positive definite")
            }
        }
    }
    data = rmvnorm(no.rows, mean = rep(0, d), sigma = sigma_star)
    p = prop.vec.bin
    q = 1 - p
    if (no.bin > 0) {
        for (i in 1:no.rows) {
            for (j in 1:no.bin) {
                if (data[i, j] <= qnorm(1 - p[j])) 
                  data[i, j] = 0
                else data[i, j] = 1
            }
        }
    }
    if (no.nor > 0) {
        temp = 1
        for (j in (no.bin + 1):d) {
            data[, j] = mean.vec.nor[temp] + (data[, j] * sqrt(var.nor[temp]))
            temp = temp + 1
        }
    }
    return(data)
}

