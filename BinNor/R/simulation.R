simulation <-
function (seed = NULL, nsim, no.rows, no.bin, no.nor, mean.vec.nor = NULL, 
    var.nor = NULL, prop.vec.bin = NULL, corr.vec = NULL, corr.mat = NULL, 
    continue.with.warning = TRUE) 
{
    d = no.bin + no.nor
    if ((nsim < 1) | (floor(nsim) != nsim)) {
        stop("Number of simulations must be an integer whose value is at least 1!\n")
    }
    if (is.null(seed)) {
        seed = runif(1, 1e+06, 9e+06)
    }
    set.seed(seed)
    d = no.bin + no.nor
    if (is.null(corr.mat)) {
        corr.mat = lower.tri.to.corr.mat(corr.vec, d)
    }
    sigma.star = compute.sigma.star(no.bin, no.nor, prop.vec.bin, 
        corr.vec, corr.mat)
    if (sigma.star$PD == FALSE) {
        if (continue.with.warning == TRUE) {
            warning("sigma_star is not positive definite.\nAlgorithm used the nearest positive definite matrix!!!!", 
                immediate. = TRUE)
        }
        else {
            stop("User has chosen to stop as the final correlation matrix is not positive definite")
        }
    }
    emp.mean = matrix(0, nsim, d)
    emp.corr = matrix(0, nsim, d^2)
    emp.var = matrix(0, nsim, no.nor)
    for (i in 1:nsim) {
        print(c(i, date()))
        mydata = jointly.generate.binary.normal(no.rows, no.bin, 
            no.nor, prop.vec.bin, mean.vec.nor, var.nor, sigma_star = sigma.star$sigma_star, 
            continue.with.warning = TRUE)
        emp.mean[i, ] = apply(mydata, 2, mean)
        if (no.nor > 0) {
            emp.var[i, ] = apply(mydata[, (no.bin + 1):d], 2, 
                var)
        }
        emp.corr[i, ] = as.vector(cor(mydata))
    }
    emp.cormat.mean = matrix(apply(emp.corr, 2, mean), d, d)
    if (sigma.star$PD == FALSE) {
        cat("==============================:\n")
        cat("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!:\n")
        warning("sigma_star is not positive definite.\nAlgorithm used the nearest positive definite matrix!!!!", 
            immediate. = TRUE)
        print(sigma.star$nonPD)
    }
    cat("==============================:\n")
    cat("Desired correlation matrix:\n")
    print(corr.mat)
    cat("Averaged correlation matrix:\n")
    print(round(emp.cormat.mean, 7))
    cat("==============================:\n")
    cat("Desired proportion and mean parameters: \n")
    print(c(prop.vec.bin, mean.vec.nor))
    cat("Averaged proportion and mean parameters: \n")
    print(apply(emp.mean, 2, mean))
    cat("==============================:\n")
    cat("Target variance of normal variate(s): \n")
    print(var.nor)
    cat("Actual variance of normal variates(s): \n")
    print(apply(emp.var, 2, mean))
    cat("==============================:\n")
    cat("Seed used:", seed, "\n")
}

