svpls <- function (k1, k2, Y, pmax = 3, fdr = 0.05){
    Y <- as.matrix(Y)
    G <- nrow(Y)
    n <- k1 + k2
    val <- function(u) return(fitModel(k1, k2, Y, u))
    model_res <- vector("list", length = pmax + 1)
    model_res <- apply(as.matrix(c(0:pmax)), 1, val)
    model_aic <- numeric()
    model_aic[1] <- model_res[[1]][[6]]
    for (i in 2:(pmax + 1)) model_aic[i] <- model_res[[i]][[11]]
    opt_model <- which(model_aic == min(model_aic))
    if (opt_model == 1) {
        beta_hat <- NULL
        PLS.imp <- matrix(0, G, n)
        gv_diff <- (1 + k1/k2) * model_res[[1]][[4]][, 1]
        mse <- model_res[[1]][[5]]
        aic.opt <- model_res[[1]][[6]]
        t <- pval.opt <- rep(0, G)
        for (i in 1:G) {
            t[i] <- gv_diff[i]/sqrt((1 - 1/G) * (1/k1 + 1/k2) * 
                mse)
            pval.opt[i] <- 2 * (1 - pt(abs(t[i]), (n-2)*G))
        }
        p.opt <- sort(pval.opt)
        sig <- 0
        gene <- index <- rep(0, G)
        for (i in 1:G) {
            if (p.opt[i] <= i * fdr/G) {
                sig <- sig + 1
                gene[i] <- order(pval.opt)[i]
                index[i] <- i
            }
        }
        index <- index[index != 0]
        sig_genes.opt <- row.names(Y)[gene[index]]
    }
    if (opt_model != 1) {
        gv_diff <- (1 + k1/k2) * model_res[[opt_model]][[4]][, 
            1]
        sc <- model_res[[opt_model]][[5]]
        beta_hat <- model_res[[opt_model]][[6]]
        est_gz1 <- model_res[[opt_model]][[7]]
        est_vz1.1 <- model_res[[opt_model]][[8]][1]
        est_vz1.2 <- model_res[[opt_model]][[8]][2]
        gv.diff_var <- ((1 + k1/k2)^2) * model_res[[opt_model]][[9]]
        mse.opt <- model_res[[opt_model]][[10]]
        aic.opt <- model_res[[opt_model]][[11]]
        R <- VZ.Z <- matrix(0, G, n)
        surr <- rep(0, n)
        for (l in 1:(opt_model - 1)) surr <- surr + beta_hat[l] * 
            sc[, l]
        for (i in 1:G) R[i, ] <- surr
        GZ.Z <- outer(est_gz1, sc[, 1], function(u, v) u * v)
        for (i in 1:G) {
            VZ.Z[i, 1:k1] <- (est_vz1.1 * sc[1:k1, 1])
            VZ.Z[i, (k1 + 1):n] <- (est_vz1.2 * sc[(k1 + 1):n, 
                1])
        }
        PLS.imp <- R + GZ.Z + VZ.Z
        t <- pval.opt <- rep(0, G)
        for (i in 1:G) {
            t[i] <- gv_diff[i]/sqrt(gv.diff_var[i] * mse.opt)
            pval.opt[i] <- 2 * (1 - pt(abs(t[i]), n*G - 3 * G - 
                pmax))
        }
        p.opt <- sort(pval.opt)
        sig <- 0
        gene <- index <- rep(0, G)
        for (i in 1:G) {
            if (p.opt[i] <= i * fdr/G) {
                sig <- sig + 1
                gene[i] <- order(pval.opt)[i]
                index[i] <- i
            }
        }
        index <- index[index != 0]
        sig_genes.opt <- row.names(Y)[gene[index]]
    }
    Y.corr <- Y - PLS.imp
    res <- list(opt.model = opt_model, PLS.imp = PLS.imp, Y.corr = Y.corr, 
        pvalues.adj = pval.opt, genes = sig_genes.opt, AIC.opt = aic.opt)
    class(res) <- c("svpls", "list", "vector")
    return(res)
}

## new summary function S3
summary.svpls <- function(object){
cat("Optimum model is \n")
print(object$opt_model)
cat("AIC value for the optimal model is \n")
print(object$AIC.opt)
cat("The PLS imputed estimate of the hidden expression heterogeneity is \n")
print(object$PLS.imp)
cat("The corrected gene expression matrix is \n")
print(object$Y.corr)
cat("The adjusted pvalues are \n")
print(object$pvalues.adj)
cat("The significant genes are \n")
print(object$genes)
}

## new print function S3
print.svpls <- function(x){
cat("The Optimal Model is: \n")
print(x$opt.model)
cat("\n")
cat("The differentially expressed genes are: \n")
print(x$genes)
cat("\n")
}