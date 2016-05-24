dunnett.test <-
function (Z = Z, select = rep(1, length(Z))) 
{
    if(sum(is.na(Z))>0){stop("No missing values allowed in Z")}
    treats <- length(Z)
    hyp.comb <- list(NULL)
    nhyp.comb <- vector(length = treats)
    for (i in 1:treats) {
        comb.dat <- combn(1:treats, i)
        hyp.comb[[i]] <- comb.dat
        nhyp.comb[i] <- ncol(hyp.comb[[i]])
        rownames(hyp.comb[[i]]) <- 1:i
        hypcol <- NULL
        for (j in 1:nhyp.comb[i]) {
            ihypcol <- paste("H", paste(hyp.comb[[i]][, j], collapse = ""), 
                sep = "")
            hypcol <- append(hypcol, ihypcol)
        }
        colnames(hyp.comb[[i]]) <- hypcol
    }
    pdunnett.test <- list(NULL)
    zscores <- list(NULL)
    int_dunnett <- function(x, z, k) {
        ((pnorm(sqrt(2) * z + x))^k) * dnorm(x)
    }
    for (i in 1:treats) {
        ptest <- NULL
        if (select[i] == 0) {
            Z[i] <- -Inf
        }
        for (j in 1:nhyp.comb[i]) {
            kselect <- sum(select[c(hyp.comb[[i]][, j])])
            Zmax <- max(Z[c(hyp.comb[[i]][, j])])
            if (kselect == 0) {
                dunnet_integral <- 0
                F_Zmax <- 1
            }
            else {
                dunnett_integral <- integrate(int_dunnett, -Inf, 
                  Inf, z = Zmax, k = kselect)
                F_Zmax <- 1 - dunnett_integral$value
            }
            ptest <- append(ptest, F_Zmax)
        }
        pdunnett.test[[i]] <- matrix(ptest, nrow = 1, ncol = length(ptest))
        colnames(pdunnett.test[[i]]) <- colnames(hyp.comb[[i]])
        rownames(pdunnett.test[[i]]) <- 1
        zscores[[i]] <- qnorm(1 - pdunnett.test[[i]])
    }
    list(pvalues = pdunnett.test, zscores = zscores, hyp.comb = hyp.comb)
}
