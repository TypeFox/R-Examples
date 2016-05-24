combn.test <-
function (stage1, stage2, weight = 0.5, method = "invnorm") 
{
    meths <- c("invnorm", "fisher")
    imeths <- as.integer(match(method, meths, -1))
    if (imeths < 1) {
        stop("Invalid combination test type")
    }
    p1 <- stage1$pvalues
    p2 <- stage2$pvalues
    hyp.comb <- stage1$hyp.comb
    treats <- length(p1)
    if (method == "invnorm") {
        w1 <- sqrt(weight)
        w2 <- sqrt(1 - w1 * w1)
    }
    else {
        w1 <- 1
        w2 <- 1
    }
    zscores <- p1
    for (i in 1:treats) {
        for (j in 1:length(p1[[i]])) {
            zscores[[i]][j] <- NA
            if (method == "invnorm") {
                zscores[[i]][j] <- w1 * qnorm(1 - p1[[i]][j]) + 
                  w2 * qnorm(1 - p2[[i]][j])
            }
            if (method == "fisher") {
                zscores[[i]][j] <- qnorm(pchisq(-2 * log(p1[[i]][j] * 
                  p2[[i]][j]), 4))
            }
        }
    }
    list(method = method, zscores = zscores, hyp.comb = hyp.comb, 
        weights = c(w1, w2))
}
