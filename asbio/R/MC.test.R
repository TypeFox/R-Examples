MC.test <- function (Y, X, perm = 1000, alternative = "not.equal", var.equal = TRUE, paired = FALSE) 
{
    initial <- t.test(Y ~ X, var.equal = TRUE)$statistic
    perm.results <- matrix(ncol = 1, nrow = perm)
    for (i in 1:perm) {
        x1 <- sample(Y, replace = TRUE)
        perm.results[i] <- t.test(x1 ~ X, var.equal = var.equal, paired = paired)$statistic
    }
    if (alternative == "less") {
        num <- length(perm.results[perm.results <= initial]) + 
            1
    }
    if (alternative == "greater") {
        num <- length(perm.results[perm.results >= initial]) + 
            1
    }
    if (alternative == "not.equal") {
        pos.init <- abs(initial)
        num <- length(perm.results[perm.results >= pos.init] + 
            1)
    }
    p.value <- num/perm
    result <- list(observed.test.statistic = initial, no_of_permutations_exceeding_observed_value = num - 
        1, p.value = p.value, alternative.hypothesis = alternative)
    result
}
