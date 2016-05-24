## Density.

dmovMF <-
function(x, theta, alpha = 1, log = FALSE)
{
    ## Be nice to users.
    theta <- rbind(theta)
    k <- max(nrow(theta), length(alpha))
    theta <- theta[rep_len(seq_len(nrow(theta)), k), , drop = FALSE]
    alpha <- rep_len(alpha, k)
    alpha <- alpha / sum(alpha)

    y <- ldmovMF(x, theta, alpha)
    if(!log)
        y <- exp(y)
    y
}

## Log-density

ldmovMF <-
function(x, theta, alpha)
{
    ## Internal use: no recycling.
    p <- ncol(x) / 2
    ## So s = d / 2 - 1 = p - 1.
    kappa <- row_norms(theta)
    G <- cadd(skmeans:::g_tcrossprod(x, theta),
              log(alpha) - lH(kappa, p - 1))
    log_row_sums(G)
}
