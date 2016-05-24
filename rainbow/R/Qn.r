Qn <- function (x)
{
    n <- length(x)
    diffs <- outer(x, x, "-")
    diffs <- diffs[!lower.tri(diffs, diag = TRUE)]
    qn <- 2.2219 * quantile(abs(diffs), 0.25)
    if (n == 2)
        dn <- 0.399
    else if (n == 3)
        dn <- 0.994
    else if (n == 4)
        dn <- 0.512
    else if (n == 5)
        dn <- 0.844
    else if (n == 6)
        dn <- 0.611
    else if (n == 7)
        dn <- 0.857
    else if (n == 8)
        dn <- 0.669
    else if (n == 9)
        dn <- 0.872
    else if (n %% 2 == 1)
        dn <- n / (n + 1.4)
    else dn <- n / (n + 3.8)
    return(dn * qn)
}
