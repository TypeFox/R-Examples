varip2 <-
function (i, p, ll, S, P) 
{
    storewrap <- NULL
    covAA <- covAB <- covBB <- 0
    for (r in 0:(2^(i - 1) - 1)) {
        covIansA <- covIwrap(S = S, m = 1, n = 1 + r, ll = ll, 
            storewrap = storewrap, P = P)
        storewrap <- covIansA$storewrap
        covAA <- covAA + (2^(i - 1) - r) * covIansA$ans * (2 - 
            (r == 0))
    }
    for (r in (1 - 2^(i - 1)):(2^(i - 1) - 1)) {
        covIansB <- covIwrap(S = S, m = 1, n = 1 + r + 2^(i - 
            1), ll = ll,  storewrap = storewrap, P = P)
        storewrap <- covIansB$storewrap
        covAB <- covAB + (2^(i - 1) - abs(r)) * covIansB$ans
    }
    covBB <- covAA
    all <- covAA - 2 * covAB + covBB
    all <- 2^(-i) * all
    list(covAA = covAA, covAB = covAB, covBB = covBB, ans = all)
}
