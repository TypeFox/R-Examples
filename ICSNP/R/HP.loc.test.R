HP.loc.test <- function(X, mu = NULL, score = "rank", angles="tyler", method = "approximation", n.perm = 1000, na.action = na.fail)
    {
    DNAME <- deparse(substitute(X))
    k <- dim(X)[2]
    if (is.null(mu)) mu <- rep(0, k) else 
    if (length(mu) != k) 
    stop("length of 'mu' must equal the number of columns of 'X'")
    X <- na.action(X)
    if (!all(sapply(X, is.numeric))) 
    stop("'X' must be numeric")
    X <- as.matrix(X)
    method <- match.arg(method,c("approximation", "permutation"))
    score <- match.arg(score, c("sign", "rank", "normal"))
    angles <- match.arg(angles, c("tyler", "interdirection"))
    ALTERNATIVE <- "two.sided"
    NVAL <- paste("c(", paste(mu, collapse = ","), ")", sep = "")
    names(NVAL) <- "location"
    
    res1 <- switch(angles,{
                tyler = HP.loc.tyler.test(X=X, mu=mu, score=score, method=method, n.perm=n.perm)}
                ,
                interdirection ={ stop("The test based on interdirections is not implemted yet")}
                )
    
    res <- c(res1, list(data.name = DNAME, alternative = ALTERNATIVE, null.value = NVAL))
    class(res) <- "htest"
    return(res)
    
    }

HP.loc.tyler.test <-
function (X, mu = NULL, score = "rank", method = "approximation", 
    n.perm = 1000) 
{
    k <- dim(X)[2]
    TylerH0 <- tyler.shape(X, location = mu)
    mahaH0 <- mahalanobis(X, center = mu, cov = TylerH0)
    ranksH0 <- rank(mahaH0)
    signsH0 <- spatial.sign(X, center = mu, shape = TylerH0)
    n <- length(ranksH0)
    scoresH0 <- switch(score, sign = 1, rank = ranksH0, normal = sqrt(qchisq(ranksH0/(n + 
        1), k)))
    signsH0w <- scoresH0 * signsH0
    switch(score, sign = {
        DIV <- k/n
        METHOD <- "TYLER ANGLES SIGN TEST"
        S.name <- "Q.S"
    }, rank = {
        DIV <- 3 * k/(n * (n + 1)^2)
        METHOD <- "TYLER ANGLES RANK TEST"
        S.name <- "Q.W"
    }, normal = {
        DIV <- 1/n
        METHOD <- "TYLER ANGLES VAN DER WAERDEN TEST"
        S.name <- "Q.N"
    })
    SUMS <- tcrossprod(signsH0w)
    Q1 <- sum(SUMS)
    STATISTIC <- DIV * Q1
    switch(method, approximation = {
        PVAL <- 1 - pchisq(STATISTIC, k)
        PARAMETER <- k
        names(PARAMETER) <- "df"
    }, permutation = {
        rep.func <- function(index, signs, n) {
            signs2 <- index*signs
            SUMS.2 <- tcrossprod(signs2)
            sum(SUMS.2)
        }
        Q.simu <- replicate(n.perm, rep.func(sample(c(1, -1), 
            n, replace = TRUE), signsH0w, n))
        PVAL <- mean(Q.simu > Q1)
        PARAMETER <- n.perm
        names(PARAMETER) <- "permutations"
    })
    names(STATISTIC) <- S.name
    res <- list(statistic = STATISTIC, p.value = PVAL, method = METHOD, 
        parameter = PARAMETER)
    return(res)
}
