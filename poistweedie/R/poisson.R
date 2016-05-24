poisson <- 
function (x, n, p, lambda1, lower.tail = TRUE, log.p = FALSE, 
    fonction = "PROBABILITE") 
{
    lambda <- double(length(lambda1))
    mu <- double(length(lambda1))
    y <- double(length(x))
    if (lambda1 <= 0) {
        print(lambda1)
        stop("Erreur : la valeur de lambda1 doit etre strictement sup\'{e}rieure a 0")
    }
    if (p < exp(20)) {
        print(p)
        stop("Erreur : la valeur de p erreur doit sup\'{e}rieure a exp(10)")
    }
    lambda <- sqrt(lambda1)
    theta0 <- -(lambda)
    mu <- lambda
    if (fonction == "PROBABILITE") {
        y <- dpoistweedie(x, p, mu, lambda, theta0, log.p)
    }
    if (fonction == "FONCTION DE REPARTITION") {
        y <- ppoistweedie(x, p, mu, lambda, theta0, lower.tail, 
            log.p)
    }
    if (fonction == "FONCTION QUANTILE") {
        if (is.element(FALSE, 0 <= x && x <= 1) == TRUE) {
            stop("Erreur : x doit etre compris entre 0 et 1.\n")
        }
        y <- qpoistweedie(x, p, mu, lambda, theta0, lower.tail, 
            log.p)
    }
    if (fonction == "FONCTION RANDOM") {
        if (n < 0) {
            print(n)
            stop("Erreur : la valeur de n doit etre un entier naturel sup\'{e}rieure a 0 ")
        }
        y <- double(n)
        y <- rpoistweedie(n, p, mu, lambda, theta0)
    }
    if (fonction == "FONCTION VARIANCE") {
        theta0 = -(max(mu) + 1)
        y <- varpt(mu, p, theta0)
    }
    y
}