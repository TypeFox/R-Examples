poistweedie <-
function (x, n, p, mu, lambda, theta0, lower.tail = TRUE, log.p = FALSE, 
    fonction = "PROBABILITE") 
{
    y <- double(length(x))
    if (p < 1) {
        print(p)
        stop("Erreur : la valeur de p erreur doit sup\'{e}rieure a 1")
    }
    if (mu <= 0) {
        print(mu)
        stop("Erreur : la valeur de mu doit etre strictement sup\'{e}rieure a 0")
    }
    if (lambda <= 0) {
        print(lambda)
        stop("Erreur : la valeur de lambda doit etre strictement sup\'{e}rieure a 0")
    }
    if (theta0 > 0) {
        print(theta0)
        stop("Erreur : la valeur de theta0 doit etre inf\'{e}rieure ou \'{e}gale a 0")
    }
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
