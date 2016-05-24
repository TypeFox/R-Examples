rpoistweedie <-
function (n, p, mu, lambda, theta0)
{
    rpt <- double(n)
    ppt <- double(n)
    if (n < 0) {
        print(n)
        stop("Erreur:n n'est pas un entier naturel")
    }
    if (mu < 0) {
        print(mu)
        stop("Erreur : la valeur de mu doit etre sup\'{e}rieure a 0 ")
    }
    if (lambda <= 0) {
        print(lambda)
        stop("Erreur : la valeur de lambda doit etre strictement sup\'{e}rieure a 0")
    }
    if (theta0 > 0) {
        print(theta0)
        stop("Erreur : la valeur de theta0 doit etre inf\'{e}rieure ou \'{e}gale a 0")
    }
    ppt <- runif(n)
    rpt <- qpoistweedie(ppt, p, mu, lambda, theta0, lower.tail = TRUE,
        log.p = FALSE)
    rpt
}

