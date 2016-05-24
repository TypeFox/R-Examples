dpoistweedie <-
function (y, p, mu, lambda, theta0, log) 
{
    m <- numeric(1)
    n <- numeric(1)
    m <- length(y)
    n <- max(y)
    denspt <- double((n + 1))
    dpt <- double(m)
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
    if (log == TRUE || log == "T") {
        denspt <- densitept2(p, n, mu, lambda, theta0)
    }
    else {
        denspt <- densitept1(p, n, mu, lambda, theta0)
    }
    for (i in 1:length(y)) {
        dpt[i] <- denspt[(y[i] + 1)]
    }
    dpt
}

