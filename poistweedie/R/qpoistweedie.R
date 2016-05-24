qpoistweedie <-
function (p1, p, mu, lambda, theta0, lower.tail, log.p) 
{
    m <- numeric(1)
    n <- numeric(1)
    k <- numeric(1)
    x <- numeric(1)
    trouve <- numeric(1)
    max <- numeric(1)
    m <- length(p1)
    n <- max(p1)
    qpt <- double(m)
    perreur <- double(m)
    ppt <- double(1000)
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
    j <- 1
    for (i in 1:length(p1)) {
        if (p1[i] < 0 || p1[i] > 1) {
            perreur[j] <- i
            j <- j + 1
        }
    }
    if (j > 1) {
        print(perreur)
        stop("Erreur:certainnes valeurs  ne sont pas comprises entre 0 et 1 ")
    }
    max <- 0
    for (i in 1:length(p1)) {
        if (0 < p1[i] && p1[i] <= 1 && p1[i] > max) {
            max <- p1[i]
        }
    }
    k <- 0
    i <- 1
    trouve = 0
    while (trouve == 0) {
        ppt[i] <- ppoistweedie(k, p, mu, lambda, theta0, lower.tail, 
            log.p)
        if ((ppt[i]) >= max) {
            x <- k
            trouve <- 1
        }
        else {
            k <- k + 1
            i <- i + 1
        }
    }
    for (i in 1:length(p1)) {
        if ((p1[i]) == 0) {
            qpt[i] <- 0
        }
        else {
            if (0 < p1[i] && p1[i] <= 1) {
                k <- 0
                j <- 1
                trouve = 0
                while (trouve == 0) {
                  if ((ppt[j]) > p1[i]) {
                    qpt[i] <- k
                    trouve <- 1
                  }
                  else {
                    k <- k + 1
                    j <- j + 1
                  }
                }
            }
        }
    }
    qpt
}

