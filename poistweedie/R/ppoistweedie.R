ppoistweedie <-
function(q, p, mu, lambda, theta0, lower.tail, log.p) 
{
    m <- numeric(1)
    n <- numeric(1)
    m <- length(q)
    n <- max(q)
    ppt <- double(m)
    dpt <- double((n + 1))
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
    for (i in 1:length(q)) {
        dpt <- dpoistweedie(0:q[i], p, mu, lambda, theta0, log = FALSE)
        if (log.p == TRUE || log.p == "T") {
            if (lower.tail == FALSE || lower.tail == "F") {
                ppt[i] <- log(1 - sum(dpt))
            }
            else {
                ppt[i] <- log(sum(dpt))
            }
        }
        else {
            if (lower.tail == FALSE || lower.tail == "F") {
                ppt[i] <- 1 - sum(dpt)
            }
            else {
                ppt[i] <- sum(dpt)
            }
        }
    }
    ppt
}

