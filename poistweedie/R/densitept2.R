densitept2 <-
function (p, n, mu, lambda, theta0) 
{
    denspt <- double((n + 1))
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
    if (p == 1) {
        denspt <- dpt1Log(p, n, mu, lambda, theta0)
    }
    else {
        if (p == 2) {
            denspt <- dpt2Log(p, n, mu, lambda, theta0)
        }
        else {
            if (1 < p) {
                denspt <- dptpLog(p, n, mu, lambda, theta0)
            }
            if (p < 1) {
                print(p)
                stop("Erreur : la valeur de p doit etre sup\'{e}rieure a 1")
            }
        }                   
    }
    denspt
}

