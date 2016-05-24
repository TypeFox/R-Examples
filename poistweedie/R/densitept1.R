densitept1 <-
function (p, n, mu, lambda, theta0) 
{
    denspt <- double((n + 1))
    if (mu < 0) {
        print(mu)
        stop("Erreur  la valeur de mu doit etre superieure a 0 ")
    }
    if (lambda <= 0) {
        print(lambda)
        stop("Erreur  la valeur de lambda doit etre strictement superieure a 0")
    }
    if (theta0 > 0) {
        print(theta0)
        stop("Erreur  la valeur de theta0 doit etre inferieure ou egale a 0")
    }
    if (p == 1) {
        denspt <- dpt1(p, n, mu, lambda, theta0)
    }
    else {
        if (p == 2) {
            denspt <- dpt2(p, n, mu, lambda, theta0)
        }
        else {
            if (1 < p) {
                denspt <- dptp(p, n, mu, lambda, theta0)
            }
            if (p < 1) {
				print(p)
				stop("Erreur  la valeur de p doit etre superieurea 1")
               
            }
        }
    }
    denspt
}

