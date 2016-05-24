varpt <-
function (mu, p, theta0) 
{
    perreur <- double(length(mu))
    omega <- double(length(mu))
    variance <- double(length(mu))
    if (p < 0) {
        print(p)
        stop("Erreur : la valeur de p doit etre sup\'{e}rieure a 1.\n")
    }
    if (theta0 > 0) {
        print(theta0)
        stop("Erreur : la valeur de theta0  doit etre inf\'{e}rieure ou \'{e}gale a 0.\n")
    }
    j <- 1
    for (i in 1:length(mu)) {
        if (mu[i] < 0) {
            perreur[j] <- i
            j <- j + 1
        }
    }
    if (j > 1) {
 
  print(perreur)
        stop("Erreur:les valeurs des entr\'{e}es ci-dessus ne sont pas sup\'{e}rieures a 0.\n ")
    }
    for (i in 1:length(mu)) {  
        omega[i] <- omega(p,mu[i],theta0)
        variance[i] <- mu[i] + exp(2 * omega[i] + p * (log(mu[i]) - 
                       omega[i]))
    }
    variance
}

