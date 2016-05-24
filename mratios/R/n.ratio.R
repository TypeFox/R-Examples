"n.ratio" <-
function (m, rho, Power, CV0, rho.star, alpha, Min.power = TRUE) 
{

    Z1.beta <- qnorm(Power)
    if(m > 1){
      correla <- (rho^2)/sqrt((rho^2 + 1) * (rho^2 + 1))
      R.matrix <- matrix(rep(correla, m * m), nrow = m)
      diag(R.matrix) <- rep(1, m)
        Equi.coord1 <- qmvnorm(1 - alpha, c(0, 10), corr = R.matrix, 
        tail = "lower.tail", abseps=1e-05)$quantile
      if (Min.power) {
          n.balance <- ((Equi.coord1 + Z1.beta)^2) * ((1 + rho^2)/((rho.star - 
             rho)^2)) * (CV0^2)
          n.balance <- ceiling(n.balance)
          cat(" ", "\n")
          cat(" Number of observations per treatment = ", n.balance, 
              "\n")
          cat("         Total number of observations = ", (m + 
              1) * n.balance, "\n")
          cat("                                      ", "\n")
      }
      else {
          Equi.coord2 <- qmvnorm(Power, c(0, 10), corr = R.matrix, 
              tail = "lower.tail", abseps=1e-05)$quantile
          n.balance <- ((Equi.coord1 + Equi.coord2)^2) * ((1 + 
              rho^2)/((rho.star - rho)^2)) * (CV0^2)
          n.balance <- ceiling(n.balance)
          cat(" ", "\n")
          cat(" Number of observations per treatment = ", n.balance, 
              "\n")
          cat("         Total number of observations = ", (m + 
              1) * n.balance, "\n")
          cat("                                      ", "\n")
      }

      } 

    else{   #  m = 1
      Z1.alpha <- qnorm(1-alpha)
      n.balance <- ((Z1.alpha + Z1.beta)^2) * ((1 + rho^2)/((rho.star - 
             rho)^2)) * (CV0^2)
      n.balance <- ceiling(n.balance)
      cat(" ", "\n")
      cat(" Number of observations per treatment = ", n.balance, 
              "\n")
      cat("         Total number of observations = ", (m + 
              1) * n.balance, "\n")
      cat("                                      ", "\n")
      }
}

