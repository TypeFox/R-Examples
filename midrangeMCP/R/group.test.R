# Function for to show the results of the SNKM and TM tests
group.test  <- function(Teste){
  n      <- nrow(Teste)
  letras <- c(1:1000000)
  w      <- Teste[order(Teste[, 1], decreasing = TRUE), ]
  M      <- rep("", n)
  k      <- 1
  aux    <- cbind(w[, -1])
  x      <- apply(aux, 1, sum)

  if (dim(aux)[2] == 1 & sum(aux) == n) {
    aux3 <- aux
       M <- rep("g1", n)
  } else {
    if (any(x == 0)) {
      for (z in 1:n) {
        s <- rep(0, n)
        if (x[z] == 0) {
          s[z] <- 1
          aux  <- cbind(aux, s)
        }
      }
    }

    pos1    <- aux[1, ] == 1
    poscol1 <- which(pos1 == T)
    aux1    <- cbind(aux[ , -poscol1])
    numcol  <- ncol(aux1)
    ncont   <- -1

    for (j in 1:numcol) {
      k <- 0
      for (i in 1:n) {
        if (aux1[i, j] ==0 ) {
          k <- k + 1
        }
        if (aux1[i+1, j] == 1) break
      }
      ncont <- cbind(ncont,k)
    }

    aux2 <- aux1[, order(ncont[-1])]
    aux3 <- cbind(aux[, poscol1], aux2)
       a <- aux3[, 1]

    for (i in 1:n) {
      if (a[i] == 1) {
        M[i] <- paste("g", letras[1], sep="")
      } else {
        M[i] <- M[i]
      }
    }

    for (j in 2:(numcol+1)) {
      for (i in 1:n) {
        if (aux3[i, j] == 1) {
          M[i] <- paste(M[i], letras[j], sep = "g")
        } else {
          M[i] <- M[i]
        }
      }
    }
  }
  return(data.frame(Means = w[, 1], Groups = M))
}

# Function for to show the results of the SKM and SKR tests
group.test2 <- function(Teste) {
  n <- nrow(Teste)
  ordertest <- Teste[order(Teste[, 1], decreasing = TRUE), ]
  letras <- c(1:1000000)
  M <- rep("", n)

  if (all(ordertest[, 2] == 1)) {
    M <- rep("g1", n)
  } else {
    M[1] <- "g1"
    j <- 2

    for (i in 2:n) {
      if (ordertest[i, 2] == ordertest[i-1, 2]) {
        M[i] <- M[i-1]
      } else {
        M[i] <- paste(M[i], letras[j], sep = "g")
           j <- j + 1
      }
    }
  }
  return(data.frame(Means = ordertest[, 1], Groups = M))
}
