order.stat.SNK <-
function (treatment, means, minimum)
{
    # minimum tera que ser agora um vetor para ter n-1 DMS associadas para calcular o SNK
    n <- length(means)
    z <- data.frame(treatment, means)
    w <- z[order(z[, 2], decreasing = TRUE), ]
    M <- rep("", n)
    k <- 1
    k1 <- 0
    j <- 1
    i <- 1
    r <- 1  # iniciando o meu indice da DMS
    cambio <- n
    cambio1 <- 0
    chequeo = 0
    M[1] <- letters[k]
    while (j < n) {
        chequeo <- chequeo + 1
        if (chequeo > n)
            break
        for (i in j:n) {
            # r para a distancia entre os tratamentos sendo minimum um vetor de n-1 DMS
            if (abs(j-i) == 0) {r<-1} else {r<-abs(j-i)}
            s <- abs(w[i, 2] - w[j, 2]) <= minimum[r]
            if (s) {
                if (lastC(M[i]) != letters[k])
                  M[i] <- paste(M[i], letters[k], sep = "")
            }
            else {
                k <- k + 1
                cambio <- i
                cambio1 <- 0
                ja <- j
                for (jj in cambio:n) M[jj] <- paste(M[jj], " ",
                  sep = "")
                M[cambio] <- paste(M[cambio], letters[k], sep = "")
                for (v in ja:cambio) {
                  if (abs(v-cambio) == 0) {r<-1} else {r<-abs(v-cambio)}
                  if (abs(w[v, 2] - w[cambio, 2]) > minimum[r]) {
                    j <- j + 1
                    cambio1 <- 1
                  }
                  else break
                }
                break
            }
        }
        if (cambio1 == 0)
            j <- j + 1
    }
    w <- data.frame(w, stat = M)
    trt <- as.character(w$treatment)
    means <- as.numeric(w$means)
    for (i in 1:n) {
        cat(M[i], "\t", trt[i], "\t    ", means[i], "\n")
    }
    output <- data.frame(trt, means, M)
    return(output)
}
