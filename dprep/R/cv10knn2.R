cv10knn2 <-
function (data, kvec) 
{
    n <- dim(data)[1]
    p <- dim(data)[2]
    salida <- matrix(0, 1, 10)
    azar <- data[rank(runif(n)), ]
    azar[, p] <- as.factor(azar[, p])
    parti <- floor(n/10)
    for (j in 1:10) {
        cc <- ((j - 1) * parti + 1):(j * parti)
        if (j == 10) {
            cc <- ((j - 1) * parti + 1):n
        }
        datap <- azar[cc, ]
        datat <- azar[-cc, ]
        tempo <- class::knn(as.matrix(datat[, 1:p - 1]), as.matrix(datap[, 
            1:p - 1]), datat[, p], kvec)
        salida[j] <- sum(tempo != datap[, p])
    }
    ECV1 <- n - sum(salida)
    return(ECV1)
}
