cv10mlp <-
function (data, units, decay = 0, maxwts = 1000, maxit = 100, 
    repet) 
{
    #require(nnet)
    n <- dim(data)[1]
    p <- dim(data)[2]
    ecv <- rep(0, repet)
    for (kk in 1:repet) {
        azar <- data[rank(runif(n)), ]
        azar[, p] <- as.factor(azar[, p])
        parti <- floor(n/10)
        salida <- rep(0, 10)
        for (j in 1:10) {
            cc <- ((j - 1) * parti + 1):(j * parti)
            if (j == 10) {
                cc <- ((j - 1) * parti + 1):n
            }
            datap <- azar[cc, ]
            datat <- azar[-cc, ]
            clasest = nnet::class.ind(azar[-cc, p])
            tempo <- nnet::nnet(as.matrix(datat[, 1:(p - 1)]), clasest, 
                entropy = TRUE, size = units, decay = decay, 
                MaxNWts = maxwts, maxit = maxit)
            tempo1 = predict(tempo, datap)
            pd = max.col(tempo1)
            salida[j] <- sum(pd != datap[, p])
        }
        ecv[kk] <- sum(salida)/n
    }
    cat("The misclassification errors of each repetition are:", 
        "\n")
    print(ecv)
    cat("The mean misclassifcation error is", "\n")
    ECV1 = mean(ecv)
    ECV1
}
