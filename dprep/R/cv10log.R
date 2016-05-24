cv10log <-
function (data, repet, maxwts = 2500) 
{
#    require(nnet)
    n <- dim(data)[1]
    p <- dim(data)[2]
    nombres <- colnames(data)
    f1 <- as.formula(paste(nombres[p], ".", sep = "~"))
    print(f1)
    ECV <- rep(0, repet)
    for (i in 1:repet) {
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
            tempo = nnet::multinom(f1, data = datat, MaxNWts = 2500)
            tempo1 = predict(tempo, datap)
            salida[j] <- sum(tempo1 != as.numeric(datap[, p]))
        }
        ECV[i] <- sum(salida)/n
    }
    cat("The error estimations in each repetition are:\n")
    print(ECV)
    ECV1 <- mean(ECV)
    cat("The mean error estimation by cross-validation using alll the repetititons is: \n")
    ECV1
}
