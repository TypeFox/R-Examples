crossval <-
function (data, nparts = 10, method = c("lda", "knn", "rpart","logistic","naiveBayes"), 
    kvec = 5, maxwts=2500, repet) 
{
    if (method == "lda") {
        #require(MASS)
        n <- dim(data)[1]
        p <- dim(data)[2]
        errorcv <- rep(0, repet)
        for (i in 1:repet) {
            salida <- matrix(0, 1, nparts)
            azar <- data[rank(runif(n)), ]
            parti <- floor(n/nparts)
            for (j in 1:nparts) {
                cc <- ((j - 1) * parti + 1):(j * parti)
                if (j == nparts) {
                  cc <- ((j - 1) * parti + 1):n
                }
                datap <- azar[cc, ]
                datat <- azar[-cc, ]
                tempo <- MASS::lda(as.matrix(datat[, 1:p - 1]), datat[, 
                  p])
                tempo1 <- predict(tempo, as.matrix(datap[, 1:p - 
                  1]))$class
                salida[j] <- sum(as.numeric(tempo1) != as.numeric(datap[, 
                  p]))
            }
            errorcv[i] <- sum(salida)/n
        }
        errorp = mean(errorcv)
    }
    if (method == "rpart") {
        #require(rpart)
        datos = as.data.frame(data)
        n <- dim(datos)[1]
        p <- dim(datos)[2]
        errorcv <- rep(0, repet)
        for (kk in 1:repet) {
            nombres <- colnames(datos)
            f1 <- as.formula(paste(nombres[p], ".", sep = "~"))
            salida <- matrix(0, 1, nparts)
            azar <- datos[rank(runif(n)), ]
            azar[, p] <- as.factor(azar[, p])
            parti <- floor(n/nparts)
            for (j in 1:nparts) {
                cc <- ((j - 1) * parti + 1):(j * parti)
                if (j == nparts) {
                  cc <- ((j - 1) * parti + 1):n
                }
                datap <- azar[cc, ]
                datat <- azar[-cc, ]
                arbol <- rpart::rpart(f1, data = datat, method = "class")
                pd1 <- predict(arbol, datap)
                pd2 = max.col(pd1)
                salida[j] <- sum(pd2 != as.numeric(datap[, p]))
            }
            errorcv[kk] <- sum(salida)/n
        }
        errorp = mean(errorcv)
    }
    if (method == "knn") {
        #require(class)
        n <- dim(data)[1]
        p <- dim(data)[2]
        errorcv <- rep(0, repet)
        for (kk in 1:repet) {
            azar <- data[rank(runif(n)), ]
            azar[, p] <- as.factor(azar[, p])
            parti <- floor(n/nparts)
            salida <- rep(0, nparts)
            for (j in 1:nparts) {
                cc <- ((j - 1) * parti + 1):(j * parti)
                if (j == nparts) {
                  cc <- ((j - 1) * parti + 1):n
                }
                datap <- azar[cc, ]
                datat <- azar[-cc, ]
                tempo <- class::knn(as.matrix(datat[, 1:p - 1]), as.matrix(datap[, 
                  1:p - 1]), datat[, p], kvec)
                salida[j] <- sum(as.numeric(tempo) != as.numeric(datap[, 
                  p]))
            }
            errorcv[kk] <- sum(salida)/n
        }
        errorp = mean(errorcv)
    }
    if(method=="logistic")
    {#require(nnet)
    n <- dim(data)[1]
    p <- dim(data)[2]
    nombres <- colnames(data)
    f1 <- as.formula(paste(nombres[p], ".", sep = "~"))
    print(f1)
    errorcv <- rep(0, repet)
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
        errorcv[i] <- sum(salida)/n
    }
    errorp = mean(errorcv)
    }
    if(method=="naiveBayes")
    {#require(e1071)
    n <- dim(data)[1]
    p <- dim(data)[2]
    nombres<-colnames(data)
    f1<-as.formula(paste(nombres[p],".",sep="~"))
    #print(f1)
    errorcv <- rep(0, repet)
    for(i in 1:repet) {
        salida <- matrix(0, 1, 10)
        azar <- data[rank(runif(n)),  ]
        parti <- floor(n/10)
        for(j in 1:10) 
        {
            cc <- ((j - 1) * parti + 1):(j * parti
            )
            if(j == 10) {
                cc <- ((j - 1) * parti + 1):
                    n
            }
            datap <- azar[cc,  ]
            datat <- azar[ - cc,  ]
            tempo =e1071::naiveBayes(f1,data=datat)
            tempo1=predict(tempo,datap[,-p],type="raw")
            tempo1=max.col(tempo1)
            salida[j] <- sum(tempo1 != as.numeric(datap[, p]))
        }
        errorcv[i] <- sum(salida)/n
    }
    errorp <- mean(errorcv)
    }
    cat("The error estimation in each repetition are:\n")
    print(errorcv)
    cat("The mean error estimation by cross-validation using all the repetititons is: \n")
    errorp
    errorp
}
