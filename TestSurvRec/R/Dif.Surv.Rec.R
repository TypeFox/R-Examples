Dif.Surv.Rec <-
function (XX, type = "LRrec", alfa = 0, beta = 0, gamma = 0, 
    eta = 0) 
{
    XL <- XX
    x <- factor(XL$group)
    Factores <- c(levels(x))
    Nivelesdefactores <- matrix(x)
    fit1 <- survfitr(Survr(id, time, event) ~ as.factor(group), 
        data = XL, type = "pe")
    fit2 <- survfitr(Survr(id, time, event) ~ 1, data = XL, type = "pe")
    failed <- matrix(fit2$failed)
    censored <- matrix(fit2$censored)
    time <- matrix(fit2$time)
    n.event <- matrix(fit2$n.event)
    AtRisk <- matrix(fit2$AtRisk)
    m <- fit2$m
    k <- matrix(m)
    n <- ncol(matrix(fit2$time))
    m <- nrow(matrix(fit2$time))
    n1 <- ncol(matrix(fit1[[Nivelesdefactores[1, 1]]]$time))
    m1 <- nrow(matrix(fit1[[Nivelesdefactores[1, 1]]]$time))
    I <- matrix(1, m, 1)
    X <- matrix(fit2$time)
    Y <- matrix(fit1[[Nivelesdefactores[1, 1]]]$time)
    Z <- matrix(fit1[[Nivelesdefactores[1, 1]]]$censored)
    Tiempo1 <- matrix(0, m, 1)
    n.eventG1 <- matrix(0, m, 1)
    NEG1 <- matrix(fit1[[Nivelesdefactores[1, 1]]]$n.event)
    Ariesgo <- matrix(0, m, 1)
    cuenta <- matrix(0, m, 1)
    Ncensura <- matrix(0, m, 1)
    AriesgoG1 <- matrix(0, m, 1)
    mm <- nrow(X)
    tiempocensura <- matrix(fit1[[Nivelesdefactores[1, 1]]]$censored)
    for (z in 1:m) {
        for (zz in 1:m1) {
            if (X[z, 1] == Y[zz, 1]) {
                Tiempo1[z, 1] <- Y[zz, 1]
                n.eventG1[z, 1] <- NEG1[zz, 1]
                zz <- m1
            }
            else {
                Tiempo1[z, 1] <- X[z, 1]
            }
        }
    }
    AriesgoG1tiempo1 <- t(fit1[[Nivelesdefactores[1, 1]]]$AtRisk) %*% 
        matrix(1, nrow(fit1[[Nivelesdefactores[1, 1]]]$AtRisk), 
            1)
    N <- fit1[[Nivelesdefactores[1, 1]]]$n + t(fit1[[Nivelesdefactores[1, 
        1]]]$m) %*% matrix(1, nrow(fit1[[Nivelesdefactores[1, 
        1]]]$m), 1)
    for (z in 1:m) {
        cuenta[z, 1] <- 0
        for (zz in 1:nrow(tiempocensura)) {
            if (tiempocensura[zz, 1] <= Tiempo1[z, 1]) 
                cuenta[z, 1] <- cuenta[z, 1] + 1
        }
    }
    for (z in 1:m - 1) {
        Ncensura[z, 1] <- -cuenta[z, 1] + cuenta[z + 1, 1]
    }
    Ariesgot0 <- N
    Ariesgo[1, 1] <- N - cuenta[1, 1]
    AriesgoG1tiempoG1 <- t(fit1[[Nivelesdefactores[1, 1]]]$AtRisk) %*% 
        matrix(1, nrow(fit1[[Nivelesdefactores[1, 1]]]$AtRisk), 
            1)
    for (z in 1:m) {
        for (zz in 1:m1) {
            if (X[z, 1] == Y[zz, 1]) 
                Ariesgo[z, 1] <- AriesgoG1tiempoG1[zz, 1]
        }
    }
    for (z in 2:m) {
        if (Ariesgo[z, 1] == 0) 
            Ariesgo[z, 1] <- Ariesgo[z - 1, 1] - Ncensura[z - 
                1, 1] - n.eventG1[z - 1, 1]
    }
    for (z in 1:m) {
        AriesgoG1[z, 1] <- Ariesgo[z, 1]
    }
    AriesgoTotal <- t(fit2$AtRisk) %*% matrix(1, nrow(fit2$AtRisk), 
        1)
    NTotal <- fit2$n + t(fit2$m) %*% matrix(1, nrow(fit2$m), 
        1)
    G1 <- AriesgoG1 * (AriesgoTotal - AriesgoG1) * n.event * 
        (AriesgoTotal - n.event)
    Varianzan.eventG1 <- G1/((AriesgoTotal * AriesgoTotal) * 
        (AriesgoTotal - matrix(1, nrow(n.event), 1)))
    for (jj in 1:nrow(Varianzan.eventG1)) if (Varianzan.eventG1[jj, 
        1] == "NaN") 
        Varianzan.eventG1[jj, 1] <- 0
    ValoresesperadoG1 <- AriesgoG1 * n.event * (1/AriesgoTotal)
    Diferencias <- n.eventG1 - ValoresesperadoG1
    GruporeferenciaGT <- cbind(Tiempo1, AriesgoTotal, n.event)
    GruporeferenciaG1 <- cbind(Tiempo1, AriesgoG1, n.eventG1)
    Gruporeferencia <- cbind(Tiempo1, ValoresesperadoG1, Varianzan.eventG1)
    GruporeferenciaGG <- cbind(GruporeferenciaGT, AriesgoG1, 
        n.eventG1)
    supervivencia <- matrix(fit2$survfunc)
    nnn <- nrow(supervivencia)
    Supervivencia <- supervivencia
    switch(type, LRrec = {
        Pesos <- matrix(1, nrow(n.event), 1)
        Numerador <- t(Pesos * Diferencias) %*% matrix(1, nrow(n.event), 
            1)
        Denominador <- t((Pesos * Pesos) * Varianzan.eventG1) %*% 
            matrix(1, nrow(n.event), 1)
        EstadisticoWilcoxonChi.cuadrado <- Numerador * Numerador/Denominador
        p.valorW <- 1 - pchisq(EstadisticoWilcoxonChi.cuadrado, 
            df = 1)
        Nomb.Est <- matrix(c("LRrec"))
        Estadisticos <- matrix(c(EstadisticoWilcoxonChi.cuadrado))
        p.valores <- round(matrix(c(p.valorW)), 7)
        tabla <- data.frame(Nomb.Est, Chi.cuadrado = Estadisticos, 
            p.valor = p.valores)
        print(tabla)
    }, Grec = {
        PesosGehan <- GruporeferenciaGG[, 2]
        NumeradorGehan <- t(PesosGehan) %*% Diferencias
        DenominadorGehan <- t((PesosGehan * PesosGehan) * Varianzan.eventG1) %*% 
            matrix(1, nrow(n.event), 1)
        EstadisticoGehan <- NumeradorGehan * NumeradorGehan/DenominadorGehan
        p.valorGehan <- 1 - pchisq(EstadisticoGehan, df = 1)
        Nomb.Est <- matrix(c("Grec"))
        Estadisticos <- matrix(c(EstadisticoGehan))
        p.valores <- round(matrix(c(p.valorGehan)), 7)
        tabla <- data.frame(Nomb.Est, Chi.cuadrado = Estadisticos, 
            p.valor = p.valores)
        print(tabla)
    }, TWrec = {
        PesosTaroneWare <- sqrt(GruporeferenciaGG[, 2])
        NumeradorTaroneWare <- t(PesosTaroneWare) %*% Diferencias
        DenominadorTaroneWare <- t((PesosTaroneWare * PesosTaroneWare) * 
            Varianzan.eventG1) %*% matrix(1, nrow(n.event), 1)
        EstadisticoTaroneWare <- NumeradorTaroneWare * NumeradorTaroneWare/DenominadorTaroneWare
        p.valorTaroneWare <- 1 - pchisq(EstadisticoTaroneWare, 
            df = 1)
        Nomb.Est <- matrix(c("TWrec"))
        Estadisticos <- matrix(c(EstadisticoTaroneWare))
        p.valores <- round(matrix(c(p.valorTaroneWare)), 7)
        tabla <- data.frame(Nomb.Est, Chi.cuadrado = Estadisticos, 
            p.valor = p.valores)
        print(tabla)
    }, PPrec = {
        PesosPetoPeto <- supervivencia
        NumeradorPetoPeto <- t(PesosPetoPeto) %*% Diferencias
        DenominadorPetoPeto <- t((PesosPetoPeto * PesosPetoPeto) * 
            Varianzan.eventG1) %*% matrix(1, nrow(n.event), 1)
        EstadisticoPetoPeto <- NumeradorPetoPeto * NumeradorPetoPeto/DenominadorPetoPeto
        p.valorPetoPeto <- 1 - pchisq(EstadisticoPetoPeto, df = 1)
        Nomb.Est <- matrix(c("PPrec"))
        Estadisticos <- matrix(c(EstadisticoPetoPeto))
        p.valores <- round(matrix(c(p.valorPetoPeto)), 7)
        tabla <- data.frame(Nomb.Est, Chi.cuadrado = Estadisticos, 
            p.valor = p.valores)
        print(tabla)
    }, PMrec = {
        supervivencia <- matrix(fit2$survfunc)
        nnn <- nrow(supervivencia)
        factor <- (GruporeferenciaGG[, 2]/(GruporeferenciaGG[, 
            2] + matrix(1, nnn, 1)))
        Supervivencia <- factor * supervivencia
        PesosPetoMod <- Supervivencia
        NumeradorPetoMod <- t(PesosPetoMod) %*% Diferencias
        DenominadorPetoMod <- t((PesosPetoMod * PesosPetoMod) * 
            Varianzan.eventG1) %*% matrix(1, nrow(n.event), 1)
        EstadisticoPetoMod <- NumeradorPetoMod * NumeradorPetoMod/DenominadorPetoMod
        p.valorPetoMod <- 1 - pchisq(EstadisticoPetoMod, df = 1)
        Nomb.Est <- matrix(c("PMrec"))
        Estadisticos <- matrix(c(EstadisticoPetoMod))
        p.valores <- round(matrix(c(p.valorPetoMod)), 7)
        tabla <- data.frame(Nomb.Est, Chi.cuadrado = Estadisticos, 
            p.valor = p.valores)
        print(tabla)
    }, PPrrec = {
        supervivencia <- matrix(fit2$survfunc)
        nnn <- nrow(supervivencia)
        Supervivencia <- supervivencia
        supe <- Supervivencia
        PesosPetoPrentice <- matrix(c(1, supe[1:(nnn - 1), 1]), 
            nnn, 1)
        NumeradorPetoPrentice <- t(PesosPetoPrentice) %*% Diferencias
        DenominadorPetoPrentice <- t((PesosPetoPrentice * PesosPetoPrentice) * 
            Varianzan.eventG1) %*% matrix(1, nrow(n.event), 1)
        EstadisticoPetoPrentice <- NumeradorPetoPrentice * NumeradorPetoPrentice/DenominadorPetoPrentice
        p.valorPetoPrentice <- 1 - pchisq(EstadisticoPetoPrentice, 
            df = 1)
        Nomb.Est <- matrix(c("PPrrec"))
        Estadisticos <- matrix(c(EstadisticoPetoPrentice))
        p.valores <- round(matrix(c(p.valorPetoPrentice)), 7)
        tabla <- data.frame(Nomb.Est, Chi.cuadrado = Estadisticos, 
            p.valor = p.valores)
        print(tabla)
    }, HFrec = {
        PesosPetoPeto <- Supervivencia
        PesosHF <- (PesosPetoPeto)^(alfa) * (matrix(1, nrow(n.event)) - 
            PesosPetoPeto)^(beta)
        NumeradorHF <- t(PesosHF) %*% Diferencias
        DenominadorHF <- t((PesosHF * PesosHF) * Varianzan.eventG1) %*% 
            matrix(1, nrow(n.event), 1)
        EstadisticoHF <- NumeradorHF * NumeradorHF/DenominadorHF
        p.valorHF <- 1 - pchisq(EstadisticoHF, df = 1)
        Nomb.Est <- matrix(c("HFrec"))
        Estadisticos <- matrix(c(EstadisticoHF))
        p.valores <- round(matrix(c(p.valorHF)), 7)
        tabla <- data.frame(Nomb.Est, Chi.cuadrado = Estadisticos, 
            p.valor = p.valores)
        print(tabla)
    }, CMrec = {
        PesosPetoPeto <- Supervivencia
        PesosGehan <- GruporeferenciaGG[, 2]
        PesosCM <- (PesosPetoPeto)^(alfa) * (matrix(1, nrow(n.event)) - 
            PesosPetoPeto)^(beta) * (PesosGehan)^(gamma) * (PesosGehan + 
            matrix(1, nrow(n.event), 1))^(eta)
        NumeradorCM <- t(PesosCM) %*% Diferencias
        DenominadorCM <- t((PesosCM * PesosCM) * Varianzan.eventG1) %*% 
            matrix(1, nrow(n.event), 1)
        EstadisticoCMabcd <- NumeradorCM * NumeradorCM/DenominadorCM
        p.valorCMabcd <- 1 - pchisq(EstadisticoCMabcd, df = 1)
        Nomb.Est <- matrix(c("CMrec"))
        Estadisticos <- matrix(c(EstadisticoCMabcd))
        p.valores <- round(matrix(c(p.valorCMabcd)), 7)
        tabla <- data.frame(Nomb.Est, Chi.cuadrado = Estadisticos, 
            p.valor = p.valores)
        print(tabla)
    }, Mrec = {
        PesosGehan <- GruporeferenciaGG[, 2]
        PesosCarlos <- n.event
        NumeradorCarlos <- t(PesosCarlos) %*% Diferencias
        DenominadorCarlos <- t((PesosCarlos * PesosCarlos) * 
            Varianzan.eventG1) %*% matrix(1, nrow(n.event), 1)
        EstadisticoCarlos <- NumeradorCarlos * NumeradorCarlos/DenominadorCarlos
        p.valorCarlos <- 1 - pchisq(EstadisticoCarlos, df = 1)
        Nomb.Est <- matrix(c("Mrec"))
        Estadisticos <- matrix(c(EstadisticoCarlos))
        p.valores <- round(matrix(c(p.valorCarlos)), 7)
        tabla <- data.frame(Nomb.Est, Chi.cuadrado = Estadisticos, 
            p.valor = p.valores)
        print(tabla)
    }, all = {
        Pesos <- matrix(1, nrow(n.event), 1)
        Diferencias <- n.eventG1 - ValoresesperadoG1
        Numerador <- t(Pesos * Diferencias) %*% matrix(1, nrow(n.event), 
            1)
        Denominador <- t((Pesos * Pesos) * Varianzan.eventG1) %*% 
            matrix(1, nrow(n.event), 1)
        EstadisticoWilcoxonChi.cuadrado <- Numerador * Numerador/Denominador
        p.valorW <- 1 - pchisq(EstadisticoWilcoxonChi.cuadrado, 
            df = 1)
        PesosGehan <- GruporeferenciaGG[, 2]
        NumeradorGehan <- t(PesosGehan) %*% Diferencias
        DenominadorGehan <- t((PesosGehan * PesosGehan) * Varianzan.eventG1) %*% 
            matrix(1, nrow(n.event), 1)
        EstadisticoGehan <- NumeradorGehan * NumeradorGehan/DenominadorGehan
        p.valorGehan <- 1 - pchisq(EstadisticoGehan, df = 1)
        PesosTaroneWare <- sqrt(GruporeferenciaGG[, 2])
        NumeradorTaroneWare <- t(PesosTaroneWare) %*% Diferencias
        DenominadorTaroneWare <- t((PesosTaroneWare * PesosTaroneWare) * 
            Varianzan.eventG1) %*% matrix(1, nrow(n.event), 1)
        EstadisticoTaroneWare <- NumeradorTaroneWare * NumeradorTaroneWare/DenominadorTaroneWare
        p.valorTaroneWare <- 1 - pchisq(EstadisticoTaroneWare, 
            df = 1)
        supervivencia <- matrix(fit2$survfunc)
        nnn <- nrow(supervivencia)
        Supervivencia <- supervivencia
        supe <- Supervivencia
        PesosPetoPeto <- Supervivencia
        NumeradorPetoPeto <- t(PesosPetoPeto) %*% Diferencias
        DenominadorPetoPeto <- t((PesosPetoPeto * PesosPetoPeto) * 
            Varianzan.eventG1) %*% matrix(1, nrow(n.event), 1)
        EstadisticoPetoPeto <- NumeradorPetoPeto * NumeradorPetoPeto/DenominadorPetoPeto
        p.valorPetoPeto <- 1 - pchisq(EstadisticoPetoPeto, df = 1)
        supervivencia <- matrix(fit2$survfunc)
        nnn <- nrow(supervivencia)
        Supervivencia <- supervivencia
        supe <- Supervivencia
        factor <- (GruporeferenciaGG[, 2]/(GruporeferenciaGG[, 
            2] + matrix(1, nnn, 1)))
        Supervivencia <- factor * supervivencia
        PesosPetoMod <- Supervivencia
        NumeradorPetoMod <- t(PesosPetoMod) %*% Diferencias
        DenominadorPetoMod <- t((PesosPetoMod * PesosPetoMod) * 
            Varianzan.eventG1) %*% matrix(1, nrow(n.event), 1)
        EstadisticoPetoMod <- NumeradorPetoMod * NumeradorPetoMod/DenominadorPetoMod
        p.valorPetoMod <- 1 - pchisq(EstadisticoPetoMod, df = 1)
        supervivencia <- matrix(fit2$survfunc)
        nnn <- nrow(supervivencia)
        Supervivencia <- supervivencia
        supe <- Supervivencia
        PesosPetoPrentice <- matrix(c(1, supe[1:(nnn - 1), 1]), 
            nnn, 1)
        NumeradorPetoPrentice <- t(PesosPetoPrentice) %*% Diferencias
        DenominadorPetoPrentice <- t((PesosPetoPrentice * PesosPetoPrentice) * 
            Varianzan.eventG1) %*% matrix(1, nrow(n.event), 1)
        EstadisticoPetoPrentice <- NumeradorPetoPrentice * NumeradorPetoPrentice/DenominadorPetoPrentice
        p.valorPetoPrentice <- 1 - pchisq(EstadisticoPetoPrentice, 
            df = 1)
        PesosHF <- (PesosPetoPeto)^(alfa) * (matrix(1, nrow(n.event)) - 
            PesosPetoPeto)^(beta)
        NumeradorHF <- t(PesosHF) %*% Diferencias
        DenominadorHF <- t((PesosHF * PesosHF) * Varianzan.eventG1) %*% 
            matrix(1, nrow(n.event), 1)
        EstadisticoHF <- NumeradorHF * NumeradorHF/DenominadorHF
        p.valorHF <- 1 - pchisq(EstadisticoHF, df = 1)
        PesosCM <- (PesosPetoPeto)^(alfa) * (matrix(1, nrow(n.event)) - 
            PesosPetoPeto)^(beta) * (PesosGehan)^(gamma) * (PesosGehan + 
            matrix(1, nrow(n.event), 1))^(eta)
        NumeradorCM <- t(PesosCM) %*% Diferencias
        DenominadorCM <- t((PesosCM * PesosCM) * Varianzan.eventG1) %*% 
            matrix(1, nrow(n.event), 1)
        EstadisticoCM <- NumeradorCM * NumeradorCM/DenominadorCM
        p.valorCM <- 1 - pchisq(EstadisticoCM, df = 1)
        PesosCarlos <- n.event
        NumeradorCarlos <- t(PesosCarlos) %*% Diferencias
        DenominadorCarlos <- t((PesosCarlos * PesosCarlos) * 
            Varianzan.eventG1) %*% matrix(1, nrow(n.event), 1)
        EstadisticoCarlos <- NumeradorCarlos * NumeradorCarlos/DenominadorCarlos
        p.valorCarlos <- 1 - pchisq(EstadisticoCarlos, df = 1)
        Nomb.Est <- matrix(c("LRrec   ", "Grec    ", "TWrec   ", 
            "PPrec   ", "PMrec   ", "PPrrec  ", "HFrec   ", "CMrec   ", 
            "Mrec    "))
        Estadisticos <- matrix(c(EstadisticoWilcoxonChi.cuadrado, 
            EstadisticoGehan, EstadisticoTaroneWare, EstadisticoPetoPeto, 
            EstadisticoPetoMod, EstadisticoPetoPrentice, EstadisticoHF, 
            EstadisticoCM, EstadisticoCarlos))
        p.valores <- round(matrix(c(p.valorW, p.valorGehan, p.valorTaroneWare, 
            p.valorPetoPeto, p.valorPetoMod, p.valorPetoPrentice, 
            p.valorHF, p.valorCM, p.valorCarlos)), 7)
        tabla <- data.frame(Nomb.Est, Chi.square = Estadisticos, 
            p.value = p.valores)
        print(tabla)
    }, ` ` = print("Arguments of Dif.Surv.Rec"))
}
