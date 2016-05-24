`dixon` <-
function (datos, nsim = 99, fortran=TRUE) 
{
    datos <- as.matrix(datos) #19/03/2012: in response to nndistG changes in splancs [storage.mode(pts) <- "double"]
    info = mNNinfoc(xy = datos[, 1:2], label = datos[, 3], fortran=fortran)
    datos.test = mNNtest(info)
    Ni = rowSums(info$ON)
    S = log10((info$ON/(Ni - info$ON))/(info$EN/(Ni - info$EN)))
    ON = info$ON
    EN = info$EN
    Z = datos.test$Z
    Z.obs = Z
    pZas = 2 * (ifelse(Z >= 0, 1 - pnorm(Z), pnorm(Z)))
    C = datos.test$C[1]
    C.obs = C
    Ci = datos.test$Ci[, 1]
    Ci.obs = Ci
    pCas = datos.test$C[2]
    pCias = datos.test$Ci[, 2]
    pZr = NULL
    pCr = NULL
    pCir = NULL
    if (nsim > 0) {
        for (i in 1:nsim) {
            progressreport(i,nsim)
            datos[, 3] = sample(datos[, 3])
            info = mNNinfoc(xy = datos[, 1:2], label = datos[,3],fortran=fortran)
#           info = mNNinfob(xy = datos[, 1:2], label = datos[,3])
            datos.test = mNNtest(info)
            Z = cbind(Z, datos.test$Z)
            C = c(C, datos.test$C[1])
            Ci = cbind(Ci, datos.test$Ci[, 1])
        }
        pZr = apply(Z, 1, p2colasr)
        pCr = 1 - rank(C)[1]/(length(C))
        pCir = apply(Ci, 1, function(x) 1 - rank(x)[1]/(length(x)))
    }
    St = as.data.frame(as.table(round(S, 2)))
    ONt = as.data.frame(as.table(ON))
    ENt = as.data.frame(as.table(round(EN, 2)))
    Zt = round(datos.test$Z, 2)
    round(pZas, 4)
    tableZ = cbind(ONt[order(ONt[, 1]), ], ENt[order(ENt[, 1]), 
        3], St[order(St[, 1]), 3], round(Z.obs, 2), round(pZas, 
        4))
    names(tableZ) = c("From", "To", "    Obs.Count", "    Exp. Count", 
        "S ", "Z ", "  p-val.as")
    if (length(pZr) != 0) {
        tableZ = cbind(tableZ, round(pZr, 4))
        names(tableZ) = c(names(tableZ)[-8], "  p-val.rnd")
    }
    rownames(tableZ) = NULL
    k = length(unique(datos[, 3]))
    df = c(k * (k - 1), rep(k - 1, k))
    nombres.test = c("Overall segregation", paste("From ", dimnames(EN)[[1]], 
        "          "))
    tablaC = data.frame(cbind(df, round(c(C.obs, Ci.obs), 2), 
        round(c(pCas, pCias), 4)))
    row.names(tablaC) = nombres.test
    names(tablaC) = c("  df ", "Chi-sq", "P.asymp")
    if (length(pCir) != 0) {
        tablaC = cbind(tablaC, round(c(pCr, pCir), 4))
        names(tablaC) = c(names(tablaC)[-4], "  P.rand")
    }
    return(list(ON = ON, EN = EN, Z = Z.obs, S = S, pZas = pZas, 
        pZr = pZr, C = C.obs, Ci = Ci.obs, pCas = pCas, pCias = pCias, 
        pCr = pCr, pCir = pCir, tablaZ = tableZ, tablaC = tablaC))
}
