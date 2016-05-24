fInclin <- function(compI, angGen, iS = 2, alb = 0.2, horizBright = TRUE, HCPV = FALSE){
    ##compI es class='G0'
    ##angGen es 'zoo', resultado de fTheta

    ##Preparo argumentos
    stopifnot(iS %in% 1:4)
    ang <- as.data.frame(coredata(angGen))
    Beta <- ang$Beta
    Alfa <- ang$Alfa
    cosTheta <- ang$cosTheta
  
    comp <- as.data.frameI(compI, complete=TRUE)
    aman <- comp$aman
    B0 <- comp$B0
    Bo0 <- comp$Bo0
    D0 <- comp$D0
    G0 <- comp$G0
    cosThzS <- comp$cosThzS
    is.na(cosThzS) <- !aman
    
    ##Método N.Martin para suciedad e incidencia no perpendicular
    Suc <- rbind(c(1, 0.17, -0.069),
                 c(0.98,.2,-0.054),
                 c(0.97,0.21,-0.049),
                 c(0.92,0.27,-0.023))
    FTb <- (exp(-cosTheta/Suc[iS,2]) - exp(-1/Suc[iS,2]))/(1 - exp(-1/Suc[iS,2]))
    FTd <- exp(-1/Suc[iS,2] * (4/(3*pi) * (sin(Beta) + (pi - Beta - sin(Beta))/(1 + cos(Beta))) +
                               Suc[iS,3] * (sin(Beta) + (pi - Beta - sin(Beta))/(1 + cos(Beta)))^2))
    FTr <- exp(-1/Suc[iS,2] * (4/(3*pi) * (sin(Beta) + (Beta - sin(Beta))/(1 - cos(Beta))) +
                               Suc[iS,3] * (sin(Beta) + (Beta - sin(Beta))/(1 - cos(Beta)))^2))
    
    ##Metodo Hay and Davies para tratamiento difusa
    B <- B0 * cosTheta/cosThzS * (cosThzS>0.007) #El factor cosThzS>0.007 hace falta para eliminar resultados erroneos cerca del amanecer
    k1 <- B0/(Bo0)
    Di <- D0 * (1-k1) * (1+cos(Beta))/2
    if (horizBright) Di <- Di * (1+sqrt(B0/G0) * sin(Beta/2)^3)
    Dc <- D0 * k1 * cosTheta/cosThzS * (cosThzS>0.007)
    R <- alb * G0 * (1-cos(Beta))/2
    D <- (Di + Dc)
    ##Irradiancia extraterrestre en el plano inclinado
    Bo <- Bo0 * cosTheta/cosThzS * (cosThzS>0.007) 
    ##Irradiancia directa normal (DNI)
    Bn <- B0/cosThzS
    ##Suma de componentes
    G <- B + D + R
    Ref <- R * Suc[iS,1] * (1-FTr) * (!HCPV)
    Ref[is.nan(FTr)] <- 0     #Cuando cos(Beta)=1, FTr=NaN. Anulo Ref.
    Dief <- Di * Suc[iS,1] * (1 - FTd) * (!HCPV)
    Dcef <- Dc * Suc[iS,1] * (1 - FTb) * (!HCPV)
    Def <- Dief + Dcef
    Bef <- B * Suc[iS,1] * (1 - FTb)
    Gef <- Bef + Def + Ref

    ##Resultado
    result <- zoo(data.frame(Bo, Bn,
                             G, D, Di, Dc, B, R,
                             FTb, FTd, FTr,
                             Dief, Dcef, Gef, Def, Bef, Ref), 
                  order.by = indexI(compI))
    ## Use 0 instead of NA for irradiance values
    result[!aman,] <- 0
    result
}
