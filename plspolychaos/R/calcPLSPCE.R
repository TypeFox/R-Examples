###############################################
# Compute PLS-PCE sensivity indexes
###############################################

calcPLSPCE <- function(pce, nc = 2) {
    #############################################
    ## regression PLS (idem sivipm package)
    ############################################
    Y <- pce@.Data[, ncol(pce@.Data), drop=FALSE]
    # XM: inputs without Y but with the constant term
    XM <- pce@.Data[, -ncol(pce@.Data)]
    # data.exp: inputs without Y and without the constant term
    data.exp <- pce@.Data[, -c(1, ncol(pce@.Data))]
    
    regpls <- regpls2(Y, data.exp, nc)
    # Ecrire les R2
    R2pourcent <- (regpls$R2y * 100)/sum(regpls$R2y)
    R2pourcentcum <- cumsum(R2pourcent)
    R2 <- matrix(c(regpls$R2y, R2pourcent, R2pourcentcum), ncol = 3, dimnames = list(rownames(regpls$R2y), 
        c("R2", "%R2", "%R2cumulated")))

    ## cat('\nExplanation level per component, percentage and cumulated
    ## percentage\n')
    ## print(R2)
    
    # Ecrire les Q2
    Q2 <- matrix(c(regpls$Q2, regpls$Q2cum[, ncol(regpls$Q2cum)]), ncol = 2, dimnames = list(rownames(regpls$Q2), 
        c("Q2", "Q2cum")))

    # Il ne doit pas y avoir de Q2<0
    Q2[Q2[,"Q2"]<0, "Q2"] <- 0
    q2zero <- which(Q2[,"Q2"]==0) # indice des zeros
    if (length(q2zero) == 0) {
        stop("\nIncrease the number of components")
    }
    q2zero <- min(q2zero) # indice du 1er zero
    # Les Q2 ne doivent pas remonter apres le 1er zero
    Q2[q2zero:nrow(Q2), "Q2"] <- Q2[q2zero, "Q2"]
    
    
    ## cat('\nExplanation-prediction level per component and cumulated\n')
    ## print(Q2)
    
    # Calculer RMSEP
    rmsep <- sqrt(regpls$PRESS/length(Y))
    colnames(rmsep) <- "rmsep"
    ## cat('\nrmsep\n'); print(rmsep)
    
    y.hat <- regpls$y.hat
    
    
    
    #############################################
    ## determiner la composante optimale
    ## c'est l'indice du dernier Q2 non nul
    ############################################
    ncopt <- q2zero - 1
    if (ncopt == 0) {
        stop("\nInternal error: all the Q2 are zeros !!!")
    }
    
#    cat("\nNumber of the optimal component: ", ncopt, "\n")
    #############################################
    ## Afficher les resultats correspondants a la composante optimale
    ############################################
    R2opt <- matrix(R2[ncopt, ], ncol = ncol(R2), dimnames = list(rownames(R2)[ncopt], 
        colnames(R2)))
    
    ## cat('\nExplanation level of the optimal component\n'); print(R2opt)
    
    Q2opt <- matrix(Q2[ncopt, ], ncol = ncol(Q2), dimnames = list(rownames(Q2)[ncopt], 
        colnames(Q2)))
    ## cat('\nExplanation-prediction level of the optimal component\n'); print(Q2opt)
    
    rmsepopt <- matrix(rmsep[ncopt, ], ncol = ncol(rmsep), dimnames = list(rownames(rmsep)[ncopt], 
        colnames(rmsep)))
    ## print(rmsepopt)

    #############################################
    ## On calcule les betas naturels de toutes les composantes
    #############################################
    
    # calcul des betaNat
##     betaNat <- matrix(NA, nrow=(nrow(regpls$ret$beta)+1), ncol=nc)
##     for (hcur in 1:nc) {
##       betaNath <- calcbetaNat(regpls$ret$beta[, hcur, drop=FALSE], Y, data.exp, regpls$mu.x, regpls$sd.x,regpls$sd.y)
##       betaNath <- c(betaNath$betaNat0, betaNath$betaNat)
##       betaNat[, hcur] <-betaNath
##     }

    oneBeta <- function(beta, Y, data.exp, mu.x, sd.x, sd.y) {
      # fonction pour apply
      betaNath <- calcbetaNat(as.matrix(beta), Y, data.exp, mu.x, sd.x, sd.y)
      return(c(betaNath$betaNat0, betaNath$betaNat))
    }

    betaNat <- apply(regpls$ret$beta, 2,  oneBeta,
               Y, data.exp, regpls$mu.x, regpls$sd.x,regpls$sd.y)
    
      
    #############################################
    ## labeller les monomes dans les beta
    ## pour les identifier
    #############################################
    rownames(betaNat) <- rownames(pce@STRUC)
    betaCR <- regpls$ret$beta
    rownames(betaCR) <- rownames(pce@STRUC)[-1]
    
    
##     #############################################
##     ## On calcule les SI (idem chaosbasics) a partir
##      du betaNat[, ncopt]
##     #############################################
##     

    indicesopt <- calcPLSPCESI(betaNat[,ncopt, drop=FALSE], data.exp, XM, Y,
                         regpls$mu.x,  regpls$sd.x,regpls$sd.y,
                     pce)

# cat('\nINDICES PLS-PCE \n'); print(indicesopt)
#############################################
# les SI en pourcentage
#############################################
    indicespercent <- apply(indicesopt, 2, function(X) {
        (X * 100)/sum(X)
    })
    
    # cat('\n%INDICES PLS-PCE \n') print(indicespercent)
    
    
    #############################################
    ## retour
    #############################################
    retour <- new("PLSPCE", rmsep = rmsep, y.hat = y.hat,
                  COEF = betaNat, R2 = R2, 
        Q2 = Q2, ncopt = ncopt, indexes = indicesopt, indexes.percent = indicespercent, betaCR=betaCR,
         STRUC = pce@STRUC)

    return(retour)
    
    
} 
##################################################
# Compute the PLS-PCE sensivity indexes (idem chaosbasics)
# for the optimal component
##################################################
calcPLSPCESI <- function(betaNat, data.exp, XM, Y,
                         mu.x,  sd.x,sd.y, pce) {
  # betaNat:  centered-reducted beta of the optimal
  # (nmono +1)
  
    
    indices <- indexes(XM, pce@nvx, betaNat, pce@STRUC@.Data)$indexes

  return(indices)
}
