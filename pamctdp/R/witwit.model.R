# -------------------------------------------------------
# función para with-with con un modelo C o B
# modificación de
# witwit.coa de ade4
# Campo E. Pardo - Julio 3/05
# Septiembre 21 se agrega factor de homotecia en salida
# PARAMETROS
# weight:   "coa" pesos del ACS, 
#       "mfa" pesos del AFM, 
#       "mfac" pesos del AFM solo para columnas,
#       "mfar" pesos del AFM solo para filas,
#       "wclo" pesos dentro de las nubes
#       "mfawclo" pesos iniciales dentro de las nubes y 
#               luego pesos AFM
# model: "C" modelo del ACI, "B" modelo intra-bloques 
# Mayo 19/08 se modifica para el cálculo de X se fijan las marginales de F
# los pesos que entran solo se usan en las matrices de métrica y pesos
# modificaciones para procedimiento iterativo en opción "afm"
# entran nuevos parámetros eps = 1e-15; iter=100
# Actualizada octubre 8/2013
#---------------------------------------------------------
"witwit.model" <- function (dudi, row.blocks, col.blocks,
        pfil=dudi$lw,pcol=dudi$cw,model="C",weight="coa",
        scannf = TRUE, nf = 2,eps=1e-15,iter=100) 
{
    if (!inherits(dudi, "coa")) stop("Object of class coa expected")
    model <- toupper(model) 
    if (model != "C" & model != "B") 
       {    cat("\n Model ", model, " is not a valid option, model C  is used \n")
        model <- "C"
       }      
    if (weight != "coa" & weight != "mfa" &  weight != "mfac" & weight != "mfar") 
       {    cat("\n Weight ", weight, " is not a valid option, weight coa is used \n")
        weight <- "coa"
       }     
   
    lig <- nrow(dudi$tab) # row number
    col <- ncol(dudi$tab) # column number
    row.fac <- rep(1:length(row.blocks),row.blocks) # row factor
    col.fac <- rep(1:length(col.blocks),col.blocks) # column factor
    if (length(col.fac)!=col) stop ("Non convenient col.fac")
    if (length(row.fac)!=lig) stop ("Non convenient row.fac")
    tabinit <- as.matrix(eval(as.list(dudi$call)$df, sys.frame(0))) # contingency table
    F <- tabinit/sum(tabinit) # table F
    #---------------------------------------------
    f.kljmat <- rowsum(F,row.fac, reorder = FALSE)[row.fac,]
    f..l.vec <- tapply(dudi$lw,row.fac,sum)[row.fac]
    f..l.vec <- as.numeric(f..l.vec)
    f.k.j_f..l.vec <- dudi$lw/f..l.vec
    AL <-   f.kljmat*f.k.j_f..l.vec
    # -----------------------------
    fi.ljmat <- rowsum(t(F),col.fac, reorder = FALSE)[col.fac,]
    f...jvec <- tapply(dudi$cw,col.fac,sum)[col.fac]
    f...jvec <- as.numeric(f...jvec)
    f.k.j_f...jvec <- dudi$cw/f...jvec
    AJ <- t(fi.ljmat*f.k.j_f...jvec)
    # -----------------------------------------
    f..ljmat <- rowsum(F,row.fac, reorder = FALSE)
    f..ljmat <- t(rowsum(t(f..ljmat),col.fac, reorder = FALSE))
    f..ljmat <-f..ljmat[row.fac,col.fac]
    AG <- f..ljmat*f.k.j_f..l.vec
    AG <- t(t(AG)*f.k.j_f...jvec)
    C <- AL+AJ-AG
    #C[C < 0] <- 0
    row.names(C) <- row.names(dudi$tab)
    # model B
    B <- f.kljmat*t(fi.ljmat)/f..ljmat
    row.names(B) <- row.names(dudi$tab)
    # subtabla de ceros -> B = 0
    B[B=="NaN"] <- 0
    if (model=="C") X <- F - C 
    if (model=="B") X <- F - B
    X <- X/dudi$lw  #X/pfil
    X <- t(t(X)/dudi$cw) #t(t(X)/pcol)
    X <- data.frame(X)
    #------------------------------------------------------------------
    # función para reponderación con inverso del primer valor propio
    # entra X, blocks, pesos y salen nuevos pesos. Bloques columna
    reweight <- function(X,blocks,pfil,pcol) { 
        X <- data.frame(X)
        # ACP bloques columna
        nblo <- length(blocks)
        col.fac <- rep(1:nblo,blocks)
        eigj <- numeric(nblo)
        for (j in 1:nblo) {
               Xj <- X[, col.fac == j]
               wj <- pcol[col.fac==j]
               eigj[j] <- as.dudi(Xj,col.w=wj,row.w=pfil,scannf=FALSE,nf=2,
                    call = match.call(),type="pca")$eig[1]
 
        }
        rew <- list(wcol=pcol/eigj[col.fac],eigen=eigj)
        return(rew)
    }# FIN de función reweight    
    #------------------------------------------------------------------
        sepeig.col <- NULL
        sepeig.row <- NULL
    
    if(weight=="coa") hom <- c(1,1) 
    #------------------------------------------------------------------
    if(weight=="mfa") {
    pcolin <- pcol 
        pfilin <- pfil
    for (it in 1:iter){
            rewcol <- reweight(X,col.blocks,pfil,pcol)
            pcol <- rewcol$wcol
            hom <- sum(pcol)
            pcol <- pcol/sum(pcol)
            rewfil <- reweight(t(X),row.blocks,pcol,pfil) #pcolin,pfil)
            pfil <- rewfil$wcol
            hom <- c(sum(pfil),hom)
            pfil <- pfil/sum(pfil)
        sum2difcol <- t(pcol-pcolin) %*% (pcol-pcolin)
        sum2diffil <- t(pfil-pfilin) %*% (pfil-pfilin)
                if (sum2difcol <= eps & sum2diffil <= eps) break
        pcolin <- pcol 
        pfilin <- pfil
    }
        cat("\n AFM iter = ",it,"; Max iter = ",iter,"\n")
        sepeig.col <- rewcol$eigen
        sepeig.row <- rewfil$eigen
     } 
    #------------------------------------------------------------------
    if(weight=="mfac") {
     #     pcolin <- pcol
        pcol <- reweight(X,col.blocks,pfil,pcol)$wcol
        hom <- sum(pcol)
        pcol <- pcol/hom
        hom <- c(1,hom)

 
     #  Xprov <- t(t(X) * pcolin / pcol)
     #      pfil <- reweight(t(X),row.blocks,pcol,pfil)$wcol
      }   
    #------------------------------------------------------------------
        if(weight=="mfar") {
    #      pfilin <- pfil
        pfil <- reweight(t(X),row.blocks,pcol,pfil)$wcol
        hom <- sum(pfil)
        pfil <- pfil/hom
        hom <- c(hom,1) 
    #   X <- X * pfilin / pfil
    #       pcol <- reweight(X,col.blocks,pfil,pcol)$wcol
       }   
    #------------------------------------------------------------------
      
  #    X <- data.frame(X)    
      ww <- as.dudi(X, pcol, pfil, scannf = scannf, nf = nf, 
        call = match.call(), type = "wwmodel")
   
    wr <- ww$li*ww$li*f.k.j_f..l.vec
    wr <- rowsum(as.matrix(wr),row.fac, reorder = FALSE)
    cha <- names(row.blocks)
    if (is.null(cha)) cha <- as.character(1:length(row.blocks))
    wr <- data.frame(wr)
    names(wr) <- names(ww$li)
    row.names(wr) <- cha
    ww$lbvar <- wr
 #   ww$lbw <- tapply(dudi$lw,row.fac,sum)
 ww$lbw <- tapply(pfil,row.fac,sum)

 #   wr <- ww$co*ww$co*wcvec
    wr <- ww$co*ww$co*f.k.j_f...jvec
    wr <- rowsum(as.matrix(wr),col.fac, reorder = FALSE)
    cha <- names(col.blocks)
    if (is.null(cha)) cha <- as.character(1:length(col.blocks))
    wr <- data.frame(wr)
    names(wr) <- names(ww$co)
    row.names(wr) <- cha
    ww$cbvar <- wr
 #   ww$cbw <- tapply(dudi$cw,col.fac,sum)
    ww$cbw <- tapply(pcol,col.fac,sum)
    ww$hom <- hom
    ww$rbl <- row.blocks
    ww$cbl <- col.blocks
    ww$sepeig.col <- sepeig.col
    ww$sepeig.row <- sepeig.row 
    if(!is.null(sepeig.col)) names(ww$sepeig.col) <- rownames(ww$cbvar)
    if(!is.null(sepeig.row)) names(ww$sepeig.row) <- rownames(ww$lbvar)
    if (model=="C") ww$model <- C
    if (model=="B") ww$model <- B
    class(ww) <- c("wwmodel","witwit", "coa", "dudi")
   return(ww)
   }
   #===========================
   # function print.wwmodel
   # modification of print.dudi
   "print.wwmodel" <- function(x, ...) 
   {
      if (!inherits(x, "wwmodel")) stop("non convenient data")
    cat("CA respect to a wiwit model\n")
    cat("class: ")
    cat(class(x))
    cat("\n$call: ")
    print(x$call)
    cat("\n$hom:",x$hom,"homotecia (rows,columns)")
    cat("\n$nf:", x$nf, "axis-components saved")
    cat("\n Subespace S = R^",x$nf,"\n")  
 
    cat("\n$rank: ")
    cat(x$rank)
    cat("\neigen values: ")
    l0 <- length(x$eig)
    cat(signif(x$eig, 4)[1:(min(5, l0))])
    if (l0 > 5) 
        cat(" ...\n")
    else cat("\n")
    
    sumry <- array("", c(7, 4), list(1:7, c("vector", "length", 
        "mode", "content")))
    sumry[1, ] <- c("$cw", length(x$cw), mode(x$cw), "column weights")
    sumry[2, ] <- c("$lw", length(x$lw), mode(x$lw), "row weights")
    sumry[3, ] <- c("$eig", length(x$eig), mode(x$eig), "eigen values")
    sumry[4, ] <- c("$lbw", length(x$lbw), mode(x$lbw), "row-band weights")
    sumry[5, ] <- c("$cbw", length(x$cbw), mode(x$cbw), "column-band weights")
    sumry[6, ] <- c("$rbl", length(x$rbl), mode(x$rbl), "row-band row numbers")
    sumry[7, ] <- c("$cbl", length(x$cbl), mode(x$cbl), "column-band column numbers")

    class(sumry) <- "table"
    cat("\n")
    print(sumry)
    cat("\n")
    sumry <- array("", c(8, 4), list(1:8, c("data.frame", "nrow", 
        "ncol", "content")))
    sumry[1, ] <- c("$tab", nrow(x$tab), ncol(x$tab), "modified array")
    sumry[2, ] <- c("$li", nrow(x$li), ncol(x$li), "row coordinates")
    sumry[3, ] <- c("$l1", nrow(x$l1), ncol(x$l1), "row normed scores")
    sumry[4, ] <- c("$co", nrow(x$co), ncol(x$co), "column coordinates")
    sumry[5, ] <- c("$c1", nrow(x$c1), ncol(x$c1), "column normed scores")
    sumry[6, ] <- c("$lbvar", nrow(x$lbvar), ncol(x$lbvar), "inertia row-bands")
    sumry[7, ] <- c("$cbvar", nrow(x$cbvar), ncol(x$cbvar), "inertia column-bands")
    sumry[8, ] <- c("$model", nrow(x$model), ncol(x$model), "used model")
  
    class(sumry) <- "table"
    print(sumry)
}
  
   #========================
   "summary.wwmodel" <- function (object, ...) {
    if (!inherits(object, "wwmodel")) 
        stop("For 'wwmodel' object")
    cat("Internal models correspondence analysis\n")
    cat("class: ")
    cat(class(object))
    cat("\n$call: ")
    print(object$call)
    cat(object$nf, "axis-components saved")
    cat("\neigen values: ")
    l0 <- length(object$eig)
    cat(signif(object$eig, 4)[1:(min(5, l0))])
    if (l0 > 5) 
        cat(" ...\n")
    else cat("\n")
    cat("\n")
    cat("Eigen value decomposition among row blocks\n")
    nf <- object$nf
    nrb <- nrow(object$lbvar)
    aa <- as.matrix(object$lbvar)
    sumry <- array("", c(nrb + 1, nf + 1), list(c(row.names(object$lbvar), 
        "mean"), c(names(object$lbvar), "weights")))
    sumry[(1:nrb), (1:nf)] <- round(aa, digits = 4)
    sumry[(1:nrb), (nf + 1)] <- round(object$lbw, digits = 4)
    sumry[(nrb + 1), (1:nf)] <- round(object$eig[1:nf], digits = 4)
    class(sumry) <- "table"
    print(sumry)
    cat("\n")
    sumry <- array("", c(nrb + 1, nf), list(c(row.names(object$lbvar), 
        "sum"), names(object$lbvar)))
    aa <- object$lbvar * object$lbw
    aa <- 1000 * t(t(aa)/object$eig[1:nf])
    sumry[(1:nrb), (1:nf)] <- round(aa, digits = 0)
    sumry[(nrb + 1), (1:nf)] <- rep(1000, nf)
    class(sumry) <- "table"
    print(sumry)
    cat("\n")
    cat("Eigen value decomposition among column blocks\n")
    nrb <- nrow(object$cbvar)
    aa <- as.matrix(object$cbvar)
    sumry <- array("", c(nrb + 1, nf + 1), list(c(row.names(object$cbvar), 
        "mean"), c(names(object$cbvar), "weights")))
    sumry[(1:nrb), (1:nf)] <- round(aa, digits = 4)
    sumry[(1:nrb), (nf + 1)] <- round(object$cbw, digits = 4)
    sumry[(nrb + 1), (1:nf)] <- round(object$eig[1:nf], digits = 4)
    class(sumry) <- "table"
    print(sumry)
    cat("\n")
    sumry <- array("", c(nrb + 1, nf), list(c(row.names(object$cbvar), 
        "sum"), names(object$cbvar)))
    aa <- object$cbvar * object$cbw
    aa <- 1000 * t(t(aa)/object$eig[1:nf])
    sumry[(1:nrb), (1:nf)] <- round(aa, digits = 0)
    sumry[(nrb + 1), (1:nf)] <- rep(1000, nf)
    class(sumry) <- "table"
    print(sumry)
    cat("\n")
}
