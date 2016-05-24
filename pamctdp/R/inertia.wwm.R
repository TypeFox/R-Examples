#----------------------------------------------------------------------------------
# Contribuciones de las subnubes a la Inercia en wwmodel
#----------------------------------------------------------------------------------
# Función inertia.wwm. Campo Elías Pardo
# Ayudas ACI (no parciales) Creada Octubre 04, modificada Feb-23-05 y agosto 05
# organiza las ayudas dadas en witwit.coa
# ENTRA
# ACww salida de witwit.model; 
# xax=eje horizontal, yax = eje vertical. Default es 1 y 2
# SALE
# porcentaje de inercia de la subnube fila(columna) al eje
# calidad de representación de la subnube,
# contribución a la inercia de los ejes y total
#tabla con dos ejes
# LAS INERCIAS ESTAN MULTIPLICADAS POR ti EN LA TABLA
# Ago 08/05 se incluye inercia de subtablas l,j
# Se actualiza a marzo 7/07 con cambio de nombre
# dec número de decimales en tablas $fil y $col
# ago 15 2010 se incluye parametro para multiplicar inercia ti antes
# estaba fijo en 1000
#---------------------------------------------------------------------------------
inertia.wwm <- function(ACww,xax=1,yax=2,dec=1,ti=1000)
{
	rbl <- ACww$rbl
      cbl <- ACww$cbl
    J <- nrow(ACww$cbvar); L <- nrow(ACww$lbvar)
    nf <- ACww$nf
   # inercia total subnubes (no parciales)
    X <- as.matrix(ACww$tab)
    normK <- apply(X*X*ACww$lw,2,sum);normI <- apply(t(X*X)*ACww$cw,2,sum)
    cbl.fac <- rep(1:J,cbl);rbl.fac <- rep(1:L,rbl)
    # columnas
    inJ <- tapply(normK*ACww$cw,cbl.fac,sum);names(inJ) <- names(cbl)
    cinJ <- inJ/sum(inJ)*100
    #filas
    inL <- tapply(normI*ACww$lw,rbl.fac,sum);names(inL) <- names(rbl)
    cinL <- inL/sum(inL)*100
    # calidad o contribucion relativa
    Unonf <- t(rep(1,nf))
    # tabla para dos ejes
    icol <- NULL
    icol <- cbind(icol,inJ*ti,cinJ) 
    icol <- cbind(icol,((ACww$cbvar * (ACww$cbw %*% t(rep(1,nf))))[xax])*ti)
    
    inertia <- NULL
    inertia$cloud.col <- cinJ
    inertia$cloud.row <- cinL
    inertia$rel.col <-   ((ACww$cbvar * (ACww$cbw %*% t(rep(1,nf))))/
                                                    (inJ %*% Unonf))*100
    inertia$rel.row <-   ((ACww$lbvar * (ACww$lbw %*% t(rep(1,nf))))/
                                                    (inL %*% Unonf))*100
    inertia$abs.col <- ((ACww$cbvar * (ACww$cbw %*% t(rep(1,nf))))/
                                    (rep(1,J) %*% t(ACww$eig[1:nf])))*100
    inertia$abs.row <- ((ACww$lbvar * (ACww$lbw %*% t(rep(1,nf))))/
                                    (rep(1,L) %*% t(ACww$eig[1:nf])))*100
    # tabla para dos ejes COLUMNAS
    icol <- NULL
    # I total
    icol <- cbind(icol,inJ*ti,cinJ); colnames(icol) <- c("Inertia(j)","Contr(j)") 
    # comp 1
    icol <- cbind(icol,((ACww$cbvar * (ACww$cbw %*% t(rep(1,nf))))[xax])*ti)
    colnames(icol)[3] <- paste("Inertia",xax,sep="")
    icol <- cbind(icol,inertia$abs.col[,xax])
    colnames(icol)[4] <- paste("Contr",xax,sep="") 
    icol <- cbind(icol,inertia$rel.col[,xax])
    colnames(icol)[5] <- paste("Qual",xax,sep="") 
    # comp 2
     icol <- cbind(icol,((ACww$cbvar * (ACww$cbw %*% t(rep(1,nf))))[yax])*ti)
    colnames(icol)[6] <- paste("Inertia",yax,sep="")
    icol <- cbind(icol,inertia$abs.col[,yax])
    colnames(icol)[7] <- paste("Contr",yax,sep="") 
    icol <- cbind(icol,inertia$rel.col[,yax])
    colnames(icol)[8] <- paste("Qual",yax,sep="") 
    # calidad plano
    icol <- cbind(icol,inertia$rel.col[,xax]+inertia$rel.col[,yax])
    colnames(icol)[9] <- paste("QualPlane",xax,"-",yax,sep="") 
    # peso
    icol <- cbind(icol,ACww$cbw*100)
    colnames(icol)[10] <- "Weight(j)" 
    # sumas 
    suma <- colSums(icol);suma[5] <- NA;suma[8:9] <- NA
    icol <- rbind(icol,suma)
    rownames(icol)[(J+1)] <- "Total" 
    # tabla para dos ejes FILAS
    ifil <- NULL
    # I total
    ifil <- cbind(ifil,inL*ti,cinL); colnames(ifil) <- c("Inertia(l)","Contr(l)") 
    # comp 1
    ifil <- cbind(ifil,((ACww$lbvar * (ACww$lbw %*% t(rep(1,nf))))[xax])*ti)
    colnames(ifil)[3] <- paste("Inertia",xax,sep="")
    ifil <- cbind(ifil,inertia$abs.row[,xax])
    colnames(ifil)[4] <- paste("Contr",xax,sep="") 
    ifil <- cbind(ifil,inertia$rel.row[,xax])
    colnames(ifil)[5] <- paste("Qual",xax,sep="") 
    # comp 2
     ifil <- cbind(ifil,((ACww$lbvar * (ACww$lbw %*% t(rep(1,nf))))[yax])*ti)
    colnames(ifil)[6] <- paste("Inertia",yax,sep="")
    ifil <- cbind(ifil,inertia$abs.row[,yax])
    colnames(ifil)[7] <- paste("Contr",yax,sep="") 
    ifil <- cbind(ifil,inertia$rel.row[,yax])
    colnames(ifil)[8] <- paste("Qual",yax,sep="") 
    # calidad plano
    ifil <- cbind(ifil,inertia$rel.row[,xax]+inertia$rel.row[,yax])
    colnames(ifil)[9] <- paste("QualPlane",xax,"-",yax,sep="") 
    # peso
    ifil <- cbind(ifil,ACww$lbw*100)
    colnames(ifil)[10] <- "Weight(l)" 
    # sumas de filas
    suma <- colSums(ifil);suma[5] <- NA;suma[8:9] <- NA
    ifil <- rbind(ifil,suma)
    rownames(ifil)[(L+1)] <- "Total" 
    inertia$col <- round(icol,dec)
    inertia$row <- round(ifil,dec)
    # correlación con primer eje de subtablas
        
        # ACP bloques columna
         X <- data.frame(ACww$tab)
    nblo <- length(cbl)
    col.fac <- rep(1:nblo,cbl)
    coraxisJ <- matrix(0,nblo,2)
    for (j in 1:nblo) {
               Xj <- X[, col.fac == j]
               wj <- ACww$cw[col.fac==j]
               coor <- as.dudi(Xj,col.w=wj,row.w=ACww$lw,scannf=FALSE,nf=2,call = match.call(),type="pca")$li[,1]
             coraxisJ[j,1] <- cov.wt(cbind(ACww$li[,xax],coor), wt = ACww$lw, cor = TRUE, center = FALSE)$cor[1,2] 
              coraxisJ[j,2] <-  cov.wt(cbind(ACww$li[,yax],coor), wt = ACww$lw, cor = TRUE, center = FALSE)$cor[1,2] 
    }
            # ACP bloques fila
             X <- data.frame(t(ACww$tab))
    nblo <- length(rbl)
    col.fac <- rep(1:nblo,rbl)
    coraxisL <- matrix(0,nblo,2)
    for (j in 1:nblo) {
               Xj <- X[, col.fac == j]
               wj <- ACww$lw[col.fac==j]
               coor <- as.dudi(Xj,col.w=wj,row.w=ACww$cw,scannf=FALSE,nf=2,call = match.call(),
                        type="pca")$li[,1]
               coraxisL[j,1] <- cov.wt(cbind(ACww$co[,xax],coor), wt = ACww$cw, cor = TRUE, 
                                center = FALSE)$cor[1,2] 
               coraxisL[j,2] <- cov.wt(cbind(ACww$co[,yax],coor), wt = ACww$cw, cor = TRUE, 
                                center = FALSE)$cor[1,2] 
    }
    rownames(coraxisJ) <- names(cbl)
    rownames(coraxisL) <- names(rbl)
          colnames(coraxisJ) <- c(paste("Axis",xax,sep=""),paste("Axis",yax,sep="")) 
          colnames(coraxisL) <- colnames(coraxisJ) 
    
    inertia$coraxisJ <- coraxisJ
    inertia$coraxisL <- coraxisL
     # distancia al cuadrado al origen de filas
     X <- as.matrix(ACww$tab)
   dis.row <- (X*X) %*% ACww$cw  
   rownames(dis.row) <- rownames(X) 
   inertia$dis.row <- dis.row
 # columnas 
   dis.col <- t(ACww$lw %*% (X*X))
   rownames(dis.col) <- colnames(X)
   inertia$dis.col <- dis.col
   class(inertia) <- c("wwinertia")
 # coordenadas para los grupos =
 # suma de las contribuciones (no en porcentaje) a la inercia del eje
 # de los elemntos de la banda
 # incorporado agosto 26 de 2010
 hL <- ACww$hom[1];hJ<-ACww$hom[2]
 inertia$mfa.inL <- (hL*cbind(ifil[3],ifil[6])/ti)[1:L,]
 colnames(inertia$mfa.inL) <- colnames(ifil)[c(3,6)]     
 inertia$mfa.inJ <- (hL*cbind(icol[3],icol[6])/ti)[1:J,]
 colnames(inertia$mfa.inJ) <- colnames(icol)[c(3,6)]     
return(inertia)
}
#-------------------------------------------------------------------------------------------
# methods print 
"print.wwinertia" <- function (x, ...) {
    if (!inherits(x, "wwinertia")) 
        stop("non convenient data")

    cat("Inertia of Subclouds of ICA ")
    cat("class: ")
    cat(class(x))
    cat("\n$call: ")
    print(x$call)
       sumry <- array("", c(4, 4), list(1:4, c("vector", "length", 
        "mode", "content")))
    sumry[1, ] <- c("$dis.row", length(x$dis.row), mode(x$dis.row), "square distances of rows to the origin")
    sumry[2, ] <- c("$dis.col", length(x$dis.col), mode(x$dis.col), "square distances of colums to the origin")
    sumry[3, ] <- c("$cloud.col", length(x$cloud.col), mode(x$cloud.col), "column band weights")
    sumry[4, ] <- c("$cloud.row", length(x$cloud.row),  mode(x$cloud.row), "row band weights")

 
    class(sumry) <- "table"
    cat("\n")
    print(sumry)
      cat("\n")

    sumry <- array("", c(6, 4), list(1:6, c("data.frame", "nrow","ncol", "content")))

    sumry[1, ] <- c("$rel.col", nrow(x$rel.col),  ncol(x$rel.col), "column band qualities of the representation")
    sumry[2, ] <- c("$rel.row", nrow(x$rel.row),  ncol(x$rel.row), "row band qualities of the representation")
    sumry[3, ] <- c("$abs.col", nrow(x$abs.col),  ncol(x$abs.col), "column band absolute contributions")
    sumry[4, ] <- c("$abs.row", nrow(x$abs.row),  ncol(x$abs.row), "row band absolute contributions")
    sumry[5, ] <- c("$coraxisJ", nrow(x$coraxisJ), ncol(x$coraxisJ), "correlations between the axis separate and global analysis for column bands")
    sumry[6, ] <- c("$coraxisL", nrow(x$coraxisL), ncol(x$coraxisL), "correlations between the axis separate and global analysis for row bands")
 
    class(sumry) <- "table"
    cat("\n")
    print(sumry)
    cat("\n")


     cat("\n$col: table of inertia for the column clouds (",nrow(x$col),",",ncol(x$col),")")
     cat("\n$row: table of inertia for the row clouds (",nrow(x$row),",",ncol(x$row),")\n")
 
 
 }



