#--------------------------------------------------------
# función para adicionar indicadores para coordenadas parciales
# modificación de ktab.util.addfactor de ade4
########### ktab.util.addfactor<- ########### 
# utilitaire utilisé dans les ktab
# ajoute les composantes TL TC et T4
# x est un ktab presque achevé
# value lista que contiene rbl,cbl,nr,nc
# et le nombre de lignes
# on récupère avec le nombre de tableaux, le nombre de variables par tableaux
#############################################################
#---------------------------------------------------------
"wwm.util.addfactor<-" <- function (x,value) {
#rbl,cbl,nr,nc
    rbl <- value[[1]]
    cbl <- value[[2]]
    nr <- value[[3]]
    nc <- value[[4]]
    nblor <- length(rbl)
    nbloc <- length(cbl)
    w <- cbind.data.frame(gl(nbloc, nr), factor(rep(1:nr, 
        nbloc)))
    names(w) <- c("j", "i")
    x$ji <- w
    w <- cbind.data.frame(gl(nblor, nc), factor(rep(1:nc, 
        nblor)))
    names(w) <- c("l", "k")
    x$lk <- w
      x
}



# ----------------------------------------------------------------
# Función partial.wwm
# cálculo de coordenadas parciales de las filas y columnas
# entra ACww "wwmodel", cbl (rbl) número de columnas (fila) por banda
# sale "list": coor coordenadas parciales y ayudas a la interpretacion
# Modificaciín dilatación según peso de la banda Julio 2005
# variable logica dil para emplificar o no
# ACTUALIZACIÓN: Marzo junio 6 de 2010
# modificada noviembre 5 de 1010
# se retoma la dilatacion por el inverso del peso de la banda
# ----------------------------------------------------------------
partial.wwm <- function(ACww,dil=TRUE){ 
    # control de entrada
        #dil# if(dil) Xp <- Xp/ACww$cbw[j] #dilatación
    if (!inherits(ACww, "wwmodel")) stop("Object of class wwmodel expected")
    # parámetros iniciales
    rbl <- ACww$rbl
    cbl=ACww$cbl
    if(dil) cat("\n Coordinates were amplified by weight cloud inverses \n")###nov5
    nf <- ACww$nf
    L <- ncol(t(rbl))   # número de bloques fila
    J <- ncol(t(cbl))   # número de bloques columna
    I <- nrow(ACww$tab) # número de filas 
    K <- ncol(ACww$tab) # número de columnas
    unol <- matrix(1,L,1)
    unoj <- matrix(1,J,1)
    unonf <- matrix(1,nf,1)
    fhl <- factor(rep(1:K,L))   # factor para columnas parciales 
    fh <- factor(rep(1:I,J))    # factor para filas parciales
    cbl.fac <- rep(1:J,cbl)     # factor para bloque columna
    rbl.fac <- rep(1:L,rbl)     # factor para bloques fila
    fj <- as.factor(rep(1:J,(unoj*I)))
    fl <- as.factor(rep(1:L,(unol*K)))
    # matriz para proyección de filas parciales
    # matriz X~ Xpar
    Xpar <- NULL
    cin <- 1; cfin <- 0
    for (j in 1:J) {
        Xp <- matrix(0,I,K)
        cfin <- cfin + cbl[j]
        Xp[,cin:cfin] <- as.matrix(ACww$tab[,cin:cfin])
        if(dil) Xp <- Xp/ACww$cbw[j] ###dilatación
        rownames(Xp) <- paste(rownames(ACww$tab),j,sep="")
        colnames(Xp) <- colnames(ACww$tab)
        Xpar <- rbind(Xpar,Xp)
        cin <- cin + cbl[j]        
        }
    # fin calculo Xpar    
    # matriz para proyección de columnas parciales
    # matriz X~ Xparc
    Xparc <- NULL
    cin <- 1; cfin <- 0
    for (l in 1:L) {
        Xpc <- matrix(0,K,I)
        cfin <- cfin + rbl[l]
        Xpc[,cin:cfin] <- as.matrix(t(ACww$tab)[,cin:cfin])
        if(dil) Xpc <- Xpc/ACww$lbw[l] ###dilatación
        rownames(Xpc) <- paste(colnames(ACww$tab),l,sep="")
        colnames(Xpc) <- rownames(ACww$tab)
        Xparc <- rbind(Xparc,Xpc)
        cin <- cin + rbl[l]        
        }
    # fin calculo Xparc    
    #inercia parcial fila espacio completo
    ipar <- numeric(I*J)
    M <- diag(ACww$cw)
    for (ip in 1:(I*J)) ipar[ip] <- t(Xpar[ip,]) %*% M %*%  Xpar[ip,]
    ipar <- as.numeric((unoj %x% ACww$lw) * ipar)
    names(ipar) <- rownames(Xpar)
    #inercia subtablas l,j ============================
    inLJ <- matrix(0,L,J)
    colnames(inLJ) <- names(cbl)
    rownames(inLJ) <- names(rbl)
    for(j in 1:J) {
      ini <- 1 + (j-1)*I
      fin <- j*I
      inLJ[,j] <- tapply(ipar[ini:fin],rbl.fac,sum) 
    }
    # inercia parcial columna espacio completo
    iparc <- numeric(K*L)
    D <- diag(ACww$lw)
    for (ipc in 1:(K*L)) iparc[ipc] <- t(Xparc[ipc,]) %*% D %*%  Xparc[ipc,]
    iparc <- as.numeric((unol %x% ACww$cw) * iparc)
    names(iparc) <- rownames(Xparc)
    # contribución inercia de la nube j en el espacio completo
    X <- as.matrix(ACww$tab)
    norm <- apply(X*X*ACww$lw,2,sum)
    Ij <- tapply(norm*ACww$cw,cbl.fac,sum);names(Ij) <- names(cbl)
    # contribución inercia de la nube l en el espacio completo
    X <- t(ACww$tab)
    norm <- apply(X*X*ACww$cw,2,sum)
    Il <- tapply(norm*ACww$lw,rbl.fac,sum);names(Il) <- names(rbl)
    # coordenadas parciales fila 
    coor <- NULL
    cin <- 1; cfin <- 0
    for (j in 1:J){
        cfin <- cfin + cbl[j]
        Xj <- as.matrix(ACww$tab[,cin:cfin])
        Mj <- diag(ACww$cw[cin:cfin])
        Uj <- as.matrix(ACww$c1[cin:cfin,])
        proy <- Xj %*% Mj %*% Uj
        if(dil) proy <- proy/ACww$cbw[j] ###dilatación
       
        rownames(proy) <- paste(rownames(proy),j,sep="-")
        coor <- rbind(coor,proy)
        cin <- cin + cbl[j]
    }  
    colnames(coor) <- names(ACww$li)  
    row.coor <- coor   # sale 1
    # coordenadas parciales columna
    coorc <- NULL
    cin <- 1; cfin <- 0
    for (l in 1:L){
        cfin <- cfin + rbl[l]
        Xl <- t(ACww$tab)[,cin:cfin]
        Dl <- diag(ACww$lw[cin:cfin])
        Vl <- as.matrix(ACww$l1[cin:cfin,])
        proyc <- Xl %*% Dl %*% Vl
        if(dil) proyc <- proyc/ACww$lbw[l] ###dilatación
       
        rownames(proyc) <- paste(rownames(proyc),l,sep="-")
        coorc <- rbind(coorc,proyc)
        cin <- cin + rbl[l]
    }  
    colnames(coorc) <- names(ACww$co)  
    col.coor <- coorc   # sale 
    # inercia filas parciales en el eje s
    ipars <- coor * coor * ((unoj %x% ACww$lw) %*% t(unonf))
    row.rel <- ipars / (ipar %*% t(unonf)) * 100  # sale
    #=========================================================
        #inercia subtablas l,j sobre eje 1
    in1LJ <- matrix(0,L,J)
    colnames(in1LJ) <- names(cbl)
    rownames(in1LJ) <- names(rbl)
    for(j in 1:J) {
      ini <- 1 + (j-1)*I
      fin <- j*I
      in1LJ[,j] <- tapply(ipars[ini:fin,1],rbl.fac,sum) 
    }

    
    imed <- unoj %x% as.matrix(ACww$li) # media con dilatación
    if(!dil) imed <- imed/J ###media sin dilatación
    # inercia columnas parciales en el eje s
    iparsc <- coorc * coorc * ((unol %x% ACww$cw) %*% t(unonf))
    col.rel <- iparsc / (iparc %*% t(unonf)) * 100  # sale
    imedc <- unol %x% as.matrix(ACww$co) #media columnas parciales con dilatación
    if(!dil) imedc <- imedc/L ###media sin dilatación
    # contribucion fila parcial inercia intra eje s
    clintra <- ((coor-imed)^2)*(unoj %x% ACww$lw %x% t(unonf))
    if(dil) clintra <- clintra * ACww$cbw[fj] %x% t(unonf) ###con dilatación
     
    row.cwit <- clintra   # sale
    # contribución fila parcial inercia intra subespacio S
    row.cwitS <- apply(row.cwit,1,sum) # sale
    # contribucion columna parcial inercia intra eje s
    clintrac <- ((coorc-imedc)^2)*(unol %x% ACww$cw %x% t(unonf))
    if(dil) clintrac * ACww$lbw[fl] %x% t(unonf) ###con dilatación
    
    
    col.cwit <- clintrac   # sale
    # contribución columna parcial inercia intra subespacio S
    col.cwitS <- apply(col.cwit,1,sum) #sale
    # inercia intra punto fila eje s
    lintra <- NULL
    for (s in 1:ACww$nf) {
        lintra <- cbind(lintra,tapply(clintra[,s],fh,sum))
        }
    rownames(lintra) <- rownames(ACww$tab)
    colnames(lintra) <- colnames(row.cwit)    
    row.wit <- lintra  # sale
    # inercia intra punto fila subespacio S
    row.witS <- apply(row.wit,1,sum) #sale
    # inercia intra punto columna eje s
    lintrac <- NULL
    for (s in 1:ACww$nf) {
        lintrac <- cbind(lintrac,tapply(clintrac[,s],fhl,sum))
        }
    rownames(lintrac) <- colnames(ACww$tab)
    colnames(lintrac) <- colnames(col.cwit)    
    col.wit <- lintrac  # sale
    # inercia intra punto columna subespacio S
    col.witS <- apply(col.wit,1,sum) #sale
    # J nubes y L nubes
    # inercia total de j sobre eje s
    Ijs <- matrix(0,J,ACww$nf)
    for (nf in 1:ACww$nf) Ijs[,nf] <- tapply(ipars[,nf],fj,sum)
    if(dil) Ijs <- Ijs*ACww$cbw %*% t(unonf) ### efecto dilatación
   
    # calidad representación nube j
    quaj <- Ijs/(Ij %*% t(unonf) )*100 # cal repr nubes par  
    colnames(quaj) <- colnames(ACww$li)
    
    # inercia total de l sobre eje s
    Ils <- matrix(0,L,ACww$nf)
    for (nf in 1:ACww$nf) Ils[,nf] <- tapply(iparsc[,nf],fl,sum)
    if(dil) Ils <- Ils * ACww$lbw %*% t(unonf) ###efecto dilatación 
    
    # calidad representación nube l
    qual <- Ils/(Il %*% t(unonf) )*100 # cal repr nubes par      
    colnames(qual) <- colnames(ACww$co)
    # % entre nubes j
    betj <- ACww$eig[1:ACww$nf]/apply(Ijs,2,sum)*100###/J 
    if(!dil) betj <- betj/J # sin dilatación
    # similaridad entre nubes en el plano
    betjS <- sum(ACww$eig[1:nf])/sum(sum(Ijs[,1:nf]))*100###/J
    if(!dil) betjS <- betjS/J # sin dilatación
      names(betj) <- colnames(ACww$li)
    # % entre nubes l
    betl <- ACww$eig[1:ACww$nf]/apply(Ils,2,sum)*100###/L
    if(!dil) betl <- betl/L # sin dilatación
    names(betl) <- colnames(ACww$co)
       # similaridad entre nubes en el plano
    betlS <- sum(ACww$eig[1:nf])/sum(sum(Ils[,1:nf]))*100###/L
    if(!dil) betlS <- betlS/L # sin dilatación
    #colnames(quaj) <- colnames(ACww$li)
    partial <- NULL
    partial$dil <- dil # dilatación T/F
    partial$nf <- ACww$nf # ejes retenidos
    partial$lw <- ACww$lw  # peso filas
    partial$cw <- ACww$cw  # peso columnas 
    if (!dil) partial$row.coor <- row.coor*J   # coordenadas parciales fila
    if (dil) partial$row.coor <- row.coor
    if (!dil)  partial$col.coor <- col.coor*L # coordenadas parciales columna
    if (dil) partial$col.coor <- col.coor
    partial$row.rel <- row.rel # calidad representacion fila
    partial$col.rel <- col.rel # calidad representacion columna
    partial$row.cwit <- row.cwit # cont. inercia intra parcial fila
    partial$row.cwitS <- row.cwitS # cont. inercia intra parcial fila subespacio S
    partial$col.cwit <- col.cwit # cont. inercia intra parcial columna
    partial$col.cwitS <- col.cwitS # cont. inercia intra parcial columna subespacio S
    partial$row.wit <- row.wit    # inercia intra de filas eje s
    partial$row.witS <- row.witS    # inercia intra de filas subespacio S
    partial$col.wit <- col.wit    # inercia intra de columnas eje s
    partial$col.witS <- col.witS    # inercia intra de columnas subespacio S
    partial$quaj <- quaj   # calidad repr. nubes j
    partial$qual <- qual   # calidad repr. nubes l
    partial$betj <- betj          # % entre j/total eje s
    partial$betjS <- betjS
    partial$betl <- betl          # % entre l/total eje s
    partial$betlS <- betlS
    partial$inLJ <- inLJ          # inercia subtablas l,j
#    partial$incLJ <- incLJ
    # indicadores de filas y columnas parciales
    wwm.util.addfactor(partial) <- list(rbl,cbl,length(ACww$lw),length(ACww$cw)) 
    rownames(partial$ji) <- rownames(partial$row.coor)
    rownames(partial$lk) <- rownames(partial$col.coor)
    
# correlaciones canónicas
    nf <- ACww$nf
    #bandas-columna J
    cancorj <- matrix(NA,J*nf,nf)
    rnames <- NULL
    for (j in 1:J) {
        CP <- row.coor[partial$ji[,1]==j,]
       norm <- sqrt(diag(t(CP) %*% diag(ACww$lw) %*% CP))  
  #    cat("\n dimensiones ", dim(t(CP)),dim(ACww$lw),dim(ACww$l1),dim(diag(1/norm)), "\n")
 
      cancorj[((j-1)*nf+1):(j*nf),] <- t(CP) %*% diag(ACww$lw) %*% as.matrix(ACww$l1) %*% diag(1/norm)
     rnames <- c(rnames,paste("F",1:nf,rep(rownames(quaj)[j],nf),sep=""))
    }
    colnames(cancorj) <- colnames(row.coor)
    rownames(cancorj) <-rnames
    partial$cancorj <- cancorj
    #bandas fila L
    cancorl <- matrix(NA,L*nf,nf)
    rnames <- NULL
    for (l in 1:L) {
        CPC <- col.coor[partial$lk[,1]==l,]
       norm <- sqrt(diag(t(CPC) %*% diag(ACww$cw) %*% CPC))  
   #   cat("\n dimensiones ", dim(t(CP)),dim(ACww$lw),dim(ACww$l1),dim(diag(1/norm)), "\n")
 
      cancorl[((l-1)*nf+1):(l*nf),] <- t(CPC) %*% diag(ACww$cw) %*% as.matrix(ACww$c1) %*% diag(1/norm)
     rnames <- c(rnames,paste("F",1:nf,rep(rownames(qual)[l],nf),sep=""))
    }
    colnames(cancorl) <- colnames(col.coor)
    rownames(cancorl) <-rnames
    partial$cancorl <- cancorl
# fin correlaciones canónicas    
     
    class(partial) <- c("parwwm")
    partial$call <-  match.call()
    return(partial)   
} # fin función partial.wwm 

# methods print
"print.parwwm" <- function (x, ...) {
    if (!inherits(x, "parwwm")) 
        stop("non convenient data")

    cat("Partial coordinates on witwit.model")
    cat("class: ")
    cat(class(x))
    cat("\n$call: ")
    print(x$call)
    if(x$dil)  cat("\n Partial coordinates are dilated by number of bands \n")
        
    cat("\n$nf:", x$nf, "axis-components saved  ")
    cat(" Subespace S = R^",x$nf,"\n")  
    cat("\n$betjS:",x$betjS, "band-column global quality on S") 
    cat("\n$betlS:",x$betlS, "band-row global quality on S \n") 
       sumry <- array("", c(8, 4), list(1:8, c("vector", "length", 
        "mode", "content")))
    sumry[1, ] <- c("$cw", length(x$cw), mode(x$cw), "column weights")
    sumry[2, ] <- c("$lw", length(x$lw), mode(x$lw), "row weights")
    sumry[3, ] <- c("$betl", length(x$betl), mode(x$betl), "band-row global quality")
    sumry[4, ] <- c("$betj", length(x$betj), mode(x$betj), "band-column global quality")
    sumry[5, ] <- c("$row.witS", length(x$row.witS), mode(x$row.witS), "row within-inertia on S")
    sumry[6, ] <- c("$col.witS", length(x$col.witS), mode(x$col.witS), "column within-inertia on S")
    sumry[7, ] <- c("$row.cwitS", length(x$row.cwitS), mode(x$row.cwitS), "partial row within-inertia on S")
    sumry[8, ] <- c("$col.cwitS", length(x$col.cwitS), mode(x$col.cwitS), "partial column within-inertia on S")
  
    class(sumry) <- "table"
    cat("\n")
    print(sumry)
      cat("\n")
    sumry <- array("", c(15, 4), list(1:15, c("data.frame", "nrow", 
        "ncol", "content")))
    sumry[1, ] <- c("$row.coor", nrow(x$row.coor), ncol(x$row.coor), "partial row coordinates")
    sumry[2, ] <- c("$col.coor", nrow(x$col.coor), ncol(x$col.coor), "partial column coordinates")
    sumry[3, ] <- c("$row.rel", nrow(x$row.rel), ncol(x$row.rel), "partial row quality")
    sumry[4, ] <- c("$col.rel", nrow(x$col.rel), ncol(x$col.rel), "partial column quality")
    sumry[5, ] <- c("$row.wit", nrow(x$row.wit), ncol(x$row.wit), "row within-inertia")
    sumry[6, ] <- c("$col.wit", nrow(x$col.wit), ncol(x$col.wit), "column within-inertia")
    sumry[7, ] <- c("$row.cwit", nrow(x$row.cwit), ncol(x$row.cwit), "partial-row within-inertia")
    sumry[8, ] <- c("$col.cwit", nrow(x$col.cwit), ncol(x$col.cwit), "partial-column within-inertia")
    sumry[9, ] <- c("$qual", nrow(x$qual), ncol(x$qual), "representtion quality cloud l")
    sumry[10, ] <- c("$quaj", nrow(x$quaj), ncol(x$quaj), "representtion quality cloud j")
    sumry[11, ] <- c("$inLJ", nrow(x$inLJ), ncol(x$inLJ), "inertia blocks (l,j)")
    sumry[12, ] <- c("$ji", nrow(x$ji), ncol(x$ji), "partial row indicators")
    sumry[13, ] <- c("$lk", nrow(x$lk), ncol(x$lk), "partial column indicators")
     sumry[14, ] <- c("$cancorj", nrow(x$cancorj), ncol(x$cancorj), "band-column cannonical correlations")
     sumry[15, ] <- c("$cancorl", nrow(x$cancorj), ncol(x$cancorl), "row-column cannonical correlations")

    class(sumry) <- "table"
    print(sumry)
 }



          

          
