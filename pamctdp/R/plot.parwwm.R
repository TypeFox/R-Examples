# ---------------------------------------------------------------------------------------
# función para graficar iterativamente planos factoriales de un doble análisis intra
# algunas idedeas tomadas de plot.GPApartial de FactoMineR
# x objeto de clase parwwm
# parámetros similares a plot.dudi (planfac), se tomamparte del código
# programdo por Campo Elías Pardo
# julio 28 de 2009
# ultima correccion julio 10/2010, fijar limites de la grafica cuando se
# llama, para evitar derroche de espacio
#---------------------------------------------------------------------------------------
plot.parwwm <- function(x,xy=c(1,2),graph="rows",namesg=NULL,
                        xlim=NULL,ylim=NULL,main=NULL,
                        rotx=FALSE,roty=FALSE,roweti=row.names(dudi$li),
                        coleti=row.names(dudi$co),axislabel=TRUE,asp=1,grid=TRUE,
                        col.row="black",col.col="black",cex=0.8,cex.row=0.8,cex.col=0.8,
                        cframe=1.2,cex.global=1,
                        col.own= c("darkred","darkgreen" ,"darkblue", "darkmagenta",
                                   "red","darkorange","green3",palette()),...)
{
    parwwm <- x
    if (!inherits(parwwm, "parwwm")) stop("non convenient data")
    x <- xy[1]
    y <- xy[2]
    if (rotx) rotx <- -1 else rotx <- 1
    if (roty) roty <- -1 else roty <- 1
    cex <- cex*cex.global
    cex.lab <- 0.8*cex.global
    cex.axis <- 0.8*cex.global
    cex.main <- 0.8*cex.global 
    cex.row <- cex.row*cex.global
    cex.col <- cex.col*cex.global
    # para grafica de filas
    if (graph=="rows") {
      Trow <- TRUE
      Tcol <- FALSE
      colp <- col.row
      cexp <- cex.row
    }
    # para grafica de columnas
    else {
      Trow <-FALSE
      Tcol <- TRUE
      colp <- col.col
      cexp <- cex.col
    }
    xlimi <- xlim
    ylimi <- ylim
    dudi <- eval(as.list(parwwm$call)$ACww) # dudi
    # valores propios y porcentajes para ejes
    eigx <- dudi$eig[x]
    peigx <- round(eigx/sum(dudi$eig)*100,1)
    eigx <- round(eigx,4)
    eigy <- dudi$eig[y]
    peigy <- round(eigy/sum(dudi$eig)*100,1)
    eigy <- round(eigy,4)                    

    if (Trow) {
      tab <- dudi$li[,xy]
      tab[,1]<- tab[,1]*rotx; tab[,2]<- tab[,2]*roty
      nrt <- nrow(tab)
      # coordenadas parciales con puntos homologos seguidos
      copar <- parwwm$row.coor[order(parwwm$ji[,2]),xy]
      copar[,1]<-copar[,1]*rotx; copar[,2]<-copar[,2]*roty
      ng <- ncol(parwwm$inLJ)
      if(is.null(namesg)) namesg <- colnames(parwwm$inLJ)
      eti <- roweti
     }
    if (Tcol) {
      tab <- dudi$co[,xy]
      tab[,1]<- tab[,1]*rotx; tab[,2]<- tab[,2]*roty
      nrt <- nrow(tab)
      # coordenadas parciales con puntos homologos seguidos
      copar <- parwwm$col.coor[order(parwwm$lk[,2]),xy]
      copar[,1]<-copar[,1]*rotx; copar[,2]<-copar[,2]*roty
      ng <- nrow(parwwm$inLJ)
       if(is.null(namesg)) namesg <- rownames(parwwm$inLJ)
      eti <- coleti
     }
    # subtablas para puntos seleccionados 
    if (nrow(tab)!=length(eti))
    {
      sico <- data.frame(ide=rownames(tab),si=FALSE)
      sico[,1]<-as.character(sico[,1])
      for (neti in 1:length(eti))
        {
        for (rowi in 1:nrow(tab))
           {
           if(sico[rowi,1]==eti[neti]) sico[rowi,2] <- TRUE
           }
       } 
      fhsi <- factor(rep(sico[,2],each=ng))
      copar  <- subset(copar,fhsi==TRUE)
      tab <- tab[eti,]
      all.point<- FALSE
    }# fin if   
   
    # para identificacion de puntos parciales a graficar
    fh <- factor(rep(1:nrow(tab),each=ng))
    draw.partial <- rep(FALSE,nrow(tab)) 
    disto <- matrix(0,nrow(tab),1)
    inicia <- TRUE
    fin <- FALSE
    pos <- list(x=0,y=0) 
    point.haut <- 1
    fh1 <-NULL
    # dev.new()
    # ciclo hasta clic arriba
    while (!fin){
      if(!inicia) pos <- locator(n=1) # lee coordenadas de clic en la grafica 
        if (pos$y> point.haut) fin <- TRUE
      else{
        if(is.null(xlimi)) xlim <- c(min(tab[,1]),max(tab[,1]))
        else xlim <- xlimi
        if(is.null(ylimi)) ylim <- c(min(tab[,2]),max(tab[,2]))
        else ylim <- ylimi
        # inicia salto primera vez
        if(!inicia)
        {
          for (i in 1:nrow(tab)) disto[i] <- (tab[i,1]-pos$x)^2+(tab[i,2]-pos$y)^2
          draw.partial[order(disto)[1]] <- !draw.partial[order(disto)[1]]
          partial <- rownames(tab)[draw.partial]
    
          if (length(partial)==0) partial <- NULL
          # para identificar filas homólogas a graficar  
          if (!is.null(draw.partial))
             fil <-  order(draw.partial,decreasing=TRUE)[1:sum(draw.partial)]
          # coordenadas parciales
          fh1 <-NULL
          for (nel in 1:length(fil)) fh1 <- rbind(fh1,copar[fh==fil[nel],])
          # cambiar límites de la grafica si no esta fija en el llamado
          if(is.null(xlimi)) if (min(fh1[,1]) < xlim[1]) xlim[1] <- min(fh1[,1])
          if(is.null(xlimi)) if (max(fh1[,1]) > xlim[2]) xlim[2] <- max(fh1[,1])
          if(is.null(ylimi)) if (min(fh1[,2]) < ylim[1]) ylim[1] <- min(fh1[,2])
          if(is.null(ylimi)) if (ylim[2] < max(fh1[,2])) ylim[2] <- max(fh1[,2])
          if(is.null(xlimi)) if(xlim[1]<0) xlim[1] <- round(xlim[1]-0.05,1) else xlim[1] <-  round(xlim[1]+0.05,1)
          if(is.null(xlimi)) if(xlim[2]<0) xlim[2] <- round(xlim[2]-0.05,1) else xlim[2] <-  round(xlim[2]+0.05,1)
          if(is.null(ylimi)) if(ylim[1]<0) ylim[1] <- round(ylim[1]-0.05,1) else ylim[1] <-  round(ylim[1]+0.05,1)
          if(is.null(ylimi)) if(ylim[2]<0) ylim[2] <- round(ylim[2]-0.05,1) else ylim[2] <-  round(ylim[2]+0.05,1)
         
        } # fin salto primera vez
        #
        #dev.off()
        plot.default(0, 0,xlim=xlim, ylim=ylim,main = main,
                   xlab = paste("Factor ",x,": ",eigx," (",peigx,"%)",sep=""), 
                   ylab = paste("Factor ",y,": ",eigy," (",peigy,"%)",sep=""), 
                   col = "white",cex=cex,asp=asp,
                   cex.lab=cex.lab,cex.axis=cex.axis,cex.main=cex.main,las=1)
        if (grid) sutil.grid(cex,FALSE)
        abline(h = 0, v = 0, lty = 2)#,col="darkgrey") # ejes por el centro
        points(cbind(tab[,1],tab[,2]),pch = 20, col = colp, cex = cexp)
        exy <- tab
    #    exy[,1] <- rotx*exy[,1] 
    #    exy[,2] <- roty*exy[,2]
        exyB <- subset(exy,abs(exy[,2])>abs(exy[,1]) & exy[,2] < 0) 
        if (nrow(exyB)>0) 
               text(x=exyB[,1],y=exyB[,2],
                labels=rownames(exyB),col=colp,pos=1,cex=cexp)
        exyL <- subset(exy,abs(exy[,2])<abs(exy[,1]) & exy[,1] < 0) 
        if (nrow(exyL)>0) 
               text(x=exyL[,1],y=exyL[,2],
                labels=rownames(exyL),col=colp,pos=2,cex=cexp)
        exyA <- subset(exy,abs(exy[,2])>abs(exy[,1]) & exy[,2] > 0) 
        if (nrow(exyA)>0) 
               text(x=exyA[,1],y=exyA[,2],
                labels=rownames(exyA),col=colp,pos=3,cex=cexp)
        exyR <- subset(exy,abs(exy[,2])<abs(exy[,1]) & exy[,1] > 0)
        if (nrow(exyR)>0) 
               text(x=exyR[,1],y=exyR[,2],
                labels=rownames(exyR),col=colp,pos=4,cex=cexp)
       
       
        point.haut <- ylim[2]
        if(!inicia) {
          segments(rep(tab[fil,1],each=ng),rep(tab[fil,2],each=ng),
                 fh1[,1],fh1[,2],col=col.own[1:ng],lty=c(1:ng))
           text(fh1,namesg,col=col.own[1:ng],cex=cexp*0.8,pos=1)
        }
        else  inicia <- FALSE
      } # fin else
    } # fin while
    #return(list(global=tab,partial=fh1,xlim=xlim,ylim=ylim,draw=draw.partial))
} # fin plot.parwwm   
#-------------------------------------------------    
