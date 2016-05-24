#---------------------------------------------------------------------------------------------
# gráfica de un plano factorial
# Campo Elías Pardo
# Octubre 11-06, modificada enero 17/07
# requiere ade4
#PARAMETROS (entre paréntesis valores por defecto)
#   dudi es un objeto dudi (coa,acm,fca)
#   x,y ejes a graficar (1,2)
#   roweti: filas a etiquetar (todas)
#   coleti: columnas a etiquetar (todas)
#   main: título de la gráfica (NULL)
#   axislabel:
#   col.row: color para las filas (black)
#     col.col: color para las colulmnas (blue)
#     cex (0.8)
#   cex.row: escala para etiquetas de filas (0.8)
#   cex.col: escala para etiquetas de columnas (0.8)
#     all.point: cierto para graficar todos los puntos aunque no estén etiquetados (TRUE)
#     Trow: cierto para graficar filas (TRUE)
#   Tcol: cierto para graficar columnas (TRUE)
#   cframe: aumento de los límites de la gráfica (1.2)
#   ucal: umbral (%) de calidad de representación (0), se etiquetan puntos por encima
#     del umbral en el plano
#   cex.global: factor de escala para todas las etiquetas
#     infaxes: lugar para imprimir información de ejes: "out","in","no" ("out")
#---------------------------------------------------------------------------------------------
plot.dudi <- function(x,ex=1,ey=2,xlim=NULL,ylim=NULL,main=NULL,rotx=FALSE,roty=FALSE,roweti=row.names(dudi$li),
                        coleti=row.names(dudi$co),axislabel=TRUE,
                        col.row="black",col.col="blue",cex=0.8,cex.row=0.8,cex.col=0.8,
                        all.point=TRUE,Trow=TRUE,Tcol=TRUE,cframe=1.2,ucal=0,
                cex.global=1,infaxes="out",...)
{   
   dudi <- x
   if (!inherits(dudi, "dudi"))  stop("Object of class  'dudi' expectes")
   x <- ex
   y <- ey
# rotación de ejes
if (rotx) rotx=-1 else rotx=1
if (roty) roty=-1 else roty=1
# selección de puntos por umbral de calidad de representación en el plano
if (ucal>0){
    cosfil <- inertia.dudi(dudi,TRUE)$row.rel
    coscol <- inertia.dudi(dudi,,TRUE)$col.rel
      roweti <- row.names(subset(dudi$li,(abs(cosfil[,x])+abs(cosfil[,y]))>ucal*100))
      coleti <- row.names(subset(dudi$co,(abs(coscol[,x])+abs(coscol[,y]))>ucal*100))   
    }
       
    eigx <- dudi$eig[x]
    peigx <- round(eigx/sum(dudi$eig)*100,1)
    eigx <- round(eigx,4)
    eigy <- dudi$eig[y]
    peigy <- round(eigy/sum(dudi$eig)*100,1)
    eigy <- round(eigy,4)                    
    if (is.null(xlim)) xlim <- c(min(min(rotx*dudi$li[,x],rotx*min(dudi$co[,x]))),
                max(max(rotx*dudi$li[,x],max(rotx*dudi$co[,x]))))
    if (is.null(ylim)) ylim <- c(min(min(roty*dudi$li[,y],min(roty*dudi$co[,y]))),
                max(max(roty*dudi$li[,y],max(roty*dudi$co[,y])))) 
    xlim <- xlim*cframe
    ylim <- ylim*cframe      
    cex <- cex*cex.global
    cex.lab <- 0.8*cex.global
    cex.axis <- 0.8*cex.global
    cex.main <- 0.8*cex.global 
    cex.row <- cex.row*cex.global
    cex.col <- cex.col*cex.global

    # estilo ade4
    if (infaxes != "out"){
        opar <- par(mar = par("mar")) # tomado de s.label de ade4
            on.exit(par(opar))      # quita los márgenes
            par(mar = c(0.1, 0.1, 0.1, 0.1)) # externos

        plot.default(0, 0, type = "n", asp = 1, xlab = "", ylab = "", 
        xaxt = "n", yaxt = "n", xlim = xlim, ylim = ylim, xaxs = "i", 
        yaxs = "i", frame.plot = TRUE)
          sutil.grid(cex)
 
          scatterutil.sub(main, cex)
        if (infaxes=="in"){
        text(xlim[2],ylim[1],adj=c(1,0),paste("Factor ",x,": ",eigx," (",peigx,"%)",sep=""),cex=cex) 
            text(xlim[1],ylim[2],adj=c(0,1),paste("Factor ",y,": ",eigy," (",peigy,"%)",sep=""),cex=cex)
        }
    }
    # estilo normal 
    if (infaxes=="out"){ 
    plot(0, 0, main = main, xlab = paste("Factor ",x,": ",eigx," (",peigx,"%)",sep=""), 
               ylab = paste("Factor ",y,": ",eigy," (",peigy,"%)",sep=""), 
               xlim = xlim, ylim = ylim, col = "white", asp=1, cex=cex,
               cex.lab=cex.lab,cex.axis=cex.axis,cex.main=cex.main,las=1)

          sutil.grid(cex,FALSE)

    }
    abline(h = 0, v = 0, lty = 2)#,col="darkgrey")
    if(all.point){                                                                      
        if(Trow) points(cbind(rotx*dudi$li[,x],roty*dudi$li[,y]), 
                        pch = 20, col = col.row, cex = cex.row)
        if (Tcol) points(cbind(rotx*dudi$co[,x],roty*dudi$co[,y]), 
                        pch = 17, col = col.col, cex = cex.col)
    } else {
        if(Trow) points(rotx*dudi$li[roweti,x],roty*dudi$li[roweti,y], 
                        pch = 20, col = col.row, cex = cex.row)
        if(Tcol) points(rotx*dudi$co[coleti,x],roty*dudi$co[coleti,y], 
                        pch = 17, col = col.col, cex = cex.col)
    }
if(Trow) {
    exy <- subset(dudi$li[roweti,],select=c(x,y)) 
    exy[,1] <- rotx*exy[,1] 
    exy[,2] <- roty*exy[,2]
    exyB <- subset(exy,abs(exy[,2])>abs(exy[,1]) & exy[,2] < 0) 
    if (nrow(exyB)>0) 
        text(x=exyB[,1],y=exyB[,2],
                labels=rownames(exyB),col=col.row,pos=1,cex=cex.row)
    exyL <- subset(exy,abs(exy[,2])<abs(exy[,1]) & exy[,1] < 0) 
    if (nrow(exyL)>0) 
        text(x=exyL[,1],y=exyL[,2],
                labels=rownames(exyL),col=col.row,pos=2,cex=cex.row)
    exyA <- subset(exy,abs(exy[,2])>abs(exy[,1]) & exy[,2] > 0) 
    if (nrow(exyA)>0) 
        text(x=exyA[,1],y=exyA[,2],
                labels=rownames(exyA),col=col.row,pos=3,cex=cex.row)
    exyR <- subset(exy,abs(exy[,2])<abs(exy[,1]) & exy[,1] > 0)
    if (nrow(exyR)>0) 
        text(x=exyR[,1],y=exyR[,2],
                labels=rownames(exyR),col=col.row,pos=4,cex=cex.row)
}
    # columnas
if(Tcol) {
    fxy <- subset(dudi$co[coleti,],select=c(x,y))
    fxy[,1] <- rotx*fxy[,1]
    fxy[,2] <- roty*fxy[,2]
    fxyB <- subset(fxy,abs(fxy[,2])>=abs(fxy[,1]) & fxy[,2] <= 0)
    if (nrow(fxyB)>0) 
        text(x=fxyB[,1],y=fxyB[,2],
                labels=rownames(fxyB),col=col.col,pos=1,cex=cex.col)
    fxyL <- subset(fxy,abs(fxy[,2])<=abs(fxy[,1]) & fxy[,1] <= 0)
    if (nrow(fxyL)>0) 
        text(x=fxyL[,1],y=fxyL[,2],
                labels=rownames(fxyL),col=col.col,pos=2,cex=cex.col)
    fxyA <- subset(fxy,abs(fxy[,2])>=abs(fxy[,1]) & fxy[,2] >= 0)
    if (nrow(fxyA)>0) 
        text(x=fxyA[,1],y=fxyA[,2],
                labels=rownames(fxyA),col=col.col,pos=3,cex=cex.col)
    fxyR <- subset(fxy,abs(fxy[,2])<=abs(fxy[,1]) & fxy[,1] >= 0)
    if (nrow(fxyR)>0) 
        text(x=fxyR[,1],y=fxyR[,2],
                labels=rownames(fxyR),col=col.col,pos=4,cex=cex.col)
}
 }
#------------------fin de plot.dudi antiguo planfac--------------------
# grilla tomada de ade4
"sutil.grid" <- function (cgrid,scale=TRUE) {
    col <- "lightgray"
    lty <- 1
    xaxp <- par("xaxp")
    ax <- (xaxp[2] - xaxp[1])/xaxp[3]
    yaxp <- par("yaxp")
    ay <- (yaxp[2] - yaxp[1])/yaxp[3]
    a <- min(ax, ay)
    v0 <- seq(xaxp[1], xaxp[2], by = a)
    h0 <- seq(yaxp[1], yaxp[2], by = a)
    abline(v = v0, col = col, lty = lty)
    abline(h = h0, col = col, lty = lty)
    if (cgrid <= 0) 
        return(invisible())
    cha <- paste(" d = ", a, " ", sep = "")
    cex0 <- par("cex") * cgrid
    xh <- strwidth(cha, cex = cex0)
    yh <- strheight(cha, cex = cex0) * 5/3
    x1 <- par("usr")[2]
    y1 <- par("usr")[4]
#    rect(x1 - xh, y1 - yh, x1 + xh, y1 + yh, col = "white", border = 0)
    if (scale) text(x1 - xh/2, y1 - yh/2, cha, cex = cex0)
}

