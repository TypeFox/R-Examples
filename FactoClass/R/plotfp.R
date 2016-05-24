#---------------------------------------------------------------------------------------------
# gráfica de un plano factorial
# Campo Elias Pardo
# A partir de plot.dudi
# para graficar plano factorial a partir de las coordenadas
# informaci adidional ingresada directamente
#   x,y ejes a graficar (1,2)
#   eti: puntos a etiquetar (todas)
#   main: título de la gráfica (NULL)
#   axislabel:
#   col.row: color para los puntos (black)
#     cex (0.8)
#   cex.row: escala para etiquetas de los puntos (0.8)
#     all.point: cierto para graficar todos los puntos aunque no estén etiquetados (TRUE)
#   cframe: aumento de los límites de la gráfica (1.2)
#   cal: calidad de representación de los puntos
#   ucal: umbral (%) de calidad de representación (0), se etiquetan puntos por encima
#     del umbral en el plano
#   cex.global: factor de escala para todas las etiquetas
#     infaxes: lugar para imprimir información de ejes: "out","in","no" ("out")
#---------------------------------------------------------------------------------------------
plotfp <- function(co,x=1,y=2,eig=NULL,cal=NULL,ucal=0,xlim=NULL,ylim=NULL,main=NULL,
                   rotx=FALSE,roty=FALSE,eti=row.names(co),
                        axislabel=TRUE,col.row="black",cex=0.8,cex.row=0.8,
                        all.point=TRUE,cframe=1.2,cex.global=1,infaxes="out",asp=1)
{
if (!is.null(eig))
  {
    eigx <- eig[x]
    peigx <- round(eigx/sum(eig)*100,1)
    eigx <- round(eigx,4)
    eigy <- eig[y]
    peigy <-round(eigy/sum(eig)*100,1)
    eigy <- round(eigy,4)                    
  } 
# rotación de ejes
if (rotx) rotx=-1 else rotx=1
if (roty) roty=-1 else roty=1
# selección de puntos por umbral de calidad de representación en el plano
if (ucal>0) eti <- row.names(subset(co,(abs(cal[,x])+abs(cal[,y]))>ucal*100))

 
    if (is.null(xlim)) xlim <- c(min(c(rotx*co[,x],0)),max(rotx*co[,x]))
    if (is.null(ylim)) ylim <- c(min(c(rotx*co[,y],0)),max(rotx*co[,y]))
    xlim <- xlim*cframe
    ylim <- ylim*cframe      
    cex <- cex*cex.global
    cex.lab <- 0.8*cex.global
    cex.axis <- 0.8*cex.global
    cex.main <- 0.8*cex.global 
    cex.row <- cex.row*cex.global
    xlabel <- paste("Factor ",x,": ",sep="")
    if (!is.null(eig)) xlabel <- paste(xlabel,eigx," (",peigx,"%)",sep="")
    ylabel <- paste("Factor ",y,": ",sep="")
    if (!is.null(eig)) ylabel <- paste(ylabel,eigy," (",peigy,"%)",sep="")
 
    # estilo ade4
    if (infaxes != "out")
      {
        opar <- par(mar = par("mar")) # tomado de s.label de ade4
            on.exit(par(opar))      # quita los márgenes
            par(mar = c(0.1, 0.1, 0.1, 0.1)) # externos

        plot.default(0, 0, type = "n", asp = asp, xlab = "", ylab = "", 
        xaxt = "n", yaxt = "n", xlim = xlim, ylim = ylim, xaxs = "i", 
        yaxs = "i", frame.plot = TRUE)
          sutil.grid(cex)
 
          scatterutil.sub(main, cex)
        if (infaxes=="in")
	  {
	    text(xlim[2],ylim[1],adj=c(1,0),xlabel,cex=cex) 
            text(xlim[1],ylim[2],adj=c(0,1),ylabel,cex=cex)
	  }
	}
    # estilo normal 
    if (infaxes=="out")
      { 
	plot(0, 0, main = main, xlab = xlabel,ylab = ylabel, 
               xlim = xlim, ylim = ylim, col = "white", asp=asp, cex=cex,
               cex.lab=cex.lab,cex.axis=cex.axis,cex.main=cex.main,las=1)

          sutil.grid(cex,FALSE)

      }
    abline(h = 0, v = 0, lty = 2)#,col="darkgrey")
    if(all.point)
      {                                                                      
        points(cbind(rotx*co[,x],roty*co[,y]), 
                        pch = 20, col = col.row, cex = cex.row)
      
      } else 
      {
        points(rotx*co[eti,x],roty*co[eti,y], 
                        pch = 20, col = col.row, cex = cex.row)
      }

    exy <- subset(co[eti,],select=c(x,y)) 
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
#------------------fin de plotfp---------------------------------------------------------------
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

