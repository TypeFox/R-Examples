# plot trees as multivariate graphics
#
setClass("branch",representation(LR="numeric",
w="numeric",h="numeric",El="numeric",LeafL="branch",LeafR="branch",Bole="numeric"),
prototype(LR=0,w=0,h=1,El=0,LeafL=NULL,LeafR=NULL,Bole=2))

setMethod("show", "branch",function(object){
    cat("Tree with ",length(object@El)," Elements \n")
               # Anzahl der Element
    if(object@LR==0) cat("shows to the left \n")
               # Richtung des Astes
    else cat("shows to the rigth \n")
    cat("Bole:  ",object@Bole," \n")
               # Ast oder Stamm
    cat("angle: ",object@w," \n") # Winkel
    cat("heigh: ",object@h," \n") # Höhe
    if((object@Bole==2)|(object@Bole==3))
        ("Tree is part of the Bole")
    if(class(object@LeafL)!="NULL") {
               # gleiche Anzeig des linkes Astes
        print("Left Branch")
        cat("Tree with ",length(object@LeafL@El)," Elements \n")
        if(object@LeafL@LR==0) cat("shows to the left \n")
        else cat("shows to the rigth \n")
        cat("Bole:  ",object@LeafL@Bole," \n")
        cat("angle: ",object@LeafL@w," \n")
        cat("heigh: ",object@LeafL@h," \n")
        if((object@LeafL@Bole==2)|(object@LeafL@Bole==3))
            ("Tree is part of the Bole")
    }
    if(class(object@LeafR)!="NULL") {
               # gleiche Anzeig des rechten Astes
        print("Left Branch")
        cat("Tree with ",length(object@LeafR@El)," Elements \n")
        if(object@LeafR@LR==0) cat("shows to the left \n")
        else cat("shows to the rigth \n")
        cat("Bole:  ",object@LeafR@Bole," \n")
        cat("angle: ",object@LeafR@w," \n")
        cat("heigh: ",object@LeafR@h," \n")
        if((object@LeafR@Bole==2)|(object@LeafR@Bole==3))
            ("Tree is part of the Bole")
    }    
    if(class(object@LeafL)=="NULL") ("Tree is Leaf") # einelementiger Ast?
})

setMethod("plot", "branch",function(a,x,y,len=1,lh=1,leg=NULL,...){
          # a Branch, x x-Koord, y y-Koord, len generelle Vergrößerung
          # len Zoomfaktor, lh Multiplikator der Höhe
    points=matrix(0,nrow=4,ncol=2)
          # Matrix mit den Eckpunkten des Astes
    if(a@LR==0) { # für linken Ast
        points[1,1]=x
        points[1,2]=y
        points[2,1]=points[1,1]-lh*len*a@h*sin(a@w*2*pi/360)
        points[2,2]=points[1,2]+lh*len*a@h*cos(a@w*2*pi/360)
        points[3,1]=points[2,1]+len*length(a@El)
        points[3,2]=points[2,2]
        points[4,1]=points[1,1]+len*length(a@El)
        points[4,2]=points[1,2]
        polygon(points)
    if((class(leg)!="NULL")&length(a@El)==1) {
        text(xy.coords(points[2,1],points[2,2]),colnames(leg)[a@El],cex=0.8,adj=c(0.3,-0.2))
    }
    }
    if(a@LR==1) { # für rechten Ast
        points[1,1]=x
        points[1,2]=y
        points[2,1]=points[1,1]+lh*len*a@h*sin(a@w*2*pi/360)
        points[2,2]=points[1,2]+lh*len*a@h*cos(a@w*2*pi/360)
        points[3,1]=points[2,1]+len*length(a@El)
        points[3,2]=points[2,2]
        points[4,1]=points[1,1]+len*length(a@El)
        points[4,2]=points[1,2]
        polygon(points)
    if((class(leg)!="NULL")&length(a@El)==1) {
        text(xy.coords(points[2,1],points[2,2]),colnames(leg)[a@El],cex=0.8,adj=c(0,-0.2))
    }
    }
    
         # Rekursiver Aufruf der Funktion plot
    if(class(a@LeafL)!="NULL") {
        if(a@LR==0) plot(a@LeafL,x-lh*len*a@h*sin(a@w*2*pi/360),y+lh*len*a@h*cos(a@w*2*pi/360),len,lh,leg)
        if(a@LR==1) plot(a@LeafL,x+lh*len*a@h*sin(a@w*2*pi/360),y+lh*len*a@h*cos(a@w*2*pi/360),len,lh,leg)
    }
    if(class(a@LeafR)!="NULL") {
        if(a@LR==0) plot(a@LeafR,x-lh*len*a@h*sin(a@w*2*pi/360)+len*length(a@LeafL@El),y+lh*len*a@h*cos(a@w*2*pi/360),len,lh,leg)
        if(a@LR==1) plot(a@LeafR,x+lh*len*a@h*sin(a@w*2*pi/360)+len*length(a@LeafL@El),y+lh*len*a@h*cos(a@w*2*pi/360),len,lh,leg)
    }  
})

conv<-function(x,i=1){
            #x hcl Matrix mit cutree
            #i te Stufe
    a<-new("branch") # Anlegen eines neuen Astes
    a@El<-x[length(as.matrix(x[,1])),]
            #Abspeichern der Elemente
    if(length(x[1,])==1) return(a) # Bei Astende abbruch
    repeat{ # Stufe finden wo es zur Astteilung kommt
        i=i+1
        if(prod(x[i,]==x[i,1])!=1) break 
    }
    xl<-as.matrix(x[i:(length(x)/length(x[1,])),x[i,1]==x[i,]])
    xr<-as.matrix(x[i:(length(x)/length(x[1,])),x[i,1]!=x[i,]])
            # Rekursiver Aufruf der Funktion conv
    a@LeafL<-conv(xl)
    a@LeafR<-conv(xr)
    return(a)
}

setlr<-function(x) {
    if(class(x@LeafL)=="branch") { # Überprüfung auf Klasse
        if((x@Bole==0)|(x@Bole==1)){ # Bei links von Stamm
            if(length(x@LeafL@El)<length(x@LeafR@El)) {
                y<-x@LeafL         # Wenn rechte Fortsetzung
                x@LeafL<-x@LeafR   # größer als linke
                x@LeafR<-y         # Vertauschung
            }
            x@LeafL@Bole=0  # links vom Stamm zeigt links
            x@LeafR@Bole=1  # links vom Stamm zeigt rechts
        }
        if(x@Bole==2){             # Stamm zeigt nach links
            if(length(x@LeafL@El)>length(x@LeafR@El)) {
                y<-x@LeafL         # Wenn linke Fortsetzung
                x@LeafL<-x@LeafR   # größer als rechte
                x@LeafR<-y         # Vertauschung
            }
            x@LeafL@Bole=0  # Links vom Stamm zeigt links
            x@LeafR@Bole=3  # Fortsetzung Stamm, zeigt rechts
        }
        if(x@Bole==3){             # Stamm zeigt nach rechts
            if(length(x@LeafL@El)<length(x@LeafR@El)) {
                y<-x@LeafL         # Wenn rechte Fortsetzung
                x@LeafL<-x@LeafR   # größer als linke
                x@LeafR<-y         # Vertauschung
            }
            x@LeafL@Bole=2  # Fortsetzung Stamm, zeigt links
            x@LeafR@Bole=5  # rechts vom Stamm zeigt rechts
        }
        if((x@Bole==4)|(x@Bole==5)){ # Bei rechts von Stamm
            if(length(x@LeafL@El)>length(x@LeafR@El)) {
                y<-x@LeafL         # Wenn linke Fortsetzung
                x@LeafL<-x@LeafR   # größer als rechte
                x@LeafR<-y         # Vertauschung
            }
            x@LeafL@Bole=4  # rechts vom Stamm zeigt links
            x@LeafR@Bole=5  # rechts vom Stamm zeigt rechts
        }
    x@LeafL@LR=0            # Linker branch zeigt links
    x@LeafR@LR=1            # rechter branch zeigt rechts
    # nur wenn nicht Fortsetzung Astende rekursiver Aufruf
    if(length(x@LeafL@El)>1) x@LeafL<-setlr(x@LeafL)
    if(length(x@LeafR@El)>1) x@LeafR<-setlr(x@LeafR)
    return(x)
    }
}

setw<-function(x,unm,wmax=0,wmin=0) {
         # x Branch, unm Unähnlichkeitsmatrix
         # wmax Winkel des Cluster aus allen El
         # wmin Winkel des homogensten Cluster aus 2 El 
     if((class(x@LeafL)=="branch") & (class(x@LeafR)=="branch")) {
         gX=max(unm[names(x@El),x@El[1]]) 
         # x cluster max nach complete linkage 
         # von 1 El zu allen anderen
         if(length(x@El)>1) {
              for(i in 2:length(x@El[1])) {
                  gX<-c(gX,max(unm[names(x@El),x@El[i]]))
              }
         # von den weiteren El zu den anderen
         }
         gX<-max(gX) # max Heterogenität von El in cluster
         gA<-max(unm) # max aller Merkmale
         gmin<-min(unm[unm!=0]) # Min aller Merkmale
         w<-(wmin*(log(gA+1)-log(gX+1))+wmax*(log(gX+1)-log(gmin+1)))/(log(gA+1)-log(gmin+1))
         if((x@LeafL@Bole!=2)&(x@LeafL@Bole!=3)&(x@LeafR@Bole!=2)&(x@LeafR@Bole!=3)) {
             x@LeafL@w<-w*length(x@LeafL@El)/length(x@El)
             x@LeafR@w<-w*length(x@LeafR@El)/length(x@El)
         # Wenn Ast proportionale Aufteilung nach El
         }
         if((x@LeafL@Bole==2)|(x@LeafL@Bole==3)|(x@LeafR@Bole==2)|(x@LeafR@Bole==3)) {
             x@LeafL@w<-w*length(x@LeafR@El)/length(x@El)
             x@LeafR@w<-w*length(x@LeafL@El)/length(x@El)
         # wenn stamm umgekehrt proportional
         }
         #rekursiver Aufruf
         if((class(x@LeafL@LeafL)=="branch") & (class(x@LeafL@LeafR)=="branch")) {
             x@LeafL<-setw(x@LeafL,unm,wmax,wmin)
         }
         if((class(x@LeafR@LeafL)=="branch") & (class(x@LeafR@LeafR)=="branch")) {
             x@LeafR<-setw(x@LeafR,unm,wmax,wmin)
         }
         return(x)
    }
}

seth<-function(x,a,i) {
       # x Branch, a Merkmalwerte, i Ort 
    if(class(x)=="branch") {
        x@h=sum(a[i,x@El])/length(x@El) # Mittelwert
        if((class(x@LeafL)=="branch") & (class(x@LeafR)=="branch")) {
        # bei Astfortsetzung rekursiver Aufruf
            x@LeafL<-seth(x@LeafL,a,i)
            x@LeafR<-seth(x@LeafR,a,i)
        }
    }
    return(x)
}
    
            

#tree<-function (x,locations=NULL,wmax=0,wmin=0,len=1,lh=1,legco=NULL,leglen=1, ...) 
tree <-
function(x,wmax=0,wmin=0,lh=1,
    labels = dimnames(x)[[1]], locations = NULL, nrow = NULL, ncol = NULL,
    key.loc = NULL, key.labels = dimnames(x)[[2]], key.xpd = TRUE, xlim = NULL,
    ylim = NULL, flip.labels = NULL, len=1, leglen=1, leglh=1,
    axes = FALSE, frame.plot = axes, main = NULL,
    sub = NULL, xlab = "", ylab = "", cex = 0.8, lwd = 0.25,
    lty = par("lty"), xpd = FALSE, mar = pmin(par("mar"), 1.1 +
        c(2 * axes + (xlab != ""), 2 * axes + (ylab != ""), 1,
            0)), add = FALSE, plot = TRUE, ...)

{
# draw trees as multivariate graphics
#
# x ... multivariate data in form of matrix or data frame
# wmax, wmin ... maximum and minimum angle for the leaves of the tree
# lh ...
# labels ... vector of character strings for labeling the plots
# locations ... locations for the boxes on the plot (e.g. X/Y coordinates)
# nrow, ncol ... integers giving the number of rows and columns to use when
#          'locations' is 'NULL'.  By default, 'nrow == ncol', a square
#          layout will be used.
# key.loc ... vector with x and y coordinates of the unit key.
# key.labels: vector of character strings for labeling the segments of
#         the unit key.  If omitted, the second component of
#         'dimnames(x)' is used, if available.
# key.xpd: clipping switch for the unit key (drawing and labeling), see
#          'par("xpd")'.
# xlim: vector with the range of x coordinates to plot.
# ylim: vector with the range of y coordinates to plot.
# flip.labels: logical indicating if the label locations should flip up
#          and down from diagram to diagram. Defaults to a somewhat
#          smart heuristic.
# len, leglen, leglh: multiplicative values for the space of the labels on the legend
# axes: logical flag: if 'TRUE' axes are added to the plot.
# frame.plot: logical flag: if 'TRUE', the plot region is framed.
#   main: a main title for the plot.
#    sub: a sub title for the plot.
#   xlab: a label for the x axis.
#   ylab: a label for the y axis.
#    cex: character expansion factor for the labels.
#    lwd: line width used for drawing.
#    lty: line type used for drawing.
#    xpd: logical or NA indicating if clipping should be done, see
#         'par(xpd = .)'.
#    mar: argument to 'par(mar = *)', typically choosing smaller
#         margings than by default.
#    add: logical, if 'TRUE' _add_ boxes to current plot.
#   plot: logical, if 'FALSE', nothing is plotted.
#    ...: further arguments, passed to the first call of 'plot()', see
#         'plot.default' and to 'box()' if 'frame.plot' is true.


    if (is.data.frame(x))
        x <- data.matrix(x)
    else if (!is.matrix(x))
        stop("'x' must be a matrix or a data frame")
    if (!is.numeric(x))
        stop("data in 'x' must be numeric")
    n.loc <- nrow(x)
    n.seg <- ncol(x)
    if (is.null(locations)) {
        if (is.null(nrow))
            nrow <- ceiling(if (!is.numeric(ncol)) sqrt(n.loc) else n.loc/ncol)
        if (is.null(ncol))
            ncol <- ceiling(n.loc/nrow)
        if (nrow * ncol < n.loc)
            stop("nrow * ncol <  number of observations")
        ff <- if (!is.null(labels)) 2.3
        else 2.1
        locations <- expand.grid(ff * 1:ncol, ff * nrow:1)[1:n.loc, ]
        if (!is.null(labels) && (missing(flip.labels) || !is.logical(flip.labels)))
            flip.labels <- ncol * mean(nchar(labels, type = "c")) > 30
    }
    else {
        if (is.numeric(locations) && length(locations) == 2) {
            locations <- cbind(rep.int(locations[1], n.loc),
                rep.int(locations[2], n.loc))
            if (!missing(labels) && n.loc > 1)
                warning("labels do not make sense for a single location")
            else labels <- NULL
        }
        else {
            if (is.data.frame(locations))
                locations <- data.matrix(locations)
            if (!is.matrix(locations) || ncol(locations) != 2)
                stop("'locations' must be a 2-column matrix.")
            if (n.loc != nrow(locations))
                stop("number of rows of 'locations' and 'x' must be equal.")
        }
        if (missing(flip.labels) || !is.logical(flip.labels))
            flip.labels <- FALSE
    }
    xloc <- locations[, 1] # Auslesen der Koordinaten
    yloc <- locations[, 2]

    x[is.na(x)] <- 0
    mx <- max(x <- x * len)
    if (is.null(xlim))
        xlim <- range(xloc) + c(-mx, mx)
    if (is.null(ylim))
        ylim <- range(yloc) + c(-mx, mx)
    op <- par(mar = mar, xpd = xpd)
    on.exit(par(op))
    if (!add)
        plot(0, type = "n", ..., xlim = xlim, ylim = ylim, main = main,
            sub = sub, xlab = xlab, ylab = ylab, asp = 1, axes = axes)
    if (!plot)
        return()

    corm = cor(x)          # Distanzmatrix
    unm  = 1-abs(corm)
    unm[upper.tri(unm)]<-0
    hcl  = hclust(dist(unm)) # Cluster erstellen
    bcl=matrix(0,nrow=ncol(x),ncol=ncol(x))
                           # einzelne Stufen durchlaufen
    for(i in 1:length(x[1,])) bcl[i,]=cutree(hcl,k=i)
    colnames(bcl)<-c(names(cutree(hcl,k=length(x[1,]))))
    a<-conv(bcl)           # in Objekt der Klasse "branch" umwandeln
    a<-setlr(a)            # Verzweigungsverlauf regeln
    a<-setw(a,unm,wmax=wmax,wmin=wmin) # Winkelsetzung
    dmin=apply(x,2,min)
    dmax=apply(x,2,max)
    q=(x-rep(1,nrow(x))%*%t(dmin))/(rep(1,nrow(x))%*%t(dmax-dmin))
                           # Merkmalswerte normieren
    for(i in 1:length(locations[,1])) {
                           # Ablaufen aller zu plotenden Punkte
        a<-seth(a,q,i)
        plot(a,xloc[i]-len*lh,yloc[i],len,lh)
    }

   if (!is.null(labels)) {
        y.off <- mx
        if (flip.labels)
            y.off <- y.off + cex * par("cxy")[2] * ((1:n.loc)%%2 - 0.4)
        text(xloc, yloc - y.off/2, labels, cex = cex, adj = c(0.5, 1))
    }

    if (!is.null(key.loc)) {
        par(xpd = key.xpd)

        xleg=key.loc[1]
        yleg=key.loc[2]

legh=matrix(seq(from=length(x[1,]),to=1,length=length(x[1,])),nrow=1,ncol=length(x[1,]))
        a<-seth(a,legh,1)
        plot(a,xleg,yleg,leglen,leglh,bcl)
    }
    if (frame.plot)
        box(...)
    invisible(locations)

}



