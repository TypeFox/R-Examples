bar3d.ade <-
function(x, y=NULL, data=NULL, xw=0.5, zw=1, main=NULL, xlab=NULL, ylab=NULL, zlab=NULL, xticks=NULL, yticks=NULL, zticks=NULL, col=NULL, tcol=NULL, bgcol=NULL, lcol=NULL, alpha=NULL, axes=TRUE, fgbox=TRUE, bgbox=TRUE, wall=0){
if(any(par('mfg')!=c(1,1,1,1)) & any(par('mai') < c(1.02, 0.82, 0.82, 0.42))){
maidiff<-rep(0, 4)
norm<-c(1.02, 0.82, 0.82, 0.42)
maidiff[par('mai')<norm]<-  norm[par('mai')<norm] - par('mai')[par('mai')<norm]
par(mai=par('mai')+maidiff)
}
oldpar<-par(no.readonly =TRUE)
oldpar<-oldpar[-which(names(oldpar)%in%c('usr',  'plt', 'pin',  'fin', 'fig', 'mfg', 'mfcol', 'mfrow', 'omd', 'omi', 'oma'))]
on.exit(par(oldpar))
if(wall!=5){
newmai<-rep(0, 4)
oldmai<-par('mai')
#if(oldmai[1]<1.2) newmai[1]<- 1.2 - oldmai[1]
if(oldmai[2]<1) newmai[2]<- 1 - oldmai[2]
if(oldmai[3]>0.75 & oldmai[3]<=0.82) newmai[3]<- 0.75-oldmai[3]
if(oldmai[4]<1) newmai[4]<- 1 - oldmai[4]
par(mai=(oldmai+newmai))


}
if(wall==5){
newmai<-rep(0, 4)
oldmai<-par('mai')
#if(oldmai[1]<1.2) newmai[1]<- 1.2 - oldmai[1]
if(oldmai[2]<1) newmai[2]<- 1 - oldmai[2]
if(oldmai[3]>0.55 & oldmai[3]<=0.82) newmai[3]<- 0.55-oldmai[3]
if(oldmai[4]<0.85) newmai[4]<- 0.85 - oldmai[4]
par(mai=(oldmai+newmai))
}

if(xw>1) xw<-1


####################
##############################
if(!is.character(x)){
xt<-deparse(substitute(x))
if(regexpr('~', xt)>=0){
x<-gsub('[~].*$', '', xt)
y<-gsub('^.*[~]', '', xt)
}
}
##############################


#Als Matrix oder Tabelle
if(is.table(x) | is.matrix(x))    tab<-x

# Mit Data Vectoren
if((is.numeric(x) | is.factor(x)) & (is.numeric(y) | is.factor(y))  & !is.matrix(x) & !is.table(x)) {
tab<-table(y, x)
if(is.null(xlab))  xlab<- gsub('[(]{0}[A-Za-z0-9]*[$]', '' ,deparse(substitute(x)))
if(is.null(zlab))  zlab<- gsub('[(]{0}[A-Za-z0-9]*[$]', '' ,deparse(substitute(y)))
}

# Mit Strings und Data-frame
if(is.character(x) & is.character(y) & !is.null(data) & length(x)==1){
if(is.null(xlab)) xlab<-x
if(is.null(zlab)) zlab<-y
x<-eval(parse(text=paste("data$",x)))
y<-eval(parse(text=paste("data$",y)))
tab<-table(y, x)
}
if(is.matrix(tab))  tab <- as.table(tab)
####################

################################################################################
##  add 3d box to a plot                                                      ##
a.bar3d.ade<-function(x0, x1, y0, y1, z0, z1, zdist=1, hx=NULL, hy=NULL, col = NA, border = NULL, lcol=NULL, density=NULL, angle=0, lty = 1, lwd = 1, lty2 = 1, lwd2 = 1, backside=TRUE){

################
if(is.null(lcol)) lcol<-a.coladd.ade(col, -100)
col1<- col
col2<- a.coladd.ade(col, -50)
col3<- a.coladd.ade(col, -75)
################
xrange<- diff(range(par('usr')[1:2]))
yrange<- diff(range(par('usr')[3:4]))
if(is.null(hx)) hx<-(xrange/40)*(zdist)
if(is.null(hy)) hy<-(yrange/40)*(zdist)
if(!is.null(hx)) hx<-hx*(zdist)
if(!is.null(hy)) hy<-hy*(zdist)
rx0<-hx*z0
rx1<-hx*(z1-z0)
ry0<-hy*z0
ry1<-hy*(z1-z0)

if(backside){
x_01<-c(x0,   x0,   x1,    x1)+rx1+rx0
y_01<-c(y0,   y1,   y1,    y0)+ry1+ry0
polygon(x_01,y_01,col=col,lty=lty,lwd=lwd,border=border)

x_05<-c(x0,   x0+rx1,   x1+rx1,    x1)+rx0
y_05<-c(y0,   y0+ry1,   y0+ry1,    y0)+ry0
polygon(x_05,y_05, col=col,lty=lty,lwd=lwd,border=border)

x_06<-c(x0,   x0,   x0+rx1,    x0+rx1)+rx0
y_06<-c(y0,   y1,   y1+ry1,    y0+ry1)+ry0
polygon(x_06,y_06,col=col,lty=lty,lwd=lwd,border=border)
}

x_03<-c(x0,   x0+rx1,   x1+rx1,    x1)+rx0
y_03<-c(y1,   y1+ry1,   y1+ry1,    y1)+ry0
polygon(x_03,y_03,col=col2,lty=lty,lwd=lwd,border=border)
if(!is.null(density)) polygon(x_03,y_03,col=lcol,lty=lty,lwd=lwd,border=border, density=density, angle=angle+45)

x_04<-c(x1,   x1,   x1+rx1,    x1+rx1)+rx0
y_04<-c(y0,   y1,   y1+ry1,    y0+ry1)+ry0
polygon(x_04,y_04,col=col3,lty=lty,lwd=lwd,border=border)
if(!is.null(density)) polygon(x_04,y_04,col=lcol,lty=lty,lwd=lwd,border=border, density=density, angle=angle+45)

x_02<-c(x0,   x0,   x1,    x1)+rx0
y_02<-c(y0,   y1,   y1,    y0)+ry0
polygon(x_02,y_02,col=col1,lty=lty,lwd=lwd,border=border)
if(!is.null(density)) polygon(x_02,y_02,col=lcol,lty=lty,lwd=lwd,border=border, density=density, angle=angle)
}
##____________________________________________________________________________##
################################################################################


#######################
#Colors
zn<-nrow(tab)
xn<-ncol(tab)
if(is.na(zn)) zn<-1
if(is.na(xn)) xn<-1
if(is.null(alpha) & wall==0) alpha<-1
if(is.null(alpha) & wall!=0) alpha<-0.5
if(is.null(col)) col <- a.getcol.ade(zn*xn, alpha=alpha)
if(!is.null(col) & length(col)==1)   col<-matrix(col, nrow = zn, ncol =xn)
if(!is.null(col) & length(col)==zn & length(col)!=xn)  col<-matrix(rep(col, xn), nrow = zn, ncol =xn)
if(!is.null(col) & length(col)==xn & length(col)!=zn)  col<-matrix(rep(col, zn), nrow = zn, ncol =xn, byrow = T)
if(!is.null(col) & length(col)==xn & length(col)==zn)  col<-matrix(rep(col, zn), nrow = zn, ncol =xn)
if(!is.null(col) & length(col)==xn*zn )  col<-matrix(col, nrow = zn, ncol =xn)
if(!is.null(col) & length(col)<zn) col <- rep(col, round((zn/length(col))+1))
if(!is.null(alpha)) col<- apply(col, 2, a.alpha.ade,  alpha)


#####################################
# Colors
if(is.null(tcol)  & wall==0)   tcol<-1
if(is.null(tcol)  & wall!=0)   tcol<-rgb(0.1,0.1,0.25)
if(is.null(bgcol) & wall==0)   bgcol<-1
if(is.null(bgcol) & wall!=0)   bgcol<-'#DBE0E8'
if(is.null(lcol))  lcol<- tcol
#####################################


#######################
ylim<-c(0, max(tab, na.rm=TRUE))
xlim<-c(0.5,(xn+0.5))
zlim<-c(-0.5,zn+0.5)
################################################################################



################################################################################
# Axis Labels
if(is.null(xticks) & !is.null(colnames(tab))) xlabels<-colnames(tab)
if(!is.null(xticks)) xlabels<-xticks
if(!is.null(zticks)) zlabels<-zticks
if(is.null(zticks) & !is.null(rownames(tab))) zlabels<-rownames(tab)
if(is.null(yticks))ylabels<-pretty(ylim, 8)
if(!is.null(yticks) & length(yticks)==1) ylabels<-pretty(ylim, yticks)
if(!is.null(yticks) & length(yticks)>1) ylabels<-yticks
ylabels<-ylabels[ylabels<=ylim[2] & ylabels>=ylim[1]]
yat<-ylabels


################################################################################

################################################################################
# Calculations
tieff<-(1:zn)-1
xx<-1:xn
zw<-zw*2
zx<-zw*zn * diff(range(xlim))/40
zy<-zw*zn * diff(range(ylim))/40
xlim2<-xlim
ylim2<-ylim
xlim2[2]<-xlim2[2]+zx
ylim2[2]<-ylim2[2]+zy
zss<- (tieff+0.5)*(zx/zn)+xlim[2]
x1<-c(xlim[1],xlim[1]+zx,xlim[2]+zx,xlim[2])
y1<-c(ylim[1],ylim[1]+zy,ylim[1]+zy,ylim[1])
x2<-c(xlim[1],xlim[1],xlim[1]+zx,xlim[1]+zx)
y2<-c(ylim[1],ylim[2],ylim[2]+zy,ylim[1]+zy)
x3<-c(xlim[1]+zx,xlim[1]+zx,xlim[2]+zx,xlim[2]+zx)
y3<-c(ylim[1]+zy,ylim[2]+zy,ylim[2]+zy,ylim[1]+zy)
################################################################################
################################################################################
################################################################################




################################################################################
if(wall==0){
par(col.axis=tcol)
plot(0,0, type='n', col=rgb(1,1,1,0), xlim=xlim2, ylim=ylim2, xlab='', ylab='',  main='',  axes=F)

#Rahmen
if(bgbox){
polygon(x1,y1, col=NA, lty=1, lwd=1, border=bgcol)  #Boden
polygon(x2,y2, col=NA, lty=1, lwd=1, border=bgcol)  #Linke Wand
polygon(x3,y3, col=NA, lty=1, lwd=1, border=bgcol)  #Hintere Wand
}


# Cubes
for(j in length(tieff):1){
for(i in 1:length(xx)){
a.bar3d.ade(x0=xx[i]-xw/2, x1=xx[i]+xw/2,  y0=0, y1=tab[j,i], z0=(tieff[j]), z1=(tieff[j]+1),  zdist=zw  , hx=(zx/zn)/(zw), hy=(zy/zn)/(zw), col=col[j, i], border=lcol, lcol=lcol, density=0, angle=0, lty = 1, lwd = 1, lty2 = 1, lwd2 = 1, backside=TRUE)
}
}

if(axes){
axis(1, at=xx, labels=xlabels, line=-0.9, tick=TRUE,   col=rgb(0,0,0,0), col.ticks=bgcol)
axis(2, at=yat, labels=ylabels, line=-0.95, tick=TRUE, col=rgb(0,0,0,0), col.ticks=bgcol, las=1)
segments(zss, (tieff+0.5)*((zy/zn)),  zss+(diff(range(xlim))/40),  (tieff+0.5)*((zy/zn)), col=bgcol, lty=1)
text(zss+(diff(range(xlim))/20), (tieff+0.5)*((zy/zn)), labels=zlabels, adj=0, xpd=TRUE, col=tcol, cex=0.9)
}

# Überschriften
par(xpd=TRUE)
dx<-8/par('din')[1]
dy<-8/par('din')[2]
xr<-(diff(par('usr')[1:2])/10)*dx
yr<-(diff(par('usr')[3:4])/10)*dy
text(par('usr')[1]-xr, y= mean(ylim) , labels=ylab, cex = 1, font=2, col=tcol, srt=90, xpd=TRUE)
text(mean(xlim), y= par('usr')[3]-yr ,labels=xlab, cex = 1, font=2, col=tcol, xpd=TRUE)
text(mean(zss+(diff(range(xlim))/4.5)), ylim[1]+(zy/2), labels=zlab, cex = 1, font=2, col=tcol, xpd=TRUE, srt=45, adj=0.5)
text(mean(xlim)+zx, y= par('usr')[4]+yr/1.85 ,  labels=main, cex = 1.25, font=2, col=tcol, xpd=TRUE)
par(xpd=FALSE)

# PseudoWände, Drüber
if(fgbox & bgbox){
segments(xlim[2], ylim[1], xlim[2],    ylim[2],    col=bgcol, lty=3)
segments(xlim[1], ylim[2], xlim[2],    ylim[2],    col=bgcol, lty=3)
segments(xlim[2], ylim[2], xlim[2]+zx, ylim[2]+zy, col=bgcol, lty=3)
}
}
################################################################################


################################################################################
if(wall==1){
par(col.axis=tcol)
plot(0,0, type='n', col=rgb(1,1,1,0), xlim=xlim2, ylim=ylim2, xlab='', ylab='',  main='',  axes=F)

#Rahmen
if(bgbox){
polygon(x1,y1, col=bgcol, lty=1, lwd=1, border=rgb(1,1,1), lwd=1)  #Boden
polygon(x2,y2, col=bgcol, lty=1, lwd=1, border=rgb(1,1,1), lwd=1)  #Linke Wand
polygon(x3,y3, col=bgcol, lty=1, lwd=1, border=rgb(1,1,1), lwd=1)  #Hintere Wand
segments(xlim[1], yat,  xlim[1]+zx,  yat+zy, col=rgb(1,1,1), lty=1, lwd=1)
segments(xlim[1]+zx, yat+zy,  xlim[2]+zx,  yat+zy, col=rgb(1,1,1), lty=1)
segments(xx+0.5, 0,  xx+zx+0.5,  zy, col=rgb(1,1,1), lty=3, lwd=1)
segments(xlim[1]+(tieff)*((zx/zn)), (tieff)*((zy/zn)),  xlim[2]+(tieff)*((zx/zn)),  (tieff)*((zy/zn)), col=rgb(1,1,1), lty=3, lwd=1)
}


# Cubes
for(j in length(tieff):1){
for(i in 1:length(xx)){
a.bar3d.ade(x0=xx[i]-xw/2, x1=xx[i]+xw/2,  y0=0, y1=tab[j,i], z0=(tieff[j]), z1=(tieff[j]+1),  zdist=zw  , hx=(zx/zn)/(zw), hy=(zy/zn)/(zw), col=col[j, i], border=lcol, lcol=lcol, density=0, angle=0, lty = 1, lwd = 1, lty2 = 1, lwd2 = 1, backside=TRUE)
}
}

if(axes){
axis(1, at=xx, labels=xlabels, line=-0.9, tick=TRUE,   col=rgb(0,0,0,0), col.ticks=bgcol)
axis(2, at=yat, labels=ylabels, line=-0.95, tick=TRUE, col=rgb(0,0,0,0), col.ticks=bgcol, las=1)
segments(zss, (tieff+0.5)*((zy/zn)),  zss+(diff(range(xlim))/40),  (tieff+0.5)*((zy/zn)), col=bgcol, lty=1)
text(zss+(diff(range(xlim))/20), (tieff+0.5)*((zy/zn)), labels=zlabels, adj=0, xpd=TRUE, col=tcol, cex=0.9)
}




# Überschriften
par(xpd=TRUE)
dx<-8/par('din')[1]
dy<-8/par('din')[2]
xr<-(diff(par('usr')[1:2])/10)*dx
yr<-(diff(par('usr')[3:4])/10)*dy
text(par('usr')[1]-xr, y= mean(ylim) , labels=ylab, cex = 1, font=2, col=tcol, srt=90, xpd=TRUE)
text(mean(xlim), y= par('usr')[3]-yr ,labels=xlab, cex = 1, font=2, col=tcol, xpd=TRUE)
text(mean(zss+(diff(range(xlim))/4.5)), ylim[1]+(zy/2), labels=zlab, cex = 1, font=2, col=tcol, xpd=TRUE, srt=45, adj=0.5)
text(mean(xlim)+zx, y= par('usr')[4]+yr/1.85 ,  labels=main, cex = 1.25, font=2, col=tcol, xpd=TRUE)
par(xpd=FALSE)

# PseudoWände, Drüber
if(fgbox & bgbox){
segments(xlim[2], ylim[1], xlim[2],    ylim[2],    col=rgb(1,1,1), lty=3)
segments(xlim[1], ylim[2], xlim[2],    ylim[2],    col=rgb(1,1,1), lty=3)
segments(xlim[2], ylim[2], xlim[2]+zx, ylim[2]+zy, col=rgb(1,1,1), lty=3)
}
}
################################################################################



################################################################################
if(wall==2){
par(col.axis=tcol)
plot(0,0, type='n', col=rgb(1,1,1,0), xlim=xlim2, ylim=ylim2, xlab='', ylab='',  main='',  axes=F)

#Rahmen
if(bgbox){
segments(xlim[1], yat,  xlim[1]+zx,  yat+zy, col=bgcol, lty=1, lwd=1)
segments(xlim[1]+zx, yat+zy,  xlim[2]+zx,  yat+zy, col=bgcol, lty=1)
segments(xx+0.5, 0,  xx+zx+0.5,  zy, col=bgcol, lty=2, lwd=1)
segments(xlim[1]+(tieff)*((zx/zn)), (tieff)*((zy/zn)),  xlim[2]+(tieff)*((zx/zn)),  (tieff)*((zy/zn)), col=bgcol, lty=2, lwd=1)
polygon(x1,y1, col=rgb(1,1,1, 0), lty=1, lwd=1, border=a.coladd.ade(bgcol, -75), lwd=1)  #Boden
polygon(x2,y2, col=rgb(1,1,1, 0), lty=1, lwd=1, border=a.coladd.ade(bgcol, -75), lwd=1)  #Linke Wand
polygon(x3,y3, col=rgb(1,1,1, 0), lty=1, lwd=1, border=a.coladd.ade(bgcol, -75), lwd=1)  #Hintere Wand
}


# Cubes
for(j in length(tieff):1){
for(i in 1:length(xx)){
a.bar3d.ade(x0=xx[i]-xw/2, x1=xx[i]+xw/2,  y0=0, y1=tab[j,i], z0=(tieff[j]), z1=(tieff[j]+1),  zdist=zw  , hx=(zx/zn)/(zw), hy=(zy/zn)/(zw), col=col[j, i], border=lcol, lcol=lcol, density=0, angle=0, lty = 1, lwd = 1, lty2 = 1, lwd2 = 1, backside=TRUE)
}
}

if(axes){
axis(1, at=xx, labels=xlabels, line=-0.9, tick=TRUE,   col=rgb(0,0,0,0), col.ticks=a.coladd.ade(bgcol, -75))
axis(2, at=yat, labels=ylabels, line=-0.95, tick=TRUE, col=rgb(0,0,0,0), col.ticks=a.coladd.ade(bgcol, -75), las=1)
segments(zss, (tieff+0.5)*((zy/zn)),  zss+(diff(range(xlim))/40),  (tieff+0.5)*((zy/zn)), col=a.coladd.ade(bgcol, -75), lty=1)
text(zss+(diff(range(xlim))/20), (tieff+0.5)*((zy/zn)), labels=zlabels, adj=0, xpd=TRUE, col=tcol, cex=0.9)
}




# Überschriften
par(xpd=TRUE)
dx<-8/par('din')[1]
dy<-8/par('din')[2]
xr<-(diff(par('usr')[1:2])/10)*dx
yr<-(diff(par('usr')[3:4])/10)*dy
text(par('usr')[1]-xr, y= mean(ylim) , labels=ylab, cex = 1, font=2, col=tcol, srt=90, xpd=TRUE)
text(mean(xlim), y= par('usr')[3]-yr ,labels=xlab, cex = 1, font=2, col=tcol, xpd=TRUE)
text(mean(zss+(diff(range(xlim))/4.5)), ylim[1]+(zy/2), labels=zlab, cex = 1, font=2, col=tcol, xpd=TRUE, srt=45, adj=0.5)
text(mean(xlim)+zx, y= par('usr')[4]+yr/1.85 ,  labels=main, cex = 1.25, font=2, col=tcol, xpd=TRUE)
par(xpd=FALSE)

# PseudoWände, Drüber
if(fgbox & bgbox){
segments(xlim[2], ylim[1], xlim[2],    ylim[2],    col=a.coladd.ade(bgcol, -75), lty=3)
segments(xlim[1], ylim[2], xlim[2],    ylim[2],    col=a.coladd.ade(bgcol, -75), lty=3)
segments(xlim[2], ylim[2], xlim[2]+zx, ylim[2]+zy, col=a.coladd.ade(bgcol, -75), lty=3)
}
}
################################################################################



################################################################################
if(wall==3){
par(col.axis=tcol)
plot(0,0, type='n', col=rgb(1,1,1,0), xlim=xlim2, ylim=ylim2, xlab='', ylab='',  main='',  axes=F)

#Rahmen
if(bgbox){
polygon(x1,y1, col=bgcol, lty=1, lwd=1, border=a.coladd.ade(bgcol, -50), lwd=1)  #Boden
polygon(x2,y2, col=bgcol, lty=1, lwd=1, border=a.coladd.ade(bgcol, -50), lwd=1)  #Linke Wand
polygon(x3,y3, col=bgcol, lty=1, lwd=1, border=a.coladd.ade(bgcol, -50), lwd=1)  #Hintere Wand
segments(xlim[1], yat,  xlim[1]+zx,  yat+zy, col=a.coladd.ade(bgcol, -50), lty=1, lwd=1)
segments(xlim[1]+zx, yat+zy,  xlim[2]+zx,  yat+zy, col=a.coladd.ade(bgcol, -50), lty=1)
segments(xx+0.5, 0,  xx+zx+0.5,  zy, col=a.coladd.ade(bgcol, -50), lty=2, lwd=1)
segments(xlim[1]+(tieff)*((zx/zn)), (tieff)*((zy/zn)),  xlim[2]+(tieff)*((zx/zn)),  (tieff)*((zy/zn)), col=a.coladd.ade(bgcol, -50), lty=2, lwd=1)

}


# Cubes
for(j in length(tieff):1){
for(i in 1:length(xx)){
a.bar3d.ade(x0=xx[i]-xw/2, x1=xx[i]+xw/2,  y0=0, y1=tab[j,i], z0=(tieff[j]), z1=(tieff[j]+1),  zdist=zw  , hx=(zx/zn)/(zw), hy=(zy/zn)/(zw), col=col[j, i], border=lcol, lcol=lcol, density=0, angle=0, lty = 1, lwd = 1, lty2 = 1, lwd2 = 1, backside=TRUE)
}
}

if(axes){
axis(1, at=xx, labels=xlabels, line=-0.9, tick=TRUE,   col=rgb(0,0,0,0), col.ticks=a.coladd.ade(bgcol, -50))
axis(2, at=yat, labels=ylabels, line=-0.95, tick=TRUE, col=rgb(0,0,0,0), col.ticks=a.coladd.ade(bgcol, -50), las=1)
segments(zss, (tieff+0.5)*((zy/zn)),  zss+(diff(range(xlim))/40),  (tieff+0.5)*((zy/zn)), col=a.coladd.ade(bgcol, -50), lty=1)
text(zss+(diff(range(xlim))/20), (tieff+0.5)*((zy/zn)), labels=zlabels, adj=0, xpd=TRUE, col=tcol, cex=0.9)
}




# Überschriften
par(xpd=TRUE)
dx<-8/par('din')[1]
dy<-8/par('din')[2]
xr<-(diff(par('usr')[1:2])/10)*dx
yr<-(diff(par('usr')[3:4])/10)*dy
text(par('usr')[1]-xr, y= mean(ylim) , labels=ylab, cex = 1, font=2, col=tcol, srt=90, xpd=TRUE)
text(mean(xlim), y= par('usr')[3]-yr ,labels=xlab, cex = 1, font=2, col=tcol, xpd=TRUE)
text(mean(zss+(diff(range(xlim))/4.5)), ylim[1]+(zy/2), labels=zlab, cex = 1, font=2, col=tcol, xpd=TRUE, srt=45, adj=0.5)
text(mean(xlim)+zx, y= par('usr')[4]+yr/1.85 ,  labels=main, cex = 1.25, font=2, col=tcol, xpd=TRUE)
par(xpd=FALSE)

# PseudoWände, Drüber
if(fgbox & bgbox){
segments(xlim[2], ylim[1], xlim[2],    ylim[2],    col=a.coladd.ade(bgcol, -50), lty=3)
segments(xlim[1], ylim[2], xlim[2],    ylim[2],    col=a.coladd.ade(bgcol, -50), lty=3)
segments(xlim[2], ylim[2], xlim[2]+zx, ylim[2]+zy, col=a.coladd.ade(bgcol, -50), lty=3)
}
}
################################################################################




################################################################################
if(wall==4){
par(col.axis=tcol)
plot(0,0, type='n', col=rgb(1,1,1,0), xlim=xlim2, ylim=ylim2, xlab='', ylab='',  main='',  axes=F)

#Rahmen
if(bgbox){
polygon(x1,y1, col=tcol, lty=1, lwd=1, border=rgb(1,1,1),  lwd=1)  #Boden
polygon(x2,y2, col=bgcol, lty=1, lwd=1, border=rgb(1,1,1), lwd=1)  #Linke Wand
polygon(x3,y3, col=bgcol, lty=1, lwd=1, border=rgb(1,1,1), lwd=1)  #Hintere Wand
segments(xlim[1], yat,  xlim[1]+zx,  yat+zy, col=rgb(1,1,1), lty=1, lwd=1)
segments(xlim[1]+zx, yat+zy,  xlim[2]+zx,  yat+zy, col=rgb(1,1,1), lty=1)
segments(xx+0.5, 0,  xx+zx+0.5,  zy, col=rgb(1,1,1), lty=1, lwd=1)
segments(xlim[1]+(tieff)*((zx/zn)), (tieff)*((zy/zn)),  xlim[2]+(tieff)*((zx/zn)),  (tieff)*((zy/zn)), col=rgb(1,1,1), lty=1, lwd=1)

}


# Cubes
for(j in length(tieff):1){
for(i in 1:length(xx)){
a.bar3d.ade(x0=xx[i]-xw/2, x1=xx[i]+xw/2,  y0=0, y1=tab[j,i], z0=(tieff[j]), z1=(tieff[j]+1),  zdist=zw  , hx=(zx/zn)/(zw), hy=(zy/zn)/(zw), col=col[j, i], border=lcol, lcol=lcol, density=0, angle=0, lty = 1, lwd = 1, lty2 = 1, lwd2 = 1, backside=TRUE)
}
}

if(axes){
axis(1, at=xx, labels=xlabels, line=-0.9, tick=TRUE,   col=rgb(0,0,0,0), col.ticks=a.coladd.ade(bgcol, -50))
axis(2, at=yat, labels=ylabels, line=-0.95, tick=TRUE, col=rgb(0,0,0,0), col.ticks=a.coladd.ade(bgcol, -50), las=1)
segments(zss, (tieff+0.5)*((zy/zn)),  zss+(diff(range(xlim))/40),  (tieff+0.5)*((zy/zn)), col=a.coladd.ade(bgcol, -50), lty=1)
text(zss+(diff(range(xlim))/20), (tieff+0.5)*((zy/zn)), labels=zlabels, adj=0, xpd=TRUE, col=tcol, cex=0.9)
}




# Überschriften
par(xpd=TRUE)
dx<-8/par('din')[1]
dy<-8/par('din')[2]
xr<-(diff(par('usr')[1:2])/10)*dx
yr<-(diff(par('usr')[3:4])/10)*dy
if(!is.null(ylab)) polygon( c(par('usr')[1]-xr/1.5, par('usr')[1]-xr/1.5, par('usr')[1]-xr*1.3,par('usr')[1]-xr*1.3 ), c(ylim[1],ylim[2],ylim[2],ylim[1] ), col=bgcol, border=rgb(1,1,1))
if(!is.null(xlab)) polygon( c(xlim[c(1,1,2,2)]), c(par('usr')[3]-yr*1.35,par('usr')[3]-yr/1.4,par('usr')[3]-yr/1.4,par('usr')[3]-yr*1.35 ), col=bgcol, border=rgb(1,1,1))
if(!is.null(main)) polygon( c(xlim[c(1,1,2,2)])+zx, c(ylim[2]+zy,par('usr')[4]+yr/1.4,par('usr')[4]+yr/1.4,ylim[2]+zy ), col=tcol, border=rgb(1,1,1))
text(par('usr')[1]-xr, y= mean(ylim) , labels=ylab, cex = 1, font=2, col=tcol, srt=90, xpd=TRUE)
text(mean(xlim), y= par('usr')[3]-yr ,labels=xlab, cex = 1, font=2, col=tcol, xpd=TRUE)
text(mean(zss+(diff(range(xlim))/4.5)), ylim[1]+(zy/2), labels=zlab, cex = 1, font=2, col=tcol, xpd=TRUE, srt=45, adj=0.5)
text(mean(c(xlim[c(1,1,2,2)])+zx), y= par('usr')[4]+yr/4.25 ,  labels=main, cex = 1.25, font=2, col=rgb(1,1,1), xpd=TRUE)
par(xpd=FALSE)

# PseudoWände, Drüber
if(fgbox & bgbox){
segments(xlim[2], ylim[1], xlim[2],    ylim[2],    col=rgb(1,1,1), lty=3, lwd=1)
segments(xlim[1], ylim[2], xlim[2],    ylim[2],    col=rgb(1,1,1), lty=3, lwd=1)
segments(xlim[2], ylim[2], xlim[2]+zx, ylim[2]+zy, col=rgb(1,1,1), lty=3, lwd=1)
}
}
################################################################################


################################################################################
if(wall==5){
par(col.axis=tcol)
plot(0,0, type='n', col=rgb(1,1,1,0), xlim=xlim2, ylim=ylim2, xlab='', ylab='',  main='',  axes=F)
#Rahmen
if(bgbox){
polygon(x1,y1, col=rgb(1,1,1,0), lty=1, lwd=1, border=tcol,  lwd=1)  #Boden
polygon(x2,y2, col=rgb(1,1,1,0), lty=1, lwd=1, border=tcol, lwd=1)  #Linke Wand
polygon(x3,y3, col=rgb(1,1,1,0), lty=1, lwd=1, border=tcol, lwd=1)  #Hintere Wand
segments(xx+0.5, 0,  xx+zx+0.5,  zy, col=tcol, lty=1, lwd=1)
segments(xlim[1]+(tieff)*((zx/zn)), (tieff)*((zy/zn)),  xlim[2]+(tieff)*((zx/zn)),  (tieff)*((zy/zn)), col=tcol, lty=1, lwd=1)

}


# Cubes
for(j in length(tieff):1){
for(i in 1:length(xx)){
a.bar3d.ade(x0=xx[i]-xw/2, x1=xx[i]+xw/2,  y0=0, y1=tab[j,i], z0=(tieff[j]), z1=(tieff[j]+1),  zdist=zw  , hx=(zx/zn)/(zw), hy=(zy/zn)/(zw), col=col[j, i], border=lcol, lcol=lcol, density=0, angle=0, lty = 1, lwd = 1, lty2 = 1, lwd2 = 1, backside=TRUE)
}
}

if(axes){
axis(1, at=xx, labels=xlabels, line=-0.9, tick=TRUE,   col=rgb(0,0,0,0), col.ticks=a.coladd.ade(bgcol, -50))
axis(2, at=yat, labels=ylabels, line=-0.95, tick=TRUE, col=rgb(0,0,0,0), col.ticks=a.coladd.ade(bgcol, -50), las=1)
segments(zss, (tieff+0.5)*((zy/zn)),  zss+(diff(range(xlim))/40),  (tieff+0.5)*((zy/zn)), col=a.coladd.ade(bgcol, -50), lty=1)
text(zss+(diff(range(xlim))/20), (tieff+0.5)*((zy/zn)), labels=zlabels, adj=0, xpd=TRUE, col=tcol, cex=0.9)
}




# Überschriften
par(xpd=TRUE)
dx<-8/par('din')[1]
dy<-8/par('din')[2]
xr<-(diff(par('usr')[1:2])/10)*dx
yr<-(diff(par('usr')[3:4])/10)*dy
polygon( c(xlim[c(1,1,2,2)])+zx, c(ylim[2]+yr/4,ylim[2]+yr,ylim[2]+yr,ylim[2]+yr/4  )+zy, col=rgb(1,1,1,0), border=tcol)
polygon( c(xlim[1], xlim[1], xlim[1]+zx, xlim[1]+zx), c(ylim[2]+yr/4,ylim[2]+yr,ylim[2]+yr+zy,ylim[2]+yr/4+zy  ), col=bgcol, border=tcol)
polygon( c(xlim[2],xlim[2],xlim[2]+xr/5,xlim[2]+xr/5)+zx, c(ylim[2]+yr/4,ylim[2]+yr,ylim[2]+yr,ylim[2]+yr/4  )+zy, col=bgcol, border=tcol)
polygon( c(xlim[2],xlim[2],xlim[2]+xr/5,xlim[2]+xr/5 )+zx, c(ylim[1],ylim[2],ylim[2],ylim[1] )+zy, col=bgcol, border=tcol)
polygon( c(par('usr')[1]-xr*1.5,par('usr')[1]-xr*1.5,par('usr')[1]-xr*1.3,par('usr')[1]-xr*1.3 ), c(ylim[1]-(yr),ylim[2]+(yr/4),ylim[2]+(yr/4),ylim[1]-(yr) ), col=bgcol, border=tcol)
polygon( c(par('usr')[1]-xr*1.5,par('usr')[1]-xr*1.5,xlim[1],xlim[1] ), c(ylim[2]+(yr/4),ylim[2]+yr,ylim[2]+yr,ylim[2]+(yr/4) ), col=bgcol, border=tcol)
polygon( c(par('usr')[1]-xr*1.5,par('usr')[1]-xr*1.5,xlim[1],xlim[1] ), c(ylim[1]-(yr*1.75),ylim[1]-yr,ylim[1]-yr,ylim[1]-(yr*1.75) ), col=bgcol, border=tcol)
polygon( c(xlim[1],xlim[1], xlim[2],xlim[2]),  c(ylim[1]-(yr*1.75),ylim[1]-yr,ylim[1]-yr,ylim[1]-(yr*1.75) ), col=rgb(1,1,1,0), border=tcol)
polygon( c(xlim[2],xlim[2],xlim[2]+xr/5,xlim[2]+xr/5 ), c(ylim[1]-(yr*1.75),ylim[1]-yr,ylim[1]-yr,ylim[1]-(yr*1.75) ), col=bgcol, border=tcol)

text(par('usr')[1]-xr, y= mean(ylim) , labels=ylab, cex = 1, font=2, col=tcol, srt=90, xpd=TRUE)
text(mean(xlim), y= par('usr')[3]-yr ,labels=xlab, cex = 1, font=2, col=tcol, xpd=TRUE)
text(mean(zss+(diff(range(xlim))/4.5)), ylim[1]+(zy/2), labels=zlab, cex = 1, font=2, col=tcol, xpd=TRUE, srt=45, adj=0.5)
text(mean(c(xlim[c(1,1,2,2)])+zx), y= par('usr')[4]+yr/4 ,  labels=main, cex = 1.25, font=2, col=tcol, xpd=TRUE)
par(xpd=FALSE)

# PseudoWände, Drüber
if(fgbox & bgbox){
segments(xlim[2], ylim[1], xlim[2],    ylim[2],    col=tcol, lty=3, lwd=1)
segments(xlim[1], ylim[2], xlim[2],    ylim[2],    col=tcol, lty=3, lwd=1)
segments(xlim[2], ylim[2], xlim[2]+zx, ylim[2]+zy, col=tcol, lty=3, lwd=1)
}
}
################################################################################


################################################################################
if(wall==6){
par(col.axis=tcol)
plot(0,0, type='n', col=rgb(1,1,1,0), xlim=xlim2, ylim=ylim2, xlab='', ylab='',  main='',  axes=F)

#Rahmen
if(bgbox){
polygon(x1,y1, col=bgcol, lty=1, lwd=1, border=a.coladd.ade(bgcol, -35), lwd=1)  #Boden
polygon(x2,y2, col=bgcol, lty=1, lwd=1, border=a.coladd.ade(bgcol, -35), lwd=1)  #Linke Wand
polygon(x3,y3, col=bgcol, lty=1, lwd=1, border=a.coladd.ade(bgcol, -35), lwd=1)  #Hintere Wand

segments(xlim[1], yat,  xlim[1]+zx,  yat+zy, col=a.coladd.ade(bgcol, -35), lty=1, lwd=3)
segments(xlim[1], yat,  xlim[1]+zx,  yat+zy, col=rgb(1,1,1), lty=1, lwd=1)

segments(xlim[1]+zx, yat+zy,  xlim[2]+zx,  yat+zy, col=a.coladd.ade(bgcol, -35), lty=1, lwd=3)
segments(xlim[1]+zx, yat+zy,  xlim[2]+zx,  yat+zy, col=rgb(1,1,1), lty=1, lwd=1)

segments(xlim[1]+(tieff)*((zx/zn)), (tieff)*((zy/zn)),  xlim[2]+(tieff)*((zx/zn)),  (tieff)*((zy/zn)), col=a.coladd.ade(bgcol, -35), lty=1, lwd=3)
segments(xlim[1]+(tieff)*((zx/zn)), (tieff)*((zy/zn)),  xlim[2]+(tieff)*((zx/zn)),  (tieff)*((zy/zn)), col=rgb(1,1,1), lty=1, lwd=1)

}

# Cubes
for(j in length(tieff):1){
for(i in 1:length(xx)){
a.bar3d.ade(x0=xx[i]-xw/2, x1=xx[i]+xw/2,  y0=0, y1=tab[j,i], z0=(tieff[j]), z1=(tieff[j]+1),  zdist=zw  , hx=(zx/zn)/(zw), hy=(zy/zn)/(zw), col=col[j, i], border=lcol, lcol=lcol, density=0, angle=0, lty = 1, lwd = 1, lty2 = 1, lwd2 = 1, backside=TRUE)
}
}

if(axes){
axis(1, at=xx, labels=xlabels, line=-0.9, tick=TRUE,   col=rgb(0,0,0,0), col.ticks=a.coladd.ade(bgcol, -35), lwd.ticks=3)
axis(1, at=xx, labels=xlabels, line=-0.9, tick=TRUE,   col=rgb(0,0,0,0), col.ticks=rgb(1,1,1), lwd.ticks=1)
axis(2, at=yat, labels=ylabels, line=-0.95, tick=TRUE, col=rgb(0,0,0,0), col.ticks=a.coladd.ade(bgcol, -35), las=1, lwd.ticks=3)
axis(2, at=yat, labels=ylabels, line=-0.95, tick=TRUE, col=rgb(0,0,0,0), col.ticks=rgb(1,1,1), las=1, lwd.ticks=1)
segments(zss, (tieff+0.5)*((zy/zn)),  zss+(diff(range(xlim))/40),  (tieff+0.5)*((zy/zn)), col=a.coladd.ade(bgcol, -35), lty=1, lwd=3)
segments(zss, (tieff+0.5)*((zy/zn)),  zss+(diff(range(xlim))/40),  (tieff+0.5)*((zy/zn)), col=rgb(1,1,1), lty=1)

segments(xlim[2], ylim[1],  xlim[2]+zx,  ylim[1]+zy, col=a.coladd.ade(bgcol, -35), lwd=3)
segments(xlim[2], ylim[1],  xlim[2]+zx,  ylim[1]+zy, col=rgb(1,1,1), lwd=1)

text(zss+(diff(range(xlim))/20), (tieff+0.5)*((zy/zn)), labels=zlabels, adj=0, xpd=TRUE, col=tcol, cex=0.9)
}



# Überschriften
par(xpd=TRUE)
dx<-8/par('din')[1]
dy<-8/par('din')[2]
xr<-(diff(par('usr')[1:2])/10)*dx
yr<-(diff(par('usr')[3:4])/10)*dy
text(par('usr')[1]-xr, y= mean(ylim) , labels=ylab, cex = 1, font=2, col=tcol, srt=90, xpd=TRUE)
text(mean(xlim), y= par('usr')[3]-yr ,labels=xlab, cex = 1, font=2, col=tcol, xpd=TRUE)
text(mean(zss+(diff(range(xlim))/4.5)), ylim[1]+(zy/2), labels=zlab, cex = 1, font=2, col=tcol, xpd=TRUE, srt=45, adj=0.5)
text(mean(xlim)+zx, y= par('usr')[4]+yr/1.85 ,  labels=main, cex = 1.25, font=2, col=tcol, xpd=TRUE)
par(xpd=FALSE)

# PseudoWände, Drüber
if(fgbox & bgbox){
segments(xlim[2], ylim[1], xlim[2],    ylim[2],    col=rgb(1,1,1), lty=3)
segments(xlim[1], ylim[2], xlim[2],    ylim[2],    col=rgb(1,1,1), lty=3)
segments(xlim[2], ylim[2], xlim[2]+zx, ylim[2]+zy, col=rgb(1,1,1), lty=3)
}
}
################################################################################
}
