ratio.plot.ade <-
function( M, vnames=NULL, sectext=NULL, main=NULL,xlab=NULL, ylab=NULL, legenlab=NULL,  rlab=NULL,  col=NULL, tcol=NULL,  bgcol=NULL,  lcol=NULL, r=NULL, v=c(0,1), lty=c(1,2), xticks=18,  hlines=TRUE,  legends=TRUE, logaxe=FALSE, wall=0){
if(any(par('mfg')!=c(1,1,1,1)) & any(par('mai') < c(1.02, 0.82, 0.82, 0.42))){
maidiff<-rep(0, 4)
norm<-c(1.02, 0.82, 0.82, 0.42)
maidiff[par('mai')<norm]<-  norm[par('mai')<norm] - par('mai')[par('mai')<norm]
par(mai=par('mai')+maidiff)
}
oldpar<-par(no.readonly =TRUE)
oldpar<-oldpar[-which(names(oldpar)%in%c('usr', 'plt', 'pin', 'fin', 'fig', 'mfg', 'mfcol', 'mfrow', 'omd', 'omi', 'oma'))]
on.exit(par(oldpar))



if(!is.list(M)){
ML<-NULL
ML<-as.list(ML)
ML[[1]]<-M
M<-ML
}

if(is.null(vnames)) legends=FALSE

n  <- dim(as.matrix(M[[1]]))[1]
if(!is.null(vnames)) vnames<-as.character(as.vector(vnames[1:n]))
if(is.null(vnames))  vnames<-''
N<-length(M)
if(length(legenlab)<N)   legenlab<- c(legenlab, rep('?', N-length(legenlab)))

################################################################################
# Colors
if(length(col)<N) col<-NULL
if(is.null(tcol)  & wall==0)   tcol<-1
if(is.null(tcol)  & wall!=0)   tcol<-rgb(0.1,0.1,0.25)
if(is.null(bgcol) & wall==0)   bgcol<-1
if(is.null(bgcol) & wall!=0)   bgcol<-'#DBE0E8'
if(is.null(lcol) & (wall==0 | wall==2| wall==5))              lcol<-bgcol
if(is.null(lcol) & (wall==1 | wall==4))   lcol<-rgb(1,1,1)
if(is.null(lcol) & (wall==3))   lcol<-a.coladd.ade(bgcol, -50)
if(is.null(col) & N==1) col <- tcol
if(is.null(col) & N>1)  col <- a.getcol.ade(N)

rcol=col
col2<-a.coladd.ade(col, 175)
fcol= a.coladd.ade(bgcol, -35)
fcol2=a.coladd.ade(bgcol, -100)
fcol3=bgcol
bgcol2<-a.coladd.ade(bgcol, -50)

#
################################################################################

if(is.list(M)){
for(k in 1:length(M)) {
M[[k]]<- apply(as.matrix(M[[k]]), c(1,2), as.numeric)
if(logaxe) M[[k]]<-log(M[[k]])
}
xmin <- min(unlist(M), na.rm=TRUE)
xmax <- max(unlist(M), na.rm=TRUE)
}

if(is.matrix(M)){
xmin <- min(M, na.rm=TRUE)
xmax <- max(M, na.rm=TRUE)
}

if(is.data.frame(M)){
M<-as.matrix(M)
xmin <- min(M, na.rm=TRUE)
xmax <- max(M, na.rm=TRUE)
}




if(logaxe) v=log(v[v!=0])
xrange <- xmax - xmin


if(!is.null(r)){
r<-(r*2.5)
print(xmin)
print(xmax)
xlimw <- c(xmin , (xmax+(xrange*r)))


}


if(is.null(r)){
vlength <-max(nchar(c(vnames, rlab, sectext)))
if(wall==4) r= vlength*0.025 + (sqrt((vlength^4))*0.0004)+0.18


if(wall!=4) r= vlength*0.025 + (sqrt((vlength^4))*0.0002)+0.18


xlimw <- c(xmin , (xmax+(xrange*r)))
}









schift <- (xrange)/6
if(N==1) legends<-F
if(legends)  ylimw <- c(0.5, (n+(sqrt(n)/2.25)))
if(!legends & is.null(rlab))  ylimw <- c(0.5, (n+0.5))
if(!legends & !is.null(rlab)) ylimw <- c(0.5, (n+0.75))
par(lend='square')
tud <- (diff(ylimw)/20)/n^0.5



################################################################################
################################################################################
#Walltype 0
if(wall==0){
par(col.axis=tcol)

# Plot
plot(0, 0, type='p', pch='', bg=col, main=main, cex=1, xlim=xlimw, ylim=ylimw, ylab='', xlab='', axes = FALSE, col.main=tcol, col=rgb(1,1,1,0))

if(length(xticks)==1)  ticksade<-pretty(c(xmin, xmax), n = xticks, min.n = n %/% 3)
if(length(xticks)>1 )  ticksade<-xticks

onetick<-ticksade[length(ticksade)]-ticksade[length(ticksade)-1]
lastline<-xmax+(0.15*xrange)
abline(v=lastline, col=bgcol, lwd=1)
ys<-seq(n,1)
yz<-seq(0.5, -0.5, length.out=(N+2))
if(N==1)   yz<-c(0, 0)
if(hlines & N==1)   segments(par('usr')[1], 1:n ,lastline, 1:n, col = bgcol, lty = 1, lwd = 1)
if(hlines & N>1) abline(h=(0:n)+0.5, col=bgcol, lwd=1)
if(legends) lagendram<-legend("topleft", legenlab, pch=c(22,22), col=col, pt.bg=col, horiz=TRUE, bg=rgb(1,1,1, 0), box.col=bgcol, text.col=tcol)

if(!legends) abline(v=v, lty=lty, col=lcol)
if(legends)  segments(v, par('usr')[3] , v, par('usr')[4]-lagendram$rect$h, lty=lty, col=lcol)

for(k in 1:N) {
segments(M[[k]][ ,2], ys+yz[k+1] ,M[[k]][ ,3] , ys+yz[k+1], col = col[k], lty = 1, lwd = 3)
points(M[[k]][ ,1], ys+yz[k+1], type='p', pch=22, cex=1.1, col=col[k], bg=col[k])
}


if(is.null(sectext))   text( rep(lastline+diff(par('usr')[c(1, 2)])/75 , n)  ,seq(n, 1),     labels =vnames,   pos=4 , col=tcol)
if(!is.null(sectext))  text( rep(lastline+diff(par('usr')[c(1, 2)])/75 , n)  ,seq(n, 1)+tud, labels =vnames,   pos=4 , col=tcol)
if(!is.null(sectext))  text( rep(lastline+diff(par('usr')[c(1, 2)])/75 , n)  ,seq(n, 1)-tud, labels =sectext,  pos=4 , col=tcol)

if(logaxe & length(xticks)==1)  logis<- ticksade<-pretty(c(exp(xmin), exp(xmax)), n = xticks, min.n = n %/% 3)
if(logaxe & length(xticks)> 1)  logis<- xticks
if(logaxe)  axis(1, at=log(logis), labels=logis, col=bgcol, col.ticks=bgcol, lwd.ticks=1)
if(!logaxe) axis(1, at=ticksade, col=bgcol, col.ticks=bgcol)

text( lastline+diff(par('usr')[c(1, 2)])/75  , n+0.5+tud*1.5, labels = rlab,     pos=4 , font=2, cex = 1.1, col=tcol)
mtext(xlab, side = 1, at = lastline,  adj = 0, padj=1.5, font=1, col=tcol, cex=par('cex.lab'))
mtext(ylab, line=1.5, side = 2,  font=1, col=tcol, cex=par('cex.lab'))
box(col=bgcol)
}
################################################################################
################################################################################



################################################################################
################################################################################
#Walltype 1
if(wall==1){
par(col.axis=tcol)

# Plot
plot(0, 0, type='p', pch='', bg=col, main=main, cex=1, xlim=xlimw, ylim=ylimw, ylab='', xlab='', axes = FALSE, col.main=tcol, col=rgb(1,1,1,0))
if(length(xticks)==1)  ticksade<-pretty(c(xmin, xmax), n = xticks, min.n = n %/% 3)
if(length(xticks)>1 )  ticksade<-xticks

onetick<-ticksade[length(ticksade)]-ticksade[length(ticksade)-1]
lastline<-xmax+(0.15*xrange)
polygon( c(par('usr')[c(1,1)], lastline, lastline), par('usr')[c(3,4,4,3)], col=bgcol, border=FALSE)
abline(v=v, lty=lty, col=lcol)
if(hlines & N==1)   segments(par('usr')[1], 1:n ,lastline, 1:n, col = rgb(1,1,1), lty = 1, lwd = 1)
if(hlines & N>1)  segments(par('usr')[1], (0:n)+0.5 ,lastline, (0:n)+0.5, col = rgb(1,1,1), lty = 1, lwd = 1)
if(hlines & N>1)  segments(lastline, (0:n)+0.5 ,par('usr')[2], (0:n)+0.5, col = bgcol, lty = 1, lwd = 1)

ys<-seq(n,1)
yz<-seq(0.5, -0.5, length.out=(N+2))
if(N==1)   yz<-c(0, 0)
for(k in 1:N) {
segments(M[[k]][ ,2], ys+yz[k+1] ,M[[k]][ ,3] , ys+yz[k+1], col = a.coladd.ade(col[k], -75), lty = 1, lwd = 3)
segments(M[[k]][ ,2], ys+yz[k+1] ,M[[k]][ ,3] , ys+yz[k+1], col = col[k], lty = 1, lwd = 1)
points(M[[k]][ ,1], ys+yz[k+1], type='p', pch=22, cex=1.1, col=a.coladd.ade(col[k], -75), bg=col[k])
}

if(legends) legend("topleft", legenlab, fill=col, border=a.coladd.ade(col, -75), horiz=TRUE, bg=bgcol, box.col=rgb(1,1,1), box.lwd=2, text.col=tcol)

if(is.null(sectext))   text( rep(lastline+diff(par('usr')[c(1, 2)])/75 , n)  ,seq(n, 1),     labels =vnames,   pos=4 , col=tcol)
if(!is.null(sectext))  text( rep(lastline+diff(par('usr')[c(1, 2)])/75 , n)  ,seq(n, 1)+tud, labels =vnames,   pos=4 , col=tcol)
if(!is.null(sectext))  text( rep(lastline+diff(par('usr')[c(1, 2)])/75 , n)  ,seq(n, 1)-tud, labels =sectext,  pos=4 , col=tcol)

if(logaxe & length(xticks)==1)  logis<- ticksade<-pretty(c(exp(xmin), exp(xmax)), n = xticks, min.n = n %/% 3)
if(logaxe & length(xticks)> 1)  logis<- xticks

if(logaxe)  axis(1, at=log(logis), labels=logis, col=bgcol, col.ticks=tcol, lwd.ticks=1)
if(!logaxe) axis(1, at=ticksade, col=bgcol, col.ticks=tcol)

text( lastline+diff(par('usr')[c(1, 2)])/75  , n+0.5+tud*1.5, labels = rlab,     pos=4 , font=2, cex = 1.1, col=tcol)
mtext(xlab, side = 1, at = lastline,  adj = 0, padj=1.5, font=1, col=tcol, cex=par('cex.lab'))
mtext(ylab, line=1.5, side = 2,  font=1, col=tcol, cex=par('cex.lab'))
box(col=rgb(1,1,1), lwd=1)
}
################################################################################
################################################################################



################################################################################
################################################################################
#Walltype 2
if(wall==2){
par(col.axis=tcol)

# Plot
plot(0, 0, type='p', pch='', bg=col, main=main, cex=1, xlim=xlimw, ylim=ylimw, ylab='', xlab='', axes = FALSE, col.main=tcol, col=rgb(1,1,1,0))
if(length(xticks)==1)  ticksade<-pretty(c(xmin, xmax), n = xticks, min.n = n %/% 3)
if(length(xticks)>1 )  ticksade<-xticks

onetick<-ticksade[length(ticksade)]-ticksade[length(ticksade)-1]
lastline<-xmax+(0.15*xrange)
if(hlines & N==1)   segments(par('usr')[1], 1:n ,lastline, 1:n, col = bgcol, lty = 1, lwd = 1)
if(hlines & N>1)  segments(par('usr')[1], (0:n)+0.5 ,lastline, (0:n)+0.5, col = bgcol, lty = 1, lwd = 1)
if(hlines & N>1)  segments(lastline, (0:n)+0.5 ,par('usr')[2], (0:n)+0.5, col = bgcol, lty = 1, lwd = 1)

abline(v=v, lty=lty, col=lcol)
ys<-seq(n,1)
yz<-seq(0.5, -0.5, length.out=(N+2))
if(N==1)   yz<-c(0, 0)
for(k in 1:N) {
arrows(M[[k]][ ,2], ys+yz[k+1] ,M[[k]][ ,3] , ys+yz[k+1], length = 0.0275, angle = 90, code = 3, col = a.coladd.ade(col[k], -75), lty = 1, lwd = 3)
arrows(M[[k]][ ,2], ys+yz[k+1] ,M[[k]][ ,3] , ys+yz[k+1], length = 0.0275, angle = 90, code = 3, col = col[k], lty = 1, lwd = 1)
points(M[[k]][ ,1], ys+yz[k+1], type='p', pch=22, cex=1.1, col=a.coladd.ade(col[k], -75), bg=col[k])

}

abline(v=lastline, col=a.coladd.ade(bgcol, -75), lwd=1)

if(legends) legend("topleft", legenlab, fill=col, border=a.coladd.ade(col, -75), horiz=TRUE, bg=rgb(1,1,1), box.col=a.coladd.ade(bgcol, -75), box.lwd=1, text.col=tcol, text.width=max(strwidth(legenlab,font = 2)))

if(is.null(sectext))   text( rep(lastline+diff(par('usr')[c(1, 2)])/75 , n)  ,seq(n, 1),     labels =vnames,   pos=4 , col=tcol)
if(!is.null(sectext))  text( rep(lastline+diff(par('usr')[c(1, 2)])/75 , n)  ,seq(n, 1)+tud, labels =vnames,   pos=4 , col=tcol)
if(!is.null(sectext))  text( rep(lastline+diff(par('usr')[c(1, 2)])/75 , n)  ,seq(n, 1)-tud, labels =sectext,  pos=4 , col=tcol)

if(logaxe & length(xticks)==1)  logis<- ticksade<-pretty(c(exp(xmin), exp(xmax)), n = xticks, min.n = n %/% 3)
if(logaxe & length(xticks)> 1)  logis<- xticks
if(logaxe)  axis(1, at=log(logis), labels=logis, col=bgcol, col.ticks=a.coladd.ade(bgcol, -75), lwd.ticks=1)
if(!logaxe) axis(1, at=ticksade, col=bgcol, col.ticks=a.coladd.ade(bgcol, -75))

text( lastline+diff(par('usr')[c(1, 2)])/75  , n+0.5+tud*1.5, labels = rlab,     pos=4 , font=2, cex = 1.1, col=tcol)
mtext(xlab, side = 1, at = lastline,  adj = 0, padj=1.5, font=1, col=tcol, cex=par('cex.lab'))
mtext(ylab, line=1.5, side = 2,  font=1, col=tcol, cex=par('cex.lab'))
box(col=a.coladd.ade(bgcol, -75))
}
################################################################################
################################################################################




################################################################################
################################################################################
#Walltype 3
if(wall==3){
par(col.axis=tcol)

# Plot
plot(0, 0, type='p', pch='', bg=col, main=main, cex=1, xlim=xlimw, ylim=ylimw, ylab='', xlab='', axes = FALSE, col.main=tcol, col=rgb(1,1,1,0))
if(length(xticks)==1)  ticksade<-pretty(c(xmin, xmax), n = xticks, min.n = n %/% 3)
if(length(xticks)>1 )  ticksade<-xticks

onetick<-ticksade[length(ticksade)]-ticksade[length(ticksade)-1]
lastline<-xmax+(0.15*xrange)
polygon( c(par('usr')[c(1,1)], lastline, lastline), par('usr')[c(3,4,4,3)], col=bgcol, border=FALSE)
if(hlines & N==1)   segments(par('usr')[1], 1:n ,lastline, 1:n, col = rgb(1,1,1), lty = 1, lwd = 1)
if(hlines & N>1)    segments(par('usr')[1], (0:n)+0.5 ,lastline, (0:n)+0.5, col = a.coladd.ade(bgcol, -50), lty = 1, lwd = 1)
if(hlines & N>1)    segments(lastline, (0:n)+0.5 ,par('usr')[2], (0:n)+0.5, col = a.coladd.ade(bgcol, -50), lty = 1, lwd = 1)

abline(v=v, lty=lty, col=lcol)
ys<-seq(n,1)
yz<-seq(0.5, -0.5, length.out=(N+2))
if(N==1)   yz<-c(0, 0)
for(k in 1:N) {
arrows(M[[k]][ ,2], ys+yz[k+1] ,M[[k]][ ,3] , ys+yz[k+1], length = 0.0275, angle = 90, code = 3, col = a.coladd.ade(col[k], -75), lty = 1, lwd = 3)
arrows(M[[k]][ ,2], ys+yz[k+1] ,M[[k]][ ,3] , ys+yz[k+1], length = 0.0275, angle = 90, code = 3, col = col[k], lty = 1, lwd = 1)
points(M[[k]][ ,1], ys+yz[k+1], type='p', pch=22, cex=1.1, col=a.coladd.ade(col[k], -75), bg=col[k])

}


abline(v=lastline, col=a.coladd.ade(bgcol, -75), lwd=1)


if(legends) legend("topleft", legenlab, fill=col, border=a.coladd.ade(col, -75), horiz=TRUE, bg=rgb(1,1,1), box.col=a.coladd.ade(bgcol, -75), box.lwd=1, text.col=tcol)

if(is.null(sectext))   text( rep(lastline+diff(par('usr')[c(1, 2)])/75 , n)  ,seq(n, 1),     labels =vnames,   pos=4 , col=tcol)
if(!is.null(sectext))  text( rep(lastline+diff(par('usr')[c(1, 2)])/75 , n)  ,seq(n, 1)+tud, labels =vnames,   pos=4 , col=tcol)
if(!is.null(sectext))  text( rep(lastline+diff(par('usr')[c(1, 2)])/75 , n)  ,seq(n, 1)-tud, labels =sectext,  pos=4 , col=tcol)

if(logaxe & length(xticks)==1)  logis<- ticksade<-pretty(c(exp(xmin), exp(xmax)), n = xticks, min.n = n %/% 3)
if(logaxe & length(xticks)> 1)  logis<- xticks
if(logaxe)  axis(1, at=log(logis), labels=logis, col=bgcol, col.ticks=a.coladd.ade(bgcol, -75), lwd.ticks=1)
if(!logaxe) axis(1, at=ticksade, col=bgcol, col.ticks=a.coladd.ade(bgcol, -75))

text( lastline+diff(par('usr')[c(1, 2)])/75  , n+0.5+tud*1.5, labels = rlab,     pos=4 , font=2, cex = 1.1, col=tcol)
mtext(xlab, side = 1, at = lastline,  adj = 0, padj=1.5, font=1, col=tcol, cex=par('cex.lab'))
mtext(ylab, line=1.5, side = 2,  font=1, col=tcol, cex=par('cex.lab'))
box(col=a.coladd.ade(bgcol, -75))
}
################################################################################
################################################################################


################################################################################
################################################################################
#Walltype 4
if(wall==4){
par(col.axis=tcol)
par(font=2)
# Plot
plot(0, 0, type='p', pch='', bg=col, main=NULL, cex=1, xlim=xlimw, ylim=ylimw, ylab='', xlab='', axes = FALSE, col.main=tcol, col=rgb(1,1,1,0))
if(length(xticks)==1)  ticksade<-pretty(c(xmin, xmax), n = xticks, min.n = n %/% 3)
if(length(xticks)>1 )  ticksade<-xticks

onetick<-ticksade[length(ticksade)]-ticksade[length(ticksade)-1]
lastline<-xmax+(0.15*xrange)
polygon( c(par('usr')[c(1,1)], lastline, lastline), par('usr')[c(3,4,4,3)], col=bgcol, border=FALSE)
polygon( c(lastline, lastline, par('usr')[c(2,2)]), par('usr')[c(3,4,4,3)], col=tcol, border=rgb(1,1,1))
if(hlines & N==1)   segments(par('usr')[1], 1:n ,lastline, 1:n, col = rgb(1,1,1), lty = 1, lwd = 1)
if(hlines & N>1)    segments(par('usr')[1], (0:n)+0.5 ,lastline, (0:n)+0.5, col = rgb(1,1,1), lty = 1, lwd = 1)
if(hlines & N>1)    segments(lastline, (0:n)+0.5 ,par('usr')[2], (0:n)+0.5, col = rgb(1,1,1), lty = 1, lwd = 1)

par(xpd=TRUE)
dx<-7/par('din')[1]
dy<-7/par('din')[2]
xr<-(diff(par('usr')[1:2])/11)*dx
yr<-(diff(par('usr')[3:4])/10)*dy

polygon(a.glc(side=c(2,2,4,4), line=c(0,0,0,0)), a.glc(side=3, line=c(0, 2.75,  2.75, 0)), col=tcol, border=rgb(1,1,1))
if(ylab!='' & ylab!=' ') polygon( a.glc(side=2, line=c(2, 2, 0, 0)), a.glc(side=c(1, 3, 3, 1), line=0), col=bgcol, border=rgb(1,1,1))
text(a.glc(side=0), a.glc(side=3, line=1),    labels=main, cex = 1.25, font=2, col=rgb(1,1,1), adj=c(0.5,0))
text(a.glc(side=2, line=0.75), a.glc(side=5),  labels=ylab, cex = 1.1,  font=2,  col=tcol, adj=c(0.5,0), srt=90)

par(xpd=FALSE)


abline(v=v, lty=lty, col=lcol)
ys<-seq(n,1)
yz<-seq(0.5, -0.5, length.out=(N+2))
if(N==1)   yz<-c(0, 0)
for(k in 1:N) {
arrows(M[[k]][ ,2], ys+yz[k+1] ,M[[k]][ ,3] , ys+yz[k+1], length = 0.0275, angle = 90, code = 3, col = a.coladd.ade(col[k], -75), lty = 1, lwd = 3)
arrows(M[[k]][ ,2], ys+yz[k+1] ,M[[k]][ ,3] , ys+yz[k+1], length = 0.0275, angle = 90, code = 3, col = col[k], lty = 1, lwd = 1)
points(M[[k]][ ,1], ys+yz[k+1], type='p', pch=22, cex=1.1, col=a.coladd.ade(col[k], -75), bg=col[k])

}

if(legends) legend("topleft", legenlab, fill=col, border=rgb(1,1,1), horiz=TRUE, bg=tcol, box.col=rgb(1,1,1), box.lwd=1, text.col=rgb(1,1,1), text.width=max(strwidth(legenlab,font = 2)))

if(is.null(sectext))   text( rep(lastline+diff(par('usr')[c(1, 2)])/75 , n)  ,seq(n, 1),     labels =vnames,   pos=4 , col=rgb(1,1,1))
if(!is.null(sectext))  text( rep(lastline+diff(par('usr')[c(1, 2)])/75 , n)  ,seq(n, 1)+tud, labels =vnames,   pos=4 , col=rgb(1,1,1))
if(!is.null(sectext))  text( rep(lastline+diff(par('usr')[c(1, 2)])/75 , n)  ,seq(n, 1)-tud, labels =sectext,  pos=4 , col=rgb(1,1,1))

if(logaxe & length(xticks)==1)  logis<- ticksade<-pretty(c(exp(xmin), exp(xmax)), n = xticks, min.n = n %/% 3)
if(logaxe & length(xticks)> 1)  logis<- xticks
if(logaxe)  axis(1, at=log(logis), labels=logis, col=bgcol, col.ticks=tcol, lwd.ticks=1)
if(!logaxe) axis(1, at=ticksade, col=bgcol, col.ticks=tcol)

text( lastline+diff(par('usr')[c(1, 2)])/75  , n+0.5+tud*1.5, labels = rlab,     pos=4 , font=2,  col=rgb(1,1,1))
mtext(xlab, side = 1, at = lastline,  adj = 0, padj=1.5, font=1, col=tcol)

box(col=rgb(1,1,1))
}
################################################################################
################################################################################


################################################################################
################################################################################
#Walltype 5
if(wall==5){
par(col.axis=tcol)
par(col.lab=tcol)
par(col.main=tcol)
newmai<-rep(0, 4)
oldmai<-par('mai')

if(oldmai[2]>0.80 & oldmai[2]<=0.82) newmai[2]<- 0.8-oldmai[3]
if(oldmai[3]>0.75 & oldmai[3]<=0.82) newmai[3]<- 0.75-oldmai[3]
if(oldmai[4]>0.25 & oldmai[4]<=0.42) newmai[4]<- 0.25-oldmai[4]
par(mai=(oldmai+newmai))

# Plot
plot(0, 0, type='p', pch='', bg=col, main=NULL, cex=1, xlim=xlimw, ylim=ylimw, ylab='', xlab='', axes = FALSE, col.main=tcol, col=rgb(1,1,1,0))
if(length(xticks)==1)  ticksade<-pretty(c(xmin, xmax), n = xticks, min.n = n %/% 3)
if(length(xticks)>1 )  ticksade<-xticks
onetick<-ticksade[length(ticksade)]-ticksade[length(ticksade)-1]
lastline<-xmax+(0.15*xrange)


par(xpd=TRUE)
dx<-7/par('din')[1]
dy<-7/par('din')[2]
xr<-(diff(par('usr')[1:2])/10)*dx
yr<-(diff(par('usr')[3:4])/10)*dy

polygon(a.glc(side=2, line=c(3.25, 3.25, 0, 0)), a.glc(side=3, line=c(0.6, 3, 3, 0.6)), col=bgcol,        border=tcol)
polygon(a.glc(side=c(2,2,4,4), line=c(0,0,0,0)), a.glc(side=3, line=c(0.6, 3, 3, 0.6)), col=rgb(1,1,1,0), border=tcol)
polygon(a.glc(side=4, line=c(0, 0 ,0.6, 0.6)),   a.glc(side=3, line=c(0.6, 3, 3, 0.6)), col=bgcol,        border=tcol)
polygon(a.glc(side=2, line=c(3.25, 3.25 ,2.65, 2.65)),  a.glc(side=c(1,3,3,1), line=c(2.6, 0.6, 0.6, 2.6)), col=bgcol,  border=tcol)
polygon(a.glc(side=4, line=c(0, 0 ,0.6, 0.6)), a.glc(side=c(1, 3, 3, 1), line=0), col=bgcol, border=tcol)
polygon(a.glc(side=2, line=c(3.25, 3.25, 0, 0)), a.glc(side=1, line=c(2.6, 4.5, 4.5, 2.6)), col=bgcol, border=tcol)
polygon(a.glc(side=c(2, 2, 4, 4), line=0), a.glc(side=1, line=c(2.6, 4.5, 4.5, 2.6)), col=rgb(1,1,1,0), border=tcol)
polygon(a.glc(side=4, line=c(0, 0, 0.6, 0.6)), a.glc(side=1, line=c(2.6, 4.5, 4.5, 2.6)), col=bgcol, border=tcol)
text(a.glc(side=0), a.glc(side=3, line=1.5),  labels=main, cex = 1.25, font=2, col=tcol, adj=c(0.5,0))
text(a.glc(side=0), a.glc(side=1, line=3.75), labels=xlab, cex = 1.1,  font=2, col=tcol, adj=c(0.5,0))
text(a.glc(side=2, line=1.5), a.glc(side=5), labels=ylab, cex = 1.1,   font=2,  col=tcol, adj=c(0.5,0), srt=90)
par(xpd=FALSE)

if(hlines & N==1)   segments(par('usr')[1], 1:n ,lastline, 1:n, col = bgcol, lty = 1, lwd = 1)
if(hlines & N>1)    segments(par('usr')[1], (0:n)+0.5 ,lastline, (0:n)+0.5, col = bgcol, lty = 1, lwd = 1)
if(hlines & N>1)    segments(lastline, (0:n)+0.5 ,par('usr')[2], (0:n)+0.5, col = bgcol, lty = 1, lwd = 1)

abline(v=v, lty=lty, col=lcol)
ys<-seq(n,1)
yz<-seq(0.5, -0.5, length.out=(N+2))
if(N==1)   yz<-c(0, 0)
for(k in 1:N) {
arrows(M[[k]][ ,2], ys+yz[k+1] ,M[[k]][ ,3] , ys+yz[k+1], length = 0.0275, angle = 90, code = 3, col = a.coladd.ade(col[k], -75), lty = 1, lwd = 3)
arrows(M[[k]][ ,2], ys+yz[k+1] ,M[[k]][ ,3] , ys+yz[k+1], length = 0.0275, angle = 90, code = 3, col = col[k], lty = 1, lwd = 1)
points(M[[k]][ ,1], ys+yz[k+1], type='p', pch=22, cex=1.1, col=a.coladd.ade(col[k], -75), bg=col[k])

}


abline(v=lastline, col=tcol, lwd=1)


if(legends) legend("topleft", legenlab, fill=col, border=a.coladd.ade(col, -75), horiz=TRUE, bg=rgb(1,1,1), box.col=tcol, box.lwd=1, text.col=tcol)

if(is.null(sectext))   text( rep(lastline+diff(par('usr')[c(1, 2)])/75 , n)  ,seq(n, 1),     labels =vnames,   pos=4 , col=tcol)
if(!is.null(sectext))  text( rep(lastline+diff(par('usr')[c(1, 2)])/75 , n)  ,seq(n, 1)+tud, labels =vnames,   pos=4 , col=tcol)
if(!is.null(sectext))  text( rep(lastline+diff(par('usr')[c(1, 2)])/75 , n)  ,seq(n, 1)-tud, labels =sectext,  pos=4 , col=tcol)

if(logaxe & length(xticks)==1)  logis<- ticksade<-pretty(c(exp(xmin), exp(xmax)), n = xticks, min.n = n %/% 3)
if(logaxe & length(xticks)> 1)  logis<- xticks
if(logaxe)  axis(1, at=log(logis), labels=logis, col=bgcol, col.ticks=tcol, lwd.ticks=1)
if(!logaxe) axis(1, at=ticksade, col=bgcol, col.ticks=tcol)

text( lastline+diff(par('usr')[c(1, 2)])/75  , n+0.5+tud*1.5, labels = rlab,     pos=4 , font=2, cex = 1.1, col=tcol)


box(col=tcol)
}
################################################################################
################################################################################


################################################################################
################################################################################
#Walltype 6
if(wall==6){
par(col.axis=tcol)

# Plot
plot(0, 0, type='p', pch='', bg=col, main=main, cex=1, xlim=xlimw, ylim=ylimw, ylab='', xlab='', axes = FALSE, col.main=tcol, col=rgb(1,1,1,0))
if(length(xticks)==1)  ticksade<-pretty(c(xmin, xmax), n = xticks, min.n = n %/% 3)
if(length(xticks)>1 )  ticksade<-xticks

onetick<-ticksade[length(ticksade)]-ticksade[length(ticksade)-1]
lastline<-xmax+(0.15*xrange)
polygon( c(par('usr')[c(1,1)], lastline, lastline), par('usr')[c(3,4,4,3)], col=bgcol, border=NA)
abline(v=v, lty=lty, col=lcol)
if(hlines & N==1) segments(par('usr')[1], 1:n ,lastline, 1:n, col = a.coladd.ade(bgcol, -35), lty = 1, lwd = 3)
if(hlines & N==1) segments(par('usr')[1], 1:n ,lastline, 1:n, col = rgb(1,1,1), lty = 1, lwd = 1)
if(hlines & N>1)  segments(par('usr')[1], (0:n)+0.5 ,par('usr')[2], (0:n)+0.5, col = a.coladd.ade(bgcol, -35), lty = 1, lwd = 3)
if(hlines & N>1)  segments(par('usr')[1], (0:n)+0.5 ,par('usr')[2], (0:n)+0.5, col = rgb(1,1,1), lty = 1, lwd = 1)


ys<-seq(n,1)
yz<-seq(0.5, -0.5, length.out=(N+2))
if(N==1)   yz<-c(0, 0)
for(k in 1:N) {
segments(M[[k]][ ,2], ys+yz[k+1] ,M[[k]][ ,3] , ys+yz[k+1], col = a.coladd.ade(col[k], -75), lty = 1, lwd = 3)
segments(M[[k]][ ,2], ys+yz[k+1] ,M[[k]][ ,3] , ys+yz[k+1], col = col[k], lty = 1, lwd = 1)
points(M[[k]][ ,1], ys+yz[k+1], type='p', pch=22, cex=1.1, col=a.coladd.ade(col[k], -75), bg=col[k])
}

if(legends) legend("topleft", legenlab, fill=col, border=a.coladd.ade(col, -75), horiz=TRUE, bg=bgcol,        box.col=rgb(1,1,1), box.lwd=3, text.col=tcol)
if(legends) legend("topleft", legenlab, fill=col, border=a.coladd.ade(col, -75), horiz=TRUE, bg=rgb(1,1,1,0), box.col=a.coladd.ade(bgcol, -35), box.lwd=1, text.col=tcol)


if(is.null(sectext))   text( rep(lastline+diff(par('usr')[c(1, 2)])/75 , n)  ,seq(n, 1),     labels =vnames,   pos=4 , col=tcol)
if(!is.null(sectext))  text( rep(lastline+diff(par('usr')[c(1, 2)])/75 , n)  ,seq(n, 1)+tud, labels =vnames,   pos=4 , col=tcol)
if(!is.null(sectext))  text( rep(lastline+diff(par('usr')[c(1, 2)])/75 , n)  ,seq(n, 1)-tud, labels =sectext,  pos=4 , col=tcol)

if(logaxe & length(xticks)==1)  logis<- ticksade<-pretty(c(exp(xmin), exp(xmax)), n = xticks, min.n = n %/% 3)
if(logaxe & length(xticks)> 1)  logis<- xticks

if(logaxe)  axis(1, at=log(logis), labels=logis, col=bgcol, col.ticks=a.coladd.ade(bgcol, -35), lwd.ticks=3)
if(logaxe)  axis(1, at=log(logis), labels=logis, col=bgcol, col.ticks=rgb(1,1,1), lwd.ticks=1)

if(!logaxe) axis(1, at=ticksade, col=bgcol, col.ticks=a.coladd.ade(bgcol, -35), lwd.ticks=3)
if(!logaxe) axis(1, at=ticksade, col=bgcol, col.ticks=rgb(1,1,1), lwd.ticks=1)


text( lastline+diff(par('usr')[c(1, 2)])/75  , n+0.5+tud*1.5, labels = rlab,     pos=4 , font=2, cex = 1.1, col=tcol)
mtext(xlab, side = 1, at = lastline,  adj = 0, padj=1.5, font=1, col=tcol, cex=par('cex.lab'))
mtext(ylab, line=1.5, side = 2,  font=1, col=tcol, cex=par('cex.lab'))
box(col=rgb(1,1,1), lwd=3)
box(col=a.coladd.ade(bgcol, -35), lwd=1)
}
################################################################################
################################################################################


}
