bar.plot.wtd <-
function(x,  y=NULL, z=NULL, w=NULL, data=NULL, vnames.x=NULL, vnames.y=NULL, vnames.z=NULL, btext=NULL, cutz=F, zperc=NULL, b=NULL, b2=0.5, v=NULL, h=NULL, gradient=FALSE, xlab='', ylab='', main='',  ylim=NULL, yticks=NULL, col=NULL, tcol=NULL,  bgcol=NULL, lcol=NULL, alpha=NULL, beside=TRUE, legendon='topright', wall=0, lhoriz=NULL, prozent=FALSE, ploc=0,  form='r', border=TRUE, density = NULL, angle = NULL, density2 = NULL, angle2 = NULL, fill=NULL, lwd=1, lty=1,  blwd=1, blty=1){
if(any(par('mfg')!=c(1,1,1,1)) & any(par('mai') < c(1.02, 0.82, 0.82, 0.42))){
maidiff<-rep(0, 4)
norm<-c(1.02, 0.82, 0.82, 0.42)
maidiff[par('mai')<norm]<-  norm[par('mai')<norm] - par('mai')[par('mai')<norm]
par(mai=par('mai')+maidiff)
}
oldpar<-par(no.readonly =TRUE)
oldpar<-oldpar[-which(names(oldpar)%in%c('usr',  'plt', 'pin',   'fin', 'fig', 'mfg', 'mfcol', 'mfrow', 'omd', 'omi', 'oma'))]
on.exit(par(oldpar))
xlim=NULL
side<-1
#library(Hmisc)


##############################
if(!is.character(x)){
xt<-deparse(substitute(x))
if(regexpr('~', xt)>=0){
x<-gsub('[~].*$', '', xt)
xpart<-gsub('^.*[~]', '', xt)
if(nchar(gsub('[^+]', '', xpart))==0) y<-gsub('[+].*$', '', xpart)
if(nchar(gsub('[^+]', '', xpart))==1){
z<-gsub('^.*[+]', '', xpart)
y<-gsub('[+].*$', '', gsub(paste(x,'[+]', sep=''), '', xpart))
}
if(nchar(gsub('[^+]', '', xpart))>1)  stop('too many factors!')
}}
##############################



if(!is.null(data) & is.character(x) & length(x)==1){
x<-eval(parse(text=paste("data$",x, sep='')))
if(!is.null(y)) y<-eval(parse(text=paste("data$",y, sep='')))
if(!is.null(z)) z<-eval(parse(text=paste("data$",z, sep='')))
if(!is.null(w)) w<-eval(parse(text=paste("data$",w, sep='')))

if(is.null(z)) cutz<-F
}

################################################################
################################
if(!is.table(x) & length(x)>1 & !is.matrix(x)){
if(is.null(y)  & is.null(z))  a.tab<-table(x)
if(!is.null(y) & is.null(z))  a.tab<-table(x, y)
if(is.null(y) & !is.null(z)){
a.tab<-table(x, z)

if(is.numeric(cutz)){
all.tab<- a.tab
a.tab<- a.tab[,cutz]

if(zperc=='overall')    a.tab<- (a.tab/sum(all.tab))*100
if(zperc=='zells')       a.tab<- (a.tab/rowSums(all.tab))*100

a.tab<- as.table(a.tab)
}
}


if(!is.null(y) & !is.null(z)){
a.tab<-table(x, y, z)
if(is.numeric(cutz)){
all.tab<- a.tab
a.tab<- a.tab[,,cutz]
sum.tab<- NULL
for(k in 1:nlevels(z)){
if(!is.null(sum.tab)) sum.tab<- sum.tab + all.tab[,,k]
if(is.null(sum.tab))  sum.tab<- all.tab[,,k]
}


if(zperc=='overall')    a.tab<- (a.tab/sum(all.tab))*100
if(zperc=='rows')       a.tab<- (a.tab/rowSums(sum.tab))*100
if(zperc=='cols')       a.tab<- t(t(a.tab)/colSums(sum.tab))*100
if(zperc=='zells')      a.tab<- (a.tab/sum.tab)*100
}

}
##############

# Gewichtung

if(!is.null(w) & is.numeric(w)){

if(is.null(y)  &  is.null(z)){
a.tab <- as.table(tapply(w, INDEX=list(x), FUN=sum, na.rm=T, simplify = T))
}

if(is.null(y) & !is.null(z)){
a.tab <- as.table(tapply(w, INDEX=list(x, z), FUN=sum, na.rm=T, simplify = T))
if(is.numeric(cutz)){
all.tab<- a.tab
a.tab<- a.tab[,cutz]

if(zperc=='overall')    a.tab<- (a.tab/sum(all.tab))*100
if(zperc=='zells')       a.tab<- (a.tab/rowSums(all.tab))*100

a.tab<- as.table(a.tab)
}

}


if(!is.null(y) & is.null(z)){
a.tab <- as.table(tapply(w, INDEX=list(x, y), FUN=sum, na.rm=T, simplify = T))
}

if(!is.null(y) & !is.null(z)){
a.tab <- as.table(tapply(w, INDEX=list(x, y, z), FUN=sum, na.rm=T, simplify = T))
if(is.numeric(cutz)){
all.tab<- a.tab
a.tab<- a.tab[,,cutz]
sum.tab<- NULL
for(k in 1:nlevels(z)){
if(!is.null(sum.tab)) sum.tab<- sum.tab + all.tab[,,k]
if(is.null(sum.tab))  sum.tab<- all.tab[,,k]
}


if(zperc=='overall') a.tab<- (a.tab/sum(all.tab))*100
if(zperc=='rows')    a.tab<- (a.tab/rowSums(sum.tab))*100
if(zperc=='cols')    a.tab<- t(t(a.tab)/colSums(sum.tab))*100
if(zperc=='zells')    a.tab<- (a.tab/sum.tab)*100
}
}
}
}
################################
################################################################

if(is.table(x)) a.tab<-x


a.n1<- dim(a.tab)[1]
a.n2<- dim(a.tab)[2]
a.n3<- dim(a.tab)[3]



if(is.null(vnames.x)) vnames1<- unlist(dimnames(a.tab)[1])
if(is.null(vnames.y)) vnames2<- unlist(dimnames(a.tab)[2])
if(is.null(vnames.z)) vnames3<- unlist(dimnames(a.tab)[3])

if(!is.null(vnames.x)) vnames1<- vnames.x
if(!is.null(vnames.y)) vnames2<- vnames.y
if(!is.null(vnames.z)) vnames3<- vnames.z




if(form=='z' | form=='Zylinder' | form=='zylinder' | form=='Zyl' | form=='zyl'){
oldgradient<-gradient
gradient<-'x'
form<-'z'
}
if(form=='z'){
if(length(dim(a.tab))==3){
gradient <- oldgradient
form<-'r'
}
}

legendon2<-NULL
if(length(legendon)==2){
legendon2<-legendon[2]
legendon<-legendon[1]
}

if(legendon=='none')  legendon2<-'none'

if(is.null(legendon2)){
legendon2<-'topright'
if(legendon=='topleft')   legendon2<-'topright'
if(legendon=='topright')  legendon2<-'topleft'
}




###############
# Colors
if(is.null(col) & length(dim(a.tab))<3)  col<-a.getcol.ade(a.n1, type='p')
if(is.null(col) & beside  & length(dim(a.tab))==3)  col<-a.getcol.ade(a.n1, type='p')
if(is.null(col) & !beside & length(dim(a.tab))==3)  col<-a.getcol.ade(a.n3, type='p')
col2<-col
if(is.null(tcol)  & wall==0)   tcol<-1
if(is.null(tcol)  & wall!=0)   tcol<-rgb(0.1,0.1,0.25)
if(is.null(bgcol) & wall==0)   bgcol<-1
if(is.null(bgcol) & wall!=0)   bgcol<-'#DBE0E8'
if(is.null(lcol) & (wall==0 | wall==2 | wall==5))  lcol<-bgcol
if(is.null(lcol) & (wall==1 | wall==6 | wall==4))  lcol<-rgb(1,1,1)

if(length(density)==1 & length(angle)>1) density<-rep(density, length(angle))
if(length(density2)==1 & length(angle2)>1) density<-rep(density2, length(angle2))
if(length(angle)==1 & length(density)>1) angle<-rep(angle, length(density))
if(length(angle2)==1 & length(density2)>1) angle2<-rep(angle2, length(density2))

###############



################################
# set defaults
if(length(dim(a.tab))==1){
if(is.null(lhoriz)) lhoriz<-T
}
if(length(dim(a.tab))==2){
if(is.null(lhoriz)) lhoriz<-T
}
if(length(dim(a.tab))==3){
if(is.null(alpha)) alpha<-0.25
fill<-col
if(form=='c') col<-a.alpha.ade(col, alpha)
if(form=='r') col<-a.coladd.ade(col, -100)
if(is.null(lhoriz)) lhoriz<-F
if(is.null(fill)){
}

if(is.null(density) & is.null(angle) & is.null(density2) & is.null(angle2) & beside){



if(a.n3<=4){
density  <- c(0, 20, 20, 20)
angle    <- c(0, 45, 45, 0)
density2 <- c(0, 0,  20, 20)
angle2   <- c(0, 0,  135, 90)
}

if(a.n3>4){
density  <- c(0, 20, 20,  20,  20,  20, 20,  20)
angle    <- c(0, 45,  0, 135,  45,   0, 25,  70)
density2 <- c(0,  0,  0,   0,  20,  20, 20,  20)
angle2   <- c(0,  0,  0,   0,  135, 90, 155, 110)
}

}

if(is.null(density) & is.null(angle) & is.null(density2) & is.null(angle2) & !beside){
density  <- c(20, 20, 20, 20, 20)
angle    <- c(45, 135, 0, 45, 0)
density2 <- c(0,  0,   0, 20, 20)
angle2   <- c(0,  0,   0, 135, 90)
}

}

if(beside){
if(is.null(b)) b<-1
}


if(!beside){
if(is.null(b)) b<-0.5
}

if(wall==5){
newmai<-rep(0, 4)
oldmai<-par('mai')
#if(oldmai[1]<1.2) newmai[1]<- 1.2 - oldmai[1]
if(oldmai[2]<1) newmai[2]<- 1 - oldmai[2]
if(oldmai[3]>0.75 & oldmai[3]<=0.82) newmai[3]<- 0.75-oldmai[3]
if(oldmai[4]>0.25 & oldmai[4]<=0.42) newmai[4]<- 0.25-oldmai[4]
par(mai=(oldmai+newmai))
}



################################################################################
if(length(dim(a.tab))==1){
if(beside){
xlim <- c(0.5, dim(a.tab)+0.5)
if(is.null(ylim)) ylim <- c(0, max(a.tab, na.rm = TRUE))
ylim[2]<-ylim[2]+diff(ylim)/20

if(side==2){
xlim2<-xlim
xlim<-ylim
ylim<-xlim2
}

xran<-1:dim(a.tab)
if(!is.null(btext))  ylim[2]<-ylim[2]+(diff(ylim)/10)
plot(0,0, xlim=xlim, ylim=ylim, col=rgb(0,0,0,0), xlab='', ylab='', main='', axes=FALSE)
}

################################################################################

if(!beside){
xlim <- c(0.5, 1.5)
if(is.null(ylim)) ylim <- c(0, sum(a.tab, na.rm = TRUE))
ylim[2]<-ylim[2]+diff(ylim)/20

if(side==2){
xlim2<-xlim
xlim<-ylim
ylim<-xlim2
}

xran<-1
if(!is.null(btext))  ylim[2]<-ylim[2]+(diff(ylim)/10)
plot(0,0, xlim=xlim, ylim=ylim, col=rgb(0,0,0,0), xlab='', ylab='', main='', axes=FALSE)
}


}
################################################################################



################################################################################
if(length(dim(a.tab))==2){

if(beside){
xlim <- c(0.5, (ncol(a.tab)+0.5))
if(is.null(ylim)) ylim <- c(0, max(a.tab, na.rm = TRUE))
ylim[2]<-ylim[2]+diff(ylim)/10

if(side==2){
xlim2<-xlim
xlim<-ylim
ylim<-xlim2
}

xran<-1:ncol(a.tab)
if(!is.null(btext))  ylim[2]<-ylim[2]+(diff(ylim)/10)
plot(0,0, xlim=xlim, ylim=ylim, col=rgb(0,0,0,0), xlab='', ylab='', main='', axes=FALSE)

}



if(!beside){
xlim <- c(0.5, (ncol(a.tab)+0.5))
if(is.null(ylim)) ylim <- c(0, max(colSums(a.tab, na.rm = TRUE, dims = 1), na.rm=TRUE))
ylim[2]<-ylim[2]+diff(ylim)/10

if(side==2){
xlim2<-xlim
xlim<-ylim
ylim<-xlim2
}

xran<-1:ncol(a.tab)
if(!is.null(btext))  ylim[2]<-ylim[2]+(diff(ylim)/10)
plot(0,0, xlim=xlim, ylim=ylim, col=rgb(0,0,0,0), xlab='', ylab='', main='', axes=FALSE)
}


}
################################################################################





################################################################################
if(length(dim(a.tab))==3){

if(beside){

xran<-1:ncol(a.tab)
xlim <- c(0.5, (ncol(a.tab)+0.5))
maxy<-0
for(i in 1:dim(a.tab)[3]){
maxy<- maxy+a.tab[, ,i]
}

if(is.null(ylim)) ylim <- c(0, max(maxy, na.rm = TRUE))
ylim[2]<-ylim[2]+diff(ylim)/10

if(side==2){
xlim2<-xlim
xlim<-ylim
ylim<-xlim2
}

plot(0,0, xlim=xlim, ylim=ylim, col=rgb(0,0,0,0), xlab='', ylab='', main='', axes=FALSE)

}


if(!beside){

xran<-1:ncol(a.tab)
xlim <- c(0.5, (ncol(a.tab)+0.5))
maxy<-0
for(i in 1:dim(a.tab)[3]){
maxy<- maxy+a.tab[, ,i]
}

if(is.null(ylim)) ylim <- c(0, max(maxy, na.rm = TRUE))
ylim[2]<-ylim[2]+diff(ylim)/10

if(side==2){
xlim2<-xlim
xlim<-ylim
ylim<-xlim2
}

plot(0,0, xlim=xlim, ylim=ylim, col=rgb(0,0,0,0), xlab='', ylab='', main='', axes=FALSE)

}



}
################################################################################



################################################################################
#Walltype 0
if(wall==0){
par(col.axis=tcol)

bcol<-bgcol
cumat<-(a.tab/2)
for(i in 2:a.n1){
if(length(dim(a.tab))==1) cumat[i]<-(a.tab[i]/2)+cumat[i-1] + a.tab[i-1]/2

if(length(dim(a.tab))==2) cumat[i, ]<-(a.tab[i, ]/2)+cumat[i-1, ] + a.tab[i-1, ]/2
}

if(!is.na(a.n2))  axis(side, at=xran, labels = vnames2, tick=FALSE, col=bgcol, lwd.ticks=1, lwd=0)
if(!is.na(a.n3) & !beside) axis(side, at=xran, labels = vnames2, tick=FALSE, col=bgcol, lwd.ticks=1, lwd=0)
if(is.na(a.n2) &  beside)  axis(side, at=xran, labels = vnames1, tick=FALSE, col=bgcol, lwd.ticks=1, lwd=0)
if(is.na(a.n2) & is.na(a.n3) & !beside & side==1) axis(4, at=cumat, labels = vnames1, tick=FALSE, col=bgcol, lwd.ticks=1, lwd=0, line=-1)
if(is.na(a.n2) & is.na(a.n3) & !beside & side==2) axis(3, at=cumat, labels = vnames1, tick=FALSE, col=bgcol, lwd.ticks=1, lwd=0, line=-0.8)



if(!is.null(yticks)) axis((3-side), at=pretty(ylim, n=yticks), labels =TRUE,  tick=TRUE, col=bgcol, lwd.ticks=1, lwd=0)
if(is.null(yticks))  axis((3-side), col=bgcol, lwd.ticks=1, lwd=0)

abline(v=v, h=h, col=lcol, lwd=1, lty=lty)

text(x=a.glc(0), y=a.glc(3, 2),   labels=main, cex=par('cex.main'), col=tcol, adj=c(0.5, 0.5), xpd=TRUE, font=par("font.main"))
text(x=a.glc(0), y=a.glc(1, 3.5), labels=xlab, cex=par('cex.lab'),  col=tcol, adj=c(0.5, 0.5), xpd=TRUE, font=par("font.lab"))
text(x=a.glc(2, 3), y=a.glc(5), labels=ylab, cex=par('cex.lab'),  col=tcol, adj=c(0.5, 0.5), xpd=TRUE, font=par("font.lab"), srt=90)



if(!is.na(a.n2) & beside & legendon!='none')  legend(legendon,  legend=vnames1, border=bgcol, fill = col2, box.lwd=1, box.col=bgcol, text.col=tcol, bg=rgb(1,1,1, 0), horiz=lhoriz, text.width=max(strwidth(vnames1,font = 2)))
if(!is.na(a.n3) & beside & legendon2!='none')  legend(legendon2, legend=vnames3, border=bgcol, col = col2,  angle = angle,  density=density,  box.lwd=1, box.col=bgcol, horiz=lhoriz, text.col=rgb(1,1,1,0), bg=rgb(1,1,1,0), text.width=max(strwidth(vnames3,font = 2)))
if(!is.na(a.n3) & beside & legendon2!='none')  legend(legendon2, legend=vnames3, border=bgcol, col = col2,  angle = angle2, density=density2, box.lwd=1, box.col=bgcol, horiz=lhoriz, text.col=tcol, bg=rgb(1,1,1, 0), text.width=max(strwidth(vnames3,font = 2)))
if(!is.na(a.n2) & is.na(a.n3) & !beside & legendon!='none') legend(legendon,  legend=vnames1, border=bgcol, fill = col2, box.lwd=1, box.col=bgcol, text.col=tcol, bg=rgb(1,1,1, 0), horiz=lhoriz, text.width=max(strwidth(vnames1,font = 2)))
if(!is.na(a.n3) & !beside & legendon!='none')               legend(legendon,  legend=vnames3, border=bgcol, fill = col2, box.lwd=1, box.col=bgcol, text.col=tcol, bg=rgb(1,1,1, 0), horiz=lhoriz, text.width=max(strwidth(vnames3,font = 2)))
if(!is.na(a.n3) & !beside & legendon2!='none') legend(legendon2, legend=vnames1, border=bgcol, col = col2,  angle = angle,  density=density,  box.lwd=1, box.col=bgcol, horiz=lhoriz, text.col=rgb(1,1,1,0), bg=rgb(1,1,1, 0), text.width=max(strwidth(vnames1,font = 2)))
if(!is.na(a.n3) & !beside & legendon2!='none') legend(legendon2, legend=vnames1, border=bgcol, col = col2,  angle = angle2, density=density2, box.lwd=1, box.col=bgcol, horiz=lhoriz, text.col=tcol, bg=rgb(1,1,1, 0), text.width=max(strwidth(vnames1,font = 2)))


box(col=bgcol, lwd=1)
}
################################################################################



################################################################################
#Walltype 1
if(wall==1){
par(col.axis=tcol)
polygon( par('usr')[c(1,1,2,2)], par('usr')[c(3,4,4,3)], col=bgcol)
bcol<-tcol

cumat<-(a.tab/2)
for(i in 2:a.n1){
if(length(dim(a.tab))==1) cumat[i]<-(a.tab[i]/2)+cumat[i-1] + a.tab[i-1]/2

if(length(dim(a.tab))==2) cumat[i, ]<-(a.tab[i, ]/2)+cumat[i-1, ] + a.tab[i-1, ]/2
}

if(!is.na(a.n2))  xa<-axis(side, at=xran, labels = vnames2, tick=FALSE, col=bgcol, lwd.ticks=1, lwd=0)
if(!is.na(a.n3) & !beside) xa<-axis(side, at=xran, labels = vnames2, tick=FALSE, col=bgcol, lwd.ticks=1, lwd=0)
if(is.na(a.n2) &  beside)  xa<-axis(side, at=xran, labels = vnames1, tick=FALSE, col=bgcol, lwd.ticks=1, lwd=0)
if(is.na(a.n2) & is.na(a.n3) & !beside & side==1) axis(4, at=cumat, labels = vnames1, tick=FALSE, col=bgcol, lwd.ticks=1, lwd=0, line=-1)
if(is.na(a.n2) & is.na(a.n3) & !beside & side==2) axis(3, at=cumat, labels = vnames1, tick=FALSE, col=bgcol, lwd.ticks=1, lwd=0, line=-0.8)


if(!is.null(yticks)) ya<-axis((3-side), at=pretty(ylim, n=yticks), labels =TRUE,  tick=TRUE, col=bgcol, lwd.ticks=1, lwd=0)
if(is.null(yticks))  ya<-axis((3-side), col=bgcol, lwd.ticks=1, lwd=0)

abline(v=c((xran-0.5), xran[length(xran)]+0.5), h=ya, col=rgb(1,1,1), lwd=1, lty=1)

abline(v=v, h=h, col=lcol, lwd=1, lty=lty)

text(x=a.glc(0), y=a.glc(3, 2),   labels=main, cex=par('cex.main'), col=tcol, adj=c(0.5, 0.5), xpd=TRUE, font=par("font.main"))
text(x=a.glc(0), y=a.glc(1, 3.5), labels=xlab, cex=par('cex.lab'),  col=tcol, adj=c(0.5, 0.5), xpd=TRUE, font=par("font.lab"))
text(x=a.glc(2, 3), y=a.glc(5), labels=ylab, cex=par('cex.lab'),  col=tcol, adj=c(0.5, 0.5), xpd=TRUE, font=par("font.lab"), srt=90)



if(!is.na(a.n2) & beside & legendon!='none')  legend(legendon,  legend=vnames1, border=tcol, fill = col2, box.lwd=2, box.col=rgb(1,1,1), text.col=tcol, bg=bgcol, horiz=lhoriz, text.width=max(strwidth(vnames1,font = 2)))
if(!is.na(a.n2) & is.na(a.n3) & !beside & legendon!='none') legend(legendon,  legend=vnames1, border=tcol, fill = col2, box.lwd=2, box.col=rgb(1,1,1), text.col=tcol, bg=bgcol, horiz=lhoriz, text.width=max(strwidth(vnames1,font = 2)))
if(!is.na(a.n3) & !beside & legendon!='none') legend(legendon,  legend=vnames3, border=tcol, fill = col2, box.lwd=2, box.col=rgb(1,1,1), text.col=tcol, bg=bgcol, horiz=lhoriz, text.width=max(strwidth(vnames3,font = 2)))
if(!is.na(a.n3) & beside & legendon2!='none')  legend(legendon2, legend=vnames3, border=tcol, col = col2,  angle = angle,  density=density,  box.lwd=2, box.col=rgb(1,1,1), horiz=lhoriz, text.col=rgb(1,1,1,0), bg=bgcol, text.width=max(strwidth(vnames3,font = 2)))
if(!is.na(a.n3) & beside & legendon2!='none')  legend(legendon2, legend=vnames3, border=tcol, col = col2,  angle = angle2, density=density2, box.lwd=2, box.col=rgb(1,1,1), horiz=lhoriz, text.col=tcol, bg=rgb(1,1,1, 0), text.width=max(strwidth(vnames3,font = 2)))
if(!is.na(a.n3) & !beside & legendon2!='none') legend(legendon2, legend=vnames1, border=tcol, col = col2,  angle = angle,  density=density,  box.lwd=2, box.col=rgb(1,1,1), horiz=lhoriz, text.col=rgb(1,1,1,0), bg=bgcol, text.width=max(strwidth(vnames1,font = 2)))
if(!is.na(a.n3) & !beside & legendon2!='none') legend(legendon2, legend=vnames1, border=tcol, col = col2,  angle = angle2, density=density2, box.lwd=2, box.col=rgb(1,1,1), horiz=lhoriz, text.col=tcol, bg=rgb(1,1,1, 0), text.width=max(strwidth(vnames1,font = 2)))


box(col=rgb(1,1,1), lwd=2)
}
################################################################################




################################################################################
#Walltype 2
if(wall==2){
par(col.axis=tcol)
bcol<-a.coladd.ade(bgcol, -75)

cumat<-(a.tab/2)
for(i in 2:a.n1){
if(length(dim(a.tab))==1) cumat[i]<-(a.tab[i]/2)+cumat[i-1] + a.tab[i-1]/2

if(length(dim(a.tab))==2) cumat[i, ]<-(a.tab[i, ]/2)+cumat[i-1, ] + a.tab[i-1, ]/2
}

if(!is.na(a.n2))  xa<-axis(side, at=xran, labels = vnames2, tick=FALSE, col=bgcol, lwd.ticks=1, lwd=0)
if(!is.na(a.n3) & !beside) xa<-axis(side, at=xran, labels = vnames2, tick=FALSE, col=bgcol, lwd.ticks=1, lwd=0)
if(is.na(a.n2) &  beside)  xa<-axis(side, at=xran, labels = vnames1, tick=FALSE, col=bgcol, lwd.ticks=1, lwd=0)
if(is.na(a.n2) & is.na(a.n3) & !beside & side==1) axis(4, at=cumat, labels = vnames1, tick=FALSE, col=bgcol, lwd.ticks=1, lwd=0, line=-1)
if(is.na(a.n2) & is.na(a.n3) & !beside & side==2) axis(3, at=cumat, labels = vnames1, tick=FALSE, col=bgcol, lwd.ticks=1, lwd=0, line=-0.8)


if(!is.null(yticks)) ya<-axis((3-side), at=pretty(ylim, n=yticks), labels =TRUE,  tick=TRUE, col=bgcol, lwd.ticks=1, lwd=0)
if(is.null(yticks))  ya<-axis((3-side), col=bgcol, lwd.ticks=1, lwd=0)

abline(v=c((xran-0.5), xran[length(xran)]+0.5), h=ya, col=bgcol, lwd=1, lty=1)


abline(v=v, h=h, col=lcol, lwd=1, lty=lty)

text(x=a.glc(0), y=a.glc(3, 2),   labels=main, cex=par('cex.main'), col=tcol, adj=c(0.5, 0.5), xpd=TRUE, font=par("font.main"))
text(x=a.glc(0), y=a.glc(1, 3.5), labels=xlab, cex=par('cex.lab'),  col=tcol, adj=c(0.5, 0.5), xpd=TRUE, font=par("font.lab"))
text(x=a.glc(2, 3), y=a.glc(5), labels=ylab, cex=par('cex.lab'),  col=tcol, adj=c(0.5, 0.5), xpd=TRUE, font=par("font.lab"), srt=90)


if(!is.na(a.n2) & beside & legendon!='none')  legend(legendon,  legend=vnames1, border=a.coladd.ade(bgcol, -75), fill = col2, box.lwd=1, box.col=a.coladd.ade(bgcol, -75), text.col=tcol, bg=rgb(1,1,1), horiz=lhoriz, text.width=max(strwidth(vnames1,font = 2)))
if(!is.na(a.n2) & is.na(a.n3) & !beside & legendon!='none') legend(legendon,  legend=vnames1, border=a.coladd.ade(bgcol, -75), fill = col2, box.lwd=1, box.col=a.coladd.ade(bgcol, -75), text.col=tcol, bg=rgb(1,1,1), horiz=lhoriz, text.width=max(strwidth(vnames1,font = 2)))
if(!is.na(a.n3) & !beside & legendon!='none') legend(legendon,  legend=vnames3, border=a.coladd.ade(bgcol, -75), fill = col2, box.lwd=1, box.col=a.coladd.ade(bgcol, -75), text.col=tcol, bg=rgb(1,1,1), horiz=lhoriz, text.width=max(strwidth(vnames3,font = 2)))
if(!is.na(a.n3) & beside & legendon2!='none')  legend(legendon2, legend=vnames3, border=a.coladd.ade(bgcol, -75), col = col2,  angle = angle,  density=density,  box.lwd=1, box.col=a.coladd.ade(bgcol, -75), horiz=lhoriz, text.col=rgb(1,1,1,0), bg=rgb(1,1,1), text.width=max(strwidth(vnames3,font = 2)))
if(!is.na(a.n3) & beside & legendon2!='none')  legend(legendon2, legend=vnames3, border=a.coladd.ade(bgcol, -75), col = col2,  angle = angle2, density=density2, box.lwd=1, box.col=a.coladd.ade(bgcol, -75), horiz=lhoriz, text.col=tcol, bg=rgb(1,1,1, 0), text.width=max(strwidth(vnames3,font = 2)))
if(!is.na(a.n3) & !beside & legendon2!='none') legend(legendon2, legend=vnames1, border=a.coladd.ade(bgcol, -75), col = col2,  angle = angle,  density=density,  box.lwd=1, box.col=a.coladd.ade(bgcol, -75), horiz=lhoriz, text.col=rgb(1,1,1,0), bg=rgb(1,1,1), text.width=max(strwidth(vnames1,font = 2)))
if(!is.na(a.n3) & !beside & legendon2!='none') legend(legendon2, legend=vnames1, border=a.coladd.ade(bgcol, -75), col = col2,  angle = angle2, density=density2, box.lwd=1, box.col=a.coladd.ade(bgcol, -75), horiz=lhoriz, text.col=tcol, bg=rgb(1,1,1, 0), text.width=max(strwidth(vnames1,font = 2)))


box(col=a.coladd.ade(bgcol, -75))
}
################################################################################



################################################################################
#Walltype 3
if(wall==3){
par(col.axis=tcol)
polygon( par('usr')[c(1,1,2,2)], par('usr')[c(3,4,4,3)], col=bgcol)

bcol<-a.coladd.ade(bgcol, -50)

cumat<-(a.tab/2)
for(i in 2:a.n1){
if(length(dim(a.tab))==1) cumat[i]<-(a.tab[i]/2)+cumat[i-1] + a.tab[i-1]/2

if(length(dim(a.tab))==2) cumat[i, ]<-(a.tab[i, ]/2)+cumat[i-1, ] + a.tab[i-1, ]/2
}

if(!is.na(a.n2) )  xa<-axis(side, at=xran, labels = vnames2, tick=FALSE, col=bgcol, lwd.ticks=1, lwd=0)
if(!is.na(a.n3) & !beside) xa<-axis(side, at=xran, labels = vnames2, tick=FALSE, col=bgcol, lwd.ticks=1, lwd=0)
if(is.na(a.n2) &  beside)  xa<-axis(side, at=xran, labels = vnames1, tick=FALSE, col=bgcol, lwd.ticks=1, lwd=0)
if(is.na(a.n2) & is.na(a.n3) & !beside & side==1) axis(4, at=cumat, labels = vnames1, tick=FALSE, col=bgcol, lwd.ticks=1, lwd=0, line=-1)
if(is.na(a.n2) & is.na(a.n3) & !beside & side==2) axis(3, at=cumat, labels = vnames1, tick=FALSE, col=bgcol, lwd.ticks=1, lwd=0, line=-0.8)


if(!is.null(yticks)) ya<-axis((3-side), at=pretty(ylim, n=yticks), labels =TRUE,  tick=TRUE, col=bgcol, lwd.ticks=1, lwd=0)
if(is.null(yticks))  ya<-axis((3-side), col=bgcol, lwd.ticks=1, lwd=0)

abline(v=c((xran-0.5), xran[length(xran)]+0.5), h=ya, col=a.coladd.ade(bgcol, -50), lwd=1, lty=1)

abline(v=v, h=h, col=lcol, lwd=1, lty=lty)

text(x=a.glc(0), y=a.glc(3, 2),   labels=main, cex=par('cex.main'), col=tcol, adj=c(0.5, 0.5), xpd=TRUE, font=par("font.main"))
text(x=a.glc(0), y=a.glc(1, 3.5), labels=xlab, cex=par('cex.lab'),  col=tcol, adj=c(0.5, 0.5), xpd=TRUE, font=par("font.lab"))
text(x=a.glc(2, 3), y=a.glc(5), labels=ylab, cex=par('cex.lab'),  col=tcol, adj=c(0.5, 0.5), xpd=TRUE, font=par("font.lab"), srt=90)

if(!is.na(a.n2) & beside & legendon!='none')  legend(legendon,  legend=vnames1, border=a.coladd.ade(bgcol, -75), fill = col2, box.lwd=1, box.col=a.coladd.ade(bgcol, -75), text.col=tcol, bg=rgb(1,1,1), horiz=lhoriz, text.width=max(strwidth(vnames1,font = 2)))
if(!is.na(a.n2) & is.na(a.n3) & !beside & legendon!='none') legend(legendon,  legend=vnames1, border=a.coladd.ade(bgcol, -75), fill = col2, box.lwd=1, box.col=a.coladd.ade(bgcol, -75), text.col=tcol, bg=rgb(1,1,1), horiz=lhoriz, text.width=max(strwidth(vnames1,font = 2)))
if(!is.na(a.n3) & !beside & legendon!='none') legend(legendon,  legend=vnames3, border=a.coladd.ade(bgcol, -75), fill = col2, box.lwd=1, box.col=a.coladd.ade(bgcol, -75), text.col=tcol, bg=rgb(1,1,1), horiz=lhoriz, text.width=max(strwidth(vnames3,font = 2)))
if(!is.na(a.n3) & beside & legendon2!='none')  legend(legendon2, legend=vnames3, border=a.coladd.ade(bgcol, -75), col = col2,  angle = angle,  density=density,  box.lwd=1, box.col=a.coladd.ade(bgcol, -75), horiz=lhoriz, text.col=rgb(1,1,1,0), bg=rgb(1,1,1), text.width=max(strwidth(vnames3,font = 2)))
if(!is.na(a.n3) & beside & legendon2!='none')  legend(legendon2, legend=vnames3, border=a.coladd.ade(bgcol, -75), col = col2,  angle = angle2, density=density2, box.lwd=1, box.col=a.coladd.ade(bgcol, -75), horiz=lhoriz, text.col=tcol, bg=rgb(1,1,1, 0), text.width=max(strwidth(vnames3,font = 2)))
if(!is.na(a.n3) & !beside & legendon2!='none') legend(legendon2, legend=vnames1, border=a.coladd.ade(bgcol, -75), col = col2,  angle = angle,  density=density,  box.lwd=1, box.col=a.coladd.ade(bgcol, -75), horiz=lhoriz, text.col=rgb(1,1,1,0), bg=rgb(1,1,1), text.width=max(strwidth(vnames1,font = 2)))
if(!is.na(a.n3) & !beside & legendon2!='none') legend(legendon2, legend=vnames1, border=a.coladd.ade(bgcol, -75), col = col2,  angle = angle2, density=density2, box.lwd=1, box.col=a.coladd.ade(bgcol, -75), horiz=lhoriz, text.col=tcol, bg=rgb(1,1,1, 0), text.width=max(strwidth(vnames1,font = 2)))


box(col=a.coladd.ade(bgcol, -75))
}
################################################################################



################################################################################
#Walltype 4
if(wall==4){
par(col.axis=tcol)
par(font=2)
polygon( par('usr')[c(1,1,2,2)], par('usr')[c(3,4,4,3)], col=bgcol)

bcol<-rgb(1,1,1)

cumat<-(a.tab/2)
for(i in 2:a.n1){
if(length(dim(a.tab))==1) cumat[i]<-(a.tab[i]/2)+cumat[i-1] + a.tab[i-1]/2

if(length(dim(a.tab))==2) cumat[i, ]<-(a.tab[i, ]/2)+cumat[i-1, ] + a.tab[i-1, ]/2
}


if(!is.na(a.n2))  xa<-axis(side, at=xran, labels = vnames2, tick=FALSE, col=bgcol, lwd.ticks=1, lwd=0)
if(!is.na(a.n3) & !beside) xa<-axis(side, at=xran, labels = vnames2, tick=FALSE, col=bgcol, lwd.ticks=1, lwd=0)
if(is.na(a.n2) &  beside)  xa<-axis(side, at=xran, labels = vnames1, tick=FALSE, col=bgcol, lwd.ticks=1, lwd=0)
if(is.na(a.n2) & is.na(a.n3) & !beside & side==1) axis(4, at=cumat, labels = vnames1, tick=FALSE, col=bgcol, lwd.ticks=1, lwd=0, line=-1)
if(is.na(a.n2) & is.na(a.n3) & !beside & side==2) axis(3, at=cumat, labels = vnames1, tick=FALSE, col=bgcol, lwd.ticks=1, lwd=0, line=-0.8)


if(!is.null(yticks)) ya<-axis((3-side), at=pretty(ylim, n=yticks), labels =TRUE,  tick=TRUE, col=bgcol, lwd.ticks=1, lwd=0)
if(is.null(yticks))  ya<-axis((3-side), col=bgcol, lwd.ticks=1, lwd=0)

abline(v=c((xran-0.5), xran[length(xran)]+0.5), h=ya, col=rgb(1,1,1), lwd=1, lty=1)


abline(v=v, h=h, col=lcol, lwd=1, lty=lty)
par(xpd=TRUE)
polygon(a.glc(side=c(2,2,4,4), line=c(0,0,0,0)), a.glc(side=3, line=c(0, 2.75,  2.75, 0)), col=tcol, border=rgb(1,1,1))
if(ylab!='' & ylab!=' ') polygon( a.glc(side=2, line=c(3.5, 3.5, 2, 2)), a.glc(side=c(1, 3, 3, 1), line=0), col=bgcol, border=rgb(1,1,1))
if(xlab!='' & xlab!=' ') polygon( a.glc(side=c(2, 2, 4, 4), line=0),     a.glc(side=1, line=c(4, 2.5, 2.5, 4)), col=bgcol, border=rgb(1,1,1))
text(x=a.glc(0), y=a.glc(3, 1.5),   labels=main, cex=par('cex.main'), col=rgb(1,1,1), adj=c(0.5, 0.5), xpd=TRUE, font=par("font.main"))
text(x=a.glc(0), y=a.glc(1, 3.2), labels=xlab, cex=par('cex.lab'),  col=tcol, adj=c(0.5, 0.5), xpd=TRUE, font=2)
text(x=a.glc(2, 2.75), y=a.glc(5), labels=ylab, cex=par('cex.lab'),  col=tcol, adj=c(0.5, 0.5), xpd=TRUE, font=2, srt=90)


par(xpd=FALSE)


if(!is.na(a.n2) & beside & legendon!='none')  legend(legendon,  legend=vnames1, border=rgb(1,1,1), fill = col2, box.lwd=1, box.col=rgb(1,1,1), text.col=rgb(1,1,1), bg=tcol, horiz=lhoriz, text.width=max(strwidth(vnames1,font = 2)))
if(!is.na(a.n2) & is.na(a.n3) & !beside & legendon!='none') legend(legendon,  legend=vnames1, border=rgb(1,1,1), fill = col2, box.lwd=1, box.col=rgb(1,1,1), text.col=rgb(1,1,1), bg=tcol, horiz=lhoriz, text.width=max(strwidth(vnames1,font = 2)))
if(!is.na(a.n3) & !beside & legendon!='none') legend(legendon,  legend=vnames3, border=rgb(1,1,1), fill = col2, box.lwd=1, box.col=rgb(1,1,1), text.col=rgb(1,1,1), bg=tcol, horiz=lhoriz, text.width=max(strwidth(vnames3,font = 2)))
if(!is.na(a.n3) & beside & legendon2!='none')  legend(legendon2, legend=vnames3, border=rgb(1,1,1), col = rgb(1,1,1),  fill = rgb(1,1,1),angle = angle,  density=density,  box.lwd=1, box.col=rgb(1,1,1), horiz=lhoriz, text.col=rgb(1,1,1,0), bg=tcol, text.width=max(strwidth(vnames3,font = 2)))
if(!is.na(a.n3) & beside & legendon2!='none')  legend(legendon2, legend=vnames3, border=rgb(1,1,1), col = rgb(1,1,1),  fill = rgb(1,1,1),angle = angle2, density=density2, box.lwd=1, box.col=rgb(1,1,1), horiz=lhoriz, text.col=rgb(1,1,1),   bg=rgb(1,1,1, 0), text.width=max(strwidth(vnames3,font = 2)))
if(!is.na(a.n3) & !beside & legendon2!='none') legend(legendon2, legend=vnames1, border=rgb(1,1,1), col = rgb(1,1,1),  fill = rgb(1,1,1),angle = angle,  density=density,  box.lwd=1, box.col=rgb(1,1,1), horiz=lhoriz, text.col=rgb(1,1,1,0), bg=tcol, text.width=max(strwidth(vnames1,font = 2)))
if(!is.na(a.n3) & !beside & legendon2!='none') legend(legendon2, legend=vnames1, border=rgb(1,1,1), col = rgb(1,1,1),  fill = rgb(1,1,1),angle = angle2, density=density2, box.lwd=1, box.col=rgb(1,1,1), horiz=lhoriz, text.col=rgb(1,1,1),   bg=rgb(1,1,1, 0), text.width=max(strwidth(vnames1,font = 2)))


box(col=rgb(1,1,1))
}
################################################################################


################################################################################
#Walltype 5
if(wall==5){
par(col.axis=tcol)
par(col.lab=tcol)
par(col.main=tcol)
par(font=2)
bcol<-a.coladd.ade(bgcol, -75)




cumat<-(a.tab/2)
for(i in 2:a.n1){
if(length(dim(a.tab))==1) cumat[i]<-(a.tab[i]/2)+cumat[i-1] + a.tab[i-1]/2
if(length(dim(a.tab))==2) cumat[i, ]<-(a.tab[i, ]/2)+cumat[i-1, ] + a.tab[i-1, ]/2
}


if(!is.na(a.n2))  xa<-axis(side, at=xran, labels = vnames2, tick=FALSE, col=bgcol, lwd.ticks=1, lwd=0)
if(!is.na(a.n3) & !beside) xa<-axis(side, at=xran, labels = vnames2, tick=FALSE, col=bgcol, lwd.ticks=1, lwd=0)
if(is.na(a.n2) &  beside)  xa<-axis(side, at=xran, labels = vnames1, tick=FALSE, col=bgcol, lwd.ticks=1, lwd=0)
if(is.na(a.n2) & is.na(a.n3) & !beside & side==1) axis(4, at=cumat, labels = vnames1, tick=FALSE, col=bgcol, lwd.ticks=1, lwd=0, line=-1)
if(is.na(a.n2) & is.na(a.n3) & !beside & side==2) axis(3, at=cumat, labels = vnames1, tick=FALSE, col=bgcol, lwd.ticks=1, lwd=0, line=-0.8)

if(!is.null(yticks)) ya<-axis((3-side), at=pretty(ylim, n=yticks), labels =TRUE,  tick=TRUE, col=tcol, lwd.ticks=1, lwd=0)
if(is.null(yticks))  ya<-axis((3-side), col=tcol, lwd.ticks=1, lwd=0)


abline(v=v, h=h, col=lcol, lwd=1, lty=lty)

par(xpd=TRUE)
polygon(a.glc(side=2, line=c(4.25, 4.25, 0, 0)), a.glc(side=3, line=c(0.6, 3, 3, 0.6)), col=bgcol,        border=tcol)
polygon(a.glc(side=c(2,2,4,4), line=c(0,0,0,0)), a.glc(side=3, line=c(0.6, 3, 3, 0.6)), col=rgb(1,1,1,0), border=tcol)
polygon(a.glc(side=4, line=c(0, 0 ,0.6, 0.6)),   a.glc(side=3, line=c(0.6, 3, 3, 0.6)), col=bgcol,        border=tcol)
polygon(a.glc(side=2, line=c(4.25, 4.25 ,3.65, 3.65)),  a.glc(side=c(1,3,3,1), line=c(2.6, 0.6, 0.6, 2.6)), col=bgcol,  border=tcol)
if(!(is.na(a.n2) & is.na(a.n3) & !beside)) polygon(a.glc(side=4, line=c(0, 0 ,0.6, 0.6)), a.glc(side=c(1, 3, 3, 1), line=0), col=bgcol, border=tcol)
polygon(a.glc(side=2, line=c(4.25, 4.25, 0, 0)), a.glc(side=1, line=c(2.6, 4.5, 4.5, 2.6)), col=bgcol, border=tcol)
polygon(a.glc(side=c(2, 2, 4, 4), line=0), a.glc(side=1, line=c(2.6, 4.5, 4.5, 2.6)), col=rgb(1,1,1,0), border=tcol)
polygon(a.glc(side=4, line=c(0, 0, 0.6, 0.6)), a.glc(side=1, line=c(2.6, 4.5, 4.5, 2.6)), col=bgcol, border=tcol)
text(a.glc(side=0), a.glc(side=3, line=1.5),  labels=main, cex = 1.25, font=2, col=tcol, adj=c(0.5,0))
text(a.glc(side=0), a.glc(side=1, line=3.75), labels=xlab, cex = 1.1,  font=2, col=tcol, adj=c(0.5,0))
text(a.glc(side=2, line=2.5), a.glc(side=5), labels=ylab, cex = 1.1,   font=2,  col=tcol, adj=c(0.5,0), srt=90)
par(xpd=FALSE)



if(!is.na(a.n2) & beside & legendon!='none')  legend(legendon,  legend=vnames1, border=tcol, fill = col2, box.lwd=1, box.col=tcol, text.col=tcol, bg=rgb(1,1,1), horiz=lhoriz, text.width=max(strwidth(vnames1,font = 2)))
if(!is.na(a.n2) & is.na(a.n3)  & !beside & legendon!='none') legend(legendon,  legend=vnames1, border=tcol, fill = col2, box.lwd=1, box.col=tcol, text.col=tcol, bg=rgb(1,1,1), horiz=lhoriz, text.width=max(strwidth(vnames1,font = 2)))
if(!is.na(a.n3) & !beside & legendon!='none') legend(legendon,  legend=vnames3, border=tcol, fill = col2, box.lwd=1, box.col=tcol, text.col=tcol, bg=rgb(1,1,1), horiz=lhoriz, text.width=max(strwidth(vnames3,font = 2)))
if(!is.na(a.n3) & beside & legendon2!='none')  legend(legendon2, legend=vnames3, border=tcol, col = col2,  angle = angle,  density=density,  box.lwd=1, box.col=tcol, horiz=lhoriz, text.col=rgb(1,1,1,0), bg=rgb(1,1,1), text.width=max(strwidth(vnames3,font = 2)))
if(!is.na(a.n3) & beside & legendon2!='none')  legend(legendon2, legend=vnames3, border=tcol, col = col2,  angle = angle2, density=density2, box.lwd=1, box.col=tcol, horiz=lhoriz, text.col=tcol, bg=rgb(1,1,1, 0), text.width=max(strwidth(vnames3,font = 2)))
if(!is.na(a.n3) & !beside & legendon2!='none') legend(legendon2, legend=vnames1, border=tcol, col = col2,  angle = angle,  density=density,  box.lwd=1, box.col=tcol, horiz=lhoriz, text.col=rgb(1,1,1,0), bg=rgb(1,1,1), text.width=max(strwidth(vnames1,font = 2)))
if(!is.na(a.n3) & !beside & legendon2!='none') legend(legendon2, legend=vnames1, border=tcol, col = col2,  angle = angle2, density=density2, box.lwd=1, box.col=tcol, horiz=lhoriz, text.col=tcol, bg=rgb(1,1,1, 0), text.width=max(strwidth(vnames1,font = 2)))


box(col=tcol)
}
################################################################################


################################################################################
#Walltype 6
if(wall==6){
par(col.axis=tcol)
polygon( par('usr')[c(1,1,2,2)], par('usr')[c(3,4,4,3)], col=bgcol)
bcol<-tcol

cumat<-(a.tab/2)
for(i in 2:a.n1){
if(length(dim(a.tab))==1) cumat[i]<-(a.tab[i]/2)+cumat[i-1] + a.tab[i-1]/2

if(length(dim(a.tab))==2) cumat[i, ]<-(a.tab[i, ]/2)+cumat[i-1, ] + a.tab[i-1, ]/2
}

if(!is.na(a.n2))  xa<-axis(side, at=xran, labels = vnames2, tick=FALSE, col=bgcol, lwd.ticks=3, lwd=0, col.ticks=a.coladd.ade(bgcol, -35))
if(!is.na(a.n2))  xa<-axis(side, at=xran, labels = vnames2, tick=FALSE, col=bgcol, lwd.ticks=1, lwd=0, col.ticks=rgb(1,1,1))
if(!is.na(a.n3) & !beside) xa<-axis(side, at=xran, labels = vnames2, tick=FALSE, col=bgcol, lwd.ticks=3, lwd=0, col.ticks=a.coladd.ade(bgcol, -35))
if(!is.na(a.n3) & !beside) xa<-axis(side, at=xran, labels = vnames2, tick=FALSE, col=bgcol, lwd.ticks=1, lwd=0, col.ticks=rgb(1,1,1))
if(is.na(a.n2) &  beside)  xa<-axis(side, at=xran, labels = vnames1, tick=FALSE, col=bgcol, lwd.ticks=3, lwd=0, col.ticks=a.coladd.ade(bgcol, -35))
if(is.na(a.n2) &  beside)  xa<-axis(side, at=xran, labels = vnames1, tick=FALSE, col=bgcol, lwd.ticks=1, lwd=0, col.ticks=rgb(1,1,1))
if(is.na(a.n2) & is.na(a.n3) & !beside & side==1) axis(4, at=cumat, labels = vnames1, tick=FALSE, col=bgcol, lwd.ticks=3, lwd=0, line=-1, col.ticks=a.coladd.ade(bgcol, -35))
if(is.na(a.n2) & is.na(a.n3) & !beside & side==1) axis(4, at=cumat, labels = vnames1, tick=FALSE, col=bgcol, lwd.ticks=1, lwd=0, line=-1, col.ticks=rgb(1,1,1))
if(is.na(a.n2) & is.na(a.n3) & !beside & side==2) axis(3, at=cumat, labels = vnames1, tick=FALSE, col=bgcol, lwd.ticks=3, lwd=0, line=-0.8, col.ticks=a.coladd.ade(bgcol, -35))
if(is.na(a.n2) & is.na(a.n3) & !beside & side==2) axis(3, at=cumat, labels = vnames1, tick=FALSE, col=bgcol, lwd.ticks=1, lwd=0, line=-0.8, col.ticks=rgb(1,1,1))

if(!is.null(yticks)) ya<-axis((3-side), at=pretty(ylim, n=yticks), labels =TRUE,  tick=TRUE, col=bgcol, lwd.ticks=1, lwd=0)
if(is.null(yticks))  ya<-axis((3-side), col=bgcol, lwd.ticks=1, lwd=0)

abline(v=c((xran-0.5), xran[length(xran)]+0.5), h=ya, col=a.coladd.ade(bgcol, -35), lwd=3, lty=1)
abline(v=c((xran-0.5), xran[length(xran)]+0.5), h=ya, col=rgb(1,1,1), lwd=1, lty=1)

abline(v=v, h=h, col=lcol, lwd=1, lty=lty)

text(x=a.glc(0), y=a.glc(3, 2),   labels=main, cex=par('cex.main'), col=tcol, adj=c(0.5, 0.5), xpd=TRUE, font=par("font.main"))
text(x=a.glc(0), y=a.glc(1, 3.5), labels=xlab, cex=par('cex.lab'),  col=tcol, adj=c(0.5, 0.5), xpd=TRUE, font=par("font.lab"))
text(x=a.glc(2, 3), y=a.glc(5), labels=ylab, cex=par('cex.lab'),  col=tcol, adj=c(0.5, 0.5), xpd=TRUE, font=par("font.lab"), srt=90)


if(!is.na(a.n2) & beside & legendon!='none')                legend(legendon,  legend=vnames1, border=a.coladd.ade(bgcol, -35), fill = col2, box.lwd=3, box.col=rgb(1,1,1), text.col=tcol, bg=bgcol, horiz=lhoriz, text.width=max(strwidth(vnames1,font = 2)))
if(!is.na(a.n2) & beside & legendon!='none')                legend(legendon,  legend=vnames1, border=a.coladd.ade(bgcol, -35), fill = col2, box.lwd=1, box.col=a.coladd.ade(bgcol, -35), text.col=tcol, bg=rgb(0,0,0,0), horiz=lhoriz, text.width=max(strwidth(vnames1,font = 2)))
if(!is.na(a.n2) & is.na(a.n3) & !beside & legendon!='none') legend(legendon,  legend=vnames1, border=a.coladd.ade(bgcol, -35), fill = col2, box.lwd=3, box.col=rgb(1,1,1), text.col=tcol, bg=bgcol, horiz=lhoriz, text.width=max(strwidth(vnames1,font = 2)))
if(!is.na(a.n2) & is.na(a.n3) & !beside & legendon!='none') legend(legendon,  legend=vnames1, border=a.coladd.ade(bgcol, -35), fill = col2, box.lwd=1, box.col=a.coladd.ade(bgcol, -35), text.col=tcol, bg=rgb(0,0,0,0), horiz=lhoriz, text.width=max(strwidth(vnames1,font = 2)))
if(!is.na(a.n3) & !beside & legendon!='none')               legend(legendon,  legend=vnames3, border=a.coladd.ade(bgcol, -35), fill = col2, box.lwd=3, box.col=rgb(1,1,1), text.col=tcol, bg=bgcol, horiz=lhoriz, text.width=max(strwidth(vnames3,font = 2)))
if(!is.na(a.n3) & !beside & legendon!='none')               legend(legendon,  legend=vnames3, border=a.coladd.ade(bgcol, -35), fill = col2, box.lwd=1, box.col=a.coladd.ade(bgcol, -35), text.col=tcol, bg=rgb(0,0,0,0), horiz=lhoriz, text.width=max(strwidth(vnames3,font = 2)))
if(!is.na(a.n3) & beside & legendon2!='none')                legend(legendon2, legend=vnames3, border=a.coladd.ade(bgcol, -35), col = col2,  angle = angle,  density=density,  box.lwd=3, box.col=rgb(1,1,1), horiz=lhoriz, text.col=rgb(1,1,1,0), bg=bgcol, text.width=max(strwidth(vnames3,font = 2)))
if(!is.na(a.n3) & beside & legendon2!='none')                legend(legendon2, legend=vnames3, border=a.coladd.ade(bgcol, -35), col = col2,  angle = angle,  density=density,  box.lwd=1, box.col=a.coladd.ade(bgcol, -35), horiz=lhoriz, text.col=rgb(1,1,1,0), bg=rgb(0,0,0,0), text.width=max(strwidth(vnames3,font = 2)))
if(!is.na(a.n3) & beside & legendon2!='none')                legend(legendon2, legend=vnames3, border=a.coladd.ade(bgcol, -35), col = col2,  angle = angle2, density=density2, box.lwd=3, box.col=rgb(1,1,1), horiz=lhoriz, text.col=tcol, bg=rgb(1,1,1, 0), text.width=max(strwidth(vnames3,font = 2)))
if(!is.na(a.n3) & beside & legendon2!='none')                legend(legendon2, legend=vnames3, border=a.coladd.ade(bgcol, -35), col = col2,  angle = angle2, density=density2, box.lwd=1, box.col=a.coladd.ade(bgcol, -35), horiz=lhoriz, text.col=tcol, bg=rgb(1,1,1, 0), text.width=max(strwidth(vnames3,font = 2)))
if(!is.na(a.n3) & !beside & legendon2!='none')               legend(legendon2, legend=vnames1, border=a.coladd.ade(bgcol, -35), col = col2,  angle = angle,  density=density,  box.lwd=3, box.col=rgb(1,1,1),               horiz=lhoriz, text.col=rgb(1,1,1,0), bg=bgcol, text.width=max(strwidth(vnames1,font = 2)))
if(!is.na(a.n3) & !beside & legendon2!='none')               legend(legendon2, legend=vnames1, border=a.coladd.ade(bgcol, -35), col = col2,  angle = angle,  density=density,  box.lwd=1, box.col=a.coladd.ade(bgcol, -35), horiz=lhoriz, text.col=rgb(1,1,1,0), bg=rgb(0,0,0,0), text.width=max(strwidth(vnames1,font = 2)))
if(!is.na(a.n3) & !beside & legendon2!='none')               legend(legendon2, legend=vnames1, border=a.coladd.ade(bgcol, -35), col = col2,  angle = angle2, density=density2, box.lwd=3, box.col=rgb(1,1,1), horiz=lhoriz, text.col=tcol, bg=rgb(1,1,1, 0), text.width=max(strwidth(vnames1,font = 2)))
if(!is.na(a.n3) & !beside & legendon2!='none')               legend(legendon2, legend=vnames1, border=a.coladd.ade(bgcol, -35), col = col2,  angle = angle2, density=density2, box.lwd=1, box.col=a.coladd.ade(bgcol, -35), horiz=lhoriz, text.col=tcol, bg=rgb(1,1,1, 0), text.width=max(strwidth(vnames1,font = 2)))

box(lwd=3, col=rgb(1,1,1))
box(lwd=1, col=a.coladd.ade(bgcol, -35))
}
################################################################################



################################################################################################################################################################
################################################################################
################################################################################
if(length(dim(a.tab))==1){




if(!is.numeric(cutz) | is.null(zperc))     TB<-paste(round_n.ade(as.vector((a.tab/sum(a.tab))*100) , 1), '%', sep='' )
if(is.numeric(cutz) &  !is.null(zperc)) if(zperc=='overall')   TB<-paste(round_n.ade(as.vector((all.tab[,cutz]/sum(all.tab))*100) , 1), '%', sep='' )
if(is.numeric(cutz) &  !is.null(zperc)) if(zperc=='zells')      TB<-paste(round_n.ade(as.vector((all.tab[,cutz]/rowSums(all.tab))*100) , 1), '%', sep='' )



if(beside){
a.form.ade(x = xran, y = ylim[1], h=a.tab, b = b, b2 = b2, col=col, side=side, nslices=100, gradient=gradient, border=border, bcol=bcol, form=form, density = density, angle = angle, density2 = density2, angle2 = angle2, fill=fill, lwd=lwd, lty=lty,  blwd=blwd, blty=blty, a.tab=a.tab)


################################
if(ploc==0) pfactor<- 0.5
if(ploc==0) padj <-   0.5

if(ploc==1) pfactor<-0
if(ploc==1) padj <- -0.5

if(ploc==2) pfactor<-1
if(ploc==2) padj <- -0.5

if(ploc==3) pfactor<-1
if(ploc==3) padj <-  1.5

if(ploc==4) pfactor<-0
if(ploc==4) padj <-  1.1


if(is.logical(prozent)){ if(prozent){
if(wall!=4) text(x=xran, y=as.vector(a.tab*pfactor), labels=TB, col=tcol, cex=1, font=2,   adj = c(0.5, padj))
if(wall==4) text(x=xran, y=as.vector(a.tab*pfactor), labels=TB, col=rgb(1,1,1), cex=1, font=2,   adj = c(0.5, padj))
}}
if(is.character(prozent) | is.numeric(prozent)){
if(wall!=4) text(x=xran, y=as.vector(a.tab*pfactor), labels=prozent, col=tcol, cex=1, font=2,   adj = c(0.5, padj))
if(wall==4) text(x=xran, y=as.vector(a.tab*pfactor), labels=prozent, col=rgb(1,1,1), cex=1, font=2,   adj = c(0.5, padj))
}

################################


################################
# Text lines
if(!is.null(btext)){
if(is.logical(btext)){
if(!btext)  btext<-NULL
if(btext){
btext<-NULL
btext<- paste('p:', format_p.ade(chisq.test(a.tab)$p.value , 4) )

}
}
if(!is.logical(btext)){
if(wall==0 | wall==5) wcol<-bgcol
if(wall==1) wcol<-rgb(1,1,1)
if(wall==2) wcol<-a.coladd.ade(bgcol, -75)
if(wall==3) wcol<-a.coladd.ade(bgcol, -50)
if(wall==4) wcol<-rgb(1,1,1)
if(wall==6) wcol<-a.coladd.ade(bgcol, -35)
if(wall==0 | wall==2 | wall==5) wlwd<-1
if(wall==1 | wall==3 | wall==4 | wall==6) wlwd<-2
yr<- a.glc(1, 0)-a.glc(1, 1.5)
segments(min(xran),    max(a.tab)+(yr/3) ,   max(xran),   max(a.tab)+(yr/3) , col=wcol, lwd=wlwd)
segments(xran,    max(a.tab)+(yr/3) ,   xran,   (a.tab)+(yr/8) , col=wcol, lwd=wlwd)
text(mean(xran), max(a.tab)+(yr*0.75) , labels=btext, col=tcol)
}
}
################################
}

################################################################################

if(!beside){

if(!is.numeric(cutz) | is.null(zperc))     TB<-paste(round_n.ade(as.vector((a.tab/sum(a.tab))*100) , 1), '%', sep='' )
if(is.numeric(cutz) &  !is.null(zperc)) if(zperc=='overall')   TB<-paste(round_n.ade(as.vector((all.tab[,cutz]/sum(all.tab))*100) , 1), '%', sep='' )
if(is.numeric(cutz) &  !is.null(zperc)) if(zperc=='zells')      TB<-paste(round_n.ade(as.vector((all.tab[,cutz]/rowSums(all.tab))*100) , 1), '%', sep='' )



hoch<-rep(ylim[1], dim(a.tab))
for(k in 1:dim(a.tab)){
if(k> 1) a.form.ade(x = xran, y = hoch[k], h=a.tab[k], b = b, b2 = b2, col=col[k], side=side, nslices=100, gradient=gradient, border=border, bcol=bcol, form=form, density = density, angle = angle, density2 = density2, angle2 = angle2, fill=fill[k], blwd=blwd, blty=blty, lwd=lwd, lty=lty, a.tab=a.tab)
if(k==1) a.form.ade(x = xran, y = hoch[k], h=a.tab[k], b = b, b2 = b2, col=col[k], side=side, nslices=100, gradient=gradient, border=border, bcol=bcol, form=form, density = density, angle = angle, density2 = density2, angle2 = angle2, fill=fill[k], blwd=blwd, blty=blty, lwd=lwd, lty=lty, a.tab=a.tab)
hoch<-hoch+a.tab[k]

################################
if(is.logical(prozent)){ if(prozent){
if(wall!=4) text(x=1, y=cumat[k], labels=TB[k], col=tcol, cex=1, font=2)
if(wall==4) text(x=1, y=cumat[k], labels=TB[k], col=rgb(1,1,1), cex=1, font=2)
}}

if(is.character(prozent) | is.numeric(prozent)){
if(wall!=4) text(x=1, y=cumat[k], labels=prozent[k], col=tcol, cex=1, font=2)
if(wall==4) text(x=1, y=cumat[k], labels=prozent[k], col=rgb(1,1,1), cex=1, font=2)
}
################################

}

################################
# Text lines
if(!is.null(btext)){
if(is.logical(btext)){
if(!btext)  btext<-NULL
if(btext){
btext<-NULL
btext<- paste('p:', format_p.ade(chisq.test(a.tab)$p.value , 4) )
}
}
if(!is.logical(btext)){
yr<- a.glc(1, 0)-a.glc(1, 1.5)

brt=4
if(form=='z') brt=1.5
yr<- a.glc(1, 0)-a.glc(1, 3)
text(xran, hoch+(yr/brt) , labels=btext, col=tcol)
}
}
################################

}


################################################################################
}
################################################################################



################################################################################
if(length(dim(a.tab))==2){


if(is.numeric(cutz) & !is.null(z)){
sum.tab<- NULL
for(k in 1:nlevels(z)){
if(!is.null(sum.tab)) sum.tab<- sum.tab + all.tab[,,k]
if(is.null(sum.tab))  sum.tab<- all.tab[,,k]
}
}




if(!is.numeric(cutz) | is.null(zperc))     TB<-(a.tab/sum(a.tab))*100
if(is.numeric(cutz) &  !is.null(zperc)) if(zperc=='overall')   TB<-(all.tab[,,cutz]/sum(all.tab))*100
if(is.numeric(cutz) &  !is.null(zperc)) if(zperc=='rows')      TB<-(all.tab[,,cutz]/rowSums(sum.tab))*100
if(is.numeric(cutz) &  !is.null(zperc)) if(zperc=='cols')      TB<-t(t(all.tab[,,cutz])/colSums(sum.tab))*100
if(is.numeric(cutz) &  !is.null(zperc)) if(zperc=='zells')     TB<-(all.tab[,,cutz]/sum.tab)*100


if(beside){
b<-(b/nrow(a.tab))*0.75
yrat<-seq(-nrow(a.tab)+1, nrow(a.tab)-1, length.out=nrow(a.tab))*(b/2)
for(k in 1:nrow(a.tab)){
a.form.ade(x = xran+yrat[k], y = ylim[1], h=a.tab[k, ], b = b, b2 = b2, col=col[k], side=side, nslices=100, gradient=gradient, border=border, bcol=bcol, form=form, density = density[k], angle = angle[k], density2 = density2[k], angle2 = angle2[k], blwd=blwd, blty=blwd, lwd=lwd, lty=lwd, fill=fill[k], a.tab=a.tab)

################################
if(ploc==0) pfactor<- 0.5
if(ploc==0) padj <-   0.5

if(ploc==1) pfactor<-0
if(ploc==1) padj <- -0.5

if(ploc==2) pfactor<-1
if(ploc==2) padj <- -0.5

if(ploc==3) pfactor<-1
if(ploc==3) padj <-  1.5

if(ploc==4) pfactor<-0
if(ploc==4) padj <-  1.1

if(is.logical(prozent)){ if(prozent){
if(wall!=4) text(x=xran+yrat[k], y=(a.tab*pfactor)[k,], labels=paste(round_n.ade(TB[k,], 1), '%', sep='' ), col=tcol, cex=1, font=2,   adj = c(0.5, padj))
if(wall==4) text(x=xran+yrat[k], y=(a.tab*pfactor)[k,], labels=paste(round_n.ade(TB[k,], 1), '%', sep='' ), col=rgb(1,1,1), cex=1, font=2,   adj = c(0.5, padj))
}}


if(is.character(prozent) | is.numeric(prozent)){
if(wall!=4) text(x=xran+yrat[k], y=(a.tab*pfactor)[k,], labels=prozent[k,], col=tcol, cex=1, font=2,   adj = c(0.5, padj))
if(wall==4) text(x=xran+yrat[k], y=(a.tab*pfactor)[k,], labels=prozent[k,], col=rgb(1,1,1), cex=1, font=2,   adj = c(0.5, padj))
}

################################
}


################################
# Text lines
if(!is.null(btext)){
if(is.logical(btext)){
if(!btext)  btext<-NULL
if(btext){
btext<-NULL
for(i in 1:ncol(a.tab)){
btext[i]<- paste('p:', format_p.ade(chisq.test(a.tab[ ,i])$p.value , 4) )
}
}
}
if(!is.logical(btext)){

if(wall==0) wcol<-bgcol
if(wall==1) wcol<-rgb(1,1,1)
if(wall==2) wcol<-a.coladd.ade(bgcol, -75)
if(wall==3) wcol<-a.coladd.ade(bgcol, -50)
if(wall==4) wcol<-rgb(1,1,1)
if(wall==5) wcol<-tcol
if(wall==6) wcol<-a.coladd.ade(bgcol, -35)
 
 
if(wall==0 | wall==2 | wall==5) wlwd<-1
if(wall==1 | wall==3 | wall==4 | wall==6) wlwd<-2



if(length(btext)==1)  btext<-rep(btext, ncol(a.tab))
for(j in 1:(a.n2)){
yr<- a.glc(1, 0)-a.glc(1, 1.5)
segments(xran[j]+min(yrat),    max(a.tab[,j])+(yr/3) ,   xran[j]+max(yrat),   max(a.tab[, j])+(yr/3) , col=wcol, lwd=wlwd)
segments(xran[j]+yrat,    max(a.tab[,j])+(yr/3) ,   xran[j]+yrat,   (a.tab[, j])+(yr/8) , col=wcol, lwd=wlwd)
text(xran[j], max(a.tab[, j])+(yr*0.85) , labels=btext[j], col=tcol, adj=c(0.5, 0.5))
}
}
}
################################
}


if(!beside){
if(is.numeric(cutz) & !is.null(z)){
sum.tab<- NULL
for(k in 1:nlevels(z)){
if(!is.null(sum.tab)) sum.tab<- sum.tab + all.tab[,,k]
if(is.null(sum.tab))  sum.tab<- all.tab[,,k]
}
}
if(!is.numeric(cutz) | is.null(zperc))  TB<-(a.tab/sum(a.tab))*100
if(is.numeric(cutz) &  !is.null(zperc)) if(zperc=='overall')   TB<-(all.tab[,,cutz]/sum(all.tab))*100
if(is.numeric(cutz) &  !is.null(zperc)) if(zperc=='rows')      TB<-(all.tab[,,cutz]/rowSums(sum.tab))*100
if(is.numeric(cutz) &  !is.null(zperc)) if(zperc=='cols')      TB<-t(t(all.tab[,,cutz])/colSums(sum.tab))*100
if(is.numeric(cutz) &  !is.null(zperc)) if(zperc=='zells')      TB<-(all.tab[,,cutz]/sum.tab)*100



hoch<-rep(ylim[1], ncol(a.tab))
for(k in 1:nrow(a.tab)){
if(k> 1) a.form.ade(x = xran, y = hoch, h=a.tab[k, ], b = b, b2 = b2, col=col[k], side=side, nslices=100, gradient=gradient, border=border, bcol=bcol, form=form, density = density[k], angle = angle[k], density2 = density2[k], angle2 = angle2[k], fill=fill[k], blwd=blwd, blty=blty, lwd=lwd, lty=lty, a.tab=a.tab)
if(k==1) a.form.ade(x = xran, y = hoch, h=a.tab[k, ], b = b, b2 = b2, col=col[k], side=side, nslices=100, gradient=gradient, border=border, bcol=bcol, form=form, density = density[k], angle = angle[k], density2 = density2[k], angle2 = angle2[k], fill=fill[k], blwd=blwd, blty=blty, lwd=lwd, lty=lty, a.tab=a.tab)
hoch<-hoch+a.tab[k, ]

################################
if(is.logical(prozent)){ if(prozent){
if(wall!=4) text(x=xran, y=cumat[k, ], labels=paste(round_n.ade(TB[k, ], 1), '%', sep='' ), col=tcol, cex=1, font=2)
if(wall==4) text(x=xran, y=cumat[k, ], labels=paste(round_n.ade(TB[k, ], 1), '%', sep='' ), col=rgb(1,1,1), cex=1, font=2)
}}


if(is.character(prozent) | is.numeric(prozent)){
if(wall!=4) text(x=xran, y=cumat[k, ], labels=prozent[k, ], col=tcol, cex=1, font=2)
if(wall==4) text(x=xran, y=cumat[k, ], labels=prozent[k, ], col=rgb(1,1,1), cex=1, font=2)
}
################################
}

################################
# Text lines
if(!is.null(btext)){
if(is.logical(btext)){
if(!btext)  btext<-NULL
if(btext){
btext<-NULL
for(i in 1:ncol(a.tab)){
btext[i]<- paste('p:', format_p.ade(chisq.test(a.tab[ ,i])$p.value , 4) )
}
}
}
if(!is.logical(btext)){
if(length(btext)==1)  btext<-rep(btext, nrow(a.tab))
brt=4
if(form=='z') brt=1.5
yr<- a.glc(1, 0)-a.glc(1, 3)
text(xran, hoch+(yr/brt) , labels=btext, col=tcol)
}
}
################################

}

}
################################################################################





################################################################################
if(length(dim(a.tab))==3){
TB<-(a.tab/sum(a.tab))*100

if(beside){
b<-(b/nrow(a.tab))*0.75
yrat<-seq(-nrow(a.tab)+1, nrow(a.tab)-1, length.out=nrow(a.tab))*(b/2)
colrun<-1
for(i in 1:nrow(a.tab)){
hoch<-rep(ylim[1], dim(a.tab)[2])
for(k in 1:dim(a.tab)[3]){

if(k> 1) a.form.ade(x = xran+yrat[i], y = hoch, h=a.tab[i, ,k], b = b, b2 = b2, col=col[i], side=side, gradient=gradient, nslices=100, border=border, bcol=bcol, form=form, angle=angle[k], density=density[k], angle2=angle2[k], density2=density2[k], blwd=blwd, blty=blty, lwd=lwd, lty=lty, fill=fill[i], a.tab=a.tab)
if(k==1) a.form.ade(x = xran+yrat[i], y = hoch, h=a.tab[i, ,k], b = b, b2 = b2, col=col[i], side=side, gradient=gradient, nslices=100, border=border, bcol=bcol, form=form, angle=angle[k], density=density[k], angle2=angle2[k], density2=density2[k], blwd=blwd, blty=blty, lwd=lwd, lty=lty, fill=fill[i], a.tab=a.tab)
hoch<-hoch+a.tab[i, ,k]
colrun<-colrun+1
}
}





}


if(!beside){
b<-(b/nrow(a.tab))*0.75
yrat<-seq(-nrow(a.tab)+1, nrow(a.tab)-1, length.out=nrow(a.tab))*(b/2)

colrun<-1
for(i in 1:nrow(a.tab)){
hoch<-rep(ylim[1], dim(a.tab)[2])
for(k in 1:dim(a.tab)[3]){
if(k> 1) a.form.ade(x = xran+yrat[i], y = hoch, h=a.tab[i, ,k], b = b, b2 = b2, col=col[k], side=side, gradient=gradient, nslices=100, border=border, bcol=bcol, form=form, angle=angle[i], density=density[i], angle2=angle2[i], density2=density2[i], blwd=blwd, blty=blty, lwd=lwd, lty=lty, fill=fill[k], a.tab=a.tab)
if(k==1) a.form.ade(x = xran+yrat[i], y = hoch, h=a.tab[i, ,k], b = b, b2 = b2, col=col[k], side=side, gradient=gradient, nslices=100, border=border, bcol=bcol, form=form, angle=angle[i], density=density[i], angle2=angle2[i], density2=density2[i], blwd=blwd, blty=blty, lwd=lwd, lty=lty, fill=fill[k], a.tab=a.tab)
hoch<-hoch+a.tab[i, ,k]
colrun<-colrun+1
}
}

}

}
################################################################################




}
