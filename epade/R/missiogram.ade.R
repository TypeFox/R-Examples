missiogram.ade <-
function(vars=NULL,  vnames=NULL, data=NULL , ints=50,  nvars=50, xlab='ID', ylab='Variables', main='Missing Value Plot', ylab2='N. Missings', col=NULL, tcol=NULL,  bgcol=NULL, wall=0){
if(any(par('mfg')!=c(1,1,1,1)) & any(par('mai') < c(1.02, 0.82, 0.82, 0.42))){
maidiff<-rep(0, 4)
norm<-c(1.02, 0.82, 0.82, 0.42)
maidiff[par('mai')<norm]<-  norm[par('mai')<norm] - par('mai')[par('mai')<norm]
par(mai=par('mai')+maidiff)
}
oldpar<-par(no.readonly =TRUE)
oldpar<-oldpar[-which(names(oldpar)%in%c('usr', 'plt', 'pin', 'fin', 'fig', 'mfg', 'mfcol', 'mfrow', 'omd', 'omi', 'oma'))]
on.exit(par(oldpar))



dname<- deparse(substitute(data))
if(is.data.frame(vars)){
data<-vars
vars<-NULL
}
if(is.null(vars))   nvars<-min(c(nvars, dim(data)[2]) )
if(is.null(vars))   {
vars<-names(data)[1:nvars]
vars<-paste("'", vars, "'", sep='')
}
if(is.null(vnames)) vnames <- vars
if(!is.data.frame(data))  stop("(data) must be a data.frame!")
mycex<-0.9
if(length(vars)>30) mycex<-0.85
if(length(vars)>40) mycex<-0.75
if(length(vars)>=50) mycex<-0.65
if(length(vars)>60) mycex<-0.55
if(length(vars)>70) mycex<-0.45
if(length(vars)>80) mycex<-0.35

cbig<-0.6
if(ints>99)  cbig<-0.7
if(ints>199)  cbig<-1
if(ints<75)  cbig<-0.5
if(ints==50)  cbig<-0.45
if(ints<50)  cbig<-0.35
if(ints<31)  cbig<-0.25
if(ints<21)  cbig<-0.15
if(ints<11)  cbig<-0.1

M<-NULL
N<-dim(data)[1]
if(ints>N)  ints<-N
nz<- round(N/ints)
nmissings<-NULL

for(i in 1:length(vars)){
v<- eval(parse(text=paste("data$",vars[i], sep='')))
nmissings<-c(nmissings, sum(is.na(v)))
p<-NULL
k<-1
while(k<N){
myend<-k+nz-1
if(myend>N) myend<-N
p<-c(p, sum(is.na(v[k:myend])))
k<- k+nz
}
M<-rbind(M, p)
}


MP<-(M-3)/nz
MP[MP>1]<-1
MP[MP<0]<-0
MP[MP>0 & MP<0.05]<-0.075
MP[M==0 | M==1 | M==2 | M==3]<-0


MP2<-MP
MP2[M!=1]<-0
MP2[M==1]<-0.25

MP3<-MP
MP3[M!=2] <-0
MP3[M==2]<-0.25

MP4<-MP
MP4[M!=3]<-0
MP4[M==3]<-0.25


cols<-dim(M)[2]
rows<-dim(M)[1]

yray<-length(vars)
xray<-ints



if(is.null(tcol)  & wall==0)   tcol<-1
if(is.null(tcol)  & wall!=0)   tcol<-rgb(0.1,0.1,0.25)
if(is.null(bgcol) & wall==0)   bgcol<-1
if(is.null(bgcol) & wall!=0)   bgcol<-'#DBE0E8'
if(is.null(col)  & wall==0)   col<-'gray25'
if(is.null(col)  & wall!=0)   col<-rgb(0.3,0.3,0.45)


################################################################################
if(wall==0){

newmai<-rep(0, 4)
oldmai<-par('mai')
if(oldmai[1]>0.95 & oldmai[1]<=1.02) newmai[1]<- 0.95-oldmai[1]
if(oldmai[2]<2.35) newmai[2]<- 2.35 - oldmai[2]
if(oldmai[3]>0.65 & oldmai[3]<=0.82) newmai[3]<- 0.65-oldmai[3]
if(oldmai[4]<1.25) newmai[4]<- 1.25 - oldmai[4]
par(mai=(oldmai+newmai))



par(col.axis=tcol)
par(cex.axis=mycex)
par(col.lab=tcol)
par(col.main=tcol)
axicol<-bgcol

plot(1,1, xlim=c(0+0.035*xray,xray-0.035*xray), ylim=c(0+0.035*yray,yray-0.035*yray), axes=F, xlab='', ylab='', main='', col=rgb(1,1,1, 0))
axis(2, at=(1:yray)-0.5, labels=vnames[length(vnames):1], las=2, col=axicol)
axis(4, at=(1:yray)-0.5, labels=nmissings[length(vnames):1], las=2, col=axicol)
axis(1, at=(0:(ints-1)+0.5), labels=(1:ints)*nz, col=axicol)

if(!is.null(xlab))  mtext(xlab,  side = 1, line = 2.5,  col = tcol, font = 1)
if(!is.null(ylab))  mtext(ylab,  side = 2, line = 10,   col = tcol, font = 1)
if(!is.null(main))  mtext(main,  side = 3, line = 1.25, col = tcol, font = 2)
if(!is.null(ylab2)) mtext(ylab2, side = 4, line = 4,    col = tcol, font = 1)


####################
# Labels
#text(x=a.glc(0), y=a.glc(3, 2.1), labels=main, cex=par('cex.main'), col=tcol, adj=c(0.5, 0.5), xpd=TRUE, font=par("font.main"))
#text(x=a.glc(0), y=a.glc(1, 3.6), labels=xlab, cex=par('cex.lab'), col=tcol, adj=c(0.5, 0.5), xpd=TRUE, font=par("font.lab"))
#text(x=a.glc(2, 2.8), y=a.glc(5), labels=ylab, cex=par('cex.lab'), col=tcol, adj=c(0.5, 0.5), xpd=TRUE, font=par("font.lab"), srt=90)
####################



box(col=axicol)

############################
for(i in 1:length(vars)) {
symbols((1:ints)-0.5, rep((yray-i), ints)+0.5 ,  rectangles=cbind(rep(1, rep(ints)), rep(1, rep(ints))), inches = FALSE, add = TRUE, fg = a.alpha.ade(col, MP[i, ]), bg = a.alpha.ade(col, MP[i, ]))

symbols((1:ints)-0.5, rep((yray-i), ints)+0.5 ,  circles=rep(cbig, rep(ints)), inches = FALSE, add = TRUE, fg = a.alpha.ade(col, MP2[i, ]), bg = a.alpha.ade(col, MP2[i, ]), xlim=c(0+0.035*xray,xray-0.035*xray), ylim=c(0+0.035*yray,yray-0.035*yray))

symbols((1:ints)-0.5, rep((yray-i), ints)+0.65 ,  circles=rep(cbig, rep(ints)), inches = FALSE, add = TRUE, fg = a.alpha.ade(col, MP3[i, ]), bg = a.alpha.ade(col, MP3[i, ]), xlim=c(0+0.035*xray,xray-0.035*xray), ylim=c(0+0.035*yray,yray-0.035*yray))
symbols((1:ints)-0.5, rep((yray-i), ints)+0.35 ,  circles=rep(cbig, rep(ints)), inches = FALSE, add = TRUE, fg = a.alpha.ade(col, MP3[i, ]), bg = a.alpha.ade(col, MP3[i, ]), xlim=c(0+0.035*xray,xray-0.035*xray), ylim=c(0+0.035*yray,yray-0.035*yray))

symbols((1:ints)-0.5, rep((yray-i), ints)+0.75 ,  circles=rep(cbig, rep(ints)), inches = FALSE, add = TRUE, fg = a.alpha.ade(col, MP4[i, ]), bg = a.alpha.ade(col, MP4[i, ]), xlim=c(0+0.035*xray,xray-0.035*xray), ylim=c(0+0.035*yray,yray-0.035*yray))
symbols((1:ints)-0.5, rep((yray-i), ints)+0.5 ,   circles=rep(cbig, rep(ints)), inches = FALSE, add = TRUE, fg = a.alpha.ade(col, MP4[i, ]), bg = a.alpha.ade(col, MP4[i, ]), xlim=c(0+0.035*xray,xray-0.035*xray), ylim=c(0+0.035*yray,yray-0.035*yray))
symbols((1:ints)-0.5, rep((yray-i), ints)+0.25 ,  circles=rep(cbig, rep(ints)), inches = FALSE, add = TRUE, fg = a.alpha.ade(col, MP4[i, ]), bg = a.alpha.ade(col, MP4[i, ]), xlim=c(0+0.035*xray,xray-0.035*xray), ylim=c(0+0.035*yray,yray-0.035*yray))


abline(h=(yray-i), col=bgcol)
}
############################

}
################################################################################


################################################################################
if(wall==1){
newmai<-rep(0, 4)
oldmai<-par('mai')
if(oldmai[1]>0.95 & oldmai[1]<=1.02) newmai[1]<- 0.95-oldmai[1]
if(oldmai[2]<2.35) newmai[2]<- 2.35 - oldmai[2]
if(oldmai[3]>0.65 & oldmai[3]<=0.82) newmai[3]<- 0.65-oldmai[3]
if(oldmai[4]<1.25) newmai[4]<- 1.25 - oldmai[4]
par(mai=(oldmai+newmai))

par(col.axis=tcol)
par(col.lab=tcol)
par(cex.axis=mycex)
par(col.main=tcol)
axicol<-tcol

plot(1,1, xlim=c(0+0.035*xray,xray-0.035*xray), ylim=c(0+0.035*yray,yray-0.035*yray), axes=F, xlab='', ylab='', main='', col=rgb(1,1,1, 0))
axis(2, at=(1:yray)-0.5, labels=vnames[length(vnames):1], las=2, col=axicol)
axis(4, at=(1:yray)-0.5, labels=nmissings[length(vnames):1], las=2, col=axicol)
axis(1, at=(0:(ints-1)+0.5), labels=(1:ints)*nz, col=axicol)
polygon( c(par('usr')[c(1,1,2,2)]), par('usr')[c(3,4,4,3)], col=bgcol, border= FALSE)

if(!is.null(xlab))  mtext(xlab,  side = 1, line = 2.5,  col = tcol, font = 1)
if(!is.null(ylab))  mtext(ylab,  side = 2, line = 10, col = tcol, font = 1)
if(!is.null(main))  mtext(main,  side = 3, line = 1.25,    col = tcol, font = 2)
if(!is.null(ylab2)) mtext(ylab2, side = 4, line = 4,    col = tcol, font = 1)

box(col=rgb(1,1,1))

############################
for(i in 1:length(vars)) {
symbols((1:ints)-0.5, rep((yray-i), ints)+0.5 ,  rectangles=cbind(rep(1, rep(ints)), rep(1, rep(ints))), inches = FALSE, add = TRUE, fg = a.alpha.ade(col, MP[i, ]), bg = a.alpha.ade(col, MP[i, ]))
symbols((1:ints)-0.5, rep((yray-i), ints)+0.5 ,  circles=rep(cbig, rep(ints)), inches = FALSE, add = TRUE, fg = a.alpha.ade(col, MP2[i, ]), bg = a.alpha.ade(col, MP2[i, ]), xlim=c(0+0.035*xray,xray-0.035*xray), ylim=c(0+0.035*yray,yray-0.035*yray))

symbols((1:ints)-0.5, rep((yray-i), ints)+0.65 ,  circles=rep(cbig, rep(ints)), inches = FALSE, add = TRUE, fg = a.alpha.ade(col, MP3[i, ]), bg = a.alpha.ade(col, MP3[i, ]), xlim=c(0+0.035*xray,xray-0.035*xray), ylim=c(0+0.035*yray,yray-0.035*yray))
symbols((1:ints)-0.5, rep((yray-i), ints)+0.35 ,  circles=rep(cbig, rep(ints)), inches = FALSE, add = TRUE, fg = a.alpha.ade(col, MP3[i, ]), bg = a.alpha.ade(col, MP3[i, ]), xlim=c(0+0.035*xray,xray-0.035*xray), ylim=c(0+0.035*yray,yray-0.035*yray))

symbols((1:ints)-0.5, rep((yray-i), ints)+0.75 ,  circles=rep(cbig, rep(ints)), inches = FALSE, add = TRUE, fg = a.alpha.ade(col, MP4[i, ]), bg = a.alpha.ade(col, MP4[i, ]), xlim=c(0+0.035*xray,xray-0.035*xray), ylim=c(0+0.035*yray,yray-0.035*yray))
symbols((1:ints)-0.5, rep((yray-i), ints)+0.5 ,   circles=rep(cbig, rep(ints)), inches = FALSE, add = TRUE, fg = a.alpha.ade(col, MP4[i, ]), bg = a.alpha.ade(col, MP4[i, ]), xlim=c(0+0.035*xray,xray-0.035*xray), ylim=c(0+0.035*yray,yray-0.035*yray))
symbols((1:ints)-0.5, rep((yray-i), ints)+0.25 ,  circles=rep(cbig, rep(ints)), inches = FALSE, add = TRUE, fg = a.alpha.ade(col, MP4[i, ]), bg = a.alpha.ade(col, MP4[i, ]), xlim=c(0+0.035*xray,xray-0.035*xray), ylim=c(0+0.035*yray,yray-0.035*yray))



abline(h=(yray-i), col=rgb(1,1,1))
}
############################

}
################################################################################


################################################################################
if(wall==2){
newmai<-rep(0, 4)
oldmai<-par('mai')
if(oldmai[1]>0.95 & oldmai[1]<=1.02) newmai[1]<- 0.95-oldmai[1]
if(oldmai[2]<2.35) newmai[2]<- 2.35 - oldmai[2]
if(oldmai[3]>0.65 & oldmai[3]<=0.82) newmai[3]<- 0.65-oldmai[3]
if(oldmai[4]<1.25) newmai[4]<- 1.25 - oldmai[4]
par(mai=(oldmai+newmai))

par(col.axis=tcol)
par(col.lab=tcol)
par(cex.axis=mycex)
par(col.main=tcol)
axicol<-a.coladd.ade(bgcol, -75)

plot(1,1, xlim=c(0+0.035*xray,xray-0.035*xray), ylim=c(0+0.035*yray,yray-0.035*yray), axes=F, xlab='', ylab='', main='', col=rgb(1,1,1, 0))
axis(2, at=(1:yray)-0.5, labels=vnames[length(vnames):1], las=2, col=axicol)
axis(4, at=(1:yray)-0.5, labels=nmissings[length(vnames):1], las=2, col=axicol)
axis(1, at=(0:(ints-1)+0.5), labels=(1:ints)*nz, col=axicol)
#polygon( c(par('usr')[c(1,1,2,2)]), par('usr')[c(3,4,4,3)], col=bgcol, border= FALSE)

if(!is.null(xlab))  mtext(xlab,  side = 1, line = 2.5,  col = tcol, font = 1)
if(!is.null(ylab))  mtext(ylab,  side = 2, line = 10, col = tcol, font = 1)
if(!is.null(main))  mtext(main,  side = 3, line = 1.25,    col = tcol, font = 2)
if(!is.null(ylab2)) mtext(ylab2, side = 4, line = 4,    col = tcol, font = 1)

box(col=a.coladd.ade(bgcol, -75))

############################
for(i in 1:length(vars)) {
symbols((1:ints)-0.5, rep((yray-i), ints)+0.5 ,  rectangles=cbind(rep(1, rep(ints)), rep(1, rep(ints))), inches = FALSE, add = TRUE, fg = a.alpha.ade(col, MP[i, ]), bg = a.alpha.ade(col, MP[i, ]))
symbols((1:ints)-0.5, rep((yray-i), ints)+0.5 ,  circles=rep(cbig, rep(ints)), inches = FALSE, add = TRUE, fg = a.alpha.ade(col, MP2[i, ]), bg = a.alpha.ade(col, MP2[i, ]), xlim=c(0+0.035*xray,xray-0.035*xray), ylim=c(0+0.035*yray,yray-0.035*yray))

symbols((1:ints)-0.5, rep((yray-i), ints)+0.65 ,  circles=rep(cbig, rep(ints)), inches = FALSE, add = TRUE, fg = a.alpha.ade(col, MP3[i, ]), bg = a.alpha.ade(col, MP3[i, ]), xlim=c(0+0.035*xray,xray-0.035*xray), ylim=c(0+0.035*yray,yray-0.035*yray))
symbols((1:ints)-0.5, rep((yray-i), ints)+0.35 ,  circles=rep(cbig, rep(ints)), inches = FALSE, add = TRUE, fg = a.alpha.ade(col, MP3[i, ]), bg = a.alpha.ade(col, MP3[i, ]), xlim=c(0+0.035*xray,xray-0.035*xray), ylim=c(0+0.035*yray,yray-0.035*yray))

symbols((1:ints)-0.5, rep((yray-i), ints)+0.75 ,  circles=rep(cbig, rep(ints)), inches = FALSE, add = TRUE, fg = a.alpha.ade(col, MP4[i, ]), bg = a.alpha.ade(col, MP4[i, ]), xlim=c(0+0.035*xray,xray-0.035*xray), ylim=c(0+0.035*yray,yray-0.035*yray))
symbols((1:ints)-0.5, rep((yray-i), ints)+0.5 ,   circles=rep(cbig, rep(ints)), inches = FALSE, add = TRUE, fg = a.alpha.ade(col, MP4[i, ]), bg = a.alpha.ade(col, MP4[i, ]), xlim=c(0+0.035*xray,xray-0.035*xray), ylim=c(0+0.035*yray,yray-0.035*yray))
symbols((1:ints)-0.5, rep((yray-i), ints)+0.25 ,  circles=rep(cbig, rep(ints)), inches = FALSE, add = TRUE, fg = a.alpha.ade(col, MP4[i, ]), bg = a.alpha.ade(col, MP4[i, ]), xlim=c(0+0.035*xray,xray-0.035*xray), ylim=c(0+0.035*yray,yray-0.035*yray))

abline(h=(yray-i), col=bgcol)
}
############################

}
################################################################################


################################################################################
if(wall==3){

newmai<-rep(0, 4)
oldmai<-par('mai')
if(oldmai[1]>0.95 & oldmai[1]<=1.02) newmai[1]<- 0.95-oldmai[1]
if(oldmai[2]<2.35) newmai[2]<- 2.35 - oldmai[2]
if(oldmai[3]>0.65 & oldmai[3]<=0.82) newmai[3]<- 0.65-oldmai[3]
if(oldmai[4]<1.25) newmai[4]<- 1.25 - oldmai[4]
par(mai=(oldmai+newmai))


par(col.axis=tcol)
par(col.lab=tcol)
par(cex.axis=mycex)
par(col.main=tcol)
axicol<-a.coladd.ade(bgcol, -75)

plot(1,1, xlim=c(0+0.035*xray,xray-0.035*xray), ylim=c(0+0.035*yray,yray-0.035*yray), axes=F, xlab='', ylab='', main='', col=rgb(1,1,1, 0))
axis(2, at=(1:yray)-0.5, labels=vnames[length(vnames):1], las=2, col=axicol)
axis(4, at=(1:yray)-0.5, labels=nmissings[length(vnames):1], las=2, col=axicol)
axis(1, at=(0:(ints-1)+0.5), labels=(1:ints)*nz, col=axicol)
polygon( c(par('usr')[c(1,1,2,2)]), par('usr')[c(3,4,4,3)], col=bgcol, border= FALSE)

if(!is.null(xlab))  mtext(xlab,  side = 1, line = 2.5,  col = tcol, font = 1)
if(!is.null(ylab))  mtext(ylab,  side = 2, line = 10, col = tcol, font = 1)
if(!is.null(main))  mtext(main,  side = 3, line = 1.25,    col = tcol, font = 2)
if(!is.null(ylab2)) mtext(ylab2, side = 4, line = 4,    col = tcol, font = 1)

box(col=a.coladd.ade(bgcol, -75))

############################
for(i in 1:length(vars)) {
symbols((1:ints)-0.5, rep((yray-i), ints)+0.5 ,  rectangles=cbind(rep(1, rep(ints)), rep(1, rep(ints))), inches = FALSE, add = TRUE, fg = a.alpha.ade(col, MP[i, ]), bg = a.alpha.ade(col, MP[i, ]))
symbols((1:ints)-0.5, rep((yray-i), ints)+0.5 ,  circles=rep(cbig, rep(ints)), inches = FALSE, add = TRUE, fg = a.alpha.ade(col, MP2[i, ]), bg = a.alpha.ade(col, MP2[i, ]), xlim=c(0+0.035*xray,xray-0.035*xray), ylim=c(0+0.035*yray,yray-0.035*yray))

symbols((1:ints)-0.5, rep((yray-i), ints)+0.65 ,  circles=rep(cbig, rep(ints)), inches = FALSE, add = TRUE, fg = a.alpha.ade(col, MP3[i, ]), bg = a.alpha.ade(col, MP3[i, ]), xlim=c(0+0.035*xray,xray-0.035*xray), ylim=c(0+0.035*yray,yray-0.035*yray))
symbols((1:ints)-0.5, rep((yray-i), ints)+0.35 ,  circles=rep(cbig, rep(ints)), inches = FALSE, add = TRUE, fg = a.alpha.ade(col, MP3[i, ]), bg = a.alpha.ade(col, MP3[i, ]), xlim=c(0+0.035*xray,xray-0.035*xray), ylim=c(0+0.035*yray,yray-0.035*yray))

symbols((1:ints)-0.5, rep((yray-i), ints)+0.75 ,  circles=rep(cbig, rep(ints)), inches = FALSE, add = TRUE, fg = a.alpha.ade(col, MP4[i, ]), bg = a.alpha.ade(col, MP4[i, ]), xlim=c(0+0.035*xray,xray-0.035*xray), ylim=c(0+0.035*yray,yray-0.035*yray))
symbols((1:ints)-0.5, rep((yray-i), ints)+0.5 ,   circles=rep(cbig, rep(ints)), inches = FALSE, add = TRUE, fg = a.alpha.ade(col, MP4[i, ]), bg = a.alpha.ade(col, MP4[i, ]), xlim=c(0+0.035*xray,xray-0.035*xray), ylim=c(0+0.035*yray,yray-0.035*yray))
symbols((1:ints)-0.5, rep((yray-i), ints)+0.25 ,  circles=rep(cbig, rep(ints)), inches = FALSE, add = TRUE, fg = a.alpha.ade(col, MP4[i, ]), bg = a.alpha.ade(col, MP4[i, ]), xlim=c(0+0.035*xray,xray-0.035*xray), ylim=c(0+0.035*yray,yray-0.035*yray))


abline(h=(yray-i), col=a.coladd.ade(bgcol, -50))
}
############################

}
################################################################################


################################################################################
if(wall==4){
newmai<-rep(0, 4)
oldmai<-par('mai')
if(oldmai[1]>0.95 & oldmai[1]<=1.02) newmai[1]<- 0.95-oldmai[1]
if(oldmai[2]<2.35) newmai[2]<- 2.35 - oldmai[2]
if(oldmai[3]>0.65 & oldmai[3]<=0.82) newmai[3]<- 0.65-oldmai[3]
if(oldmai[4]<1.25) newmai[4]<- 1.25 - oldmai[4]
par(mai=(oldmai+newmai))

par(col.axis=tcol)
par(cex.axis=mycex-0.05)
par(col.lab=tcol)
par(col.main=tcol)
axicol<-tcol

plot(1,1, xlim=c(0+0.035*xray,xray-0.035*xray), ylim=c(0+0.035*yray,yray-0.035*yray), axes=F, xlab='', ylab='', main='', col=rgb(1,1,1, 0))
polygon( c(par('usr')[c(1,1,2,2)]), par('usr')[c(3,4,4,3)], col=bgcol, border= FALSE)
par(xpd= TRUE)
polygon(a.glc(side=c(2,2,4,4), line=c(0,0,0,0)), a.glc(side=3, line=c(0, 2.75,  2.75, 0)), col=tcol, border=rgb(1,1,1))
if(ylab!='' & ylab!=' ') polygon( a.glc(side=2, line=c(11, 11, 9.5, 9.5)), a.glc(side=c(1, 3, 3, 1), line=0), col=bgcol, border=rgb(1,1,1))
if(xlab!='' & xlab!=' ') polygon( a.glc(side=c(2, 2, 4, 4), line=0),     a.glc(side=1, line=c(4, 2.5, 2.5, 4)), col=bgcol, border=rgb(1,1,1))
if(ylab2!='' & ylab2!=' ') polygon( a.glc(side=4, line=c(4, 4, 5.5, 5.5)), a.glc(side=c(1, 3, 3, 1), line=0), col=bgcol, border=rgb(1,1,1))
text(a.glc(side=0), a.glc(side=3, line=1),    labels=main, cex = 1.25, font=2, col=rgb(1,1,1), adj=c(0.5,0))
text(a.glc(side=0), a.glc(side=1, line=3.5),  labels=xlab, cex = 1.1,  font=2, col=tcol, adj=c(0.5,0))
text(a.glc(side=2, line=10), a.glc(side=5),  labels=ylab, cex = 1.1,  font=2,  col=tcol, adj=c(0.5,0), srt=90)
text(a.glc(side=4, line=5), a.glc(side=5),  labels=ylab2, cex = 1.1,  font=2,  col=tcol, adj=c(0.5,0), srt=90)

par(xpd= FALSE)
axis(2, at=(1:yray)-0.5, labels=vnames[length(vnames):1], las=2, col=axicol)
axis(4, at=(1:yray)-0.5, labels=nmissings[length(vnames):1], las=2, col=axicol)
axis(1, at=(0:(ints-1)+0.5), labels=(1:ints)*nz, col=tcol)

box(col=rgb(1,1,1))
############################
for(i in 1:length(vars)) {
symbols((1:ints)-0.5, rep((yray-i), ints)+0.5 ,  rectangles=cbind(rep(1, rep(ints)), rep(1, rep(ints))), inches = FALSE, add = TRUE, fg = a.alpha.ade(col, MP[i, ]), bg = a.alpha.ade(col, MP[i, ]))
symbols((1:ints)-0.5, rep((yray-i), ints)+0.5 ,  circles=rep(cbig, rep(ints)), inches = FALSE, add = TRUE, fg = a.alpha.ade(col, MP2[i, ]), bg = a.alpha.ade(col, MP2[i, ]), xlim=c(0+0.035*xray,xray-0.035*xray), ylim=c(0+0.035*yray,yray-0.035*yray))

symbols((1:ints)-0.5, rep((yray-i), ints)+0.65 ,  circles=rep(cbig, rep(ints)), inches = FALSE, add = TRUE, fg = a.alpha.ade(col, MP3[i, ]), bg = a.alpha.ade(col, MP3[i, ]), xlim=c(0+0.035*xray,xray-0.035*xray), ylim=c(0+0.035*yray,yray-0.035*yray))
symbols((1:ints)-0.5, rep((yray-i), ints)+0.35 ,  circles=rep(cbig, rep(ints)), inches = FALSE, add = TRUE, fg = a.alpha.ade(col, MP3[i, ]), bg = a.alpha.ade(col, MP3[i, ]), xlim=c(0+0.035*xray,xray-0.035*xray), ylim=c(0+0.035*yray,yray-0.035*yray))

symbols((1:ints)-0.5, rep((yray-i), ints)+0.75 ,  circles=rep(cbig, rep(ints)), inches = FALSE, add = TRUE, fg = a.alpha.ade(col, MP4[i, ]), bg = a.alpha.ade(col, MP4[i, ]), xlim=c(0+0.035*xray,xray-0.035*xray), ylim=c(0+0.035*yray,yray-0.035*yray))
symbols((1:ints)-0.5, rep((yray-i), ints)+0.5 ,   circles=rep(cbig, rep(ints)), inches = FALSE, add = TRUE, fg = a.alpha.ade(col, MP4[i, ]), bg = a.alpha.ade(col, MP4[i, ]), xlim=c(0+0.035*xray,xray-0.035*xray), ylim=c(0+0.035*yray,yray-0.035*yray))
symbols((1:ints)-0.5, rep((yray-i), ints)+0.25 ,  circles=rep(cbig, rep(ints)), inches = FALSE, add = TRUE, fg = a.alpha.ade(col, MP4[i, ]), bg = a.alpha.ade(col, MP4[i, ]), xlim=c(0+0.035*xray,xray-0.035*xray), ylim=c(0+0.035*yray,yray-0.035*yray))

abline(h=(yray-i), col=rgb(1,1,1))
}
############################

}
################################################################################



################################################################################
if(wall==5){
newmai<-rep(0, 4)
oldmai<-par('mai')
#if(oldmai[1]<2.35) newmai[1]<- 2.35 - oldmai[1]
if(oldmai[2]<2.4) newmai[2]<- 2.4 - oldmai[2]
if(oldmai[3]>0.75 & oldmai[3]<=0.82) newmai[3]<- 0.75-oldmai[3]
if(oldmai[4]<1.25) newmai[4]<- 1.25 - oldmai[4]
par(mai=(oldmai+newmai))


par(col.axis=tcol)
par(cex.axis=mycex-0.05)
par(col.lab=tcol)
par(col.main=tcol)
axicol<-tcol

plot(1,1, xlim=c(0+0.035*xray,xray-0.035*xray), ylim=c(0+0.035*yray,yray-0.035*yray), axes=F, xlab='', ylab='', main='', col=rgb(1,1,1, 0))
par(xpd= TRUE)
polygon(a.glc(side=2, line=c(11.25, 11.25, 0, 0)), a.glc(side=3, line=c(0.6, 3, 3, 0.6)), col=bgcol,        border=tcol)
polygon(a.glc(side=c(2,2,4,4), line=c(0,0,0,0)), a.glc(side=3, line=c(0.6, 3, 3, 0.6)), col=rgb(1,1,1,0), border=tcol)
polygon(a.glc(side=4, line=c(0, 0 ,5.5, 5.5)),   a.glc(side=3, line=c(0.6, 3, 3, 0.6)), col=bgcol,        border=tcol)
polygon(a.glc(side=2, line=c(11.25, 11.25 ,10.65, 10.65)),  a.glc(side=c(1,3,3,1), line=c(2.6, 0.6, 0.6, 2.6)), col=bgcol,  border=tcol)
polygon(a.glc(side=4, line=c(4.9, 4.9 ,5.5, 5.5)), a.glc(side=c(1, 3, 3, 1), line=c(2.6,0.6,0.6,2.6)), col=bgcol, border=tcol)
polygon(a.glc(side=2, line=c(11.25, 11.25, 0, 0)), a.glc(side=1, line=c(2.6, 4.5, 4.5, 2.6)), col=bgcol, border=tcol)
polygon(a.glc(side=c(2, 2, 4, 4), line=0), a.glc(side=1, line=c(2.6, 4.5, 4.5, 2.6)), col=rgb(1,1,1,0), border=tcol)
polygon(a.glc(side=4, line=c(0, 0, 5.5, 5.5)), a.glc(side=1, line=c(2.6, 4.5, 4.5, 2.6)), col=bgcol, border=tcol)
if(!is.null(main)) text(a.glc(side=0), a.glc(side=3, line=1.5),  labels=main, cex = 1.25, font=2, col=tcol, adj=c(0.5,0))
if(!is.null(xlab)) text(a.glc(side=0), a.glc(side=1, line=3.8), labels=xlab, cex = 1.1,  font=2, col=tcol, adj=c(0.5,0))
if(!is.null(ylab)) text(a.glc(side=2, line=9.5), a.glc(side=5), labels=ylab, cex = 1.1,   font=2,  col=tcol, adj=c(0.5,0), srt=90)
if(!is.null(ylab2)) text(a.glc(side=4, line=4.25), a.glc(side=5), labels=ylab2, cex = 1.1,   font=2,  col=tcol, adj=c(0.5,0), srt=90)


par(xpd= FALSE)
axis(2, at=(1:yray)-0.5, labels=vnames[length(vnames):1], las=2, col=axicol)
axis(4, at=(1:yray)-0.5, labels=nmissings[length(vnames):1], las=2, col=axicol)
axis(1, at=(0:(ints-1)+0.5), labels=(1:ints)*nz, col=tcol)


box(col=tcol)
############################
for(i in 1:length(vars)) {
symbols((1:ints)-0.5, rep((yray-i), ints)+0.5 ,  rectangles=cbind(rep(1, rep(ints)), rep(1, rep(ints))), inches = FALSE, add = TRUE, fg = a.alpha.ade(col, MP[i, ]), bg = a.alpha.ade(col, MP[i, ]))
symbols((1:ints)-0.5, rep((yray-i), ints)+0.5 ,  circles=rep(cbig, rep(ints)), inches = FALSE, add = TRUE, fg = a.alpha.ade(col, MP2[i, ]), bg = a.alpha.ade(col, MP2[i, ]), xlim=c(0+0.035*xray,xray-0.035*xray), ylim=c(0+0.035*yray,yray-0.035*yray))

symbols((1:ints)-0.5, rep((yray-i), ints)+0.65 ,  circles=rep(cbig, rep(ints)), inches = FALSE, add = TRUE, fg = a.alpha.ade(col, MP3[i, ]), bg = a.alpha.ade(col, MP3[i, ]), xlim=c(0+0.035*xray,xray-0.035*xray), ylim=c(0+0.035*yray,yray-0.035*yray))
symbols((1:ints)-0.5, rep((yray-i), ints)+0.35 ,  circles=rep(cbig, rep(ints)), inches = FALSE, add = TRUE, fg = a.alpha.ade(col, MP3[i, ]), bg = a.alpha.ade(col, MP3[i, ]), xlim=c(0+0.035*xray,xray-0.035*xray), ylim=c(0+0.035*yray,yray-0.035*yray))

symbols((1:ints)-0.5, rep((yray-i), ints)+0.75 ,  circles=rep(cbig, rep(ints)), inches = FALSE, add = TRUE, fg = a.alpha.ade(col, MP4[i, ]), bg = a.alpha.ade(col, MP4[i, ]), xlim=c(0+0.035*xray,xray-0.035*xray), ylim=c(0+0.035*yray,yray-0.035*yray))
symbols((1:ints)-0.5, rep((yray-i), ints)+0.5 ,   circles=rep(cbig, rep(ints)), inches = FALSE, add = TRUE, fg = a.alpha.ade(col, MP4[i, ]), bg = a.alpha.ade(col, MP4[i, ]), xlim=c(0+0.035*xray,xray-0.035*xray), ylim=c(0+0.035*yray,yray-0.035*yray))
symbols((1:ints)-0.5, rep((yray-i), ints)+0.25 ,  circles=rep(cbig, rep(ints)), inches = FALSE, add = TRUE, fg = a.alpha.ade(col, MP4[i, ]), bg = a.alpha.ade(col, MP4[i, ]), xlim=c(0+0.035*xray,xray-0.035*xray), ylim=c(0+0.035*yray,yray-0.035*yray))

abline(h=(yray-i), col=tcol)
}
############################


}
################################################################################


################################################################################
if(wall==6){
newmai<-rep(0, 4)
oldmai<-par('mai')
if(oldmai[1]>0.95 & oldmai[1]<=1.02) newmai[1]<- 0.95-oldmai[1]
if(oldmai[2]<2.35) newmai[2]<- 2.35 - oldmai[2]
if(oldmai[3]>0.65 & oldmai[3]<=0.82) newmai[3]<- 0.65-oldmai[3]
if(oldmai[4]<1.25) newmai[4]<- 1.25 - oldmai[4]
par(mai=(oldmai+newmai))

par(col.axis=tcol)
par(col.lab=tcol)
par(cex.axis=mycex)
par(col.main=tcol)
axicol<-tcol

plot(1,1, xlim=c(0+0.035*xray,xray-0.035*xray), ylim=c(0+0.035*yray,yray-0.035*yray), axes=F, xlab='', ylab='', main='', col=rgb(1,1,1, 0))
axis(2, at=(1:yray)-0.5, labels=vnames[length(vnames):1], las=2, col=axicol, col.ticks=a.coladd.ade(bgcol, -35), lwd.ticks=3)
axis(2, at=(1:yray)-0.5, labels=vnames[length(vnames):1], las=2, col=axicol, col.ticks=rgb(1,1,1), lwd.ticks=1)

axis(4, at=(1:yray)-0.5, labels=nmissings[length(vnames):1], las=2, col=axicol, col.ticks=a.coladd.ade(bgcol, -35), lwd.ticks=3)
axis(4, at=(1:yray)-0.5, labels=nmissings[length(vnames):1], las=2, col=axicol, col.ticks=rgb(1,1,1), lwd.ticks=1)

axis(1, at=(0:(ints-1)+0.5), labels=(1:ints)*nz, col=axicol, col.ticks=a.coladd.ade(bgcol, -35), lwd.ticks=3)
axis(1, at=(0:(ints-1)+0.5), labels=(1:ints)*nz, col=axicol, col.ticks=rgb(1,1,1), lwd.ticks=1)

polygon( c(par('usr')[c(1,1,2,2)]), par('usr')[c(3,4,4,3)], col=bgcol, border= FALSE)

if(!is.null(xlab))  mtext(xlab,  side = 1, line = 2.5,  col = tcol, font = 1)
if(!is.null(ylab))  mtext(ylab,  side = 2, line = 10, col = tcol, font = 1)
if(!is.null(main))  mtext(main,  side = 3, line = 1.25,    col = tcol, font = 2)
if(!is.null(ylab2)) mtext(ylab2, side = 4, line = 4,    col = tcol, font = 1)



############################
for(i in 1:length(vars)) {
symbols((1:ints)-0.5, rep((yray-i), ints)+0.5 ,  rectangles=cbind(rep(1, rep(ints)), rep(1, rep(ints))), inches = FALSE, add = TRUE, fg = a.alpha.ade(col, MP[i, ]), bg = a.alpha.ade(col, MP[i, ]))
symbols((1:ints)-0.5, rep((yray-i), ints)+0.5 ,  circles=rep(cbig, rep(ints)), inches = FALSE, add = TRUE, fg = a.alpha.ade(col, MP2[i, ]), bg = a.alpha.ade(col, MP2[i, ]), xlim=c(0+0.035*xray,xray-0.035*xray), ylim=c(0+0.035*yray,yray-0.035*yray))

symbols((1:ints)-0.5, rep((yray-i), ints)+0.65 ,  circles=rep(cbig, rep(ints)), inches = FALSE, add = TRUE, fg = a.alpha.ade(col, MP3[i, ]), bg = a.alpha.ade(col, MP3[i, ]), xlim=c(0+0.035*xray,xray-0.035*xray), ylim=c(0+0.035*yray,yray-0.035*yray))
symbols((1:ints)-0.5, rep((yray-i), ints)+0.35 ,  circles=rep(cbig, rep(ints)), inches = FALSE, add = TRUE, fg = a.alpha.ade(col, MP3[i, ]), bg = a.alpha.ade(col, MP3[i, ]), xlim=c(0+0.035*xray,xray-0.035*xray), ylim=c(0+0.035*yray,yray-0.035*yray))

symbols((1:ints)-0.5, rep((yray-i), ints)+0.75 ,  circles=rep(cbig, rep(ints)), inches = FALSE, add = TRUE, fg = a.alpha.ade(col, MP4[i, ]), bg = a.alpha.ade(col, MP4[i, ]), xlim=c(0+0.035*xray,xray-0.035*xray), ylim=c(0+0.035*yray,yray-0.035*yray))
symbols((1:ints)-0.5, rep((yray-i), ints)+0.5 ,   circles=rep(cbig, rep(ints)), inches = FALSE, add = TRUE, fg = a.alpha.ade(col, MP4[i, ]), bg = a.alpha.ade(col, MP4[i, ]), xlim=c(0+0.035*xray,xray-0.035*xray), ylim=c(0+0.035*yray,yray-0.035*yray))
symbols((1:ints)-0.5, rep((yray-i), ints)+0.25 ,  circles=rep(cbig, rep(ints)), inches = FALSE, add = TRUE, fg = a.alpha.ade(col, MP4[i, ]), bg = a.alpha.ade(col, MP4[i, ]), xlim=c(0+0.035*xray,xray-0.035*xray), ylim=c(0+0.035*yray,yray-0.035*yray))



abline(h=(yray-i), col=a.coladd.ade(bgcol, -35), lwd=3)
abline(h=(yray-i), col=rgb(1,1,1))

box(lwd=3, col=rgb(1,1,1))
box(lwd=1, col=a.coladd.ade(bgcol, -35))
}
############################

}
################################################################################


}
