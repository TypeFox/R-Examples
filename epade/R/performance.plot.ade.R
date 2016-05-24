performance.plot.ade <-
function(pred, event, data=NULL,  vnames=NULL, cutoffs=NULL, cutnames=NULL, main=NULL, xlab='cutoff', ylab='%', xlim=NULL, xticks=12, col=NULL, tcol=NULL, bgcol=NULL, lcol=NULL, alpha=NULL, nints=100, lty=NULL, lwd=2, stats=c(1,2), youden=TRUE, wall=0){
if(any(par('mfg')!=c(1,1,1,1)) & any(par('mai') < c(1.02, 0.82, 0.82, 0.42))){
maidiff<-rep(0, 4)
norm<-c(1.02, 0.82, 0.82, 0.42)
maidiff[par('mai')<norm]<-  norm[par('mai')<norm] - par('mai')[par('mai')<norm]
par(mai=par('mai')+maidiff)
}
oldpar<-par(no.readonly =TRUE)
oldpar<-oldpar[-which(names(oldpar)%in%c('usr', 'plt',   'pin',  'fin', 'fig', 'mfg', 'mfcol', 'mfrow', 'omd', 'omi', 'oma'))]
on.exit(par(oldpar))


##############################
if(!is.character(pred)){
xt<-deparse(substitute(pred))
if(regexpr('~', xt)>=0){
event<-gsub('[~].*$', '', xt)
pred<-gsub('^.*[~]', '', xt)
}}
##############################



#####################################
# without Data.Frame or with
if(is.numeric(pred) & is.numeric(event) ){
data<-NULL
data<-as.data.frame(pred)
data$mypred <-pred
data$myevent<-event
pred <-'mypred'
event<-'myevent'
}
#####################################



subd <- subset(data, !is.na(eval(parse(text=paste(pred, sep='')))))
subd <- subset(subd, !is.na(eval(parse(text=paste(event, sep='')))))
x<-eval(parse(text=paste("subd$",pred , sep='')))
y<-eval(parse(text=paste("subd$",event, sep='')))
if(nints>length(x)) nints<-  length(x)

######################################
a.diags.ade<-function(test , treat){
tab<-table(test, treat)

if(all(test==1)){
n.tp<- tab[1,2]
n.tn<- 0
n.fp<- tab[1,1]
n.fn<- 0
}

if(all(treat==1)){
n.tp<- tab[2,1]
n.tn<- 0
n.fp<- 0
n.fn<- tab[1,1]
}

if(all(test==0)){
n.tp<- 0
n.tn<- tab[1,1]
n.fp<- 0
n.fn<- tab[1,2]
}


if(all(treat==0)){
n.tp<- 0
n.tn<- tab[1,1]
n.fp<- tab[2,1]
n.fn<- 0
}


N<-sum(tab)
if(length(unique(test))==2 & length(unique(test))==2){
n.tp<- tab[2,2]
n.tn<- tab[1,1]
n.fp<- tab[2,1]
n.fn<- tab[1,2]
}
sens<- n.tp/(n.tp+n.fn)
spec<- n.tn/(n.tn+n.fp)
ppv <- n.tp/(n.tp+n.fp)
npv <- n.tn/(n.tn+n.fn)


yod<-sens+spec-1
out<-c((n.tp/N), (n.tn/N), (n.fp/N), (n.fn/N), sens, spec, ppv, npv, yod)

names(out)<-c('TP', 'TN', 'FP', 'FN', 'TPR', 'TNR', 'PPV', 'NPV', 'youden')
return(out)
}
######################################


######################################
# Colors
if(is.null(tcol)  & wall==0)   tcol<-1
if(is.null(tcol)  & wall!=0)   tcol<-rgb(0.1,0.1,0.25)
if(is.null(bgcol) & wall==0)   bgcol<-1
if(is.null(bgcol) & wall!=0)   bgcol<-'#DBE0E8'
if(is.null(col)){
if(wall==0) col <-c('palegreen4', 'palegreen2', 'goldenrod1', 'tomato1')
if(wall!=0) col <-c('green4', 'green', 'orange', 'red')
}
if(is.null(alpha)){
if(wall==0)  alpha<-1
if(wall!=0)  alpha<-0.5
}
col<-a.alpha.ade(col, alpha)
lcol2<-lcol
if(is.null(lcol) & wall!=4) lcol<- rep(tcol, 4)
if(is.null(lcol) & wall==4) lcol<- rep(tcol, 4)

if(is.null(lcol2) & wall!=4) lcol2<- rep(tcol, 4)
if(is.null(lcol2) & wall==4) lcol2<- rep(rgb(1,1,1), 4)


linesart<-lty
if(is.null(linesart)) linesart<-c(1,2,3,4)

######################################



############################
# Calcs
if(is.null(xlim))   cuts<-seq(min(x, na.rm=TRUE), max(x, na.rm=TRUE),   length.out=nints)
if(!is.null(xlim))  cuts<-seq(xlim[1], xlim[2], length.out=nints)
if(is.null(vnames)) vnames<-c('Sensitivity', 'Specificity', 'PPV', 'NPV')
if(is.null(xlim))  xlim<-  c(min(cuts, na.rm=TRUE), max(cuts, na.rm=TRUE))
xlim[2]<-xlim[2]+diff(xlim)/30

M<-NULL
for(k in 1:nints){
M<-rbind(M, a.diags.ade(as.numeric(x>cuts[k]) , y))
}
N<-dim(M)[1]

if(youden){
M2<-NULL
cuts2<-unique(x)
for(k in 1:length(cuts2)){
M2<-rbind(M2, a.diags.ade(as.numeric(x>cuts2[k]) , y))
}

yid<-which(M2[ ,which(colnames(M2)=='youden')]==max(M2[ ,which(colnames(M2)=='youden')], na.rm=TRUE) )
yid<-round(mean(yid[1], yid[length(yid)]))
if(M2[ ,which(colnames(M2)=='youden')][yid]<=0)  youden<-FALSE
}



################################################################################
#wall=0
if(wall==0){
par(col.axis=tcol)
par(col.lab=tcol)
par(col.main=tcol)


plot(0,0,  ylim=c(-0.025,1.05), col=rgb(0,0,0,0), xlim=xlim, axes=FALSE, ylab=ylab, xlab=xlab, main=main)
axis(2, at=seq(0, 1, 0.1), labels=seq(0, 1, 0.1)*100, col.ticks=bgcol, las=1)
if(is.null(xticks)) axis(1, col.ticks=bgcol)
if(!is.null(xticks)){
if(length(xticks)==1) axis(1, at=pretty(c(xlim[1], xlim[2]-diff(xlim)/15), xticks), col.ticks=bgcol)
if(length(xticks)>1)  axis(1, at=xticks  , col.ticks=bgcol)
}

polygon(c(cuts, cuts[N:1]), c(rep(0, N), M[N:1, 1]), col=col[1], border=rgb(0,0,0,0.1))
polygon(c(cuts, cuts[N:1]), c(M[, 1], M[N:1, 1]+M[N:1, 2]), col=col[2], border=rgb(0,0,0,0.1))
polygon(c(cuts, cuts[N:1]), c(M[, 1]+M[, 2], M[N:1, 1]+M[N:1, 2]+M[N:1, 3]), col=col[3], border=rgb(0,0,0,0.1))
polygon(c(cuts, cuts[N:1]), c(M[, 1]+M[, 2]+M[, 3], M[N:1, 1]+M[N:1, 2]+M[N:1, 3]+M[N:1, 4]), col=col[4], border=rgb(0,0,0,0.1))

legend(par('usr')[2], 0.95,  legend=c('FN'), horiz=FALSE, fill=col[4], xpd=TRUE, box.col=rgb(0,0,0,0), adj=c(1.5,2), xjust=0.55, text.col=tcol)
legend(par('usr')[2], 0.7,   legend=c('FP'), horiz=FALSE, fill=col[3], xpd=TRUE, box.col=rgb(0,0,0,0), adj=c(1.5,2), xjust=0.55, text.col=tcol)
legend(par('usr')[2], 0.45,  legend=c('TN'), horiz=FALSE, fill=col[2], xpd=TRUE, box.col=rgb(0,0,0,0), adj=c(1.5,2), xjust=0.55, text.col=tcol)
legend(par('usr')[2], 0.2,   legend=c('TP'), horiz=FALSE, fill=col[1], xpd=TRUE, box.col=rgb(0,0,0,0), adj=c(1.5,2), xjust=0.55, text.col=tcol)

legend('top', legend=vnames[stats], horiz=TRUE, lwd=2, lty=linesart, box.col=rgb(0,0,0,0), col=lcol, text.col=tcol )

if(youden) segments(cuts2[yid], -0.01 ,cuts2[yid], 1, lwd=lwd, col=2 )
if(youden) text(cuts2[yid] , 0, labels='Youden-Index', adj=c(0.5, 1.5), col=2)


if(!is.null(cutoffs)) segments(cutoffs, -0.01 ,cutoffs, 1, lwd=lwd, col=a.alpha.ade(tcol, 1) )
if(!is.null(cutoffs)) text(cutoffs , 0, labels=cutnames, adj=c(0.5, 1.5), col=tcol)

for(i in 1:length(stats)){
points(cuts,  M[, (stats[i]+4)], type='l', col=lcol[i], lwd=lwd, lty=linesart[i])
}

box(col=bgcol)
}
################################################################################



################################################################################
#wall=1
if(wall==1){
par(col.axis=tcol)
par(col.lab=tcol)
par(col.main=tcol)


plot(0,0,  ylim=c(-0.025,1.05), col=rgb(0,0,0,0), xlim=xlim, axes=FALSE, ylab=ylab, xlab=xlab, main=main)
polygon( c(par('usr')[c(1,1,2,2)]), par('usr')[c(3,4,4,3)], col=bgcol, border=FALSE)
if(is.null(xticks)) a1<-axis(1, col.ticks=tcol)
if(!is.null(xticks)){
if(length(xticks)==1) a1<-axis(1, at=pretty(c(xlim[1], xlim[2]-diff(xlim)/15), xticks), col.ticks=tcol)
if(length(xticks)>1)  a1<-axis(1, at=xticks  , col.ticks=tcol)
}
a2<-axis(2, at=seq(0, 1, 0.1), labels=seq(0, 1, 0.1)*100, col.ticks=tcol, las=1)
segments(par('usr')[1], a2, xlim[2]-diff(xlim)/30, a2, col=rgb(1,1,1))


polygon(c(cuts, cuts[N:1]), c(rep(0, N), M[N:1, 1]), col=col[1], border=rgb(0,0,0,0.1))
polygon(c(cuts, cuts[N:1]), c(M[, 1], M[N:1, 1]+M[N:1, 2]), col=col[2], border=rgb(0,0,0,0.1))
polygon(c(cuts, cuts[N:1]), c(M[, 1]+M[, 2], M[N:1, 1]+M[N:1, 2]+M[N:1, 3]), col=col[3], border=rgb(0,0,0,0.1))
polygon(c(cuts, cuts[N:1]), c(M[, 1]+M[, 2]+M[, 3], M[N:1, 1]+M[N:1, 2]+M[N:1, 3]+M[N:1, 4]), col=col[4], border=rgb(0,0,0,0.1))

legend(par('usr')[2], 0.95,  legend=c('FN'), horiz=FALSE, fill=col[4], xpd=FALSE, box.col=rgb(0,0,0,0), adj=c(1.5,2), xjust=0.55, text.col=tcol, bg=bgcol)
legend(par('usr')[2], 0.7,   legend=c('FP'), horiz=FALSE, fill=col[3], xpd=FALSE, box.col=rgb(0,0,0,0), adj=c(1.5,2), xjust=0.55, text.col=tcol, bg=bgcol)
legend(par('usr')[2], 0.45,  legend=c('TN'), horiz=FALSE, fill=col[2], xpd=FALSE, box.col=rgb(0,0,0,0), adj=c(1.5,2), xjust=0.55, text.col=tcol, bg=bgcol)
legend(par('usr')[2], 0.2,   legend=c('TP'), horiz=FALSE, fill=col[1], xpd=FALSE, box.col=rgb(0,0,0,0), adj=c(1.5,2), xjust=0.55, text.col=tcol, bg=bgcol)

legend('top', legend=vnames[stats], horiz=TRUE, lwd=2, lty=linesart, box.col=rgb(0,0,0,0), col=lcol, text.col=tcol, bg=bgcol)
abline(h=1, v=xlim[2]-diff(xlim)/30, col=rgb(1,1,1), lwd=2)


if(youden) segments(cuts2[yid], -0.01 ,cuts2[yid], 1, lwd=lwd, col=2 )
if(youden) text(cuts2[yid] , 0, labels='Youden-Index', adj=c(0.5, 1.5), col=2)


if(!is.null(cutoffs)) segments(cutoffs, -0.01 ,cutoffs, 1, lwd=lwd, col=a.alpha.ade(tcol, 0.5) )
if(!is.null(cutoffs)) text(cutoffs , 0, labels=cutnames, adj=c(0.5, 1.5), col=tcol)

for(i in 1:length(stats)){
points(cuts,  M[, (stats[i]+4)], type='l', col=lcol[i], lwd=lwd, lty=linesart[i])
}

box(col=rgb(1,1,1))
}
################################################################################



################################################################################
#wall=2
if(wall==2){
par(col.axis=tcol)
par(col.lab=tcol)
par(col.main=tcol)


plot(0,0,  ylim=c(-0.025,1.05), col=rgb(0,0,0,0), xlim=xlim, axes=FALSE, ylab=ylab, xlab=xlab, main=main)
if(is.null(xticks)) a1<-axis(1, col.ticks=a.coladd.ade(bgcol, -75))
if(!is.null(xticks)){
if(length(xticks)==1) a1<-axis(1, at=pretty(c(xlim[1], xlim[2]-diff(xlim)/15), xticks), col.ticks=a.coladd.ade(bgcol, -75))
if(length(xticks)>1)  a1<-axis(1, at=xticks  , col.ticks=a.coladd.ade(bgcol, -75))
}
a2<-axis(2, at=seq(0, 1, 0.1), labels=seq(0, 1, 0.1)*100, col.ticks=a.coladd.ade(bgcol, -75), las=1)
segments(par('usr')[1], a2, xlim[2]-diff(xlim)/30, a2, col=bgcol)


polygon(c(cuts, cuts[N:1]), c(rep(0, N), M[N:1, 1]), col=col[1], border=rgb(0,0,0,0.1))
polygon(c(cuts, cuts[N:1]), c(M[, 1], M[N:1, 1]+M[N:1, 2]), col=col[2], border=rgb(0,0,0,0.1))
polygon(c(cuts, cuts[N:1]), c(M[, 1]+M[, 2], M[N:1, 1]+M[N:1, 2]+M[N:1, 3]), col=col[3], border=rgb(0,0,0,0.1))
polygon(c(cuts, cuts[N:1]), c(M[, 1]+M[, 2]+M[, 3], M[N:1, 1]+M[N:1, 2]+M[N:1, 3]+M[N:1, 4]), col=col[4], border=rgb(0,0,0,0.1))

legend(par('usr')[2], 0.95,  legend=c('FN'), horiz=FALSE, fill=col[4], xpd=TRUE, box.col=rgb(0,0,0,0), adj=c(1.5,2), xjust=0.55, text.col=tcol)
legend(par('usr')[2], 0.7,   legend=c('FP'), horiz=FALSE, fill=col[3], xpd=TRUE, box.col=rgb(0,0,0,0), adj=c(1.5,2), xjust=0.55, text.col=tcol)
legend(par('usr')[2], 0.45,  legend=c('TN'), horiz=FALSE, fill=col[2], xpd=TRUE, box.col=rgb(0,0,0,0), adj=c(1.5,2), xjust=0.55, text.col=tcol)
legend(par('usr')[2], 0.2,   legend=c('TP'), horiz=FALSE, fill=col[1], xpd=TRUE, box.col=rgb(0,0,0,0), adj=c(1.5,2), xjust=0.55, text.col=tcol)

legend('top', legend=vnames[stats], horiz=TRUE, lwd=2, lty=linesart, box.col=rgb(0,0,0,0), col=lcol, text.col=tcol )
abline(h=1, v=xlim[2]-diff(xlim)/30, col=a.coladd.ade(bgcol, -75), lwd=1)


if(youden) segments(cuts2[yid], -0.01 ,cuts2[yid], 1, lwd=lwd, col=2 )
if(youden) text(cuts2[yid] , 0, labels='Youden-Index', adj=c(0.5, 1.5), col=2)


if(!is.null(cutoffs)) segments(cutoffs, -0.01 ,cutoffs, 1, lwd=lwd, col=a.alpha.ade(tcol, 0.5) )
if(!is.null(cutoffs)) text(cutoffs , 0, labels=cutnames, adj=c(0.5, 1.5), col=tcol)

for(i in 1:length(stats)){
points(cuts,  M[, (stats[i]+4)], type='l', col=lcol[i], lwd=lwd, lty=linesart[i])
}

box(col=a.coladd.ade(bgcol, -75))
}
################################################################################



################################################################################
#wall=3
if(wall==3){
par(col.axis=tcol)
par(col.lab=tcol)
par(col.main=tcol)

plot(0,0,  ylim=c(-0.025,1.05), col=rgb(0,0,0,0), xlim=xlim, axes=FALSE, ylab=ylab, xlab=xlab, main=main)
polygon( c(par('usr')[c(1,1)], xlim[2]-diff(xlim)/30, xlim[2]-diff(xlim)/30), c(par('usr')[3],1,1,par('usr')[3]), col=bgcol, border=FALSE)
rect(xlim[2]-diff(xlim)/30, 1, par('usr')[2], par('usr')[4], col=bgcol, border=FALSE)
if(is.null(xticks)) a1<-axis(1, col.ticks=a.coladd.ade(bgcol, -50))
if(!is.null(xticks)){
if(length(xticks)==1) a1<-axis(1, at=pretty(c(xlim[1], xlim[2]-diff(xlim)/15), xticks), col.ticks=a.coladd.ade(bgcol, -50))
if(length(xticks)>1)  a1<-axis(1, at=xticks  , col.ticks=a.coladd.ade(bgcol, -50))
}
a2<-axis(2, at=seq(0, 1, 0.1), labels=seq(0, 1, 0.1)*100, col.ticks=a.coladd.ade(bgcol, -50), las=1)
segments(par('usr')[1], a2, xlim[2]-diff(xlim)/30, a2, col=a.coladd.ade(bgcol, -50))

polygon(c(cuts, cuts[N:1]), c(rep(0, N), M[N:1, 1]), col=col[1], border=rgb(0,0,0,0.1))
polygon(c(cuts, cuts[N:1]), c(M[, 1], M[N:1, 1]+M[N:1, 2]), col=col[2], border=rgb(0,0,0,0.1))
polygon(c(cuts, cuts[N:1]), c(M[, 1]+M[, 2], M[N:1, 1]+M[N:1, 2]+M[N:1, 3]), col=col[3], border=rgb(0,0,0,0.1))
polygon(c(cuts, cuts[N:1]), c(M[, 1]+M[, 2]+M[, 3], M[N:1, 1]+M[N:1, 2]+M[N:1, 3]+M[N:1, 4]), col=col[4], border=rgb(0,0,0,0.1))

legend(par('usr')[2], 0.95,  legend=c('FN'), horiz=FALSE, fill=col[4], xpd=TRUE, box.col=rgb(0,0,0,0), adj=c(1.5,2), xjust=0.55, text.col=tcol)
legend(par('usr')[2], 0.7,   legend=c('FP'), horiz=FALSE, fill=col[3], xpd=TRUE, box.col=rgb(0,0,0,0), adj=c(1.5,2), xjust=0.55, text.col=tcol)
legend(par('usr')[2], 0.45,  legend=c('TN'), horiz=FALSE, fill=col[2], xpd=TRUE, box.col=rgb(0,0,0,0), adj=c(1.5,2), xjust=0.55, text.col=tcol)
legend(par('usr')[2], 0.2,   legend=c('TP'), horiz=FALSE, fill=col[1], xpd=TRUE, box.col=rgb(0,0,0,0), adj=c(1.5,2), xjust=0.55, text.col=tcol)

legend('top', legend=vnames[stats], horiz=TRUE, lwd=2, lty=linesart, box.col=rgb(0,0,0,0), col=lcol, text.col=tcol )
abline(h=1, v=xlim[2]-diff(xlim)/30, col=a.coladd.ade(bgcol, -75), lwd=1)


if(youden) segments(cuts2[yid], -0.01 ,cuts2[yid], 1, lwd=lwd, col=2 )
if(youden) text(cuts2[yid] , 0, labels='Youden-Index', adj=c(0.5, 1.5), col=2)


if(!is.null(cutoffs)) segments(cutoffs, -0.01 ,cutoffs, 1, lwd=lwd, col=a.alpha.ade(tcol, 0.5) )
if(!is.null(cutoffs)) text(cutoffs , 0, labels=cutnames, adj=c(0.5, 1.5), col=tcol)

for(i in 1:length(stats)){
points(cuts,  M[, (stats[i]+4)], type='l', col=lcol[i], lwd=lwd, lty=linesart[i])
}

box(col=a.coladd.ade(bgcol, -50))
}
################################################################################




################################################################################
#wall=4
if(wall==4){
par(col.axis=tcol)
par(col.lab=tcol)
par(col.main=tcol)
par(font=2)

#newmar<-rep(0, 4)
#oldmar<-par('mar')
#if(oldmar[2]<4.85) newmar[2]<- 4.85 - oldmar[2]
#if(oldmar[4]>1.3 & par('mar')[4]<=2.1) newmar[4]<- 1.3-oldmar[4]
#par(mar=(oldmar+newmar))

newmai<-rep(0, 4)
oldmai<-par('mai')
if(oldmai[2]<1) newmai[2]<- 1 - oldmai[2]
if(oldmai[4]>0.3 & oldmai[4]<=0.42) newmai[4]<- 0.3-oldmai[4]
par(mai=(oldmai+newmai))



plot(0,0,  ylim=c(-0.025,1.05), col=rgb(0,0,0,0), xlim=xlim, axes=FALSE, ylab='', xlab='', main='')
polygon( c(par('usr')[c(1,1)], xlim[2]-diff(xlim)/30, xlim[2]-diff(xlim)/30), c(0,1,1,0), col=bgcol, border=FALSE)
if(is.null(xticks)) a1<-axis(1, col.ticks=tcol)
if(!is.null(xticks)){
if(length(xticks)==1) a1<-axis(1, at=pretty(c(xlim[1], xlim[2]-diff(xlim)/15), xticks), col.ticks=tcol)
if(length(xticks)>1)  a1<-axis(1, at=xticks  , col.ticks=tcol)
}
a2<-axis(2, at=seq(0, 1, 0.1), labels=seq(0, 1, 0.1)*100, col.ticks=tcol, las=1)
segments(par('usr')[1], a2, xlim[2]-diff(xlim)/30, a2, col=tcol)

polygon(c(cuts, cuts[N:1]), c(rep(0, N), M[N:1, 1]), col=col[1], border=rgb(0,0,0,0.1))
polygon(c(cuts, cuts[N:1]), c(M[, 1], M[N:1, 1]+M[N:1, 2]), col=col[2], border=rgb(0,0,0,0.1))
polygon(c(cuts, cuts[N:1]), c(M[, 1]+M[, 2], M[N:1, 1]+M[N:1, 2]+M[N:1, 3]), col=col[3], border=rgb(0,0,0,0.1))
polygon(c(cuts, cuts[N:1]), c(M[, 1]+M[, 2]+M[, 3], M[N:1, 1]+M[N:1, 2]+M[N:1, 3]+M[N:1, 4]), col=col[4], border=rgb(0,0,0,0.1))

for(i in 1:length(stats)){
points(cuts,  M[, (stats[i]+4)], type='l', col=lcol[i], lwd=lwd, lty=linesart[i])
}

polygon( c(par('usr')[c(1,1)], xlim[2]-diff(xlim)/30, xlim[2]-diff(xlim)/30), c(1,par('usr')[c(4,4)],1), col=tcol, border=rgb(1,1,1))
polygon( c(xlim[2]-diff(xlim)/30, xlim[2]-diff(xlim)/30,  par('usr')[c(2,2)]), c(par('usr')[c(3,4,4,3)]), col=tcol, border=rgb(1,1,1))
polygon( c(par('usr')[c(1,1,2,2)]), c(par('usr')[3], 0,0, par('usr')[3]), col=tcol, border=rgb(1,1,1))

legend(par('usr')[2], 0.95,  legend=c('FN'), horiz=FALSE, fill=a.alpha.ade(col[4], alpha*1.5), xpd=FALSE, box.col=rgb(0,0,0,0), adj=c(1.5,2), xjust=0.55, text.col=rgb(1,1,1), border=rgb(1,1,1), bg=rgb(0,0,0,0))
legend(par('usr')[2], 0.7,   legend=c('FP'), horiz=FALSE, fill=a.alpha.ade(col[3], alpha*1.5), xpd=FALSE, box.col=rgb(0,0,0,0), adj=c(1.5,2), xjust=0.55, text.col=rgb(1,1,1), border=rgb(1,1,1), bg=rgb(0,0,0,0))
legend(par('usr')[2], 0.45,  legend=c('TN'), horiz=FALSE, fill=a.alpha.ade(col[2], alpha*1.5), xpd=FALSE, box.col=rgb(0,0,0,0), adj=c(1.5,2), xjust=0.55, text.col=rgb(1,1,1), border=rgb(1,1,1), bg=rgb(0,0,0,0))
legend(par('usr')[2], 0.2,   legend=c('TP'), horiz=FALSE, fill=a.alpha.ade(col[1], alpha*1.5), xpd=FALSE, box.col=rgb(0,0,0,0), adj=c(1.5,2), xjust=0.55, text.col=rgb(1,1,1), border=rgb(1,1,1), bg=rgb(0,0,0,0))

legend('top', legend=vnames[stats], horiz=TRUE, lwd=2, lty=linesart, box.col=rgb(0,0,0,0), col=lcol2, text.col=rgb(1,1,1) , border=rgb(1,1,1), bg=rgb(0,0,0,0))

par(xpd=TRUE)
polygon(a.glc(side=c(2,2,4,4), line=c(0,0,0,0)), a.glc(side=3, line=c(0, 2.75,  2.75, 0)), col=tcol, border=rgb(1,1,1))
if(ylab!='' & ylab!=' ') polygon( a.glc(side=2, line=c(4, 4, 2.5, 2.5)), a.glc(side=c(1, 3, 3, 1), line=0), col=bgcol, border=rgb(1,1,1))
if(xlab!='' & xlab!=' ') polygon( a.glc(side=c(2, 2, 4, 4), line=0),     a.glc(side=1, line=c(4, 2.5, 2.5, 4)), col=bgcol, border=rgb(1,1,1))
text(a.glc(side=0), a.glc(side=3, line=1),    labels=main, cex = 1.25, font=2, col=rgb(1,1,1), adj=c(0.5,0))
text(a.glc(side=0), a.glc(side=1, line=3.5), labels=xlab, cex = 1.1,  font=2, col=tcol, adj=c(0.5,0))
text(a.glc(side=2, line=3), a.glc(side=5),  labels=ylab, cex = 1.1,  font=2,  col=tcol, adj=c(0.5,0), srt=90)
par(xpd=FALSE)

if(youden) segments(cuts2[yid], -0.01 ,cuts2[yid], 1, lwd=lwd, col=2 )
if(youden) text(cuts2[yid] , 0, labels='Youden-Index', adj=c(0.5, 1.5), col=2)


if(!is.null(cutoffs)) segments(cutoffs, -0.01 ,cutoffs, 1, lwd=lwd, col=a.alpha.ade(tcol, 0.5) )
if(!is.null(cutoffs)) text(cutoffs , 0, labels=cutnames, adj=c(0.5, 1.5), col=rgb(1,1,1))


abline(h=1, v=xlim[2]-diff(xlim)/30, col=rgb(1,1,1), lwd=1)

box(col=rgb(1,1,1))
}
################################################################################



################################################################################
#wall=5
if(wall==5){
par(col.axis=tcol)
par(col.lab=tcol)
par(col.main=tcol)
#if(wall==5){
#newmar<-rep(0, 4)
#oldmar<-par('mar')
#if(oldmar[2]<4.85) newmar[2]<- 4.85 - oldmar[2]
#if(oldmar[3]>3.6 & par('mar')[3]<=4.1) newmar[3]<- 3.6-oldmar[3]
#if(oldmar[4]>1.3 & par('mar')[4]<=2.1) newmar[4]<- 1.3-oldmar[4]
#par(mar=(oldmar+newmar))
#}

newmai<-rep(0, 4)
oldmai<-par('mai')
if(oldmai[2]<1) newmai[2]<- 1 - oldmai[2]
if(oldmai[3]>0.75 & oldmai[3]<=0.82) newmai[3]<- 0.75-oldmai[3]
if(oldmai[4]>0.25 & oldmai[4]<=0.42) newmai[4]<- 0.25-oldmai[4]
par(mai=(oldmai+newmai))






plot(0,0,  ylim=c(-0.025,1.05), col=rgb(0,0,0,0), xlim=xlim, axes=FALSE, ylab='', xlab='', main='')
if(is.null(xticks)) a1<-axis(1, col.ticks=tcol)
if(!is.null(xticks)){
if(length(xticks)==1) a1<-axis(1, at=pretty(c(xlim[1], (xlim[2]-diff(xlim)/15)), xticks), col.ticks=tcol)
if(length(xticks)>1)  a1<-axis(1, at=xticks  , col.ticks=tcol)
}
a2<-axis(2, at=seq(0, 1, 0.1), labels=seq(0, 1, 0.1)*100, col.ticks=a.coladd.ade(bgcol, -75), las=1)
segments(par('usr')[1], a2, xlim[2]-diff(xlim)/30, a2, col=bgcol)


polygon(c(cuts, cuts[N:1]), c(rep(0, N), M[N:1, 1]), col=col[1], border=rgb(0,0,0,0.1))
polygon(c(cuts, cuts[N:1]), c(M[, 1], M[N:1, 1]+M[N:1, 2]), col=col[2], border=rgb(0,0,0,0.1))
polygon(c(cuts, cuts[N:1]), c(M[, 1]+M[, 2], M[N:1, 1]+M[N:1, 2]+M[N:1, 3]), col=col[3], border=rgb(0,0,0,0.1))
polygon(c(cuts, cuts[N:1]), c(M[, 1]+M[, 2]+M[, 3], M[N:1, 1]+M[N:1, 2]+M[N:1, 3]+M[N:1, 4]), col=col[4], border=rgb(0,0,0,0.1))

legend(par('usr')[2], 0.95,  legend=c('FN'), horiz=FALSE, fill=col[4], xpd=TRUE, box.col=rgb(0,0,0,0), adj=c(1.5,2), xjust=0.55, text.col=tcol)
legend(par('usr')[2], 0.7,   legend=c('FP'), horiz=FALSE, fill=col[3], xpd=TRUE, box.col=rgb(0,0,0,0), adj=c(1.5,2), xjust=0.55, text.col=tcol)
legend(par('usr')[2], 0.45,  legend=c('TN'), horiz=FALSE, fill=col[2], xpd=TRUE, box.col=rgb(0,0,0,0), adj=c(1.5,2), xjust=0.55, text.col=tcol)
legend(par('usr')[2], 0.2,   legend=c('TP'), horiz=FALSE, fill=col[1], xpd=TRUE, box.col=rgb(0,0,0,0), adj=c(1.5,2), xjust=0.55, text.col=tcol)


for(i in 1:length(stats)){
points(cuts,  M[, (stats[i]+4)], type='l', col=lcol[i], lwd=lwd, lty=linesart[i])
}


legend('top', legend=vnames[stats], horiz=TRUE, lwd=2, lty=linesart, box.col=rgb(0,0,0,0), col=lcol, text.col=tcol )
abline(h=1, v=xlim[2]-diff(xlim)/30, col=a.coladd.ade(bgcol, -75), lwd=1)


if(youden) segments(cuts2[yid], -0.01 ,cuts2[yid], 1, lwd=lwd, col=2 )
if(youden) text(cuts2[yid] , 0, labels='Youden-Index', adj=c(0.5, 1.5), col=2)


if(!is.null(cutoffs)) segments(cutoffs, -0.01 ,cutoffs, 1, lwd=lwd, col=a.alpha.ade(tcol, 0.5) )
if(!is.null(cutoffs)) text(cutoffs , 0, labels=cutnames, adj=c(0.5, 1.5), col=tcol)


par(xpd=TRUE)
polygon(a.glc(side=2, line=c(4.25, 4.25, 0, 0)), a.glc(side=3, line=c(0.6, 3, 3, 0.6)), col=bgcol,        border=tcol)
polygon(a.glc(side=c(2,2,4,4), line=c(0,0,0,0)), a.glc(side=3, line=c(0.6, 3, 3, 0.6)), col=rgb(1,1,1,0), border=tcol)
polygon(a.glc(side=4, line=c(0, 0 ,0.6, 0.6)),   a.glc(side=3, line=c(0.6, 3, 3, 0.6)), col=bgcol,        border=tcol)
polygon(a.glc(side=2, line=c(4.25, 4.25 ,3.65, 3.65)),  a.glc(side=c(1,3,3,1), line=c(2.6, 0.6, 0.6, 2.6)), col=bgcol,  border=tcol)
polygon(a.glc(side=4, line=c(0, 0 ,0.6, 0.6)), a.glc(side=c(1, 3, 3, 1), line=0), col=bgcol, border=tcol)
polygon(a.glc(side=2, line=c(4.25, 4.25, 0, 0)), a.glc(side=1, line=c(2.6, 4.5, 4.5, 2.6)), col=bgcol, border=tcol)
polygon(a.glc(side=c(2, 2, 4, 4), line=0), a.glc(side=1, line=c(2.6, 4.5, 4.5, 2.6)), col=rgb(1,1,1,0), border=tcol)
polygon(a.glc(side=4, line=c(0, 0, 0.6, 0.6)), a.glc(side=1, line=c(2.6, 4.5, 4.5, 2.6)), col=bgcol, border=tcol)
text(a.glc(side=0), a.glc(side=3, line=1.5),  labels=main, cex = 1.25, font=2, col=tcol, adj=c(0.5,0))
text(a.glc(side=0), a.glc(side=1, line=3.75), labels=xlab, cex = 1.1,  font=2, col=tcol, adj=c(0.5,0))
text(a.glc(side=2, line=2.5), a.glc(side=5), labels=ylab, cex = 1.1,   font=2,  col=tcol, adj=c(0.5,0), srt=90)
par(xpd=FALSE)
box(lwd=1, col=tcol)

}
################################################################################


################################################################################
#wall=6
if(wall==6){
par(col.axis=tcol)
par(col.lab=tcol)
par(col.main=tcol)


plot(0,0,  ylim=c(-0.025,1.05), col=rgb(0,0,0,0), xlim=xlim, axes=FALSE, ylab=ylab, xlab=xlab, main=main)
polygon( c(par('usr')[c(1,1,2,2)]), par('usr')[c(3,4,4,3)], col=bgcol, border=FALSE)
if(is.null(xticks)) a1<-axis(1, col.ticks=a.coladd.ade(bgcol, -35), lwd.ticks=3)
if(is.null(xticks)) a1<-axis(1, col.ticks=rgb(1,1,1), lwd.ticks=1)
if(!is.null(xticks)){
if(length(xticks)==1) a1<-axis(1, at=pretty(c(xlim[1], xlim[2]-diff(xlim)/15), xticks), col.ticks=a.coladd.ade(bgcol, -35), lwd.ticks=3)
if(length(xticks)==1) a1<-axis(1, at=pretty(c(xlim[1], xlim[2]-diff(xlim)/15), xticks), col.ticks=rgb(1,1,1), lwd.ticks=1)


if(length(xticks)>1)  a1<-axis(1, at=xticks  , col.ticks=a.coladd.ade(bgcol, -35), lwd.ticks=3)
if(length(xticks)>1)  a1<-axis(1, at=xticks  , col.ticks=rgb(1,1,1), lwd.ticks=1)
}
a2<-axis(2, at=seq(0, 1, 0.1), labels=seq(0, 1, 0.1)*100, col.ticks=a.coladd.ade(bgcol, -35), las=1, lwd.ticks=3)
a2<-axis(2, at=seq(0, 1, 0.1), labels=seq(0, 1, 0.1)*100, col.ticks=rgb(1,1,1), las=1, lwd.ticks=1)


segments(par('usr')[1], a2, xlim[2]-diff(xlim)/30, a2, col=a.coladd.ade(bgcol, -35), lwd=3)
segments(par('usr')[1], a2, xlim[2]-diff(xlim)/30, a2, col=rgb(1,1,1))


polygon(c(cuts, cuts[N:1]), c(rep(0, N), M[N:1, 1]), col=col[1], border=rgb(0,0,0,0.1))
polygon(c(cuts, cuts[N:1]), c(M[, 1], M[N:1, 1]+M[N:1, 2]), col=col[2], border=rgb(0,0,0,0.1))
polygon(c(cuts, cuts[N:1]), c(M[, 1]+M[, 2], M[N:1, 1]+M[N:1, 2]+M[N:1, 3]), col=col[3], border=rgb(0,0,0,0.1))
polygon(c(cuts, cuts[N:1]), c(M[, 1]+M[, 2]+M[, 3], M[N:1, 1]+M[N:1, 2]+M[N:1, 3]+M[N:1, 4]), col=col[4], border=rgb(0,0,0,0.1))

legend(par('usr')[2], 0.95,  legend=c('FN'), horiz=FALSE, border=a.coladd.ade(bgcol, -35), fill=col[4], xpd=FALSE, box.col=rgb(0,0,0,0), adj=c(1.5,2), xjust=0.55, text.col=tcol, bg=bgcol)
legend(par('usr')[2], 0.7,   legend=c('FP'), horiz=FALSE, border=a.coladd.ade(bgcol, -35), fill=col[3], xpd=FALSE, box.col=rgb(0,0,0,0), adj=c(1.5,2), xjust=0.55, text.col=tcol, bg=bgcol)
legend(par('usr')[2], 0.45,  legend=c('TN'), horiz=FALSE, border=a.coladd.ade(bgcol, -35), fill=col[2], xpd=FALSE, box.col=rgb(0,0,0,0), adj=c(1.5,2), xjust=0.55, text.col=tcol, bg=bgcol)
legend(par('usr')[2], 0.2,   legend=c('TP'), horiz=FALSE, border=a.coladd.ade(bgcol, -35), fill=col[1], xpd=FALSE, box.col=rgb(0,0,0,0), adj=c(1.5,2), xjust=0.55, text.col=tcol, bg=bgcol)

legend('top', legend=vnames[stats], horiz=TRUE, lwd=2, lty=linesart, box.col=rgb(0,0,0,0), col=lcol, text.col=tcol, bg=bgcol)
abline(h=1, v=xlim[2]-diff(xlim)/30, col=a.coladd.ade(bgcol, -35), lwd=3)
abline(h=1, v=xlim[2]-diff(xlim)/30, col=rgb(1,1,1), lwd=1)


if(youden) segments(cuts2[yid], -0.01 ,cuts2[yid], 1, lwd=lwd, col=2 )
if(youden) text(cuts2[yid] , 0, labels='Youden-Index', adj=c(0.5, 1.5), col=2)


if(!is.null(cutoffs)) segments(cutoffs, -0.01 ,cutoffs, 1, lwd=lwd, col=a.alpha.ade(tcol, 0.5) )
if(!is.null(cutoffs)) text(cutoffs , 0, labels=cutnames, adj=c(0.5, 1.5), col=tcol)

for(i in 1:length(stats)){
points(cuts,  M[, (stats[i]+4)], type='l', col=lcol[i], lwd=lwd, lty=linesart[i])
}

box(col=rgb(1,1,1), lwd=3)
box(col=a.coladd.ade(bgcol, -35), lwd=1)
}
################################################################################



}
