parallel.ade <-
function(vars,  vnames=NULL,  data=NULL, group=NULL, ylim=NULL , xlab=NULL, ylab=NULL, main=NULL, alpha = NULL, col=NULL, tcol=NULL,  bgcol=NULL, lcol=NULL, scale=FALSE, desc=FALSE, means=TRUE, legendon='top', wall=0){
if(any(par('mfg')!=c(1,1,1,1)) & any(par('mai') < c(1.02, 0.82, 0.82, 0.42))){
maidiff<-rep(0, 4)
norm<-c(1.02, 0.82, 0.82, 0.42)
maidiff[par('mai')<norm]<-  norm[par('mai')<norm] - par('mai')[par('mai')<norm]
par(mai=par('mai')+maidiff)
}
oldpar<-par(no.readonly =TRUE)
oldpar<-oldpar[-which(names(oldpar)%in%c('usr',  'plt', 'pin',   'fin', 'fig', 'mfg', 'mfcol', 'mfrow', 'omd', 'omi', 'oma'))]
on.exit(par(oldpar))


#####################################
# without Data.Frame or with
if(is.list(vars)){
x<- vars[[1]]
data<-as.data.frame(x)
newvars<-NULL
for(i in 1:length(vars)){
eval(parse(text=paste("data$var_",i,"<-vars[[i]]", sep='')))
newvars<-c(newvars, paste("var_",i, sep=''))
}
data$groupmy<-group
vars<-newvars
if(!is.null(group)) group<-'groupmy'
}
#####################################




for(i in 1:length(vars)){
data<- subset(data, !is.na(eval(parse(text=paste(vars[i], sep='')))))
}


if (is.null(vnames)) vnames <- vars
if(!is.null(group)) gg<- eval(parse(text=paste("data$",group,  sep='')))
if(!is.null(group)) ggf<-eval(parse(text=paste("data$",group,  sep='')))
if(!is.null(group)) ggf<-as.factor(ggf)
if(!is.null(group)) gg<-as.numeric(gg)




if(is.null(tcol)  & wall==0)   tcol<-1
if(is.null(tcol)  & wall!=0)   tcol<-rgb(0.1,0.1,0.25)
if(is.null(bgcol) & wall==0)   bgcol<-1
if(is.null(bgcol) & wall!=0)   bgcol<-'#DBE0E8'
if(is.null(alpha)){
alpha<-0.25
if(wall==0) alpha<-1
}

M<-NULL
for(i in 1:length(vars)){
M<-cbind(M, eval(parse(text=paste("data$",vars[i],  sep=''))))
}


runlines<-function(alcol, mcol=tcol, plot=TRUE){
################################################################################
################################################################################
# liniean


M<-NULL
for(i in 1:length(vars)){
M<-cbind(M, eval(parse(text=paste("data$",vars[i],  sep=''))))
}


rbc2 <-alcol

if(scale){
if(is.null(group))   yr<-range(M, na.rm=TRUE)
if(!is.null(group))  yr<-range(M[(!is.na(gg)), ], na.rm=TRUE)
for(i in 1:length(vars)){
if(is.null(group))  M[ , i]<-(M[, i]/(range(M[,i], na.rm=TRUE)[2]-range(M[,i], na.rm=TRUE)[1]))-(range(M[,i], na.rm=TRUE)[1]/(range(M[,i], na.rm=TRUE)[2]-range(M[,i], na.rm=TRUE)[1]))
if(!is.null(group)) M[ , i]<-(M[, i]/(range(M[(!is.na(gg)),i], na.rm=TRUE)[2]-range(M[(!is.na(gg)),i], na.rm=TRUE)[1]))-(range(M[(!is.na(gg)),i], na.rm=TRUE)[1]/(range(M[(!is.na(gg)),i], na.rm=TRUE)[2]-range(M[(!is.na(gg)),i], na.rm=TRUE)[1]))
ylim=(ylim/(ylim[2]-ylim[1]))-(ylim[1]/(ylim[2]-ylim[1]))
}
}

d<-dim(M)
x<-1:dim(M)[2]

if(is.null(group)){
if(is.null(col))   rbc<-rainbow(dim(M)[1], start=0.2, end=0, alpha=alpha)
if(!is.null(col))  rbc<-rep(col, round(d[1]/length(col))+1 )
o<-order(as.numeric(M[,1] ), decreasing = desc)
M<-M[o,]
}




if(!is.null(group)){

if(is.null(col)) col <- a.getcol.ade(nlevels(ggf))


if(!is.null(col)){
rbc<- a.alpha.ade(col, alpha)
rbc2<-a.alpha.ade(col, 1)
}

o<-order(as.numeric(gg), decreasing = F)
M<-M[o,]
gg<-gg[o]
}

if(!is.null(group)) if(min(gg, na.rm=TRUE)==0) gg<-gg+1

if(plot){
for(i in 1:(dim(M)[2]-1)){
for(j in 1:dim(M)[1]){
if(is.null(group))   segments(i,M[j,i],(i+1),M[j,(i+1)], col=rbc[j])
if(!is.null(group))  segments(i,M[j,i],(i+1),M[j,(i+1)], col=rbc[gg[j]])
}



# means linien
if(means){
# means
if(is.null(group) & is.null(lcol))  segments(i,mean(M[,i], na.rm=TRUE),(i+1),mean(M[,(i+1)], na.rm=TRUE), col=mcol, lwd=2)
if(is.null(group) & !is.null(lcol)) segments(i,mean(M[,i], na.rm=TRUE),(i+1),mean(M[,(i+1)], na.rm=TRUE), col=lcol, lwd=2)


if(!is.null(group)){


for(k in 1:length(rbc)){
MD<-cbind(M, gg)
x<-levels(factor(gg))[k]
MD<-as.data.frame(MD)
MD<-subset(MD, gg==x  )

segments(i,mean(MD[,i], na.rm=TRUE),(i+1),mean(MD[,(i+1)], na.rm=TRUE), col=a.alpha.ade(a.coladd.ade(rbc[k], -75), 1), lwd=3)

}
}

}

abline(v=i, col=alcol)
}
}


if(is.null(group))   yrange<-range(M, na.rm=TRUE)
if(!is.null(group))  yrange<-range(M[(!is.na(gg)), ], na.rm=TRUE)


if(plot) abline(v=(i+1), col=alcol)
return(list(rbc2, yrange[1],  yrange[2]))
################################################################################


}
################################################################################
################################################################################
################################################################################


################################################################################
################################################################################
# Plot Umgebung
################################################################################
if(wall==0){
par(col.axis=tcol)
par(col.lab=tcol)
par(col.main=tcol)
#######################


if(!is.null(ylim)) yrange<-ylim
xrange<-c(1, length(vars))
out<-runlines(a.coladd.ade(bgcol, -50), tcol, plot=FALSE)
rbc2<-out[[1]]
if(is.null(ylim)) yrange<-c(out[[2]], out[[3]])

########################
yrange[2]<-yrange[2]+yrange[2]/10
 

plot(0, 0 , xlim=xrange, ylim=yrange, axes=FALSE, xlab=xlab, ylab=ylab, main=main, col=rgb(1,1,1,0))

axis(2 , at=pretty(c(yrange[1],yrange[2])) , col.ticks=bgcol)

a1<-axis(1,  at=1:length(vars), labels=vnames, col.ticks=bgcol)


runlines(bgcol)

if(!is.null(group)) legend(legendon, legend=levels(ggf), text.col=tcol, box.col=bgcol,  fill=rbc2, border=bgcol,  yjust=0.425, horiz=TRUE, bg=rgb(1,1,1, 0.99))
box(col=bgcol)
}
################################################################################



################################################################################
if(wall==1){
par(col.axis=tcol)
par(col.lab=tcol)
par(col.main=tcol)
#######################

if(!is.null(ylim)) yrange=ylim
xrange<-c(1, length(vars))
out<-runlines(a.coladd.ade(bgcol, -50), tcol, plot=FALSE)
rbc2<-out[[1]]
if(is.null(ylim)) yrange<-c(out[[2]], out[[3]])

########################
yrange[2]<-yrange[2]+yrange[2]/10


plot(0, 0 , xlim=xrange, ylim=yrange, axes=FALSE, xlab=xlab, ylab=ylab, main=main, col=rgb(1,1,1,0))
polygon( c(par('usr')[c(1,1,2,2)]), par('usr')[c(3,4,4,3)], col=bgcol, border=FALSE)
a2<-axis(2 , at=pretty(c(yrange[1],yrange[2])) ,  col.ticks=bgcol)
a1<-axis(1, at=1:length(vars), labels=vnames, col.ticks=bgcol)
abline(v=a1, h=a2, lty=1, col=rgb(1,1,1), lwd=1)


runlines(rgb(1,1,1), tcol)

if(!is.null(group)) legend(legendon, legend=levels(ggf), text.col=tcol, box.col=rgb(1,1,1), fill=rbc2, border=tcol,  box.lwd=2,  horiz=TRUE, bg=bgcol)
box(col=rgb(1,1,1))
}
################################################################################



################################################################################
if(wall==2){
par(col.axis=tcol)
par(col.lab=tcol)
par(col.main=tcol)
#######################
if(!is.null(ylim)) yrange=ylim
xrange<-c(1, length(vars))
out<-runlines(a.coladd.ade(bgcol, -50), tcol, plot=FALSE)
rbc2<-out[[1]]
if(is.null(ylim))  yrange<-c(out[[2]], out[[3]])

########################
yrange[2]<-yrange[2]+yrange[2]/10


plot(0, 0 , xlim=xrange, ylim=yrange, axes=FALSE, xlab=xlab, ylab=ylab, main=main, col=rgb(1,1,1,0))


a2<-axis(2 , at=pretty(c(yrange[1],yrange[2])) ,  col.ticks=bgcol)
a1<-axis(1, at=1:length(vars), labels=vnames, col.ticks=bgcol)
abline(v=a1, h=a2, lty=1, col=bgcol, lwd=1)

runlines(bgcol, tcol)


if(!is.null(group)) legend(legendon, legend=levels(ggf), text.col=tcol, box.col=a.coladd.ade(bgcol, -75), fill=rbc2,  border=a.coladd.ade(bgcol, -75),  box.lwd=1,  horiz=TRUE, bg=rgb(1,1,1))
box(col=a.coladd.ade(bgcol, -75))
}
################################################################################



################################################################################
if(wall==3){
par(col.axis=tcol)
par(col.lab=tcol)
par(col.main=tcol)


#######################
if(!is.null(ylim)) yrange=ylim
xrange<-c(1, length(vars))

out<-runlines(a.coladd.ade(bgcol, -50), tcol, plot=FALSE)
rbc2<-out[[1]]
if(is.null(ylim))  yrange<-c(out[[2]], out[[3]])
########################
yrange[2]<-yrange[2]+yrange[2]/10

plot(0, 0 , ylim=yrange, axes=FALSE, xlab=xlab, ylab=ylab, main=main, xlim=xrange, col=rgb(1,1,1,0))
polygon( c(par('usr')[c(1,1,2,2)]), par('usr')[c(3,4,4,3)], col=bgcol, border=FALSE)
a1<-axis(1, at=1:length(vars), labels=vnames, col.ticks=bgcol)
a2<-axis(2 , at=pretty(c(yrange[1],yrange[2])) , col.ticks=bgcol)
abline(v=a1, h=a2, lty=1, col=a.coladd.ade(bgcol, -50), lwd=1)

out<-runlines(a.coladd.ade(bgcol, -50), tcol, plot=TRUE)


if(!is.null(group)) legend(legendon, legend=levels(ggf), text.col=tcol, fill=rbc2,  box.col=a.coladd.ade(bgcol, -75), border=a.coladd.ade(bgcol, -75),  box.lwd=1, horiz=TRUE, bg=rgb(1,1,1))
box(col=a.coladd.ade(bgcol, -75))
}
################################################################################


################################################################################
if(wall==4){
par(col.axis=tcol)
par(col.lab=tcol)
par(col.main=rgb(1,1,1))
par(font=2)
#######################


if(!is.null(ylim)) yrange=ylim
xrange<-c(1, length(vars))

out<-runlines(a.coladd.ade(bgcol, -50), tcol, plot=FALSE)
rbc2<-out[[1]]
if(is.null(ylim))  yrange<-c(out[[2]], out[[3]])

########################
yrange[2]<-yrange[2]+yrange[2]/10




plot(0, 0 , xlim=xrange, ylim=yrange, axes=FALSE, xlab='', ylab='', main='', col=rgb(1,1,1,0))
polygon( c(par('usr')[c(1,1,2,2)]), par('usr')[c(3,4,4,3)], col=bgcol, border=FALSE)
a2<-axis(2 , at=pretty(c(yrange[1],yrange[2])) ,  col.ticks=bgcol)
a1<-axis(1, at=1:length(vars), labels=vnames, col.ticks=bgcol)
abline(v=a1, h=a2, lty=1, col=rgb(1,1,1), lwd=1)



# Outer
par(xpd=TRUE)
polygon(a.glc(side=c(2,2,4,4), line=c(0,0,0,0)), a.glc(side=3, line=c(0, 2.75,  2.75, 0)), col=tcol, border=rgb(1,1,1))
if(!is.null(ylab)) if(ylab!='' & ylab!=' ') polygon( a.glc(side=2, line=c(3.5, 3.5, 2, 2)), a.glc(side=c(1, 3, 3, 1), line=0), col=bgcol, border=rgb(1,1,1))
if(!is.null(xlab)) if(xlab!='' & xlab!=' ') polygon( a.glc(side=c(2, 2, 4, 4), line=0),     a.glc(side=1, line=c(4, 2.5, 2.5, 4)), col=bgcol, border=rgb(1,1,1))
text(a.glc(side=0), a.glc(side=3, line=1),    labels=main, cex = 1.25, font=2, col=rgb(1,1,1), adj=c(0.5,0))
text(a.glc(side=0), a.glc(side=1, line=3.5), labels=xlab, cex = 1.1,  font=2, col=tcol, adj=c(0.5,0))
text(a.glc(side=2, line=2.5), a.glc(side=5),  labels=ylab, cex = 1.1,  font=2,  col=tcol, adj=c(0.5,0), srt=90)
par(xpd=FALSE)
box(col=rgb(1,1,1))

runlines(rgb(1,1,1), tcol)

if(!is.null(group)) legend(legendon, legend=levels(ggf), text.col=rgb(1,1,1), fill=rbc2,   box.col=rgb(1,1,1), border=rgb(1,1,1), lty=0,box.lwd=1, pt.cex=2, col=rbc2, horiz=TRUE, bg=tcol, text.width=max(strwidth(levels(ggf),font = 2)))
box(col=rgb(1,1,1))
}
################################################################################

################################################################################
if(wall==5){
par(col.axis=tcol)
par(col.lab=tcol)
par(col.main=tcol)
par(font=2)
newmai<-rep(0, 4)
oldmai<-par('mai')
if(oldmai[2]<1) newmai[2]<- 1 - oldmai[2]
if(oldmai[3]>0.75 & oldmai[3]<=0.82) newmai[3]<- 0.75-oldmai[3]
if(oldmai[4]>0.25 & oldmai[4]<=0.42) newmai[4]<- 0.25-oldmai[4]
par(mai=(oldmai+newmai))


#######################

if(!is.null(ylim)) yrange=ylim
xrange<-c(1, length(vars))

out<-runlines(a.coladd.ade(bgcol, -50), tcol, plot=FALSE)
rbc2<-out[[1]]
if(is.null(ylim)) yrange<-c(out[[2]], out[[3]])

########################
yrange[2]<-yrange[2]+yrange[2]/10


plot(0, 0 , xlim=xrange, ylim=yrange, axes=FALSE, xlab='', ylab='', main='', col=rgb(1,1,1,0))
a2<-axis(2 , at=pretty(c(yrange[1],yrange[2])) , col.ticks=bgcol)
a1<-axis(1, at=1:length(vars), labels=vnames, col.ticks=bgcol)
abline(v=a1, h=a2, lty=1, col=rgb(1,1,1), lwd=1)

# Outer
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
box(col=rgb(1,1,1))
runlines(tcol, tcol)

if(!is.null(group)) legend(legendon, legend=levels(ggf), text.col=tcol, fill=rbc2,   box.col=tcol, border=tcol, lty=0,box.lwd=1, pt.cex=2, col=rbc2,  horiz=TRUE, bg=rgb(1,1,1), text.width=max(strwidth(levels(ggf),font = 2)))
box(col=tcol)
}
################################################################################


################################################################################
if(wall==6){
par(col.axis=tcol)
par(col.lab=tcol)
par(col.main=tcol)
#######################

if(!is.null(ylim)) yrange=ylim
xrange<-c(1, length(vars))
out<-runlines(a.coladd.ade(bgcol, -50), tcol, plot=FALSE)
rbc2<-out[[1]]
if(is.null(ylim)) yrange<-c(out[[2]], out[[3]])

########################
yrange[2]<-yrange[2]+yrange[2]/10
#print(pretty(c(yrange[1],yrange[2])))

plot(0, 0 , xlim=xrange, ylim=yrange, axes=FALSE, xlab=xlab, ylab=ylab, main=main, col=rgb(1,1,1,0))
polygon( c(par('usr')[c(1,1,2,2)]), par('usr')[c(3,4,4,3)], col=bgcol, border=NA)
a2<-axis(2 , at=pretty(c(yrange[1],yrange[2])) , col.ticks=a.coladd.ade(bgcol, -35), lwd.ticks=3)
a2<-axis(2 , at=pretty(c(yrange[1],yrange[2])) , col.ticks=rgb(1,1,1), lwd.ticks=1)

a1<-axis(1, at=1:length(vars), labels=vnames, col.ticks=a.coladd.ade(bgcol, -35), lwd.ticks=3)
a1<-axis(1, at=1:length(vars), labels=vnames, col.ticks=rgb(1,1,1), lwd.ticks=1)

abline(v=a1, h=a2, lty=1, col=a.coladd.ade(bgcol, -35), lwd=3)
abline(v=a1, h=a2, lty=1, col=rgb(1,1,1), lwd=1)

runlines(rgb(1,1,1), tcol)

if(!is.null(group)) legend(legendon, legend=levels(ggf), text.col=tcol, box.col=rgb(1,1,1),               fill=rbc2, border=tcol,  box.lwd=3,  horiz=TRUE, bg=bgcol)
if(!is.null(group)) legend(legendon, legend=levels(ggf), text.col=tcol, box.col=a.coladd.ade(bgcol, -35), fill=rbc2, border=tcol,  box.lwd=1,  horiz=TRUE, bg=rgb(0,0,0,0))

box(lwd=3, col=rgb(1,1,1))
box(lwd=1, col=a.coladd.ade(bgcol, -35))

}
################################################################################

}
