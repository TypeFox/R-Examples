box.plot.wtd <-
function(x, group=NULL, group2=NULL, w=NULL, data=NULL, vnames=NULL, main=NULL, xlab=NULL, ylab=NULL, ylim=NULL, yticks=NULL, col=NULL, tcol=NULL, bgcol=NULL, lcol=NULL, pdigs=4, alpha=NULL, cex=1, cex.axis=1, lwd=2, h=NULL, lty=2, varwidth=TRUE, means=FALSE, count=TRUE, zylinder=FALSE, outlier=TRUE, wall=0, type='box'){
library(Hmisc)
if(any(par('mfg')!=c(1,1,1,1)) & any(par('mai') < c(1.02, 0.82, 0.82, 0.42))){
maidiff<-rep(0, 4)
norm<-c(1.02, 0.82, 0.82, 0.42)
maidiff[par('mai')<norm]<-  norm[par('mai')<norm] - par('mai')[par('mai')<norm]
par(mai=par('mai')+maidiff)
}
oldpar<-par(no.readonly =TRUE)
pnames<-names(oldpar)
oldpar<-oldpar[-which(names(oldpar)%in%c('usr', 'plt',   'pin', 'fin', 'fig', 'mfg', 'mfcol', 'mfrow', 'omd', 'omi', 'oma'))]
on.exit(par(oldpar))
test=FALSE

wtd.sd<- function(x, weights=NULL, normwt=FALSE, na.rm=TRUE){
sqrt(wtd.var(x, weights, normwt, na.rm))
}

##############################
if(!is.character(x)){
xt<-deparse(substitute(x))
if(regexpr('~', xt)>=0){
x<-gsub('[~].*$', '', xt)
xpart<-gsub('^.*[~]', '', xt)
group<-gsub('[+].*$', '', xpart)
if(nchar(gsub('[^+]', '', xpart))==1) group2<-gsub('^.*[+]', '', xpart)
}}
##############################


type<-tolower(type)
if(type== 'b' | type== 'box' | type=='boxplot')    type<-1
if(type== 'sd')  type<-2
if(type== '2sd' | type== '2 sd')  type<-3
if(type== 'iqr' | type== 'IQR'| type== 'median'| type== 'm')  type<-4

require(plotrix)
if(is.null(group)){
test=FALSE
varwidth=FALSE
}

if(!is.null(group) & is.character(group)) if(is.null(xlab)) xlab <-group

##############################
# Mit Data frame
ismitdata=FALSE
if(is.numeric(x)){
data<-NULL
xname<-gsub('.*[$]', '' ,deparse(substitute(x)))
}

if(is.character(x)){
ismitdata=TRUE
xname<-x
if(!is.null(group)){
group<-eval(parse(text=paste("data$",group, sep='')))
}
if(!is.null(group2)){
group2<-eval(parse(text=paste("data$",group2, sep='')))
}
x<-eval(parse(text=paste("data$",x, sep='')))
w<-eval(parse(text=paste("data$",w, sep='')))
}
if(is.null(ylab)) ylab <-xname
if(!is.null(group)) if(is.null(xlab)) xlab <-xname<-gsub('[(]{0}[A-Za-z0-9]*[$]', '' ,deparse(substitute(group)))
############################



# Errors
if(!is.null(group)) { if(!is.factor(group))    group<-as.factor(group) }
if(!is.null(group2)){ if(!is.factor(group2))   group2<-as.factor(group2) }

if(is.null(group)){  group<-as.factor(rep(1, length(x)))}
if(!is.numeric(x))     stop('x is not numeric!!')
if(!is.numeric(w))     stop('w is not numeric!!')

if(is.null(tcol)  & wall==0)   tcol<-1
if(is.null(tcol)  & wall!=0)   tcol<-rgb(0.1,0.1,0.25)
if(is.null(bgcol) & wall==0)   bgcol<-1
if(is.null(bgcol) & wall!=0)   bgcol<-'#DBE0E8'
if(is.null(col)  & wall==0)   col<-'gray25'
if(is.null(col)  & wall!=0)   col<-rgb(0.3,0.3,0.45)
if(is.null(lcol))  lcol<- tcol
if(is.null(alpha) & wall==0)  alpha<-1
if(is.null(alpha) & wall!=0)  alpha<-0.75

#############################################################################################
#######################################################################################
################################################################################
if(is.null(group2))  n.g2<- 1
g2 <-NULL
if(!is.null(group2)) n.g2<- nlevels(group2)
g <-  group
if(!is.null(group2)) g2 <- group2
n.g<- nlevels(g)


if(is.null(group) & is.null(group2))  xrange<- range(x, na.rm=TRUE)
if(!is.null(group) & is.null(group2)) xrange<- range(x[!is.na(g)], na.rm=TRUE)
if(is.null(group) & !is.null(group2))  xrange<- range(x[!is.na(g2)], na.rm=TRUE)
if(!is.null(group) & !is.null(group2))  xrange<- range(x[!is.na(g) & !is.na(g2)], na.rm=TRUE)
par(new=FALSE)
if(is.null(ylim)) ylim<-xrange

axlin<-0

count.text<-F
if(!is.logical(count) & is.character(count)){
count.text<-T
count.t<-count
count<-T
}
nxlin<-0
if(count){
nxlin<-1
oldmar<-par('mar')
newmar<-rep(0, 4)
if(par('mar')[3]<(4.6+nxlin)) newmar[3]<- (4.6+nxlin) - oldmar[3]
par(mar=(oldmar+newmar))
}

if(wall==5){
newmar<-rep(0, 4)
oldmar<-par('mar')
if(oldmar[2]<4.85) newmar[2]<- 4.85 - oldmar[2]
if(oldmar[3]<4.6) newmar[3]<- 4.6 - oldmar[3]
if(oldmar[4]>1.2 & oldmar[4]<=2.54) newmar[4]<- 1.2-oldmar[4]
par(mar=(oldmar+newmar))
}







#
##########################################


#########################################
# BOX PLOT Per Hand
par(cex=cex)
par(cex.axis=cex.axis)
par(col.main=tcol)
par(col.axis=tcol)
par(col.lab=tcol)
par(font=2)
par(fg=rgb(1,1,1))
par(lend=1)
xrat<- 1:(n.g*n.g2)
xlim<-c(1-(1.5/(n.g*n.g2)),(n.g*+n.g2)+(1.5/(n.g*n.g2)))

xpand<-NULL
if(!varwidth){
xpand<-rep(0.3, length(xrat))
}


plot(c(0,0), col=rgb(1,1,1,0), xlim=xlim, ylim=ylim, axes=FALSE, ylab='', main='', xlab='')
bstats<-NULL
bstats<-as.list(bstats)
count.n<-NULL
for(k in 1:n.g2){
if(!is.null(group2)) bstats[[k]]<-boxplot(x[which(g2==levels(g2)[k])] ~ g[which(g2==levels(g2)[k])], plot=FALSE, pch=16, cex=1.25, col=rgb(1,1,1), staplewex=c(0.25,0.25,0.25,0.25), boxwex=c(0.75,0.75,0.75,0.75), ylim=ylim, xlim=xlim, names=NA, varwidth=varwidth)
if(is.null(group2))  bstats[[k]]<-boxplot(x ~ g, plot=FALSE, pch=16, cex=1.25, col=rgb(1,1,1), staplewex=c(0.25,0.25,0.25,0.25), boxwex=c(0.75,0.75,0.75,0.75), ylim=ylim, xlim=xlim, names=NA, varwidth=varwidth)
count.n<-c(count.n, bstats[[k]]$n)
if(varwidth)  xpand<-c(xpand, sqrt(bstats[[k]]$n/sum(!is.na(x)))/2)
}
if(count.text){
if(length(count.t)==1) count.t<-rep(count.t, length(count.n))
if(length(count.t)==n.g) count.t<-rep(count.t, n.g2)
for(i in 1:length(count.n)){
count.n[i]<- gsub('[?]', count.n[i], count.t[i])
}
}

if(!is.null(group2)) xmeans<-as.vector(unlist(by(x, list(g, g2), mean,na.rm=TRUE , simplify =FALSE)))
if(is.null(group2))  xmeans<-as.vector(unlist(by(x, g, mean,na.rm=TRUE , simplify =FALSE)))


vnames2<-NULL
if(is.list(vnames)){
vnames2<-vnames[[2]]
vnames<-vnames[[1]]
}


if(!is.list(vnames)){
if(!is.null(vnames)) vnames<-rep(vnames, n.g2)
if(is.null(vnames))  vnames<-rep(levels(g), n.g2)
if(is.null(vnames2)) vnames2<-levels(g2)
}


################################################################################
a.draw.box<-function(v, w, at, expand, bstats, means, wall, zylinder=FALSE, col=1, alpha=1, type='box'){

################
# Der Boxplot
if(type==1){

vmedian<-vmeans<-v25<-v75<-vIQR<-uppW<-lowW<- NULL
outs<- list()
for(i in 1:length(v)){
vmeans<- c(vmeans, wtd.mean(v[[i]], w[[i]], na.rm=TRUE))
km <- wtd.quantile(v[[i]], w[[i]], probs=c(0.5))
k25<- wtd.quantile(v[[i]], w[[i]], probs=c(0.25))
k75<- wtd.quantile(v[[i]], w[[i]], probs=c(0.75))

uppWk <- k75+1.5*(k75-k25)
lowWk <- k25-1.5*(k75-k25)
uppWk <- min(max(v[[i]][which(v[[i]]<uppWk)], na.rm=T),  uppWk)
lowWk <- max(min(v[[i]][which(v[[i]]>lowWk)], na.rm=T),  lowWk)

outs<-v[[i]][which(v[[i]]<lowWk | v[[i]]>uppWk)]
ows<- w[[i]][which(v[[i]]<lowWk | v[[i]]>uppWk)]

if(outlier){
points(rep(at[i], length(outs)), outs, col=a.alpha.ade(a.coladd.ade(col, 50), alpha), pch=16, cex=ows)
points(rep(at[i], length(outs)), outs, col=a.alpha.ade(a.coladd.ade(col, 25), alpha), pch=1 , cex=ows)
}

vmedian<- c(vmedian, km)
v25<- c(v25, k25)
v75<- c(v75, k75)
vIQR<- c(vIQR, k75-k25)
uppW<- c(uppW, uppWk)
lowW<- c(lowW, lowWk)
}

if(zylinder){
library(plotrix)
cylindrect(at-expand, v75, at+expand, v25,  col=col, border=col)
}

if(!zylinder) {
rect(at-expand, v75, at+expand, v25,  col=a.alpha.ade(a.coladd.ade(col, 75), alpha), border=col)
}
segments(at, lowW, at, v25, col = col, lty = 1, lwd = lwd)
segments(at, v75, at,  uppW, col = col, lty = 1, lwd = lwd)
segments(at +expand/3, uppW, at-expand/3, uppW, col = col, lty = 1, lwd = lwd)
segments(at +expand/3, lowW, at-expand/3, lowW, col = col, lty = 1, lwd = lwd)
segments(at +expand,   vmedian, at-expand,  vmedian, col = col, lty = 1, lwd = lwd+1)
if(means){
points(at, vmeans, col=col, pch=15)
}
}
################


################
if(type==2){
vmeans<-NULL
vsds  <-NULL

for(i in 1:length(v)){
vmeans<- c(vmeans, wtd.mean(v[[i]], w[[i]], na.rm=TRUE))
vsds  <- c(vsds,   wtd.sd(v[[i]],   w[[i]], na.rm=TRUE))
}

arrows(at, vmeans-vsds, at, vmeans+vsds, col = col, lty = 1, lwd = lwd, angle = 90, code=3, length = 0.1)
lines(at,vmeans, col=col, lwd=2, lty=3)
points(at, vmeans, col=col, pch=16)
}
################

################
if(type==3){
vmeans<-NULL
vsds  <-NULL

for(i in 1:length(v)){
vmeans<- c(vmeans, wtd.mean(v[[i]], w[[i]], na.rm=TRUE))
vsds  <- c(vsds,   wtd.sd(v[[i]],   w[[i]], na.rm=TRUE))
}


arrows(at, vmeans-vsds*2, at, vmeans+vsds*2, col = a.coladd.ade(col, 75), lty = 1, lwd = lwd, angle = 90, code=3, length = 0.1)
arrows(at, vmeans-vsds,    at, vmeans+vsds, col = col, lty = 1, lwd = lwd, angle = 90, code=3, length = 0.05)
lines(at,vmeans, col=col, lwd=2, lty=3)
points(at, vmeans, col=col, pch=16)
}
################


################
if(type==4){
vmedian<-v25<-v75<-NULL
for(i in 1:length(v)){
vmedian<- c(vmedian, wtd.quantile(v[[i]], w[[i]], probs=c(0.5)))
v25<- c(v25, wtd.quantile(v[[i]], w[[i]], probs=c(0.25)))
v75<- c(v75, wtd.quantile(v[[i]], w[[i]], probs=c(0.75)))
}

arrows(at, v25, at, v75, col = col, lty = 1, lwd = lwd, angle = 45, code=3, length = 0.1)
lines(at,vmedian, col=col, lwd=2, lty=3)
points(at,vmedian, col=col, pch=18, cex=1.25)
}
################



}
################################################################################


################################################################################
################################################################################
#Walltype 0
if(wall==0){
if((n.g*n.g2)>1)                           axis(1, at=xrat,                         labels = vnames, tick=FALSE, col=bgcol, lwd.ticks=1, lwd=0, padj=0.5, line =axlin-0.5)
if(!is.null(yticks) & length(yticks)==1)   axis(2, at=pretty(x, n=yticks),          labels =TRUE,       tick=TRUE, col=bgcol, lwd.ticks=1, lwd=0)
if(!is.null(yticks) & length(yticks)>1 )   axis(2, at=yticks,                       labels =TRUE,       tick=TRUE, col=bgcol, lwd.ticks=1, lwd=0)
if(is.null(yticks) )                       axis(2,                                  labels =TRUE,       tick=TRUE, col=bgcol, lwd.ticks=1, lwd=0)
if(!is.null(group2))                       axis(3, at=((1:n.g2)*(n.g))+0.5-(n.g/2), labels =vnames2, tick=FALSE, col=bgcol, lwd.ticks=1, lwd=0, line=-0.75+nxlin)

if(!is.null(group2))                       abline(v=seq(n.g, n.g*(n.g2-1), n.g)+0.5,                         col=bgcol,              lwd=1)
abline(v=NULL, h=h, col=lcol, lwd=1, lty=lty)

for(k in 1:n.g2){
xids<-(1+n.g*(k-1)):(n.g*k)
v<-NULL
v<-as.list(v)

wv<-NULL
wv<-as.list(wv)

if(!is.null(g2)){
xk<-x[which(g2==levels(g2)[k])]
wk<-w[which(g2==levels(g2)[k])]
gk<-g[which(g2==levels(g2)[k])]
}

if(is.null(g2)){
xk<-x
wk<-w
gk<-g
}
for(i in 1:n.g){
v[[i]] <-xk[which(gk==levels(gk)[i])]
wv[[i]]<-wk[which(gk==levels(gk)[i])]
}

a.draw.box(v, wv, xrat[xids], xpand[xids], bstats[[k]], means, wall=wall, zylinder=zylinder, col=col, alpha=alpha, type)
}

####################
# Labels
text(x=a.glc(0), y=a.glc(3, 2+nxlin), labels=main, cex=par('cex.main'), col=tcol, adj=c(0.5, 0.5), xpd=TRUE, font=par("font.main"))
text(x=a.glc(0), y=a.glc(1, 3.5+axlin), labels=xlab, cex=par('cex.lab'), col=tcol, adj=c(0.5, 0.5), xpd=TRUE, font=par("font.lab"))
text(x=a.glc(2, 3.5), y=a.glc(5), labels=ylab, cex=par('cex.lab'), col=tcol, adj=c(0.5, 0.5), xpd=TRUE, font=par("font.lab"), srt=90)
####################

####################
# Counts
if(count){
#rect( a.glc(side=2, line=0), a.glc(side=3, line=0),a.glc(4, line=0),a.glc(side=3, line=nxlin), col=rgb(1,1,1,0), border=bgcol, lwd=1, xpd=TRUE)
text(x=xrat, y=a.glc(3, line=(nxlin/2)),  labels=count.n,  cex=(cex-0.25),  col=tcol, adj=c(0.5, 0.5), xpd=TRUE)
}
####################

if(!is.null(group2)) segments(x0=seq(n.g, n.g*(n.g2-1), n.g)+0.5, y0=par('usr')[4]+(a.glc(3,nxlin)-a.glc(3, 0)), x1=seq(n.g, n.g*(n.g2-1), n.g)+0.5, y1=a.glc(1, 0+axlin), col=bgcol,lwd=1, xpd=TRUE)
rect(a.glc(2,0),a.glc(1,0+axlin),a.glc(4,0),a.glc(3, 0+nxlin) ,  col=rgb(1,1,1,0), border=bgcol, lwd=1, xpd=TRUE)

}
################################################################################



################################################################################
#Walltype 1
if(wall==1){
polygon( par('usr')[c(1,1,2,2)], par('usr')[c(3,4,4,3)], col=bgcol)
#grid(col = rgb(1,1,1), lty = "solid", lwd = 1)
if((n.g*n.g2)>1)                           axis(1, at=xrat,                         labels = vnames, tick=FALSE, col=tcol, lwd.ticks=1, lwd=0, padj=0.5, line =axlin-0.5)
if(!is.null(yticks) & length(yticks)==1)   zt<-axis(2, at=pretty(x, n=yticks),          labels =TRUE,       tick=TRUE, col=tcol, lwd.ticks=1, lwd=0)
if(!is.null(yticks) & length(yticks)>1 )   zt<-axis(2, at=yticks,                       labels =TRUE,       tick=TRUE, col=tcol, lwd.ticks=1, lwd=0)
if(is.null(yticks) )                       zt<-axis(2,                                  labels =TRUE,       tick=TRUE, col=tcol, lwd.ticks=1, lwd=0)
if(!is.null(group2))                       axis(3, at=((1:n.g2)*(n.g))+0.5-(n.g/2), labels =vnames2, tick=FALSE, col=tcol, lwd.ticks=1, lwd=0, line=-0.75+nxlin)
abline(h=zt, col=rgb(1,1,1), lwd=1)
abline(v=NULL, h=h, col=rgb(1,1,1), lwd=1, lty=lty)
abline(v=xrat, col=rgb(1,1,1), lwd=1, lty=1)



for(k in 1:n.g2){
xids<-(1+n.g*(k-1)):(n.g*k)

v<-NULL
v<-as.list(v)

wv<-NULL
wv<-as.list(wv)

if(!is.null(g2)){
xk<-x[which(g2==levels(g2)[k])]
gk<-g[which(g2==levels(g2)[k])]
wk<-w[which(g2==levels(g2)[k])]
}
if(is.null(g2)){
xk<-x
gk<-g
wk<-w
}
for(i in 1:n.g){
v[[i]]<-xk[which(gk==levels(gk)[i])]
wv[[i]]<-wk[which(gk==levels(gk)[i])]
}
a.draw.box(v, wv, xrat[xids], xpand[xids], bstats[[k]], means, wall=wall, zylinder=zylinder, col=col, alpha=alpha, type)
}

####################
# Labels
text(x=a.glc(0), y=a.glc(3, 2+nxlin), labels=main, cex=par('cex.main'), col=tcol, adj=c(0.5, 0.5), xpd=TRUE, font=par("font.main"))
text(x=a.glc(0), y=a.glc(1, 3.5+axlin), labels=xlab, cex=par('cex.lab'), col=tcol, adj=c(0.5, 0.5), xpd=TRUE, font=par("font.lab"))
text(x=a.glc(2, 3.5), y=a.glc(5), labels=ylab, cex=par('cex.lab'), col=tcol, adj=c(0.5, 0.5), xpd=TRUE, font=par("font.lab"), srt=90)
####################

####################
# Counts
if(count){
rect( a.glc(side=2, line=0), a.glc(side=3, line=0),a.glc(4, line=0),a.glc(side=3, line=nxlin), col=bgcol, border=rgb(1,1,1), lwd=2, xpd=TRUE)
text(x=xrat, y=a.glc(3, line=(nxlin/2)),  labels=count.n,  cex=(cex-0.25),  col=tcol, adj=c(0.5, 0.5), xpd=TRUE)
}
####################

if(!is.null(group2)) segments(x0=seq(n.g, n.g*(n.g2-1), n.g)+0.5, y0=par('usr')[4]+(a.glc(3,nxlin)-a.glc(3, 0)), x1=seq(n.g, n.g*(n.g2-1), n.g)+0.5, y1=a.glc(1, 0+axlin), col=rgb(1,1,1),lwd=2, xpd=TRUE)
box( col=rgb(1,1,1), lwd=2)
}
################################################################################



################################################################################
#Walltype 2
if(wall==2){

if((n.g*n.g2)>1)                           axis(1, at=xrat,                         labels = vnames, tick=FALSE,  col=a.coladd.ade(bgcol, -75), lwd.ticks=1, lwd=0, padj=0.5, line =axlin-0.5)
if(!is.null(yticks) & length(yticks)==1)   zt<-axis(2, at=pretty(x, n=yticks),          labels =TRUE,       tick=TRUE, col=a.coladd.ade(bgcol, -75), lwd.ticks=1, lwd=0)
if(!is.null(yticks) & length(yticks)>1 )   zt<-axis(2, at=yticks,                       labels =TRUE,       tick=TRUE, col=a.coladd.ade(bgcol, -75), lwd.ticks=1, lwd=0)
if(is.null(yticks) )                       zt<-axis(2,                                  labels =TRUE,       tick=TRUE, col=a.coladd.ade(bgcol, -75), lwd.ticks=1, lwd=0)
if(!is.null(group2))                       axis(3, at=((1:n.g2)*(n.g))+0.5-(n.g/2), labels =vnames2, tick=FALSE, col=a.coladd.ade(bgcol, -75), lwd.ticks=1, lwd=0, line=-0.75+nxlin)
abline(h=zt, col=bgcol, lwd=1)
#grid(col = bgcol, lty = "solid", lwd = 1)
abline(v=NULL, h=h, col=lcol, lwd=1, lty=lty)
abline(v=xrat, col=bgcol, lwd=1, lty=1)


if(test) rect(a.glc(2,0),a.glc(1,0),a.glc(4,0),a.glc(1, 2.5) ,  col=rgb(1,1,1,0), border=a.coladd.ade(bgcol, -75), lwd=1, xpd=TRUE)

for(k in 1:n.g2){
xids<-(1+n.g*(k-1)):(n.g*k)
v<-NULL
v<-as.list(v)

wv<-NULL
wv<-as.list(wv)

if(!is.null(g2)){
xk<-x[which(g2==levels(g2)[k])]
gk<-g[which(g2==levels(g2)[k])]
wk<-w[which(g2==levels(g2)[k])]
}
if(is.null(g2)){
xk<-x
gk<-g
wk<-w
}
for(i in 1:n.g){
v[[i]]<-xk[which(gk==levels(gk)[i])]
wv[[i]]<-wk[which(gk==levels(gk)[i])]
}
a.draw.box(v, wv, xrat[xids], xpand[xids], bstats[[k]], means, wall=wall, zylinder=zylinder, col=col, alpha=alpha, type)

}

####################
# Labels
text(x=a.glc(0), y=a.glc(3, 2+nxlin), labels=main, cex=par('cex.main'), col=tcol, adj=c(0.5, 0.5), xpd=TRUE, font=par("font.main"))
text(x=a.glc(0), y=a.glc(1, 3.5+axlin), labels=xlab, cex=par('cex.lab'), col=tcol, adj=c(0.5, 0.5), xpd=TRUE, font=par("font.lab"))
text(x=a.glc(2, 3.5), y=a.glc(5), labels=ylab, cex=par('cex.lab'), col=tcol, adj=c(0.5, 0.5), xpd=TRUE, font=par("font.lab"), srt=90)
####################

####################
# Counts
if(count){
rect( a.glc(side=2, line=0), a.glc(side=3, line=0),a.glc(4, line=0),a.glc(side=3, line=nxlin), col=rgb(1,1,1), border=a.coladd.ade(bgcol, -75), lwd=1, xpd=TRUE)
text(x=xrat, y=a.glc(3, line=(nxlin/2)),  labels=count.n,  cex=(cex-0.25),  col=tcol, adj=c(0.5, 0.5), xpd=TRUE)
}
####################

if(!is.null(group2)) segments(x0=seq(n.g, n.g*(n.g2-1), n.g)+0.5, y0=par('usr')[4]+(a.glc(3,nxlin)-a.glc(3, 0)), x1=seq(n.g, n.g*(n.g2-1), n.g)+0.5, y1=a.glc(1, 0+axlin),  col=a.coladd.ade(bgcol, -75),lwd=2, xpd=TRUE)
box(col=a.coladd.ade(bgcol, -75))
}
################################################################################



################################################################################
#Walltype 3
if(wall==3){

polygon( par('usr')[c(1,1,2,2)], par('usr')[c(3,4,4,3)], col=bgcol)
#grid(col = a.coladd.ade(bgcol, -50), lty = "solid", lwd = 1)
if((n.g*n.g2)>1)                           axis(1, at=xrat,                         labels = vnames, tick=FALSE, col=a.coladd.ade(bgcol, -50), lwd.ticks=1, lwd=0, padj=0.5, line =axlin-0.5)
if(!is.null(yticks) & length(yticks)==1)   zt<-axis(2, at=pretty(x, n=yticks),          labels =TRUE,       tick=TRUE, col=a.coladd.ade(bgcol, -50), lwd.ticks=1, lwd=0)
if(!is.null(yticks) & length(yticks)>1 )   zt<-axis(2, at=yticks,                       labels =TRUE,       tick=TRUE, col=a.coladd.ade(bgcol, -50), lwd.ticks=1, lwd=0)
if(is.null(yticks) )                       zt<-axis(2,                                  labels =TRUE,       tick=TRUE, col=a.coladd.ade(bgcol, -50), lwd.ticks=1, lwd=0)
if(!is.null(group2))                       axis(3, at=((1:n.g2)*(n.g))+0.5-(n.g/2), labels =vnames2, tick=FALSE, col=a.coladd.ade(bgcol, -50), lwd.ticks=1, lwd=0, line=-0.75+nxlin)
abline(h=zt, col=a.coladd.ade(bgcol, -50), lwd=1)
abline(v=NULL, h=h, col=lcol, lwd=1, lty=lty)
abline(v=xrat, col=a.coladd.ade(bgcol, -50), lwd=1, lty=1)


if(test) rect(a.glc(2,0),a.glc(1,0),a.glc(4,0),a.glc(1, 2.5) ,  col=rgb(1,1,1,0), border=a.coladd.ade(bgcol, -50), lwd=1, xpd=TRUE)

for(k in 1:n.g2){
xids<-(1+n.g*(k-1)):(n.g*k)

v<-NULL
v<-as.list(v)

wv<-NULL
wv<-as.list(wv)

if(!is.null(g2)){
xk<-x[which(g2==levels(g2)[k])]
gk<-g[which(g2==levels(g2)[k])]
wk<-w[which(g2==levels(g2)[k])]
}
if(is.null(g2)){
xk<-x
gk<-g
wk<-w
}
for(i in 1:n.g){
v[[i]]<-xk[which(gk==levels(gk)[i])]
wv[[i]]<-wk[which(gk==levels(gk)[i])]
}
a.draw.box(v, wv, xrat[xids], xpand[xids], bstats[[k]], means, wall=wall, zylinder=zylinder, col=col, alpha=alpha, type)



}

####################
# Labels
text(x=a.glc(0), y=a.glc(3, 2+nxlin), labels=main, cex=par('cex.main'), col=tcol, adj=c(0.5, 0.5), xpd=TRUE, font=par("font.main"))
text(x=a.glc(0), y=a.glc(1, 3.5+axlin), labels=xlab, cex=par('cex.lab'), col=tcol, adj=c(0.5, 0.5), xpd=TRUE, font=par("font.lab"))
text(x=a.glc(2, 3.5), y=a.glc(5), labels=ylab, cex=par('cex.lab'), col=tcol, adj=c(0.5, 0.5), xpd=TRUE, font=par("font.lab"), srt=90)
####################

####################
# Counts
if(count){
rect( a.glc(side=2, line=0), a.glc(side=3, line=0),a.glc(4, line=0),a.glc(side=3, line=nxlin), col=rgb(1,1,1), border=a.coladd.ade(bgcol, -50), lwd=1, xpd=TRUE)
text(x=xrat, y=a.glc(3, line=(nxlin/2)),  labels=count.n,  cex=(cex-0.25),  col=tcol, adj=c(0.5, 0.5), xpd=TRUE)
}
####################

if(!is.null(group2)) segments(x0=seq(n.g, n.g*(n.g2-1), n.g)+0.5, y0=par('usr')[4]+(a.glc(3,nxlin)-a.glc(3, 0)), x1=seq(n.g, n.g*(n.g2-1), n.g)+0.5, y1=a.glc(1, 0+axlin),  col=a.coladd.ade(bgcol, -50),lwd=2, xpd=TRUE)
box(col=a.coladd.ade(bgcol, -50))
}
################################################################################



################################################################################
#Walltype 4
if(wall==4){
par(col.lab=rgb(1,1,1))
par(col.main=rgb(1,1,1))
par(font=2)
polygon( par('usr')[c(1,1,2,2)], par('usr')[c(3,4,4,3)], col=bgcol)

#grid(col = rgb(1,1,1), lty = "solid", lwd = 1)
if((n.g*n.g2)>1)                           axis(1, at=xrat,                         labels = vnames, tick=FALSE, col=tcol, lwd.ticks=1, lwd=0, padj=0.5, line =axlin-0.5)
if(!is.null(yticks) & length(yticks)==1)   zt<-axis(2, at=pretty(x, n=yticks),          labels =TRUE,       tick=TRUE, col=tcol, lwd.ticks=1, lwd=0)
if(!is.null(yticks) & length(yticks)>1 )   zt<-axis(2, at=yticks,                       labels =TRUE,       tick=TRUE, col=tcol, lwd.ticks=1, lwd=0)
if(is.null(yticks) )                       zt<-axis(2,                                  labels =TRUE,       tick=TRUE, col=tcol, lwd.ticks=1, lwd=0)
if(!is.null(group2))                       axis(3, at=((1:n.g2)*(n.g))+0.5-(n.g/2), labels =vnames2, tick=FALSE, col=tcol, lwd.ticks=1, lwd=0, line=-0.8+nxlin)
abline(h=zt, col=rgb(1,1,1), lwd=1)
abline(v=NULL, h=h, col=lcol, lwd=1, lty=lty)
abline(v=xrat, col=rgb(1,1,1), lwd=1, lty=1)


par(xpd=TRUE)
dx<-7/par('din')[1]
dy<-7/par('din')[2]
xr<-(diff(par('usr')[1:2])/10)*dx
yr<-(diff(par('usr')[3:4])/10)*dy

if(!is.null(ylab)) if(!(ylab==''  | ylab==' ')) polygon( a.glc(side=2, line=c(3.5, 3.5, 2, 2)), a.glc(side=c(1,3,3,1), line=0) ,  col=bgcol, border=rgb(1,1,1))
if(!is.null(xlab)) if(!(xlab=='' | xlab==' '))  polygon( a.glc(side=c(2,2,4,4), line=0), a.glc(side=1, line=c(3+axlin, 4.5+axlin, 4.5+axlin, 3+axlin)), col=bgcol, border=rgb(1,1,1))
if(!is.null(main)) if(!(main=='' | main==' '))  polygon( a.glc(side=c(2,2,4,4), line=0), a.glc(side=3, line=c(1.2+nxlin,3+nxlin,3+nxlin,1.2+nxlin)), col=tcol, border=rgb(1,1,1))

par(xpd=FALSE)

if(test) rect(a.glc(2,0),a.glc(1,0),a.glc(4,0),a.glc(1, 2.5) ,  col=tcol, border=rgb(1,1,1), lwd=1, xpd=TRUE)

for(k in 1:n.g2){
xids<-(1+n.g*(k-1)):(n.g*k)

v<-NULL
v<-as.list(v)

wv<-NULL
wv<-as.list(wv)

if(!is.null(g2)){
xk<-x[which(g2==levels(g2)[k])]
gk<-g[which(g2==levels(g2)[k])]
wk<-w[which(g2==levels(g2)[k])]
}
if(is.null(g2)){
xk<-x
gk<-g
wk<-w
}
for(i in 1:n.g){
v[[i]]<-xk[which(gk==levels(gk)[i])]
wv[[i]]<-wk[which(gk==levels(gk)[i])]
}
a.draw.box(v, wv, xrat[xids], xpand[xids], bstats[[k]], means, wall=wall, zylinder=zylinder, col=col, alpha=alpha, type)


}

####################
# Labels
text(x=a.glc(0), y=a.glc(3, 2.1+nxlin), labels=main, cex=par('cex.main'), col=rgb(1,1,1), adj=c(0.5, 0.5), xpd=TRUE, font=par("font.main"))
text(x=a.glc(0), y=a.glc(1, 3.75+axlin), labels=xlab, cex=par('cex.lab'), col=tcol, adj=c(0.5, 0.5), xpd=TRUE, font=2)
text(x=a.glc(2, 2.8), y=a.glc(5), labels=ylab, cex=par('cex.lab'), col=tcol, adj=c(0.5, 0.5), xpd=TRUE, font=2, srt=90)
####################

####################
# Counts
if(count){
rect( a.glc(side=2, line=0), a.glc(side=3, line=0),a.glc(4, line=0),a.glc(side=3, line=nxlin), col=tcol, border=rgb(1,1,1), lwd=1, xpd=TRUE)
text(x=xrat, y=a.glc(3, line=(nxlin/2)),  labels=count.n,  cex=(cex-0.25),  col=rgb(1,1,1), adj=c(0.5, 0.5), xpd=TRUE)
}
####################

if(!is.null(group2)) segments(x0=seq(n.g, n.g*(n.g2-1), n.g)+0.5, y0=par('usr')[4]+(a.glc(3,nxlin)-a.glc(3, 0)), x1=seq(n.g, n.g*(n.g2-1), n.g)+0.5, y1=a.glc(1, 0+axlin), col=rgb(1,1,1),lwd=2, xpd=TRUE)
box(col=rgb(1,1,1))
}
################################################################################



################################################################################
#Walltype 5
if(wall==5){
par(col.lab=tcol)
par(col.main=tcol)

if((n.g*n.g2)>1)                           axis(1, at=xrat,                         labels = vnames, tick=FALSE, col=tcol, lwd.ticks=1, lwd=0, padj=0.5, line =axlin-0.5)
if(!is.null(yticks) & length(yticks)==1)   axis(2, at=pretty(x, n=yticks),          labels =TRUE,       tick=TRUE, col=tcol, lwd.ticks=1, lwd=0)
if(!is.null(yticks) & length(yticks)>1 )   axis(2, at=yticks,                       labels =TRUE,       tick=TRUE, col=tcol, lwd.ticks=1, lwd=0)
if(is.null(yticks) )                       axis(2,                                  labels =TRUE,       tick=TRUE, col=tcol, lwd.ticks=1, lwd=0)
if(!is.null(group2))                       axis(3, at=((1:n.g2)*(n.g))+0.5-(n.g/2), labels =vnames2, tick=FALSE, col=tcol, lwd.ticks=1, lwd=0, line=-0.7+nxlin)

abline(v=NULL, h=h, col=lcol, lwd=1, lty=lty)


par(xpd=TRUE)

polygon(a.glc(side=2, line=c(4.25, 4.25, 0, 0)), a.glc(side=3, line=c(1.6, 4, 4, 1.6)+nxlin), col=bgcol,        border=tcol)
polygon(a.glc(side=c(2,2,4,4), line=c(0,0,0,0)), a.glc(side=3, line=c(1.6, 4, 4, 1.6)+nxlin), col=rgb(1,1,1,0), border=tcol)
polygon(a.glc(side=4, line=c(0, 0 ,0.6, 0.6)),   a.glc(side=3, line=c(1.6, 4, 4, 1.6)+nxlin), col=bgcol,        border=tcol)
polygon(a.glc(side=2, line=c(4.25, 4.25 ,3.65, 3.65)),  a.glc(side=c(1,3,3,1), line=c(3.6++axlin, 1.6+nxlin, 1.6+nxlin, 3.6++axlin)), col=bgcol,  border=tcol)
polygon(a.glc(side=4, line=c(0, 0 ,0.6, 0.6)), a.glc(side=c(1, 3, 3, 1), line=c(0+axlin, 0+nxlin, 0+nxlin, 0+axlin)), col=bgcol, border=tcol)
polygon(a.glc(side=2, line=c(4.25, 4.25, 0, 0)), a.glc(side=1, line=c(2.6, 4.5, 4.5, 2.6)+axlin), col=bgcol, border=tcol)
polygon(a.glc(side=c(2, 2, 4, 4), line=0), a.glc(side=1, line=c(2.6, 4.5, 4.5, 2.6)+axlin), col=rgb(1,1,1,0), border=tcol)
polygon(a.glc(side=4, line=c(0, 0, 0.6, 0.6)), a.glc(side=1, line=c(2.6, 4.5, 4.5, 2.6)+axlin), col=bgcol, border=tcol)
par(xpd=FALSE)



if(test) rect(a.glc(2,0),a.glc(1,0),a.glc(4,0),a.glc(1, 2.5) ,  col=rgb(1,1,1), border=tcol, lwd=1, xpd=TRUE)

for(k in 1:n.g2){
xids<-(1+n.g*(k-1)):(n.g*k)

v<-NULL
v<-as.list(v)

wv<-NULL
wv<-as.list(wv)

if(!is.null(g2)){
xk<-x[which(g2==levels(g2)[k])]
gk<-g[which(g2==levels(g2)[k])]
wk<-w[which(g2==levels(g2)[k])]
}
if(is.null(g2)){
xk<-x
gk<-g
wk<-w
}
for(i in 1:n.g){
v[[i]]<-xk[which(gk==levels(gk)[i])]
wv[[i]]<-wk[which(gk==levels(gk)[i])]
}
a.draw.box(v,  wv, xrat[xids], xpand[xids], bstats[[k]], means, wall=wall, zylinder=zylinder, col=col, alpha=alpha, type)


}

####################
# Labels
text(x=a.glc(0), y=a.glc(3, 3+nxlin-0.1), labels=main, cex=par('cex.main'), col=tcol, adj=c(0.5, 0.5), xpd=TRUE, font=par("font.main"))
text(x=a.glc(0), y=a.glc(1, 3.5+axlin), labels=xlab, cex=par('cex.lab'), col=tcol, adj=c(0.5, 0.5), xpd=TRUE, font=2)
text(x=a.glc(2, 3), y=a.glc(5), labels=ylab, cex=par('cex.lab'), col=tcol, adj=c(0.5, 0.5), xpd=TRUE, font=2, srt=90)
####################

####################
# Counts
if(count){
rect( a.glc(side=2, line=0), a.glc(side=3, line=0),a.glc(4, line=0),a.glc(side=3, line=nxlin), col=rgb(1,1,1), border=tcol, lwd=1, xpd=TRUE)
text(x=xrat, y=a.glc(3, line=(nxlin/2)),  labels=count.n,  cex=(cex-0.25),  col=tcol, adj=c(0.5, 0.5), xpd=TRUE)
}
####################

if(!is.null(group2)) segments(x0=seq(n.g, n.g*(n.g2-1), n.g)+0.5, y0=par('usr')[4]+(a.glc(3,nxlin)-a.glc(3, 0)), x1=seq(n.g, n.g*(n.g2-1), n.g)+0.5, y1=a.glc(1, 0+axlin),  col=tcol,lwd=2, xpd=TRUE)
box(lwd=1, col=tcol)
}
################################################################################


################################################################################
#Walltype 6
if(wall==6){
polygon(par('usr')[c(1,1,2,2)], par('usr')[c(3,4,4,3)], col=bgcol)
grid(nx=NA,col = rgb(1,1,1), lty = "solid", lwd = 1)
par(lend=0)
if((n.g*n.g2)>1)                           axis(1, at=xrat,                         labels = vnames, tick=FALSE, col=tcol, lwd.ticks=3, lwd=0, col.ticks=a.coladd.ade(bgcol, -35), padj=0.5, line =axlin-0.5)
if((n.g*n.g2)>1)                           axis(1, at=xrat,                         labels = vnames, tick=FALSE, col=tcol, lwd.ticks=1, lwd=0, col.ticks=rgb(1,1,1), padj=0.5, line =axlin-0.5)
if(!is.null(yticks) & length(yticks)==1)  a2<-  axis(2, at=pretty(x, n=yticks),     labels =TRUE,       tick=TRUE, col=tcol, lwd.ticks=3, lwd=0, col.ticks=a.coladd.ade(bgcol, -35))
if(!is.null(yticks) & length(yticks)==1)  a2<-  axis(2, at=pretty(x, n=yticks),     labels =TRUE,       tick=TRUE, col=tcol, lwd.ticks=1, lwd=0, col.ticks=rgb(1,1,1))
if(!is.null(yticks) & length(yticks)>1 )  a2<-  axis(2, at=yticks,                  labels =TRUE,       tick=TRUE, col=tcol, lwd.ticks=3, lwd=0, col.ticks=a.coladd.ade(bgcol, -35))
if(!is.null(yticks) & length(yticks)>1 )  a2<-  axis(2, at=yticks,                  labels =TRUE,       tick=TRUE, col=tcol, lwd.ticks=1, lwd=0, col.ticks=rgb(1,1,1))
if(is.null(yticks) )                      a2<- axis(2,                              labels =TRUE,       tick=TRUE, col=tcol, lwd.ticks=3, lwd=0, col.ticks=a.coladd.ade(bgcol, -35))
if(is.null(yticks) )                      a2<- axis(2,                              labels =TRUE,       tick=TRUE, col=tcol, lwd.ticks=1, lwd=0, col.ticks=rgb(1,1,1))
if(!is.null(group2))                       axis(3, at=((1:n.g2)*(n.g))+0.5-(n.g/2), labels =vnames2, tick=FALSE, col=tcol, lwd.ticks=3, lwd=0, line=-0.75+nxlin, col.ticks=a.coladd.ade(bgcol, -35))
if(!is.null(group2))                       axis(3, at=((1:n.g2)*(n.g))+0.5-(n.g/2), labels =vnames2, tick=FALSE, col=tcol, lwd.ticks=1, lwd=0, line=-0.75+nxlin, col.ticks=rgb(1,1,1))
abline(v=NULL, h=h, col=rgb(1,1,1), lwd=1, lty=lty)
abline(v=NULL, h=h, col=rgb(1,1,1), lwd=1, lty=lty)
abline(h=a2, col=a.coladd.ade(bgcol, -35), lwd=3, lty=1)
abline(h=a2, col=rgb(1,1,1), lwd=1, lty=1)


if(test) rect(a.glc(2,0),a.glc(1,0),a.glc(4,0),a.glc(1, 2.5) ,  col=bgcol, border=rgb(1,1,1), lwd=3, xpd=TRUE)
if(test) rect(a.glc(2,0),a.glc(1,0),a.glc(4,0),a.glc(1, 2.5) ,  col=rgb(0,0,0,0), border=a.coladd.ade(bgcol, -35), lwd=1, xpd=TRUE)


for(k in 1:n.g2){
xids<-(1+n.g*(k-1)):(n.g*k)

v<-NULL
v<-as.list(v)

wv<-NULL
wv<-as.list(wv)

if(!is.null(g2)){
xk<-x[which(g2==levels(g2)[k])]
gk<-g[which(g2==levels(g2)[k])]
wk<-w[which(g2==levels(g2)[k])]
}
if(is.null(g2)){
xk<-x
gk<-g
wk<-w
}
for(i in 1:n.g){
v[[i]]<-xk[which(gk==levels(gk)[i])]
wv[[i]]<-wk[which(gk==levels(gk)[i])]
}
a.draw.box(v, wv, xrat[xids], xpand[xids], bstats[[k]], means, wall=wall, zylinder=zylinder, col=col, alpha=alpha, type)


}

####################
# Labels
text(x=a.glc(0), y=a.glc(3, 2.1+nxlin), labels=main, cex=par('cex.main'), col=tcol, adj=c(0.5, 0.5), xpd=TRUE, font=par("font.main"))
text(x=a.glc(0), y=a.glc(1, 3.6+axlin), labels=xlab, cex=par('cex.lab'), col=tcol, adj=c(0.5, 0.5), xpd=TRUE, font=par("font.lab"))
text(x=a.glc(2, 2.8), y=a.glc(5), labels=ylab, cex=par('cex.lab'), col=tcol, adj=c(0.5, 0.5), xpd=TRUE, font=par("font.lab"), srt=90)
####################

####################
# Counts
if(count){
rect( a.glc(side=2, line=0), a.glc(side=3, line=0),a.glc(4, line=0),a.glc(side=3, line=nxlin), col=bgcol, border=rgb(1,1,1), lwd=3, xpd=TRUE)
rect( a.glc(side=2, line=0), a.glc(side=3, line=0),a.glc(4, line=0),a.glc(side=3, line=nxlin), col=rgb(1,1,1,0), border=a.coladd.ade(bgcol, -35), lwd=1, xpd=TRUE)
text(x=xrat, y=a.glc(3, line=(nxlin/2)),  labels=count.n,  cex=(cex-0.25),  col=tcol, adj=c(0.5, 0.5), xpd=TRUE)
}
####################

if(!is.null(group2)) segments(x0=seq(n.g, n.g*(n.g2-1), n.g)+0.5, y0=par('usr')[4]+(a.glc(3,nxlin)-a.glc(3, 0)), x1=seq(n.g, n.g*(n.g2-1), n.g)+0.5, y1=a.glc(1, 0+axlin), col=a.coladd.ade(bgcol, -35),lwd=3, xpd=TRUE)
if(!is.null(group2)) segments(x0=seq(n.g, n.g*(n.g2-1), n.g)+0.5, y0=par('usr')[4]+(a.glc(3,nxlin)-a.glc(3, 0)), x1=seq(n.g, n.g*(n.g2-1), n.g)+0.5, y1=a.glc(1, 0+axlin), col=rgb(1,1,1),lwd=1, xpd=TRUE)


rect(a.glc(2,0),a.glc(1,0+axlin),a.glc(4,0),a.glc(3, 0+nxlin) ,  col=rgb(1,1,1,0), border=a.coladd.ade(bgcol, -35), lwd=3, xpd=TRUE)
rect(a.glc(2,0),a.glc(1,0+axlin),a.glc(4,0),a.glc(3, 0+nxlin) ,  col=rgb(1,1,1,0), border=rgb(1,1,1), lwd=1, xpd=TRUE)
}
################################################################################
################################################################################
################################################################################
}
