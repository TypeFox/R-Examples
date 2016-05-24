tornado.ade <-
function(x,  group=NULL, group2=NULL, data=NULL, vnames=NULL, gnames=NULL, gnames2=NULL, breaks=6, density=NULL, angle = NULL, xlab=NULL, glab=NULL, main='', legendon='topright', xticks=NULL, col=NULL, tcol=NULL,  bgcol=NULL, lcol=NULL, alpha=NULL, r=0.05, lwd=1, lty=2, wall=0, v=NULL, h=NULL){
if(any(par('mfg')!=c(1,1,1,1)) & any(par('mai') < c(1.02, 0.82, 0.82, 0.42))){
maidiff<-rep(0, 4)
norm<-c(1.02, 0.82, 0.82, 0.42)
maidiff[par('mai')<norm]<-  norm[par('mai')<norm] - par('mai')[par('mai')<norm]
par(mai=par('mai')+maidiff)
}
oldpar<-par(no.readonly =TRUE)
oldpar<-oldpar[-which(names(oldpar)%in%c('usr', 'plt',   'pin', 'fin', 'fig', 'mfg', 'mfcol', 'mfrow', 'omd', 'omi', 'oma'))]
on.exit(par(oldpar))


##############################
xname<-NULL
if(!is.character(x) & !is.factor(x) & !is.list(x) & !is.table(x)){
xt<-deparse(substitute(x))
if(regexpr('~', xt)>=0){
x<-gsub('[~].*$', '', xt)
xpart<-gsub('^.*[~]', '', xt)
group<-gsub('[+].*$', '', xpart)
if(nchar(gsub('[^+]', '', xpart))==1) group2<-gsub('^.*[+]', '', xpart)
}}
##############################


################################################################
# without Data.Frame or with
ismitdata=FALSE
if(is.numeric(x) & !is.list(x)){
data<-NULL
xname<-gsub('[(]{0}[A-Za-z0-9]*[$]', '' ,deparse(substitute(x)))
}
if(is.character(x) & !is.list(x)){
ismitdata=TRUE
if(ismitdata){
if(!is.null(data)){ if(!is.data.frame(data))  stop("(data) must be a data.frame!") }
}
xname<-x
if(!is.null(group)){
groupname<-group
group<-eval(parse(text=paste("data$",group, sep='')))
}
if(!is.null(group2)){
groupname2<-group2
group2<-eval(parse(text=paste("data$",group2, sep='')))
group2<-as.factor(group2)
}
x<-eval(parse(text=paste("data$",x, sep='')))
}

g=1
# Errors
if(!is.null(group) & length(dim(group))==1)  {
if(!is.factor(group))   group<-as.factor(group)
if(nlevels(group )>6) stop('To many levels in group!')
g<-group
}
if(!is.numeric(x) & !is.list(x) & !is.factor(x))     stop('x is not numeric!')
################################################################

################################################################
# Colors
if(is.null(tcol)  & wall==0)   tcol<-1
if(is.null(tcol)  & wall!=0)   tcol<-rgb(0.1,0.1,0.25)
if(is.null(bgcol) & wall==0)   bgcol<-1
if(is.null(bgcol) & wall!=0)   bgcol<-'#DBE0E8'
if(is.null(alpha))  alpha<-1

if(is.null(lcol))  lcol<- tcol
if(is.null(density)) density<-rep(20, 20)
if(is.null(xlab)) xlab <-xname
if(is.null(xticks)) xticks<-5
################################################################

################################################################
a.scale<-function(x, from, xmax){
#########################################
# make single values
if(length(x)>1){
y<-x
for(i in 1:length(x)){
y[i] <- a.scale(x[i], from=from, xmax=xmax)
}
return(y)
}
#########################################
if(length(x)==1){
if(x>=0){
x<-x/xmax
x<- (x*(1-from))+from
}
if(x<0){
x<-x/xmax
x<- (x*(1-from))-from
}
return(x)
}
}
################################################################



################################################################
# x ist continuirlich
if(length(unique(x))>8 & !is.table(x) & !is.matrix(x) & !is.factor(x)){
xrange<-range(x, na.rm=TRUE)
group<-as.factor(group)
if(is.null(gnames))  gnames<-levels(group)
########################################
counts1<-NULL
counts2<-NULL
xbreaks<-NULL
if(!is.null(group2)) g2<-group2
if(is.null(group2))  g2<-as.factor(rep(1, length(x)))
if(is.null(gnames2) & !is.null(group2)) gnames2<-levels(group2)
n.gs<-nlevels(g2)
if(nlevels(group)<2) stop('need a group factor')
for(i in 1:nlevels(g2)){
x1<-hist(x[group==levels(group)[1] & g2==levels(g2)[i]], breaks=seq(min(x, na.rm=T), max(x, na.rm=T)+(max(x, na.rm=T)%%breaks), breaks), plot = FALSE)
x2<-hist(x[group==levels(group)[2] & g2==levels(g2)[i]], breaks=seq(min(x, na.rm=T), max(x, na.rm=T)+(max(x, na.rm=T)%%breaks), breaks), plot = FALSE)
counts1[[i]]<- x1$counts
counts2[[i]]<- x2$counts
if(is.null(xbreaks)) xbreaks<- x1$breaks
if(!is.null(xbreaks)) xbreaks<- c(x1$breaks[x1$breaks<min(xbreaks)] ,xbreaks, x1$breaks[x1$breaks>max(xbreaks)])
if(!is.null(xbreaks)) xbreaks<- c(x2$breaks[x2$breaks<min(xbreaks)] ,xbreaks, x2$breaks[x2$breaks>max(xbreaks)])
}

if(is.null(vnames)){
vnames<-paste('[', format_n.ade(x1$breaks[-length(x1$breaks)], 2), '-', format_n.ade(x1$breaks[-1], 2), ')', sep='')
}
nticks<- max(nchar(vnames))

yrun  <- x1$mids
xmax<-max(c(unlist(counts1), unlist(counts2)), na.rm=TRUE)
col1<-col
col2<-col
if(is.na(n.gs)) n.gs<-1
if(is.null(col) & n.gs==1){
col1 <- 'cornflowerblue'
col2 <- 'brown2'
}
if(is.null(col) & n.gs>1){
col1 <- a.getcol.ade(n.gs, type='p')
col2 <- col1
col<-col1
}
col1<- a.alpha.ade(col1, alpha=alpha)
col2<- a.alpha.ade(col2, alpha=alpha)

#################
# sortieren
l1 <- unlist(lapply(counts1, sum))
l2 <- unlist(lapply(counts2, sum))
o<-order((l1+l2), decreasing = TRUE)
counts1<-counts1[o]
counts2<-counts2[o]
#col1<-col1[o]
#col2<-col2[o]
#density<-density[o]
#angle<-angle[o]
#################
}
################################################################

################################################################
# x ist Factor
if((length(unique(x))<=8 | is.factor(x)) & !is.matrix(x) & !is.table(x) & !is.list(x)){
x<-as.factor(x)
group<-as.factor(group)
xrange<-c(0.5, nlevels(x)+0.5)
xbreaks<-seq(0.5, (nlevels(x)+0.5), 1)
if(is.null(gnames))  gnames<-levels(group)
########################################
counts1<-NULL
counts2<-NULL

if(!is.null(group2)) g2<-group2
if(is.null(group2))  g2<-as.factor(rep(1, length(x)))
if(is.null(gnames2) & !is.null(group2)) gnames2<-levels(group2)
n.gs<-nlevels(g2)

for(i in 1:nlevels(g2)){
counts1[[i]]<- table(x[group==levels(group)[1] & g2==levels(g2)[i]])
counts2[[i]]<- table(x[group==levels(group)[2] & g2==levels(g2)[i]])
}

if(is.null(vnames)){
vnames<-levels(x)
}
nticks<- max(nchar(vnames))

yrun  <- 1:nlevels(x)
xmax<-max(c(max(unlist(counts1), na.rm=TRUE), max(unlist(counts2), na.rm=TRUE)), na.rm=TRUE)
col1<-col
col2<-col
if(is.na(n.gs)) n.gs<-1
if(is.null(col) & n.gs==1){
col1 <- 'cornflowerblue'
col2 <- 'brown2'
}
if(is.null(col) & n.gs>1){
col1 <- a.getcol.ade(n.gs, type='p')
col2 <- col1
col<-col1
}
col1<- a.alpha.ade(col1, alpha=alpha)
col2<- a.alpha.ade(col2, alpha=alpha)

#################
# sortieren
l1 <- unlist(lapply(counts1, sum))
l2 <- unlist(lapply(counts2, sum))
o<-order((l1+l2), decreasing = TRUE)
counts1<-counts1[o]
counts2<-counts2[o]
#col1<-col1[o]
#col2<-col2[o]
#density<-density[o]
#angle<-angle[o]

#################
}
################################################################

################################################################
# x Table  and matrix
if((is.table(x) | is.matrix(x))){
xrange<-c(0.5, dim(x)[1]+0.5)
if(is.null(gnames))  gnames<-colnames(x)
########################################
n.gs<-dim(x)[3]
counts1<-NULL
counts2<-NULL
if(!is.na(n.gs)){
for(i in 1:n.gs){
counts1[[i]]<- x[,1,i]
counts2[[i]]<- x[,2,i]
}
}

if(is.na(n.gs)){
counts1[[1]]<- x[,1]
counts2[[1]]<- x[,2]
}


yrun  <- 1:(dim(x)[1])
if(is.null(vnames)) vnames<- rownames(x)
nticks<- max(nchar(vnames))
xmax<-max(x, na.rm=TRUE)
xbreaks<-seq(0.5, (dim(x)[1]+0.5), 1)
col1<-col
col2<-col
if(!is.null(col) & length(col)==2 & is.na(n.gs)){
col1<-col[1]
col2<-col[2]
}

if(is.na(n.gs)) n.gs<-1
if(is.null(col) & n.gs==1){
col1 <- 'cornflowerblue'
col2 <- 'brown2'
}
if(is.null(col) & n.gs>1){
col1 <- a.getcol.ade(n.gs, type='p')
col2 <- col1
col<-col1
}
col1<- a.alpha.ade(col1, alpha=alpha)
col2<- a.alpha.ade(col2, alpha=alpha)

#################
# sortieren
l1 <- unlist(lapply(counts1, sum))
l2 <- unlist(lapply(counts2, sum))
o<-order((l1+l2), decreasing = TRUE)
counts1<-counts1[o]
counts2<-counts2[o]
col1<-col1[o]
col2<-col2[o]
#density<-density[o]
#angle<-angle[o]
#################
}
################################################################

################################################################
# x List of tables
if(is.list(x)){
lmax<-NULL
xlmax<-NULL
for(k in 1:length(x)){
lmax<-max(lmax, dim(x[[k]])[1])
xlmax<-max(xlmax, max(x[[k]], na.rm=TRUE))
}
xrange<-c(0.5, lmax+0.5)

########################################
n.gs<-length(x)
counts1<-NULL
counts2<-NULL

for(i in 1:n.gs){
counts1[[i]]<- x[[i]][,1]
counts2[[i]]<- x[[i]][,2]
}

yrun  <- 1:lmax
if(is.null(vnames)) vnames<- yrun
nticks<- max(nchar(vnames))
xmax<-xlmax
xbreaks<-seq(0.5, (lmax+0.5), 1)
col1<-col
col2<-col
if(is.na(n.gs)) n.gs<-1
if(is.null(col) & n.gs==1){
col1 <- 'cornflowerblue'
col2 <- 'brown2'
}
if(is.null(col) & n.gs>1){
col1 <- a.getcol.ade(n.gs, type='p')
col2 <- col1
}
col1<- a.alpha.ade(col1, alpha=alpha)
col2<- a.alpha.ade(col2, alpha=alpha)

#################
# sortieren
l1 <- unlist(lapply(counts1, sum))
l2 <- unlist(lapply(counts2, sum))
o<-order((l1+l2), decreasing = TRUE)
counts1<-counts1[o]
counts2<-counts2[o]
col1<-col1[o]
col2<-col2[o]
if(!is.null(gnames2)) gnames2<-gnames2[o]
density<-density[o]
angle<-angle[o]
#################
}
################################################################



################################################################################
################################################################################
if(wall==0){
################################################################################
par(col.main=tcol)
par(col.axis=tcol)
par(col.lab=tcol)
par(mai=c(1.02, 0.42, 0.875, 0.42))
par(cex.axis=0.8)

         
plot(0, axes=F, col=rgb(0,0,0,0), xlim=c(-1-r, 1+r), ylim=range(xbreaks), xlab='', ylab='', main='')
breite<-max((strwidth(vnames, units = "user", cex.axis=0.8)/2.4))+0.055
llim<-a.glc(side=0)-breite
rlim<-a.glc(side=0)+breite
rect(llim,a.glc(3,0), rlim, a.glc(1,0), border=bgcol)
text(a.glc(side=0), yrun, labels=vnames, font=1,  col=tcol, cex=0.8)

segments(llim, yrun, llim+(a.glc(2, 0)-a.glc(2, 0.4)), yrun,  col=bgcol, lwd=1)
segments(rlim, yrun, rlim-(a.glc(2, 0)-a.glc(2, 0.4)), yrun,  col=bgcol, lwd=1)


xlabs<-pretty(c(0, xmax), n=xticks)
xrun<-xlabs/xmax
xrun<-(xrun*(1-rlim))+rlim
axis(1, at=xrun,  labels=xlabs, col=bgcol)
axis(1, at=-xrun, labels=xlabs, col=bgcol)



text(x=llim-0.025, y=a.glc(3,0.75), labels=gnames[1], col=tcol, adj=c(1, 0.5), xpd=TRUE, font=par("font.main"))
text(x=rlim+0.025, y=a.glc(3,0.75), labels=gnames[2], col=tcol, adj=c(0, 0.5), xpd=TRUE, font=par("font.main"))

for(i in 1:n.gs){
rect( -a.scale(counts1[[i]], rlim, xmax), xbreaks[-length(xbreaks)], llim ,                             xbreaks[-1], col=col1[i], border=a.alpha.ade(a.coladd.ade(col1[i], -75)), density=density[i], angle=angle[i])
rect( rlim,                               xbreaks[-length(xbreaks)], a.scale(counts2[[i]], rlim, xmax), xbreaks[-1], col=col2[i], border=a.alpha.ade(a.coladd.ade(col2[i], -75)), density=density[i], angle=angle[i])
}

####################
# Labels
text(x=a.glc(0), y=a.glc(3, 2.8), labels=main, cex=par('cex.main'), col=tcol, adj=c(0.5, 0.5), xpd=TRUE, font=par("font.main"))
text(x=a.glc(0), y=a.glc(1, 3.5), labels=xlab, cex=par('cex.lab'),  col=tcol, adj=c(0.5, 0.5), xpd=TRUE, font=par("font.lab"))
####################

#Legend
if(n.gs>1 & !is.null(gnames2))  legend(legendon, legend=gnames2, fill = col1,  title=glab, border=bgcol, box.lwd=1, box.col=bgcol, text.col=tcol, bg=rgb(1,1,1, 0), density=density*2, angle=angle, text.width=max(strwidth(c(gnames2, glab),font = 2)))

abline(v=a.scale(v, rlim, xmax), h=h, col=lcol, lty=lty, lwd=lwd)
box(col=bgcol)
}
################################################################################
################################################################################

################################################################################
################################################################################
if(wall==1){
################################################################################
par(col.main=tcol)
par(col.axis=tcol)
par(col.lab=tcol)
par(mai=c(1.02, 0.42, 0.875, 0.42))
par(cex.axis=0.8)

plot(0, axes=F, col=rgb(0,0,0,0), xlim=c(-1-r, 1+r), ylim=range(xbreaks), xlab='', ylab='', main='')
polygon( c(par('usr')[c(1,1,2,2)]), par('usr')[c(3,4,4,3)], col=bgcol, border=FALSE)

breite<-max((strwidth(vnames, units = "user", cex.axis=0.8)/2.4))+0.055
llim<-a.glc(side=0)-breite
rlim<-a.glc(side=0)+breite
rect(llim,a.glc(3,0), rlim, a.glc(1,0), border=rgb(1,1,1))
text(a.glc(side=0), yrun, labels=vnames, font=1,  col=tcol, cex=0.8)

segments(llim, yrun, llim+(a.glc(2, 0)-a.glc(2, 0.4)), yrun,  col=rgb(1,1,1), lwd=1)
segments(rlim, yrun, rlim-(a.glc(2, 0)-a.glc(2, 0.4)), yrun,  col=rgb(1,1,1), lwd=1)
xlabs<-pretty(c(0, xmax), n=xticks)
xrun<-xlabs/xmax
xrun<-(xrun*(1-rlim))+rlim
a1<-axis(1, at=xrun,  labels=xlabs, col=bgcol)
a2<-axis(1, at=-xrun, labels=xlabs, col=bgcol)
abline(v=c(a1, a2), lty=1, col=rgb(1,1,1), lwd=1)


text(x=llim-0.025, y=a.glc(3,0.75), labels=gnames[1], col=tcol, adj=c(1, 0.5), xpd=TRUE, font=par("font.main"))
text(x=rlim+0.025, y=a.glc(3,0.75), labels=gnames[2], col=tcol, adj=c(0, 0.5), xpd=TRUE, font=par("font.main"))

for(i in 1:n.gs){
rect( -a.scale(counts1[[i]], rlim, xmax), xbreaks[-length(xbreaks)], llim ,                             xbreaks[-1], col=col1[i], border=a.alpha.ade(a.coladd.ade(col1[i], -75)), density=density[i], angle=angle[i])
rect( rlim,                               xbreaks[-length(xbreaks)], a.scale(counts2[[i]], rlim, xmax), xbreaks[-1], col=col2[i], border=a.alpha.ade(a.coladd.ade(col2[i], -75)), density=density[i], angle=angle[i])
}

####################
# Labels
text(x=a.glc(0), y=a.glc(3, 2.8), labels=main, cex=par('cex.main'), col=tcol, adj=c(0.5, 0.5), xpd=TRUE, font=par("font.main"))
text(x=a.glc(0), y=a.glc(1, 3.5), labels=xlab, cex=par('cex.lab'),  col=tcol, adj=c(0.5, 0.5), xpd=TRUE, font=par("font.lab"))
####################

#Legend
if(n.gs>1 & !is.null(gnames2))  legend(legendon, legend=gnames2, fill = col1, title=glab, border=col1, box.col=rgb(1,1,1),  box.lwd=2, text.col=tcol, bg=bgcol, density=density*2, angle=angle,text.width=max(strwidth(c(gnames2, glab),font = 2)))
abline(v=a.scale(v, rlim, xmax), h=h, col=lcol, lty=lty, lwd=lwd)
}
################################################################################
################################################################################


################################################################################
################################################################################
if(wall==2){
################################################################################
par(col.main=tcol)
par(col.axis=tcol)
par(col.lab=tcol)
par(mai=c(1.02, 0.42, 0.875, 0.42))
par(cex.axis=0.8)

plot(0, axes=F, col=rgb(0,0,0,0), xlim=c(-1-r, 1+r), ylim=range(xbreaks), xlab='', ylab='', main='')
breite<-max((strwidth(vnames, units = "user", cex.axis=0.8)/2.4))+0.055
llim<-a.glc(side=0)-breite
rlim<-a.glc(side=0)+breite
text(a.glc(side=0), yrun, labels=vnames, font=1,  col=tcol, cex=0.8)

segments(llim, yrun, llim+(a.glc(2, 0)-a.glc(2, 0.4)), yrun,  col=a.coladd.ade(bgcol, -75), lwd=1)
segments(rlim, yrun, rlim-(a.glc(2, 0)-a.glc(2, 0.4)), yrun,  col=a.coladd.ade(bgcol, -75), lwd=1)
xlabs<-pretty(c(0, xmax), n=xticks)
xrun<-xlabs/xmax
xrun<-(xrun*(1-rlim))+rlim
a1<-axis(1, at=xrun,  labels=xlabs, col=bgcol, col.ticks = a.coladd.ade(bgcol, -75))
a2<-axis(1, at=-xrun, labels=xlabs, col=bgcol, col.ticks = a.coladd.ade(bgcol, -75))
abline(v=c(a1, a2), lty=1, col=bgcol, lwd=1)
rect(llim,a.glc(3,0), rlim, a.glc(1,0), border=a.coladd.ade(bgcol, -75))

text(x=llim-0.025, y=a.glc(3,0.75), labels=gnames[1], col=tcol, adj=c(1, 0.5), xpd=TRUE, font=par("font.main"))
text(x=rlim+0.025, y=a.glc(3,0.75), labels=gnames[2], col=tcol, adj=c(0, 0.5), xpd=TRUE, font=par("font.main"))

for(i in 1:n.gs){
rect( -a.scale(counts1[[i]], rlim, xmax), xbreaks[-length(xbreaks)], llim ,                             xbreaks[-1], col=col1[i], border=a.alpha.ade(a.coladd.ade(col1[i], -75)), density=density[i], angle=angle[i])
rect( rlim,                               xbreaks[-length(xbreaks)], a.scale(counts2[[i]], rlim, xmax), xbreaks[-1], col=col2[i], border=a.alpha.ade(a.coladd.ade(col2[i], -75)), density=density[i], angle=angle[i])
}

####################
# Labels
text(x=a.glc(0), y=a.glc(3, 2.8), labels=main, cex=par('cex.main'), col=tcol, adj=c(0.5, 0.5), xpd=TRUE, font=par("font.main"))
text(x=a.glc(0), y=a.glc(1, 3.5), labels=xlab, cex=par('cex.lab'),  col=tcol, adj=c(0.5, 0.5), xpd=TRUE, font=par("font.lab"))
####################

#Legend
if(n.gs>1 & !is.null(gnames2))  legend(legendon, legend=gnames2, fill = col1, title=glab, border=a.coladd.ade(bgcol, -75), box.col=a.coladd.ade(bgcol, -75),  box.lwd=1, density=density*2, angle=angle, text.col=tcol, bg=rgb(1,1,1), text.width=max(strwidth(c(gnames2, glab),font = 2)))

abline(v=a.scale(v, rlim, xmax), h=h, col=lcol, lty=lty, lwd=lwd)
box(col=a.coladd.ade(bgcol, -75))
}
################################################################################
################################################################################


################################################################################
################################################################################
if(wall==3){
################################################################################
par(col.main=tcol)
par(col.axis=tcol)
par(col.lab=tcol)
par(mai=c(1.02, 0.42, 0.875, 0.42))
par(cex.axis=0.8)

plot(0, axes=F, col=rgb(0,0,0,0), xlim=c(-1-r, 1+r), ylim=range(xbreaks), xlab='', ylab='', main='')
polygon( c(par('usr')[c(1,1,2,2)]), par('usr')[c(3,4,4,3)], col=bgcol, border=FALSE)
breite<-max((strwidth(vnames, units = "user", cex.axis=0.8)/2.4))+0.055
llim<-a.glc(side=0)-breite
rlim<-a.glc(side=0)+breite
rect(llim,a.glc(3,0), rlim, a.glc(1,0), col=rgb(1,1,1), border=a.coladd.ade(bgcol, -50))
text(a.glc(side=0), yrun, labels=vnames, font=1,  col=tcol, cex=0.8)

segments(llim, yrun, llim+(a.glc(2, 0)-a.glc(2, 0.4)), yrun,  col=a.coladd.ade(bgcol, -50), lwd=1)
segments(rlim, yrun, rlim-(a.glc(2, 0)-a.glc(2, 0.4)), yrun,  col=a.coladd.ade(bgcol, -50), lwd=1)
xlabs<-pretty(c(0, xmax), n=xticks)
xrun<-xlabs/xmax
xrun<-(xrun*(1-rlim))+rlim
a1<-axis(1, at=xrun,  labels=xlabs, col=bgcol, col.ticks = a.coladd.ade(bgcol, -50))
a2<-axis(1, at=-xrun, labels=xlabs, col=bgcol, col.ticks = a.coladd.ade(bgcol, -50))
abline(v=c(a1, a2), lty=1, col=a.coladd.ade(bgcol, -50), lwd=1)

text(x=llim-0.025, y=a.glc(3,0.75), labels=gnames[1], col=tcol, adj=c(1, 0.5), xpd=TRUE, font=par("font.main"))
text(x=rlim+0.025, y=a.glc(3,0.75), labels=gnames[2], col=tcol, adj=c(0, 0.5), xpd=TRUE, font=par("font.main"))

for(i in 1:n.gs){
rect( -a.scale(counts1[[i]], rlim, xmax), xbreaks[-length(xbreaks)], llim ,                             xbreaks[-1], col=col1[i], border=a.alpha.ade(a.coladd.ade(col1[i], -75)), density=density[i], angle=angle[i])
rect( rlim,                               xbreaks[-length(xbreaks)], a.scale(counts2[[i]], rlim, xmax), xbreaks[-1], col=col2[i], border=a.alpha.ade(a.coladd.ade(col2[i], -75)), density=density[i], angle=angle[i])
}

####################
# Labels
text(x=a.glc(0), y=a.glc(3, 2.8), labels=main, cex=par('cex.main'), col=tcol, adj=c(0.5, 0.5), xpd=TRUE, font=par("font.main"))
text(x=a.glc(0), y=a.glc(1, 3.5), labels=xlab, cex=par('cex.lab'),  col=tcol, adj=c(0.5, 0.5), xpd=TRUE, font=par("font.lab"))
####################

#Legend
if(n.gs>1 & !is.null(gnames2))  legend(legendon, legend=gnames2, fill = col1, title=glab, border=a.coladd.ade(bgcol, -50), box.col=a.coladd.ade(bgcol, -50),  box.lwd=1, text.col=tcol, density=density*2, angle=angle,bg=rgb(1,1,1), text.width=max(strwidth(c(gnames2, glab),font = 2)))

abline(v=a.scale(v, rlim, xmax), h=h, col=lcol, lty=lty, lwd=lwd)
box(col=a.coladd.ade(bgcol, -50))
}
################################################################################
################################################################################


################################################################################
################################################################################
if(wall==4){
################################################################################
par(col.main=tcol)
par(col.axis=tcol)
par(col.lab=tcol)
par(mai=c(1.02, 0.42, 0.875, 0.42))
par(cex.axis=0.8)
par(font=2)

plot(0, axes=F, col=rgb(0,0,0,0), xlim=c(-1-r, 1+r), ylim=range(xbreaks), xlab='', ylab='', main='')
polygon( c(par('usr')[c(1,1,2,2)]), par('usr')[c(3,4,4,3)], col=bgcol, border=FALSE)
breite<-max((strwidth(vnames, units = "user", cex.axis=0.8)/2.4))+0.055
llim<-a.glc(side=0)-breite
rlim<-a.glc(side=0)+breite
rect(llim,a.glc(3,0), rlim, a.glc(1,0), col=rgb(1,1,1), border=a.coladd.ade(bgcol, -50))
text(a.glc(side=0), yrun, labels=vnames, font=1,  col=tcol, cex=0.8)

segments(llim, yrun, llim+(a.glc(2, 0)-a.glc(2, 0.4)), yrun,  col=tcol, lwd=1)
segments(rlim, yrun, rlim-(a.glc(2, 0)-a.glc(2, 0.4)), yrun,  col=tcol, lwd=1)
xlabs<-pretty(c(0, xmax), n=xticks)
xrun<-xlabs/xmax
xrun<-(xrun*(1-rlim))+rlim
a1<-axis(1, at=xrun,  labels=xlabs, col=bgcol, col.ticks = tcol)
a2<-axis(1, at=-xrun, labels=xlabs, col=bgcol, col.ticks = tcol)
abline(v=c(a1, a2), lty=1, col=rgb(1,1,1), lwd=1)

text(x=llim-0.025, y=a.glc(3,0.75), labels=gnames[1], col=tcol, adj=c(1, 0.5), xpd=TRUE, font=par("font.main"))
text(x=rlim+0.025, y=a.glc(3,0.75), labels=gnames[2], col=tcol, adj=c(0, 0.5), xpd=TRUE, font=par("font.main"))

for(i in 1:n.gs){
rect( -a.scale(counts1[[i]], rlim, xmax), xbreaks[-length(xbreaks)], llim ,                             xbreaks[-1], col=col1[i], border=a.alpha.ade(a.coladd.ade(col1[i], -75)), density=density[i], angle=angle[i])
rect( rlim,                               xbreaks[-length(xbreaks)], a.scale(counts2[[i]], rlim, xmax), xbreaks[-1], col=col2[i], border=a.alpha.ade(a.coladd.ade(col2[i], -75)), density=density[i], angle=angle[i])
}

par(xpd=TRUE)
polygon(a.glc(side=c(2,2,4,4), line=c(0,0,0,0)), a.glc(side=3, line=c(1.2, 3.5,  3.5, 1.2)), col=tcol, border=rgb(1,1,1))
if(!is.null(xlab)) if(xlab!='' & xlab!=' ') polygon( a.glc(side=c(2, 2, 4, 4), line=0),   a.glc(side=1, line=c(4, 2.5, 2.5, 4)), col=bgcol, border=rgb(1,1,1))
par(xpd=FALSE)


####################
# Labels
text(x=a.glc(0), y=a.glc(3, 2.5), labels=main, cex=par('cex.main'), col=rgb(1,1,1), adj=c(0.5, 0.5), xpd=TRUE, font=par("font.main"))
text(x=a.glc(0), y=a.glc(1, 3.15), labels=xlab, cex=par('cex.lab'),  col=tcol, adj=c(0.5, 0.5), xpd=TRUE, font=2)
####################

#Legend
if(n.gs>1 & !is.null(gnames2))  legend(legendon, legend=gnames2, fill = a.alpha.ade(col1, 1), text.col=rgb(1,1,1), title=glab, border=rgb(1,1,1), box.col=rgb(1,1,1),  box.lwd=1, bg=tcol, density=density*2, angle=angle,text.width=max(strwidth(c(gnames2, glab),font = 2)))

abline(v=a.scale(v, rlim, xmax), h=h, col=lcol, lty=lty, lwd=lwd)
box(col=rgb(1,1,1))
}
################################################################################
################################################################################


################################################################################
################################################################################
if(wall==5){
################################################################################
par(col.main=tcol)
par(col.axis=tcol)
par(col.lab=tcol)
par(mai=c(1.02, 0.42, 0.875, 0.42))
par(cex.axis=0.8)

plot(0, axes=F, col=rgb(0,0,0,0), xlim=c(-1-r, 1+r), ylim=range(xbreaks), xlab='', ylab='', main='')
breite<-max((strwidth(vnames, units = "user", cex.axis=0.8)/2.4))+0.055
llim<-a.glc(side=0)-breite
rlim<-a.glc(side=0)+breite
text(a.glc(side=0), yrun, labels=vnames, font=1,  col=tcol, cex=0.8)

segments(llim, yrun, llim+(a.glc(2, 0)-a.glc(2, 0.4)), yrun,  col=tcol, lwd=1)
segments(rlim, yrun, rlim-(a.glc(2, 0)-a.glc(2, 0.4)), yrun,  col=tcol, lwd=1)
xlabs<-pretty(c(0, xmax), n=xticks)
xrun<-xlabs/xmax
xrun<-(xrun*(1-rlim))+rlim
a1<-axis(1, at=xrun,  labels=xlabs, col=bgcol, col.ticks = tcol)
a2<-axis(1, at=-xrun, labels=xlabs, col=bgcol, col.ticks = tcol)
rect(llim,a.glc(3,0), rlim, a.glc(1,0), border=tcol)

text(x=llim-0.025, y=a.glc(3,0.75), labels=gnames[1], col=tcol, adj=c(1, 0.5), xpd=TRUE, font=par("font.main"))
text(x=rlim+0.025, y=a.glc(3,0.75), labels=gnames[2], col=tcol, adj=c(0, 0.5), xpd=TRUE, font=par("font.main"))

for(i in 1:n.gs){
rect( -a.scale(counts1[[i]], rlim, xmax), xbreaks[-length(xbreaks)], llim ,                             xbreaks[-1], col=col1[i], border=a.alpha.ade(a.coladd.ade(col1[i], -75)), density=density[i], angle=angle[i])
rect( rlim,                               xbreaks[-length(xbreaks)], a.scale(counts2[[i]], rlim, xmax), xbreaks[-1], col=col2[i], border=a.alpha.ade(a.coladd.ade(col2[i], -75)), density=density[i], angle=angle[i])
}

par(xpd=TRUE)
polygon(a.glc(side=2, line=c(0.6, 0.6, 0, 0)), a.glc(side=3, line=c(1.6, 4, 4, 1.6)), col=bgcol,        border=tcol)
polygon(a.glc(side=c(2,2,4,4), line=c(0,0,0,0)), a.glc(side=3, line=c(1.6, 4, 4, 1.6)), col=rgb(1,1,1,0), border=tcol)
polygon(a.glc(side=4, line=c(0, 0 ,0.6, 0.6)),   a.glc(side=3, line=c(1.6, 4, 4, 1.6)), col=bgcol,        border=tcol)
polygon(a.glc(side=2, line=c(0.6, 0.6 ,0, 0)),  a.glc(side=c(1,3,3,1), line=c(0, 0, 0, 0)), col=bgcol,  border=tcol)
polygon(a.glc(side=4, line=c(0, 0 ,0.6, 0.6)), a.glc(side=c(1, 3, 3, 1), line=0), col=bgcol, border=tcol)
polygon(a.glc(side=2, line=c(0.6, 0.6, 0, 0)), a.glc(side=1, line=c(2.6, 4.5, 4.5, 2.6)), col=bgcol, border=tcol)
polygon(a.glc(side=c(2, 2, 4, 4), line=0), a.glc(side=1, line=c(2.6, 4.5, 4.5, 2.6)), col=rgb(1,1,1,0), border=tcol)
polygon(a.glc(side=4, line=c(0, 0, 0.6, 0.6)), a.glc(side=1, line=c(2.6, 4.5, 4.5, 2.6)), col=bgcol, border=tcol)
par(xpd=FALSE)

####################
# Labels
text(x=a.glc(0), y=a.glc(3, 2.8), labels=main, cex=par('cex.main'), col=tcol, adj=c(0.5, 0.5), xpd=TRUE, font=par("font.main"))
text(x=a.glc(0), y=a.glc(1, 3.5), labels=xlab, cex=par('cex.lab'),  col=tcol, adj=c(0.5, 0.5), xpd=TRUE, font=par("font.lab"))
####################

#Legend
if(n.gs>1 & !is.null(gnames2))  legend(legendon, legend=gnames2, fill = col1, title=glab, border=tcol, box.col=tcol,  box.lwd=1, text.col=tcol, bg=rgb(1,1,1), density=density*2, angle=angle, text.width=max(strwidth(c(gnames2, glab),font = 2)))

abline(v=a.scale(v, rlim, xmax), h=h, col=lcol, lty=lty, lwd=lwd)
box(col=tcol)
}
################################################################################
################################################################################


################################################################################
################################################################################
if(wall==6){
################################################################################
par(col.main=tcol)
par(col.axis=tcol)
par(col.lab=tcol)
par(mai=c(1.02, 0.42, 0.875, 0.42))
par(cex.axis=0.8)

plot(0, axes=F, col=rgb(0,0,0,0), xlim=c(-1-r, 1+r), ylim=range(xbreaks), xlab='', ylab='', main='')
polygon( c(par('usr')[c(1,1,2,2)]), par('usr')[c(3,4,4,3)], col=bgcol, border=FALSE)
breite<-max((strwidth(vnames, units = "user", cex.axis=0.8)/2.4))+0.055
llim<-a.glc(side=0)-breite
rlim<-a.glc(side=0)+breite
rect(llim,a.glc(3,0), rlim, a.glc(1,0), col=rgb(1,1,1), border=a.coladd.ade(bgcol, -50))
text(a.glc(side=0), yrun, labels=vnames, font=1,  col=tcol, cex=0.8)

segments(llim, yrun, llim+(a.glc(2, 0)-a.glc(2, 0.4)), yrun,  col=a.coladd.ade(bgcol, -35), lwd=3)
segments(llim, yrun, llim+(a.glc(2, 0)-a.glc(2, 0.4)), yrun,  col=rgb(1,1,1), lwd=1)
segments(rlim, yrun, rlim-(a.glc(2, 0)-a.glc(2, 0.4)), yrun,  col=a.coladd.ade(bgcol, -35), lwd=3)
segments(rlim, yrun, rlim-(a.glc(2, 0)-a.glc(2, 0.4)), yrun,  col=rgb(1,1,1), lwd=1)
xlabs<-pretty(c(0, xmax), n=xticks)
xrun<-xlabs/xmax
xrun<-(xrun*(1-rlim))+rlim
a1<-axis(1, at=xrun,  labels=xlabs, col=bgcol, col.ticks=a.coladd.ade(bgcol, -35), lwd.ticks=3)
a1<-axis(1, at=xrun,  labels=xlabs, col=bgcol, col.ticks=rgb(1,1,1),               lwd.ticks=1,)
a2<-axis(1, at=-xrun, labels=xlabs, col=bgcol, col.ticks=a.coladd.ade(bgcol, -35), lwd.ticks=3)
a2<-axis(1, at=-xrun, labels=xlabs, col=bgcol, col.ticks=rgb(1,1,1),               lwd.ticks=1,)

abline(v=c(a1, a2), lty=1, col=a.coladd.ade(bgcol, -35), lwd=3)
abline(v=c(a1, a2), lty=1, col=rgb(1,1,1), lwd=1)


text(x=llim-0.025, y=a.glc(3,0.75), labels=gnames[1], col=tcol, adj=c(1, 0.5), xpd=TRUE, font=par("font.main"))
text(x=rlim+0.025, y=a.glc(3,0.75), labels=gnames[2], col=tcol, adj=c(0, 0.5), xpd=TRUE, font=par("font.main"))

for(i in 1:n.gs){
rect( -a.scale(counts1[[i]], rlim, xmax), xbreaks[-length(xbreaks)], llim ,                             xbreaks[-1], col=col1[i], border=a.alpha.ade(a.coladd.ade(col1[i], -75)), density=density[i], angle=angle[i])
rect( rlim,                               xbreaks[-length(xbreaks)], a.scale(counts2[[i]], rlim, xmax), xbreaks[-1], col=col2[i], border=a.alpha.ade(a.coladd.ade(col2[i], -75)), density=density[i], angle=angle[i])
}

####################
# Labels
text(x=a.glc(0), y=a.glc(3, 2.8), labels=main, cex=par('cex.main'), col=tcol, adj=c(0.5, 0.5), xpd=TRUE, font=par("font.main"))
text(x=a.glc(0), y=a.glc(1, 3.5), labels=xlab, cex=par('cex.lab'),  col=tcol, adj=c(0.5, 0.5), xpd=TRUE, font=par("font.lab"))
####################

#Legend
if(n.gs>1 & !is.null(gnames2))  legend(legendon, legend=gnames2, fill = col1, title=glab, border=a.coladd.ade(bgcol, -50), box.col=a.coladd.ade(bgcol, -35),  box.lwd=3, density=density*2, angle=angle, text.col=tcol, bg=bgcol, text.width=max(strwidth(c(gnames2, glab),font = 2)))
if(n.gs>1 & !is.null(gnames2))  legend(legendon, legend=gnames2, fill = col1, title=glab, border=a.coladd.ade(bgcol, -50), box.col=rgb(1,1,1),  box.lwd=1, text.col=tcol, density=density*2, angle=angle, bg=rgb(1,1,1,0), text.width=max(strwidth(c(gnames2, glab),font = 2)))

abline(v=a.scale(v, rlim, xmax), h=h, col=lcol, lty=lty, lwd=lwd)
box(lwd=3, col=rgb(1,1,1))
box(lwd=1, col=a.coladd.ade(bgcol, -35))
}
################################################################################
################################################################################

}
