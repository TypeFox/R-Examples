correlogram.ade <-
function(vars1,  vnames1='noname', vars2,  vnames2='noname', prediktors=0, data=NULL , xlab=NULL, ylab=NULL, main=NULL, method="p", digits=2, pdigs=4, pvals=TRUE, bars=TRUE, col=NULL, tcol=NULL,  bgcol=NULL, wall=0){
if(any(par('mfg')!=c(1,1,1,1)) & any(par('mai') < c(1.02, 0.82, 0.82, 0.42))){
maidiff<-rep(0, 4)
norm<-c(1.02, 0.82, 0.82, 0.42)
maidiff[par('mai')<norm]<-  norm[par('mai')<norm] - par('mai')[par('mai')<norm]
par(mai=par('mai')+maidiff)
}
oldpar<-par(no.readonly =TRUE)
oldpar<-oldpar[-which(names(oldpar)%in%c('usr',  'plt', 'pin',   'fin', 'fig', 'mfg', 'mfcol', 'mfrow', 'omd', 'omi', 'oma'))]
on.exit(par(oldpar))


if (any(vnames1=='noname')) vnames1 <- vars1
if (any(vnames2=='noname')) vnames2 <- vars2
if(!is.data.frame(data))  stop("(data) must be a data.frame!")
if (length(vars1)>26) stop('vars1 to long!')
if (length(vars2)>26) stop('vars2 to long!')
x2lab <- xlab
xlab  <- main
N1 <- length(vars1)
N2 <- length(vars2)

W  <-matrix(data = NA, nrow = N1, ncol = N2*2)
WW <-matrix(data = NA, nrow = N1, ncol = N2*2)
M  <-matrix(data = NA, nrow = N1, ncol = N2)
MM <-matrix(data = NA, nrow = N1, ncol = N2)
P  <-matrix(data = NA, nrow = N1, ncol = N2) 
PP <-matrix(data = NA, nrow = N1, ncol = N2) 
rownames(W)   <- vnames1
colnames(W)  <- c(vnames2, vnames2  )

rownames(WW)  <- vnames1
colnames(WW)  <- c(vnames2, vnames2  )

rownames(M)  <- vnames1
colnames(M)  <- vnames2

rownames(MM) <- vnames1
colnames(MM) <- vnames2

rownames(P)  <- vnames1
colnames(P)  <- paste('p -', vnames2 )

rownames(PP) <- vnames1
colnames(PP) <- paste('p -', vnames2 )


partial=TRUE
if(any(prediktors==0)) partial=FALSE

if(partial){
nl<- length(prediktors)
p <- NULL
for(i in 1:nl){
top <-eval(parse(text=paste("data$",prediktors[i], sep='')))
p <- cbind(p, top)
}
}




if(is.null(tcol)  & wall==0)   tcol<-1
if(is.null(tcol)  & wall!=0)   tcol<-rgb(0.1,0.1,0.25)
if(is.null(bgcol) & wall==0)   bgcol<-1
if(is.null(bgcol) & wall!=0)   bgcol<-'#DBE0E8'
if(is.null(col)  & wall==0)   col<-'gray25'
if(is.null(col)  & wall!=0)   col<-rgb(0.3,0.3,0.45)



################################################################################
################################################################################
k=0
for (i in 1:N1) {
for (j in 1:N2) {

x<-eval(parse(text=paste("data$",vars1[i], sep='')))
y<-eval(parse(text=paste("data$",vars2[j], sep='')))


######################
##################
if(partial){
wert<-a.pcor.test(x,y,p,use="mat",method=method,na.rm=TRUE)
          
M[i, j]<- format_n.ade(as.numeric(wert[1]),digits=digits)
P[i, j]<- format_p.ade(as.numeric(wert[2]),pgits=pdigs)
pvalue <-as.numeric(wert[2])
MM[i, j]<-a.hexrange.ade(M[i, j], '#FF0000','#999999','#0000FF', rangemy=c(-0.5,0.5), abs=FALSE)
if(pvalue> 0.05)  PP[i, j]<-'#bb2222'
if(pvalue<=0.05)  PP[i, j]<-'#000000'

}
##################
######################

######################
##################
if(!partial){
wert<-cor.test(x, y, method=method,na.rm=TRUE)

M[i, j]<- format_n.ade(as.numeric(wert[4]),digits=digits)
P[i, j]<- format_p.ade(as.numeric(wert[3]),pgits=pdigs)
pvalue <-as.numeric(wert[3])
MM[i, j]<-a.hexrange.ade(M[i, j], '#FF0000','#999999','#0000FF', rangemy=c(-0.4,0.4), abs=FALSE)
if(pvalue> 0.05)  PP[i, j]<-'#bb2222'
if(pvalue<=0.05)  PP[i, j]<-'#000000'

}
##################
######################

}
}
###########################


if(pvals) {

for (i in 1:N2) {
W[ ,(i*2)-1]<- M[ , i]
W[ ,(i*2)  ]<- P[ , i]
WW[ ,(i*2)-1]<- MM[ , i]
WW[ ,(i*2)  ]<- PP[ , i]


colnames(W)[(i*2)-1] <- vnames2[i]
colnames(W)[(i*2)]   <- 'p-value'

}


}

if(!pvals) {
W<-M
WW<-MM

}
################################################################################
################################################################################

############################################################################################################################

################################################################################
################################################################################
# Plot Bereich
if(wall==0){
xray<-dim(M)[2]
yray<-dim(M)[1]
newmai<-rep(0, 4)
oldmai<-par('mai')
if(oldmai[1]>0.6 & oldmai[1]<=1.02) newmai[1]<- 0.6-oldmai[1]
if(oldmai[2]<1.95) newmai[2]<- 1.95 - oldmai[2]
if(oldmai[3]<1.85) newmai[3]<- 1.85 - oldmai[3]
if(oldmai[4]<0.6) newmai[4]<- 0.6 - oldmai[4]
par(mai=(oldmai+newmai))
par(col.axis=tcol)


plot(0,0, xlim=c(0+0.04*xray,xray-0.04*xray), ylim=c(0+0.04*yray,yray-0.04*yray), axes=FALSE, xlab='', ylab='', main='', col=rgb(1,1,1,0))
axis(2, labels =vnames1, tick=FALSE, col=tcol, lwd=0, at=(yray-0.5):0.5, las =TRUE, line=-0.5, font=2)
axis(3, labels=FALSE,  tick=FALSE, col=tcol, at=(xray-0.5):0.5, line=-0.5)
text(x=0.5:(xray-0.5) , y=rep(yray+diff(par('usr')[3:4]/30), xray),labels=vnames2,, col=tcol,  srt =45, adj=0, font=2, xpd=TRUE)

####################
# Labels
if(!is.null(x2lab)) text(x=a.glc(0), y=a.glc(3, 7.5), labels=x2lab,  cex=par('cex.main'), col=tcol, adj=c(0.5, 0.5), xpd=TRUE, font=par("font.main"))
if(!is.null(xlab))  text(x=a.glc(0), y=a.glc(1, 1.5), labels=xlab,  cex=par('cex.lab'),  col=tcol, adj=c(0.5, 0.5), xpd=TRUE, font=par("font.lab"))
if(!is.null(ylab))  text(x=a.glc(2, 8.5), y=a.glc(5), labels=ylab,  cex=par('cex.lab'),  col=tcol, adj=c(0.5, 0.5), xpd=TRUE, font=par("font.lab"), srt=90)
####################



MM[MM=='#NANAFF']<-1
if(xray<=4 & yray<=4)  par(cex=1.5)
if(xray>4  | yray> 4)  par(cex=1.2)
if(xray>5  | yray> 5)  par(cex=1.1)
if(xray>6  | yray> 6)  par(cex=0.9)
if(xray>7  | yray> 7)  par(cex=0.8)
if(xray>8  | yray> 8)  par(cex=0.7)
if(xray>9  | yray> 9)  par(cex=0.6)

withtext<-T
circles<-F

for(i in 1:yray){
if(circles) symbols(x=0.5:(xray-0.5) , y=rep(yray-i+0.5, xray), circles=(abs(as.numeric(M[i, ]))), add =TRUE,                                               bg = a.alpha.ade(MM[i, ], 0.5), fg=a.alpha.ade(MM[i, ], 0), inches =FALSE)
if(bars)    symbols(x=0.5:(xray-0.5) , y=rep(yray-i, xray)+(abs(as.numeric(M[i, ]))/2), rectangles=cbind(rep(1, xray),abs(as.numeric(M[i, ]))) , add =TRUE, bg = a.alpha.ade(MM[i, ], 0.5), fg=a.alpha.ade(MM[i, ], 0.75), inches =FALSE)

if(bars)   color <- a.coladd.ade(MM[i, ], -50 )
if(!bars)  color <- MM[i, ]

if(withtext & pvals)  text(x=0.5:(xray-0.5) , y=yray-i+0.625 ,labels=M[i, ], col=color)
if(withtext & !pvals) text(x=0.5:(xray-0.5) , y=yray-i+0.5   ,labels=M[i, ], col=color)
if(withtext & pvals)  text(x=0.5:(xray-0.5) , y=yray-i+0.37  ,labels=P[i, ], col=color)
}

abline(col=bgcol, v=0:xray, h=0:yray, lwd=1)
box(col=bgcol, lwd=1)
par(cex=1)
}
################################################################################


################################################################################
################################################################################
# Wall 1
if(wall==1){
par(col.axis=tcol)
par(col.lab=tcol)
par(col.main=tcol)

newmai<-rep(0, 4)
oldmai<-par('mai')
if(oldmai[1]>0.6 & oldmai[1]<=1.02) newmai[1]<- 0.6-oldmai[1]
if(oldmai[2]<1.95) newmai[2]<- 1.95 - oldmai[2]
if(oldmai[3]<1.85) newmai[3]<- 1.85 - oldmai[3]
if(oldmai[4]<0.6) newmai[4]<- 0.6 - oldmai[4]
par(mai=(oldmai+newmai))

xray<-dim(M)[2]
yray<-dim(M)[1]
par(col.axis=tcol)

plot(0,0, xlim=c(0+0.04*xray,xray-0.04*xray), ylim=c(0+0.04*yray,yray-0.04*yray), axes=FALSE, xlab='', ylab='', main='', col=rgb(1,1,1,0))
polygon( c(par('usr')[c(1,1,2,2)]), par('usr')[c(3,4,4,3)], col=bgcol, border=FALSE)
axis(2, labels =vnames1, tick=FALSE, col=tcol, lwd=0, at=(yray-0.5):0.5, las =TRUE, line=-0.5, font=2)
axis(3, labels=FALSE,  tick=FALSE, col=tcol, at=(xray-0.5):0.5, line=-0.5)

text(x=0.5:(xray-0.5) , y=rep(yray+diff(par('usr')[3:4]/30), xray),labels=vnames2,, col=tcol,  srt =45, adj=0, font=2, xpd=TRUE)

####################
# Labels
if(!is.null(x2lab)) text(x=a.glc(0), y=a.glc(3, 7.5), labels=x2lab,  cex=par('cex.main'), col=tcol, adj=c(0.5, 0.5), xpd=TRUE, font=par("font.main"))
if(!is.null(xlab))  text(x=a.glc(0), y=a.glc(1, 1.5), labels=xlab,  cex=par('cex.lab'),  col=tcol, adj=c(0.5, 0.5), xpd=TRUE, font=par("font.lab"))
if(!is.null(ylab))  text(x=a.glc(2, 8.5), y=a.glc(5), labels=ylab,  cex=par('cex.lab'),  col=tcol, adj=c(0.5, 0.5), xpd=TRUE, font=par("font.lab"), srt=90)
####################


MM[MM=='#NANAFF']<-1
if(xray<=4 & yray<=4)  par(cex=1.5)
if(xray>4  | yray> 4)  par(cex=1.2)
if(xray>5  | yray> 5)  par(cex=1.1)
if(xray>6  | yray> 6)  par(cex=0.9)
if(xray>7  | yray> 7)  par(cex=0.8)
if(xray>8  | yray> 8)  par(cex=0.7)
if(xray>9  | yray> 9)  par(cex=0.6)

withtext<-T
circles<-F


for(i in 1:yray){
if(circles) symbols(x=0.5:(xray-0.5) , y=rep(yray-i+0.5, xray), circles=(abs(as.numeric(M[i, ]))), add =TRUE, bg = a.alpha.ade(MM[i, ], 0.5), fg=a.alpha.ade(MM[i, ], 0), inches =FALSE)
if(bars)   symbols(x=0.5:(xray-0.5) , y=rep(yray-i, xray)+(abs(as.numeric(M[i, ]))/2), rectangles=cbind(rep(1, xray),abs(as.numeric(M[i, ]))) , add =TRUE, bg = a.alpha.ade(MM[i, ], 0.5), fg=a.alpha.ade(MM[i, ], 0.75), inches =FALSE)
if(bars)  color <- a.coladd.ade(MM[i, ], -50 )
if(!bars)  color <- MM[i, ]
if(withtext & pvals)  text(x=0.5:(xray-0.5) , y=yray-i+0.625 ,labels=as.numeric(M[i, ]), col=MM[i, ])
if(withtext & !pvals) text(x=0.5:(xray-0.5) , y=yray-i+0.5 ,labels=as.numeric(M[i, ]), col=MM[i, ])
if(withtext & pvals)  text(x=0.5:(xray-0.5) , y=yray-i+0.375 ,labels=P[i, ], col=MM[i, ])
}

abline(col=rgb(1,1,1), v=0:xray, h=0:yray, lwd=1)
par(cex=1)
box(col=rgb(1,1,1), lwd=2)
}
################################################################################



################################################################################
################################################################################
# Wall 2
if(wall==2){
par(col.axis=tcol)
par(col.lab=tcol)
par(col.main=tcol)

newmai<-rep(0, 4)
oldmai<-par('mai')
if(oldmai[1]>0.6 & oldmai[1]<=1.02) newmai[1]<- 0.6-oldmai[1]
if(oldmai[2]<1.95) newmai[2]<- 1.95 - oldmai[2]
if(oldmai[3]<1.85) newmai[3]<- 1.85 - oldmai[3]
if(oldmai[4]<0.6) newmai[4]<- 0.6 - oldmai[4]
par(mai=(oldmai+newmai))


xray<-dim(M)[2]
yray<-dim(M)[1]
par(col.axis=tcol)

plot(0,0, xlim=c(0+0.04*xray,xray-0.04*xray), ylim=c(0+0.04*yray,yray-0.04*yray), axes=FALSE, xlab='', ylab='', main='', col=rgb(1,1,1,0))
axis(2, labels =vnames1, tick=FALSE, col=tcol, lwd=0, at=(yray-0.5):0.5, las =TRUE, line=-0.5, font=2)
axis(3, labels=FALSE,  tick=FALSE, col=tcol, at=(xray-0.5):0.5, line=-0.5)

text(x=0.5:(xray-0.5) , y=rep(yray+diff(par('usr')[3:4]/30), xray),labels=vnames2,, col=tcol,  srt =45, adj=0, font=2, xpd=TRUE)
####################
# Labels
if(!is.null(x2lab)) text(x=a.glc(0), y=a.glc(3, 7.5), labels=x2lab,  cex=par('cex.main'), col=tcol, adj=c(0.5, 0.5), xpd=TRUE, font=par("font.main"))
if(!is.null(xlab))  text(x=a.glc(0), y=a.glc(1, 1.5), labels=xlab,  cex=par('cex.lab'),  col=tcol, adj=c(0.5, 0.5), xpd=TRUE, font=par("font.lab"))
if(!is.null(ylab))  text(x=a.glc(2, 8.5), y=a.glc(5), labels=ylab,  cex=par('cex.lab'),  col=tcol, adj=c(0.5, 0.5), xpd=TRUE, font=par("font.lab"), srt=90)
####################


MM[MM=='#NANAFF']<-1
if(xray<=4 & yray<=4)  par(cex=1.5)
if(xray>4  | yray> 4)  par(cex=1.2)
if(xray>5  | yray> 5)  par(cex=1.1)
if(xray>6  | yray> 6)  par(cex=0.9)
if(xray>7  | yray> 7)  par(cex=0.8)
if(xray>8  | yray> 8)  par(cex=0.7)
if(xray>9  | yray> 9)  par(cex=0.6)

withtext<-T
circles<-F


for(i in 1:yray){
if(circles) symbols(x=0.5:(xray-0.5) , y=rep(yray-i+0.5, xray), circles=(abs(as.numeric(M[i, ]))), add =TRUE, bg = a.alpha.ade(MM[i, ], 0.5), fg=a.alpha.ade(MM[i, ], 0), inches =FALSE)
if(bars)   symbols(x=0.5:(xray-0.5) , y=rep(yray-i, xray)+(abs(as.numeric(M[i, ]))/2), rectangles=cbind(rep(1, xray),abs(as.numeric(M[i, ]))) , add =TRUE, bg = a.alpha.ade(MM[i, ], 0.5), fg=a.alpha.ade(MM[i, ], 0.75), inches =FALSE)
if(bars)  color <- a.coladd.ade(MM[i, ], -50 )
if(!bars)  color <- MM[i, ]
if(withtext & pvals) text(x=0.5:(xray-0.5) , y=yray-i+0.625 ,labels=as.numeric(M[i, ]), col=MM[i, ])
if(withtext & !pvals)text(x=0.5:(xray-0.5) , y=yray-i+0.5 ,labels=as.numeric(M[i, ]), col=MM[i, ])
if(withtext & pvals) text(x=0.5:(xray-0.5) , y=yray-i+0.375 ,labels=P[i, ], col=MM[i, ])
}

abline(col=bgcol, v=0:xray, h=0:yray, lwd=1)
par(cex=1)
box(col=a.coladd.ade(bgcol, -75), lwd=1)
}
################################################################################


################################################################################
################################################################################
# Wall 3
if(wall==3){
par(col.axis=tcol)
par(col.lab=tcol)
par(col.main=tcol)
newmai<-rep(0, 4)
oldmai<-par('mai')
if(oldmai[1]>0.6 & oldmai[1]<=1.02) newmai[1]<- 0.6-oldmai[1]
if(oldmai[2]<1.95) newmai[2]<- 1.95 - oldmai[2]
if(oldmai[3]<1.85) newmai[3]<- 1.85 - oldmai[3]
if(oldmai[4]<0.6) newmai[4]<- 0.6 - oldmai[4]
par(mai=(oldmai+newmai))

xray<-dim(M)[2]
yray<-dim(M)[1]
par(col.axis=tcol)
plot(0,0, xlim=c(0+0.04*xray,xray-0.04*xray), ylim=c(0+0.04*yray,yray-0.04*yray), axes=FALSE, xlab='', ylab='', main='', col=rgb(1,1,1,0))
polygon( c(par('usr')[c(1,1,2,2)]), par('usr')[c(3,4,4,3)], col=bgcol, border=FALSE)
axis(2, labels =vnames1, tick=FALSE, col=tcol, lwd=0, at=(yray-0.5):0.5, las =TRUE, line=-0.5, font=2)
axis(3, labels=FALSE,  tick=FALSE, col=tcol, at=(xray-0.5):0.5, line=-0.5)

text(x=0.5:(xray-0.5) , y=rep(yray+diff(par('usr')[3:4]/30), xray),labels=vnames2,, col=tcol,  srt =45, adj=0, font=2, xpd=TRUE)

####################
# Labels
if(!is.null(x2lab)) text(x=a.glc(0), y=a.glc(3, 7.5), labels=x2lab,  cex=par('cex.main'), col=tcol, adj=c(0.5, 0.5), xpd=TRUE, font=par("font.main"))
if(!is.null(xlab))  text(x=a.glc(0), y=a.glc(1, 1.5), labels=xlab,  cex=par('cex.lab'),  col=tcol, adj=c(0.5, 0.5), xpd=TRUE, font=par("font.lab"))
if(!is.null(ylab))  text(x=a.glc(2, 8.5), y=a.glc(5), labels=ylab,  cex=par('cex.lab'),  col=tcol, adj=c(0.5, 0.5), xpd=TRUE, font=par("font.lab"), srt=90)
####################


MM[MM=='#NANAFF']<-1
if(xray<=4 & yray<=4)  par(cex=1.5)
if(xray>4  | yray> 4)  par(cex=1.2)
if(xray>5  | yray> 5)  par(cex=1.1)
if(xray>6  | yray> 6)  par(cex=0.9)
if(xray>7  | yray> 7)  par(cex=0.8)
if(xray>8  | yray> 8)  par(cex=0.7)
if(xray>9  | yray> 9)  par(cex=0.6)

withtext<-T
circles<-F

for(i in 1:yray){
if(circles) symbols(x=0.5:(xray-0.5) , y=rep(yray-i+0.5, xray), circles=(abs(as.numeric(M[i, ]))), add =TRUE, bg = a.alpha.ade(MM[i, ], 0.5), fg=a.alpha.ade(MM[i, ], 0), inches =FALSE)
if(bars)   symbols(x=0.5:(xray-0.5) , y=rep(yray-i, xray)+(abs(as.numeric(M[i, ]))/2), rectangles=cbind(rep(1, xray),abs(as.numeric(M[i, ]))) , add =TRUE, bg = a.alpha.ade(MM[i, ], 0.5), fg=a.alpha.ade(MM[i, ], 0.75), inches =FALSE)
if(bars)  color <- a.coladd.ade(MM[i, ], -50 )
if(!bars)  color <- MM[i, ]
if(withtext & pvals) text(x=0.5:(xray-0.5) , y=yray-i+0.625 ,labels=as.numeric(M[i, ]), col=MM[i, ])
if(withtext & !pvals)text(x=0.5:(xray-0.5) , y=yray-i+0.5 ,labels=as.numeric(M[i, ]), col=MM[i, ])
if(withtext & pvals) text(x=0.5:(xray-0.5) , y=yray-i+0.375 ,labels=P[i, ], col=MM[i, ])
}

abline(col=a.coladd.ade(bgcol, -50), v=0:xray, h=0:yray, lwd=1)
par(cex=1)
box(col=a.coladd.ade(bgcol, -75), lwd=1)
}
################################################################################


################################################################################
################################################################################
# Wall 4
if(wall==4){
par(col.axis=tcol)
par(col.lab=rgb(1,1,1))
par(col.main=rgb(1,1,1))

newmai<-rep(0, 4)
oldmai<-par('mai')
if(oldmai[1]>0.6 & oldmai[1]<=1.02) newmai[1]<- 0.6-oldmai[1]
if(oldmai[2]<1.95) newmai[2]<- 1.95 - oldmai[2]
if(oldmai[3]<1.85) newmai[3]<- 1.85 - oldmai[3]
if(oldmai[4]<0.6) newmai[4]<- 0.6 - oldmai[4]
par(mai=(oldmai+newmai))


xray<-dim(M)[2]
yray<-dim(M)[1]
par(col.axis=tcol)
plot(0,0, xlim=c(0+0.04*xray,xray-0.04*xray), ylim=c(0+0.04*yray,yray-0.04*yray), axes=FALSE, xlab='', ylab='', main='', col=rgb(1,1,1,0))
polygon( c(par('usr')[c(1,1,2,2)]), par('usr')[c(3,4,4,3)], col=bgcol, border=FALSE)
par(xpd=TRUE)
dx<-7/par('din')[1]
dy<-7/par('din')[2]
xr<-(diff(par('usr')[1:2])/10)*dx
yr<-(diff(par('usr')[3:4])/10)*dy



polygon(a.glc(side=c(2,2,4,4), line=c(0,0,2,2)), a.glc(side=1, line=c(0,2,2,0)), col=bgcol, border=rgb(1,1,1))
polygon(a.glc(side=2, line=c(6.5,6.5,0,0)),     a.glc(side=c(1,3,3,1), line=c(2, 0, 0, 2)),  col=bgcol, border=rgb(1,1,1))
polygon(a.glc(side=2, line=c(6.5,6.5,8.5,8.5)), a.glc(side=c(1,3,3,1), line=c(2, 0, 0, 2)), col=tcol, border=rgb(1,1,1))
polygon(a.glc(side=c(2,2,4,4), line=c(0,0,2,2)), a.glc(side=3, line=c(0, 6.25, 6.25, 0)), col=bgcol, border=rgb(1,1,1))
polygon(a.glc(side=4, line=c(0,0,2,2)), a.glc(side=c(1,3,3,1), line=0), col=bgcol, border=rgb(1,1,1))
polygon(a.glc(side=2, line=c(0,0,8.5,8.5)), a.glc(side=3, line=c(0, 8.25, 8.25, 0)), col=tcol, border=rgb(1,1,1))
polygon(a.glc(side=c(2,2, 4,4), line=c(0,0,2,2)), a.glc(side=3, line=c(6.25, 8.25, 8.25, 6.25)), col=tcol, border=rgb(1,1,1))
polygon( a.glc(side=2, line=c(0,0,6.5,6.5)), a.glc(side=3, line=c(0, 6.25, 6.25, 0)), col=bgcol, border=rgb(1,1,1))

axis(2, labels =vnames1, tick=FALSE, col=tcol, lwd=0, at=(yray-0.5):0.5, las =TRUE, line=-0.5, font=2)
axis(3, labels=FALSE,  tick=FALSE, col=tcol, at=(xray-0.5):0.5, line=-0.5)
text(x=0.5:(xray-0.5) , y=rep(yray+diff(par('usr')[3:4]/24), xray),labels=vnames2, col=tcol,  srt =45, adj=0, font=2)


####################
# Labels
if(!is.null(x2lab)) text(x=a.glc(0), y=a.glc(3, 7.25), labels=x2lab, cex=par('cex.main'), col=rgb(1,1,1), adj=c(0.5, 0.5), xpd=TRUE, font=2)
if(!is.null(xlab))  text(x=a.glc(0), y=a.glc(1, 1), labels=xlab,  cex=par('cex.lab'),  col=tcol,       adj=c(0.5, 0.5),  xpd=TRUE, font=2)
if(!is.null(ylab))  text(x=a.glc(2, 7.5), y=a.glc(5), labels=ylab,  cex=par('cex.lab'),  col=rgb(1,1,1), adj=c(0.5, 0.5),  xpd=TRUE, font=2, srt=90)
####################

par(xpd=FALSE)

MM[MM=='#NANAFF']<-1
if(xray<=4 & yray<=4)  par(cex=1.5)
if(xray>4  | yray> 4)  par(cex=1.2)
if(xray>5  | yray> 5)  par(cex=1.1)
if(xray>6  | yray> 6)  par(cex=0.9)
if(xray>7  | yray> 7)  par(cex=0.8)
if(xray>8  | yray> 8)  par(cex=0.7)
if(xray>9  | yray> 9)  par(cex=0.6)

withtext<-T
circles<-F


for(i in 1:yray){
if(circles) symbols(x=0.5:(xray-0.5) , y=rep(yray-i+0.5, xray), circles=(abs(as.numeric(M[i, ]))), add =TRUE, bg = a.alpha.ade(MM[i, ], 0.5), fg=a.alpha.ade(MM[i, ], 0), inches =FALSE)
if(bars)   symbols(x=0.5:(xray-0.5) , y=rep(yray-i, xray)+(abs(as.numeric(M[i, ]))/2), rectangles=cbind(rep(1, xray),abs(as.numeric(M[i, ]))) , add =TRUE, bg = a.alpha.ade(MM[i, ], 0.5), fg=a.alpha.ade(MM[i, ], 0.75), inches =FALSE)
if(bars)  color <- a.coladd.ade(MM[i, ], -50 )
if(!bars)  color <- MM[i, ]
if(withtext & pvals) text(x=0.5:(xray-0.5) , y=yray-i+0.625 ,labels=as.numeric(M[i, ]), col=MM[i, ])
if(withtext & !pvals) text(x=0.5:(xray-0.5) , y=yray-i+0.5 ,labels=as.numeric(M[i, ]), col=MM[i, ])
if(withtext & pvals) text(x=0.5:(xray-0.5) , y=yray-i+0.375 ,labels=P[i, ], col=MM[i, ])
}

abline(col=rgb(1,1,1), v=0:xray, h=0:yray, lwd=1)
par(cex=1)
box(col=rgb(1,1,1), lwd=1)
}
################################################################################



################################################################################
################################################################################
# Wall 5
if(wall==5){
par(col.axis=tcol)
par(col.lab=tcol)
par(col.main=tcol)

xray<-dim(M)[2]
yray<-dim(M)[1]

newmai<-rep(0, 4)
oldmai<-par('mai')
if(oldmai[1]>0.75 & oldmai[1]<=1.02) newmai[1]<- 0.75-oldmai[1]
if(oldmai[2]<2.05) newmai[2]<- 2.05 - oldmai[2]
if(oldmai[3]<1.95) newmai[3]<- 1.95 - oldmai[3]
if(oldmai[4]<0.65) newmai[4]<- 0.65 - oldmai[4]
par(mai=(oldmai+newmai))



plot(0,0, xlim=c(0+0.04*xray,xray-0.04*xray), ylim=c(0+0.04*yray,yray-0.04*yray), axes=FALSE, xlab='', ylab='', main='', col=rgb(1,1,1,0))
par(xpd=TRUE)
dx<-7/par('din')[1]
dy<-7/par('din')[2]
xr<-(diff(par('usr')[1:2])/10)*dx
yr<-(diff(par('usr')[3:4])/10)*dy

polygon(a.glc(side=2, line=c(9, 9, 0, 0)), a.glc(side=3, line=c(0, 8.5, 8.5, 0)), col=rgb(1,1,1,0), border=tcol)
polygon(a.glc(side=2, line=c(9.6, 9.6 ,9, 9)),  a.glc(side=c(1,3,3,1), line=c(0.6, 8.5, 8.5, 0.6)), col=bgcol,  border=tcol)
polygon(a.glc(side=2, line=c(9.6, 9.6, 0, 0)), a.glc(side=1, line=c(0.6, 3, 3, 0.6)), col=bgcol, border=tcol)
polygon(a.glc(side=c(2,2,4,4), line=c(9.6,9.6,2.6,2.6)), a.glc(side=3, line=c(8.5, 9.1, 9.1, 8.5)), col=bgcol, border=tcol)
polygon(a.glc(side=4, line=c(2, 2 ,2.6, 2.6)), a.glc(side=c(1, 3, 3, 1), line=c(0, 8.5, 8.5, 0)), col=bgcol, border=tcol)
polygon(a.glc(side=4, line=c(0, 0 ,2, 2)), a.glc(side=c(1, 3, 3, 1), line=0), col=rgb(1,1,1,0), border=tcol)
polygon(a.glc(side=c(2, 2, 4, 4), line=c(0, 0, 2, 2)), a.glc(side=1, line=c(0.6, 3, 3, 0.6)), col=rgb(1,1,1,0), border=tcol)
polygon(a.glc(side=4, line=c(2, 2, 2.6, 2.6)), a.glc(side=1, line=c(0.6, 3, 3, 0.6)), col=bgcol, border=tcol)
text(a.glc(side=0), a.glc(side=1, line=2),  labels=main, cex = 1.25, font=2,  col=tcol, adj=c(0.5,0))
text(a.glc(side=0), a.glc(side=3, line=7), labels=x2lab, cex = 1.1,  font=2,  col=tcol, adj=c(0.5,0))
text(a.glc(side=2, line=8), a.glc(side=5), labels=ylab,  cex = 1.1,   font=2,  col=tcol, adj=c(0.5,0), srt=90)



axis(2, labels =vnames1, tick=FALSE, col=tcol, lwd=0, at=(yray-0.5):0.5, las =TRUE, line=-0.5, font=2)
axis(3, labels=FALSE,  tick=FALSE, col=tcol, at=(xray-0.5):0.5, line=-0.5)
text(x=0.5:(xray-0.5) , y=rep(yray+diff(par('usr')[3:4]/24), xray),labels=vnames2, col=tcol,  srt =45, adj=0, font=2)

par(xpd=FALSE)


MM[MM=='#NANAFF']<-1
if(xray<=4 & yray<=4)  par(cex=1.5)
if(xray>4  | yray> 4)  par(cex=1.2)
if(xray>5  | yray> 5)  par(cex=1.1)
if(xray>6  | yray> 6)  par(cex=0.9)
if(xray>7  | yray> 7)  par(cex=0.8)
if(xray>8  | yray> 8)  par(cex=0.7)
if(xray>9  | yray> 9)  par(cex=0.6)

withtext<-T
circles<-F


for(i in 1:yray){
if(circles) symbols(x=0.5:(xray-0.5) , y=rep(yray-i+0.5, xray), circles=(abs(as.numeric(M[i, ]))), add =TRUE, bg = a.alpha.ade(MM[i, ], 0.5), fg=a.alpha.ade(MM[i, ], 0), inches =FALSE)
if(bars)   symbols(x=0.5:(xray-0.5) , y=rep(yray-i, xray)+(abs(as.numeric(M[i, ]))/2), rectangles=cbind(rep(1, xray),abs(as.numeric(M[i, ]))) , add =TRUE, bg = a.alpha.ade(MM[i, ], 0.5), fg=a.alpha.ade(MM[i, ], 0.75), inches =FALSE)
if(bars)  color <- a.coladd.ade(MM[i, ], -50 )
if(!bars)  color <- MM[i, ]
if(withtext & pvals) text(x=0.5:(xray-0.5) , y=yray-i+0.625 ,labels=as.numeric(M[i, ]), col=MM[i, ])
if(withtext & !pvals) text(x=0.5:(xray-0.5) , y=yray-i+0.5 ,labels=as.numeric(M[i, ]), col=MM[i, ])
if(withtext & pvals) text(x=0.5:(xray-0.5) , y=yray-i+0.375 ,labels=P[i, ], col=MM[i, ])
}

abline(col=tcol, v=0:xray, h=0:yray, lwd=1)
par(cex=1)
box(col=tcol, lwd=1)
}
################################################################################


################################################################################
################################################################################
# Wall 6
if(wall==6){
par(col.axis=tcol)
par(col.lab=tcol)
par(col.main=tcol)

newmai<-rep(0, 4)
oldmai<-par('mai')
if(oldmai[1]>0.6 & oldmai[1]<=1.02) newmai[1]<- 0.6-oldmai[1]
if(oldmai[2]<1.95) newmai[2]<- 1.95 - oldmai[2]
if(oldmai[3]<1.85) newmai[3]<- 1.85 - oldmai[3]
if(oldmai[4]<0.6) newmai[4]<- 0.6 - oldmai[4]
par(mai=(oldmai+newmai))


xray<-dim(M)[2]
yray<-dim(M)[1]
par(col.axis=tcol)
plot(0,0, xlim=c(0+0.04*xray,xray-0.04*xray), ylim=c(0+0.04*yray,yray-0.04*yray), axes=FALSE, xlab='', ylab='', main='', col=rgb(1,1,1,0))
polygon( c(par('usr')[c(1,1,2,2)]), par('usr')[c(3,4,4,3)], col=bgcol, border=FALSE)
axis(2, labels =vnames1, tick=FALSE, col=tcol, lwd=0, at=(yray-0.5):0.5, las =TRUE, line=-0.5, font=2)
axis(3, labels=FALSE,  tick=FALSE, col=tcol, at=(xray-0.5):0.5, line=-0.5)

text(x=0.5:(xray-0.5) , y=rep(yray+diff(par('usr')[3:4]/30), xray),labels=vnames2,, col=tcol,  srt =45, adj=0, font=2, xpd=TRUE)

####################
# Labels
if(!is.null(x2lab)) text(x=a.glc(0), y=a.glc(3, 7.5), labels=x2lab,  cex=par('cex.main'), col=tcol, adj=c(0.5, 0.5), xpd=TRUE, font=par("font.main"))
if(!is.null(xlab))  text(x=a.glc(0), y=a.glc(1, 1.5), labels=xlab,  cex=par('cex.lab'),  col=tcol, adj=c(0.5, 0.5), xpd=TRUE, font=par("font.lab"))
if(!is.null(ylab))  text(x=a.glc(2, 8.5), y=a.glc(5), labels=ylab,  cex=par('cex.lab'),  col=tcol, adj=c(0.5, 0.5), xpd=TRUE, font=par("font.lab"), srt=90)
####################


MM[MM=='#NANAFF']<-1
if(xray<=4 & yray<=4)  par(cex=1.5)
if(xray>4  | yray> 4)  par(cex=1.2)
if(xray>5  | yray> 5)  par(cex=1.1)
if(xray>6  | yray> 6)  par(cex=0.9)
if(xray>7  | yray> 7)  par(cex=0.8)
if(xray>8  | yray> 8)  par(cex=0.7)
if(xray>9  | yray> 9)  par(cex=0.6)

withtext<-T
circles<-F


for(i in 1:yray){
if(circles) symbols(x=0.5:(xray-0.5) , y=rep(yray-i+0.5, xray), circles=(abs(as.numeric(M[i, ]))), add =TRUE, bg = a.alpha.ade(MM[i, ], 0.5), fg=a.alpha.ade(MM[i, ], 0), inches =FALSE)
if(bars)   symbols(x=0.5:(xray-0.5) , y=rep(yray-i, xray)+(abs(as.numeric(M[i, ]))/2), rectangles=cbind(rep(1, xray),abs(as.numeric(M[i, ]))) , add =TRUE, bg = a.alpha.ade(MM[i, ], 0.5), fg=a.alpha.ade(MM[i, ], 0.75), inches =FALSE)
if(bars)  color <- a.coladd.ade(MM[i, ], -50 )
if(!bars)  color <- MM[i, ]
if(withtext & pvals)  text(x=0.5:(xray-0.5) , y=yray-i+0.625 ,labels=as.numeric(M[i, ]), col=MM[i, ])
if(withtext & !pvals) text(x=0.5:(xray-0.5) , y=yray-i+0.5 ,labels=as.numeric(M[i, ]), col=MM[i, ])
if(withtext & pvals)  text(x=0.5:(xray-0.5) , y=yray-i+0.375 ,labels=P[i, ], col=MM[i, ])
}

abline(col=a.coladd.ade(bgcol, -35), v=0:xray, h=0:yray, lwd=3)
abline(col=rgb(1,1,1), v=0:xray, h=0:yray, lwd=1)
par(cex=1)
box(lwd=3, col=rgb(1,1,1))
box(lwd=1, col=a.coladd.ade(bgcol, -35))
}
################################################################################

################################################################################
################################################################################

}
