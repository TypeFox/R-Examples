parallel.set.ade <-
function(vars,  vnames=NULL,  data=NULL, xlab='Factors', ylab='Proportion', main=NULL, col=NULL, tcol=NULL,  bgcol=NULL, lcol=NULL, alpha = NULL, cex=NULL, wall=0, horizontal=FALSE){
if(any(par('mfg')!=c(1,1,1,1)) & any(par('mai') < c(1.02, 0.82, 0.82, 0.42))){
maidiff<-rep(0, 4)
norm<-c(1.02, 0.82, 0.82, 0.42)
maidiff[par('mai')<norm]<-  norm[par('mai')<norm] - par('mai')[par('mai')<norm]
par(mai=par('mai')+maidiff)
}
oldpar<-par(no.readonly =TRUE)
oldpar<-oldpar[-which(names(oldpar)%in%c('usr', 'plt',   'pin',  'fin', 'fig', 'mfg', 'mfcol', 'mfrow', 'omd', 'omi', 'oma'))]
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
vars<-newvars
}
#####################################


if(length(vars)>6) stop('to many factors')
for(i in 1:length(vars)){
data<- subset(data, !is.na(eval(parse(text=paste(vars[i], sep='')))))
}

if (is.null(vnames)) vnames <- vars
if(is.null(tcol)  & wall==0)   tcol<-1
if(is.null(tcol)  & wall!=0)   tcol<-rgb(0.1,0.1,0.25)
if(is.null(bgcol) & wall==0)   bgcol<-1
if(is.null(bgcol) & wall!=0)   bgcol<-'#DBE0E8'
if(is.null(lcol)) lcol<-c(rgb(1,1,1), 'gray90')

if(is.null(alpha)){
alpha<-0.25
if(wall==0) alpha<-1
}


  
xrange<-c(1-0.025, length(vars)+0.025)
yrange<-c(0, 1)


######################################################
ord<-FALSE
if(length(vars)<=6){
for(i in 1:length(vars)){
ord<-TRUE
v<- eval(parse(text=paste("data$",vars[i],  sep='')))
if(is.null(cex)){
acex<-1
if(length(unique(v))>5)  acex<-0.9
if(length(unique(v))>7)  acex<-0.85
if(length(unique(v))>9)  acex<-0.75
if(length(unique(v))>13) acex<-0.6
if(length(unique(v))>16) acex<-0.5
}
if(!is.null(cex))  acex<-cex
}
}


M<-NULL
Fs<-NULL
for(i in 1:length(vars)){
v<- eval(parse(text=paste("data$",vars[i],  sep='')))
M<-cbind(M, v)
if(ord) Fs[[i]]<-levels(as.factor(v))
}
######################################################



runpoly<-function(M, alcol, mcol=tcol, plot=TRUE, horizontal=TRUE){
################################################################################
################################################################################
rbc2 <-alcol

######################################
a.prop<-function(v, from=0, to=1){
v<-v/sum(v, na.rm=TRUE)
rdiff<-diff(c(from, to))
cuts<-from
runin<-from
mitts<-NULL
for(i in 1:length(v)){
cuts<-c(cuts, (rdiff*v[i]+runin))
mitts<-c(mitts, (cuts[i]+cuts[(i+1)])/2)
runin<-runin+rdiff*v[i]
}
return(cuts)
}
######################################

# Allgemein Calcs
if(horizontal) {
xr2<-a.glc(side=4, line=0.5)-a.glc(side=4, line=0)
xr<- xr2
yr<- a.glc(side=3, line=0.575)-a.glc(side=3, line=0)
}
if(!horizontal) {
xr2<-a.glc(side=1, line=0.5)-a.glc(side=1, line=0)
xr<- xr2
yr<- a.glc(side=4, line=0.575)-a.glc(side=4, line=0)
}

Mdim<-dim(M)
x<-as.factor(M[ , 1])
y<-as.factor(M[ , 2])
if(Mdim[2]>2) z<-as.factor(M[ , 3])
if(Mdim[2]>3) w<-as.factor(M[ , 4])
if(Mdim[2]>4) v<-as.factor(M[ , 5])
if(Mdim[2]>5) u<-as.factor(M[ , 6])

color<-col
if(is.null(color)) color <- a.getcol.ade(nlevels(x))

x_only <- a.prop(table(x), from=0, to=1)
y_only <- a.prop(table(y), from=0, to=1)
if(Mdim[2]>2) z_only <- a.prop(table(z), from=0, to=1)
if(Mdim[2]>3) w_only <- a.prop(table(w), from=0, to=1)
if(Mdim[2]>4) v_only <- a.prop(table(v), from=0, to=1)
if(Mdim[2]>5) u_only <- a.prop(table(u), from=0, to=1)

################################################################################
# X  Block
for(i in 1:nlevels(x)){

# Subby 0
y_in_x <- a.prop(table(y[x==levels(x)[i]]), from=x_only[i], to=x_only[i+1])


################################################################################
# Y  Block
for(j in 1:nlevels(y)){

# Subby 1
x_in_y <- a.prop(table(x[y==levels(y)[j]]), from=y_only[j], to=y_only[j+1])

if(all(!is.nan(x_in_y)) &  all(!is.nan(y_in_x))){
xrun<-c(1+xr, 2-xr, 2-xr, 1+xr)
yrun<-c(y_in_x[j], x_in_y[i], x_in_y[i+1], y_in_x[j+1])
if(horizontal & !(y_in_x[j]==y_in_x[j+1]) & !(x_in_y[i]==x_in_y[i+1]))  polygon(xrun, yrun, col=a.alpha.ade(color[i], alpha), border=a.alpha.ade(color[i], 0.5))
if(!horizontal & !(y_in_x[j]==y_in_x[j+1]) & !(x_in_y[i]==x_in_y[i+1])) polygon(yrun, xrun, col=a.alpha.ade(color[i], alpha), border=a.alpha.ade(color[i], 0.5))
}

################################################################################
# Z block
if(Mdim[2]>2){
z_in_y <- a.prop(table(z[y==levels(y)[j] & x==levels(x)[i]]), from=x_in_y[i], to=x_in_y[i+1])
for(k in 1:nlevels(z)){

# Subby 2
x_in_z <- a.prop(table(x[z==levels(z)[k]]), from=z_only[k], to=z_only[k+1])
y_in_z <- a.prop(table(y[z==levels(z)[k] & x==levels(x)[i]]), from=x_in_z[i], to=x_in_z[i+1])

if(all(!is.nan(y_in_z)) &  all(!is.nan(z_in_y))){
xrun<-c(2+xr, 3-xr, 3-xr, 2+xr)
yrun<-c(z_in_y[k], y_in_z[j], y_in_z[j+1], z_in_y[k+1])
if(horizontal & !(z_in_y[k]==z_in_y[k+1]) & !(y_in_z[j]==y_in_z[j+1]))  polygon(xrun, yrun, col=a.alpha.ade(color[i], alpha), border=a.alpha.ade(color[i], 0.5))
if(!horizontal & !(z_in_y[k]==z_in_y[k+1]) & !(y_in_z[j]==y_in_z[j+1])) polygon(yrun, xrun, col=a.alpha.ade(color[i], alpha), border=a.alpha.ade(color[i], 0.5))

}

################################################################################
# W block
if(Mdim[2]>3){
w_in_z <- a.prop(table(w[y==levels(y)[j] & x==levels(x)[i] & z==levels(z)[k]]), from=y_in_z[j], to=y_in_z[j+1])
for(l in 1:nlevels(w)){

# Subby 3
x_in_w <- a.prop(table(x[w==levels(w)[l]]), from=w_only[l], to=w_only[l+1])
y_in_w <- a.prop(table(y[w==levels(w)[l] & x==levels(x)[i]]), from=x_in_w[i], to=x_in_w[i+1])
z_in_w <- a.prop(table(z[w==levels(w)[l] & y==levels(y)[j] & x==levels(x)[i]]), from=y_in_w[j], to=y_in_w[j+1])

if(all(!is.nan(z_in_w)) &  all(!is.nan(w_in_z))){
xrun<-c(3+xr, 4-xr, 4-xr, 3+xr)
yrun<-c(w_in_z[l], z_in_w[k], z_in_w[k+1], w_in_z[l+1])
if(horizontal & !(w_in_z[l]==w_in_z[l+1]) & !(z_in_w[k]==z_in_w[k+1])) polygon(xrun, yrun, col=a.alpha.ade(color[i], alpha), border=a.alpha.ade(color[i], 0.5))
if(!horizontal & !(w_in_z[l]==w_in_z[l+1]) & !(z_in_w[k]==z_in_w[k+1])) polygon(yrun, xrun, col=a.alpha.ade(color[i], alpha), border=a.alpha.ade(color[i], 0.5))
}

################################################################################
# V block
if(Mdim[2]>4){
v_in_w <- a.prop(table(v[y==levels(y)[j] & x==levels(x)[i] & z==levels(z)[k] & w==levels(w)[l]]), from=z_in_w[k], to=z_in_w[k+1])

for(m in 1:nlevels(v)){

# Subby 4
x_in_v <- a.prop(table(x[v==levels(v)[m]]), from=v_only[m], to=v_only[m+1])
y_in_v <- a.prop(table(y[v==levels(v)[m] & x==levels(x)[i]]), from=x_in_v[i], to=x_in_v[i+1])
z_in_v <- a.prop(table(z[v==levels(v)[m] & y==levels(y)[j] & x==levels(x)[i]]), from=y_in_v[j], to=y_in_v[j+1])
w_in_v <- a.prop(table(w[v==levels(v)[m] & y==levels(y)[j] & x==levels(x)[i] & z==levels(z)[k]]), from=z_in_v[k], to=z_in_v[k+1])


if(all(!is.nan(w_in_v)) &  all(!is.nan(v_in_w))){
xrun<-c(4+xr, 5-xr, 5-xr, 4+xr)
yrun<-c(v_in_w[m], w_in_v[l], w_in_v[l+1], v_in_w[m+1])
if(horizontal & !(v_in_w[m]==v_in_w[m+1]) & !(w_in_v[l]==w_in_v[l+1]))  polygon(xrun, yrun, col=a.alpha.ade(color[i], alpha), border=a.alpha.ade(color[i], 0.5))
if(!horizontal & !(v_in_w[m]==v_in_w[m+1]) & !(w_in_v[l]==w_in_v[l+1]))  polygon(yrun, xrun, col=a.alpha.ade(color[i], alpha), border=a.alpha.ade(color[i], 0.5))

}


################################################################################
# U block
if(Mdim[2]>5){
u_in_v <- a.prop(table(u[y==levels(y)[j] & x==levels(x)[i] & z==levels(z)[k] & w==levels(w)[l] & v==levels(v)[m]]), from=w_in_v[l], to=w_in_v[l+1])


for(n in 1:nlevels(u)){

# Subby 5
x_in_u <- a.prop(table(x[u==levels(u)[n]]), from=u_only[n], to=u_only[n+1])
y_in_u <- a.prop(table(y[u==levels(u)[n] & x==levels(x)[i]]), from=x_in_u[i], to=x_in_u[i+1])
z_in_u <- a.prop(table(z[u==levels(u)[n] & y==levels(y)[j] & x==levels(x)[i]]), from=y_in_u[j], to=y_in_u[j+1])
w_in_u <- a.prop(table(w[u==levels(u)[n] & y==levels(y)[j] & x==levels(x)[i] & z==levels(z)[k]]), from=z_in_u[k], to=z_in_u[k+1])
v_in_u <- a.prop(table(v[u==levels(u)[n] & y==levels(y)[j] & x==levels(x)[i] & z==levels(z)[k] & w==levels(w)[l]]), from=w_in_u[l], to=w_in_u[l+1])


if(all(!is.nan(v_in_u)) & all(!is.nan(u_in_v))){
xrun<-c(5+xr, 6-xr, 6-xr, 5+xr)
yrun<-c(u_in_v[n], v_in_u[m], v_in_u[m+1], u_in_v[n+1])
if(horizontal & !(u_in_v[n]==u_in_v[n+1]) & !(v_in_u[m]==v_in_u[m+1])) polygon(xrun, yrun, col=a.alpha.ade(color[i], alpha), border=a.alpha.ade(color[i], 0.5))
if(!horizontal & !(u_in_v[n]==u_in_v[n+1]) & !(v_in_u[m]==v_in_u[m+1])) polygon(yrun, xrun, col=a.alpha.ade(color[i], alpha), border=a.alpha.ade(color[i], 0.5))
}


################################################################################
}
}
}
}
}
}
}
}
}
}
################################################################################


# Dim Legens
if(horizontal){
rect(par('usr')[1], x_only[-length(x_only)],  1-xr2 , x_only[-1], col=a.alpha.ade(color, alpha), border=tcol, xpd=TRUE)

rect(1-xr2, x_only[-length(x_only)],  1+xr2 , x_only[-1], col=c(rgb(1,1,1), 'gray90'), border=color, xpd=TRUE)
rect(2-xr2, y_only[-length(y_only)],  2+xr2 , y_only[-1], col=c(rgb(1,1,1), 'gray90'), border=bgcol, xpd=TRUE)
if(Mdim[2]>2) rect(3-xr2, z_only[-length(z_only)],  3+xr2 , z_only[-1], col=lcol, border=bgcol, xpd=TRUE)
if(Mdim[2]>3) rect(4-xr2, w_only[-length(w_only)],  4+xr2 , w_only[-1], col=lcol, border=bgcol, xpd=TRUE)
if(Mdim[2]>4) rect(5-xr2, v_only[-length(v_only)],  5+xr2 , v_only[-1], col=lcol, border=bgcol, xpd=TRUE)
if(Mdim[2]>5) rect(6-xr2, u_only[-length(u_only)],  6+xr2 , u_only[-1], col=lcol, border=bgcol, xpd=TRUE)


text(x=rep(1, nlevels(x)), y=((x_only[1:length(x_only)-1]+x_only[2:length(x_only)])/2),               cex=acex, labels=Fs[[1]] ,srt =90, adj = c(0.5, 0.45), col=tcol)
text(x=rep(2, nlevels(y)), y=((y_only[1:length(y_only)-1]+y_only[2:length(y_only)])/2),               cex=acex, labels=Fs[[2]] ,srt =90, adj = c(0.5, 0.45), col=tcol)
if(Mdim[2]>2) text(x=rep(3, nlevels(z)), y=((z_only[1:length(z_only)-1]+z_only[2:length(z_only)])/2), cex=acex, labels=Fs[[3]] ,srt =90, adj = c(0.5, 0.45), col=tcol)
if(Mdim[2]>3) text(x=rep(4, nlevels(w)), y=((w_only[1:length(w_only)-1]+w_only[2:length(w_only)])/2), cex=acex, labels=Fs[[4]] ,srt =90, adj = c(0.5, 0.45), col=tcol)
if(Mdim[2]>4) text(x=rep(5, nlevels(v)), y=((v_only[1:length(v_only)-1]+v_only[2:length(v_only)])/2), cex=acex, labels=Fs[[5]] ,srt =90, adj = c(0.5, 0.45), col=tcol)
if(Mdim[2]>5) text(x=rep(6, nlevels(u)), y=((u_only[1:length(u_only)-1]+u_only[2:length(u_only)])/2), cex=acex, labels=Fs[[6]] ,srt =90, adj = c(0.5, 0.45), col=tcol)
}
###########


if(!horizontal){
rect(x_only[-length(x_only)],  par('usr')[4], x_only[-1], 1-xr2 , col=a.alpha.ade(color, alpha), border=tcol, xpd=TRUE)

              rect(x_only[-length(x_only)], 1-xr2, x_only[-1], 1+xr2 ,col=lcol, border=color, xpd=TRUE)
              rect(y_only[-length(y_only)], 2-xr2, y_only[-1], 2+xr2 ,col=lcol, border=bgcol, xpd=TRUE)
if(Mdim[2]>2) rect(z_only[-length(z_only)], 3-xr2, z_only[-1], 3+xr2 ,col=lcol, border=bgcol, xpd=TRUE)
if(Mdim[2]>3) rect(w_only[-length(w_only)], 4-xr2, w_only[-1], 4+xr2 ,col=lcol, border=bgcol, xpd=TRUE)
if(Mdim[2]>4) rect(v_only[-length(v_only)], 5-xr2, v_only[-1], 5+xr2 ,col=lcol, border=bgcol, xpd=TRUE)
if(Mdim[2]>5) rect(u_only[-length(u_only)], 6-xr2, u_only[-1], 6+xr2 ,col=lcol, border=bgcol, xpd=TRUE)


              text(y=rep(1, nlevels(x)), x=((x_only[1:length(x_only)-1]+x_only[2:length(x_only)])/2), cex=acex, labels=Fs[[1]] , adj = c(0.5, 0.45), col=tcol)
              text(y=rep(2, nlevels(y)), x=((y_only[1:length(y_only)-1]+y_only[2:length(y_only)])/2), cex=acex, labels=Fs[[2]] , adj = c(0.5, 0.45), col=tcol)
if(Mdim[2]>2) text(y=rep(3, nlevels(z)), x=((z_only[1:length(z_only)-1]+z_only[2:length(z_only)])/2), cex=acex, labels=Fs[[3]] , adj = c(0.5, 0.45), col=tcol)
if(Mdim[2]>3) text(y=rep(4, nlevels(w)), x=((w_only[1:length(w_only)-1]+w_only[2:length(w_only)])/2), cex=acex, labels=Fs[[4]] , adj = c(0.5, 0.45), col=tcol)
if(Mdim[2]>4) text(y=rep(5, nlevels(v)), x=((v_only[1:length(v_only)-1]+v_only[2:length(v_only)])/2), cex=acex, labels=Fs[[5]] , adj = c(0.5, 0.45), col=tcol)
if(Mdim[2]>5) text(y=rep(6, nlevels(u)), x=((u_only[1:length(u_only)-1]+u_only[2:length(u_only)])/2), cex=acex, labels=Fs[[6]] , adj = c(0.5, 0.45), col=tcol)
}
###########


return(list(rbc2, 0,  1))
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

if(!horizontal & (xlab!='Factors' |  ylab!='Proportion')){
save_xlab<-xlab
xlab<-ylab
ylab<-save_xlab
}

if(horizontal)  plot(0, 0 , xlim=xrange, ylim=yrange, axes=FALSE, xlab=xlab, ylab=ylab, main=main, col=rgb(1,1,1,0))
if(!horizontal) plot(0, 0 , xlim=yrange, ylim=c(xrange[2], xrange[1]),  axes=FALSE, xlab=ylab, ylab=xlab, main=main, col=rgb(1,1,1,0))
if(horizontal){
axis(2 , at=(pretty(c(yrange[1],yrange[2]), n = 10)) , labels=rep('', 11), col.ticks=bgcol, las=1)
a1<-axis(1,  at=1:length(vars), labels=vnames, col.ticks=bgcol)
}
if(!horizontal){
axis(1 , at=(pretty(c(yrange[1],yrange[2]), n = 10)) , labels=rep('', 11), col.ticks=bgcol, las=1)
a2<-axis(2,  at=1:length(vars), labels=vnames, col.ticks=bgcol)
}

runpoly(M, bgcol, horizontal=horizontal)

box(col=bgcol)
}
################################################################################



################################################################################
if(wall==1){
par(col.axis=tcol)
par(col.lab=tcol)
par(col.main=tcol)
#######################

if(!horizontal & (xlab!='Factors' |  ylab!='Proportion')){
save_xlab<-xlab
xlab<-ylab
ylab<-save_xlab
}

if(horizontal)  plot(0, 0 , xlim=xrange, ylim=yrange, axes=FALSE,                   xlab=xlab, ylab=ylab, main=main, col=rgb(1,1,1,0))
if(!horizontal) plot(0, 0 , xlim=yrange, ylim=c(xrange[2], xrange[1]),  axes=FALSE, xlab=ylab, ylab=xlab, main=main, col=rgb(1,1,1,0))
polygon( c(par('usr')[c(1,1,2,2)]), par('usr')[c(3,4,4,3)], col=bgcol, border=FALSE)
if(horizontal){
a1<-axis(2 , at=(pretty(c(yrange[1],yrange[2]), n = 10)) , labels=rep('', 11), col.ticks=bgcol, las=1)
axis(1,  at=1:length(vars), labels=vnames, col.ticks=bgcol)
abline(h=a1, lty=1, col=rgb(1,1,1), lwd=1)
}
if(!horizontal){
a1<-axis(1 , at=(pretty(c(yrange[1],yrange[2]), n = 10)) , labels=rep('', 11), col.ticks=bgcol, las=1)
axis(2,  at=1:length(vars), labels=vnames, col.ticks=bgcol)
abline(v=a1, lty=1, col=rgb(1,1,1), lwd=1)
}

runpoly(M, rgb(1,1,1), tcol, horizontal=horizontal)

box(col=bgcol, lwd=1)
}
################################################################################



################################################################################
if(wall==2){
par(col.axis=tcol)
par(col.lab=tcol)
par(col.main=tcol)


if(!horizontal & (xlab!='Factors' |  ylab!='Proportion')){
save_xlab<-xlab
xlab<-ylab
ylab<-save_xlab
}

if(horizontal)  plot(0, 0 , xlim=xrange, ylim=yrange, axes=FALSE,                   xlab=xlab, ylab=ylab, main=main, col=rgb(1,1,1,0))
if(!horizontal) plot(0, 0 , xlim=yrange, ylim=c(xrange[2], xrange[1]),  axes=FALSE, xlab=ylab, ylab=xlab, main=main, col=rgb(1,1,1,0))

if(horizontal){
a1<-axis(2 , at=(pretty(c(yrange[1],yrange[2]), n = 10)) , labels=rep('', 11), col.ticks=a.coladd.ade(bgcol, -75), las=1)
axis(1,  at=1:length(vars), labels=vnames, col.ticks=a.coladd.ade(bgcol, -75))
abline(h=a1, lty=1, col=bgcol, lwd=1)
}
if(!horizontal){
a1<-axis(1 , at=(pretty(c(yrange[1],yrange[2]), n = 10)) , labels=rep('', 11), col.ticks=a.coladd.ade(bgcol, -75), las=1)
axis(2,  at=1:length(vars), labels=vnames, col.ticks=a.coladd.ade(bgcol, -75))
abline(v=a1, lty=1, col=bgcol, lwd=1)
}

runpoly(M, rgb(1,1,1), tcol, horizontal=horizontal)


box(col=a.coladd.ade(bgcol, -75))
}
################################################################################



################################################################################
if(wall==3){
par(col.axis=tcol)
par(col.lab=tcol)
par(col.main=tcol)

if(!horizontal & (xlab!='Factors' |  ylab!='Proportion')){
save_xlab<-xlab
xlab<-ylab
ylab<-save_xlab
}

if(horizontal)  plot(0, 0 , xlim=xrange, ylim=yrange, axes=FALSE,                   xlab=xlab, ylab=ylab, main=main, col=rgb(1,1,1,0))
if(!horizontal) plot(0, 0 , xlim=yrange, ylim=c(xrange[2], xrange[1]),  axes=FALSE, xlab=ylab, ylab=xlab, main=main, col=rgb(1,1,1,0))
polygon( c(par('usr')[c(1,1,2,2)]), par('usr')[c(3,4,4,3)], col=bgcol, border=FALSE)

if(horizontal){
a1<-axis(2 , at=(pretty(c(yrange[1],yrange[2]), n = 10)) , labels=rep('', 11), col.ticks=bgcol, las=1)
axis(1,  at=1:length(vars), labels=vnames, col.ticks=bgcol)
abline(h=a1, lty=1, col=a.coladd.ade(bgcol, -50), lwd=1)
}
if(!horizontal){
a1<-axis(1 , at=(pretty(c(yrange[1],yrange[2]), n = 10)) , labels=rep('', 11), col.ticks=bgcol, las=1)
axis(2,  at=1:length(vars), labels=vnames, col.ticks=bgcol)
abline(v=a1, lty=1, col=a.coladd.ade(bgcol, -50), lwd=1)
}

runpoly(M, rgb(1,1,1), tcol, horizontal=horizontal)


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

if(!horizontal & (xlab=='Factors' |  ylab=='Proportion')){
save_xlab<-xlab
xlab<-ylab
ylab<-save_xlab
}

if(horizontal)  plot(0, 0 , xlim=xrange, ylim=yrange, axes=FALSE,                   xlab='', ylab='', main='', col=rgb(1,1,1,0))
if(!horizontal) plot(0, 0 , xlim=yrange, ylim=c(xrange[2], xrange[1]),  axes=FALSE, xlab='', ylab='', main='', col=rgb(1,1,1,0))
polygon( c(par('usr')[c(1,1,2,2)]), par('usr')[c(3,4,4,3)], col=bgcol, border=FALSE)


if(horizontal){
a1<-axis(2 , at=(pretty(c(yrange[1],yrange[2]), n = 10)) , labels=rep('', 11), col.ticks=bgcol, las=1)
axis(1,  at=1:length(vars), labels=vnames, col.ticks=bgcol)
abline(h=a1, lty=1, col=rgb(1,1,1), lwd=1)
}
if(!horizontal){
a1<-axis(1 , at=(pretty(c(yrange[1],yrange[2]), n = 10)) , labels=rep('', 11), col.ticks=bgcol, las=1)
axis(2,  at=1:length(vars), labels=vnames, col.ticks=bgcol)
abline(v=a1, lty=1, col=rgb(1,1,1), lwd=1)
}

runpoly(M, rgb(1,1,1), tcol, horizontal=horizontal)

# Outer
par(xpd=TRUE)
polygon(a.glc(side=c(2,2,4,4), line=c(0,0,0,0)), a.glc(side=3, line=c(0, 2.75,  2.75, 0)), col=tcol, border=rgb(1,1,1))
if(ylab!='' & ylab!=' ') polygon( a.glc(side=2, line=c(3.5, 3.5, 2, 2)), a.glc(side=c(1, 3, 3, 1), line=0), col=bgcol, border=rgb(1,1,1))
if(xlab!='' & xlab!=' ') polygon( a.glc(side=c(2, 2, 4, 4), line=0),     a.glc(side=1, line=c(4, 2.5, 2.5, 4)), col=bgcol, border=rgb(1,1,1))
text(a.glc(side=0), a.glc(side=3, line=1),    labels=main, cex = 1.25, font=2, col=rgb(1,1,1), adj=c(0.5,0))
text(a.glc(side=0), a.glc(side=1, line=3.5), labels=xlab, cex = 1.1,  font=2, col=tcol, adj=c(0.5,0))
text(a.glc(side=2, line=2.5), a.glc(side=5),  labels=ylab, cex = 1.1,  font=2,  col=tcol, adj=c(0.5,0), srt=90)
box(col=rgb(1,1,1))
par(xpd=FALSE)

runpoly(M, rgb(1,1,1), tcol, horizontal=horizontal)

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

if(!horizontal & (xlab=='Factors' |  ylab=='Proportion')){
save_xlab<-xlab
xlab<-ylab
ylab<-save_xlab
}

if(horizontal)  plot(0, 0 , xlim=xrange, ylim=yrange, axes=FALSE,                   xlab='', ylab='', main='', col=rgb(1,1,1,0))
if(!horizontal) plot(0, 0 , xlim=yrange, ylim=c(xrange[2], xrange[1]),  axes=FALSE, xlab='', ylab='', main='', col=rgb(1,1,1,0))

if(horizontal){
a1<-axis(2 , at=(pretty(c(yrange[1],yrange[2]), n = 10)) , labels=rep('', 11), col.ticks=bgcol, las=1)
axis(1,  at=1:length(vars), labels=vnames, col.ticks=bgcol)
}
if(!horizontal){
a1<-axis(1 , at=(pretty(c(yrange[1],yrange[2]), n = 10)) , labels=rep('', 11), col.ticks=bgcol, las=1)
axis(2,  at=1:length(vars), labels=vnames, col.ticks=bgcol)
}


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
runpoly(M, rgb(1,1,1), tcol, horizontal=horizontal)

box(col=tcol)
}
################################################################################


################################################################################
if(wall==6){
par(col.axis=tcol)
par(col.lab=tcol)
par(col.main=tcol)
#######################

if(!horizontal & (xlab!='Factors' |  ylab!='Proportion')){
save_xlab<-xlab
xlab<-ylab
ylab<-save_xlab
}

if(horizontal)  plot(0, 0 , xlim=xrange, ylim=yrange, axes=FALSE,                   xlab=xlab, ylab=ylab, main=main, col=rgb(1,1,1,0))
if(!horizontal) plot(0, 0 , xlim=yrange, ylim=c(xrange[2], xrange[1]),  axes=FALSE, xlab=xlab, ylab=ylab, main=main, col=rgb(1,1,1,0))
polygon( c(par('usr')[c(1,1,2,2)]), par('usr')[c(3,4,4,3)], col=bgcol, border=NA)

if(horizontal){
a1<-axis(2 , at=(pretty(c(yrange[1],yrange[2]), n = 10)) , labels=rep('', 11), col.ticks=a.coladd.ade(bgcol, -35), lwd.ticks=3)
a1<-axis(2 , at=(pretty(c(yrange[1],yrange[2]), n = 10)) , labels=rep('', 11), col.ticks=rgb(1,1,1), lwd.ticks=1)
axis(1,  at=1:length(vars), labels=vnames, col.ticks=a.coladd.ade(bgcol, -35), lwd.ticks=3)
axis(1,  at=1:length(vars), labels=vnames, col.ticks=rgb(1,1,1), lwd.ticks=1)
abline(h=a1, lty=1, col=a.coladd.ade(bgcol, -35), lwd=3)
abline(h=a1, lty=1, col=rgb(1,1,1), lwd=1)
}
if(!horizontal){
a1<-axis(1 , at=(pretty(c(yrange[1],yrange[2]), n = 10)) , labels=rep('', 11), col.ticks=a.coladd.ade(bgcol, -35), lwd.ticks=3)
a1<-axis(1 , at=(pretty(c(yrange[1],yrange[2]), n = 10)) , labels=rep('', 11), col.ticks=rgb(1,1,1), lwd.ticks=1)
axis(2,  at=1:length(vars), labels=vnames, col.ticks=a.coladd.ade(bgcol, -35), lwd.ticks=3)
axis(2,  at=1:length(vars), labels=vnames, col.ticks=rgb(1,1,1), lwd.ticks=1)
abline(v=a1, lty=1, col=a.coladd.ade(bgcol, -35), lwd=3)
abline(v=a1, lty=1, col=rgb(1,1,1), lwd=1)
}


runpoly(M, rgb(1,1,1), tcol, horizontal=horizontal)

box(lwd=3, col=rgb(1,1,1))
box(lwd=1, col=a.coladd.ade(bgcol, -35))

}
################################################################################

}
