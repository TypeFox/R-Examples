histogram.ade <-
function(x,  group=NULL, w=NULL, data=NULL, vnames=NULL, freq=FALSE,  breaks="Sturges", density=NULL, angle = NULL, xlab=NULL, ylab=NULL, main='', xlim=NULL, ylim=NULL, legendon='topright', xticks=NULL, col=NULL, tcol=NULL,  bgcol=NULL, lcol=NULL, alpha=NULL, lwd=1, kern=TRUE, norm=TRUE, bars=TRUE, wall=0, v=NULL, h=NULL, lty=2){
if(any(par('mfg')!=c(1,1,1,1)) & any(par('mai') < c(1.02, 0.82, 0.82, 0.42))){
maidiff<-rep(0, 4)
norm_f<-c(1.02, 0.82, 0.82, 0.42)
maidiff[par('mai')<norm_f]<-  norm_f[par('mai')<norm_f] - par('mai')[par('mai')<norm_f]
par(mai=par('mai')+maidiff)
}
oldpar<-par(no.readonly =TRUE)
oldpar<-oldpar[-which(names(oldpar)%in%c('usr', 'plt',   'pin', 'fin', 'fig', 'mfg', 'mfcol', 'mfrow', 'omd', 'omi', 'oma'))]
on.exit(par(oldpar))



##############################
if(!is.character(x)){
xt<-deparse(substitute(x))
if(regexpr('~', xt)>=0){
x<-gsub('[~].*$', '', xt)
xpart<-gsub('^.*[~]', '', xt)
group<-gsub('[+].*$', '', xpart)
}}
##############################


#####################################
a.wtd.mean<-function(x, weights = NULL, normwt = "ignored", na.rm = TRUE)
{
    if (!length(weights))
        return(mean(x, na.rm = na.rm))
    if (na.rm) {
        s <- !is.na(x + weights)
        x <- x[s]
        weights <- weights[s]
    }
    sum(weights * x)/sum(weights)
}
#####################################

#####################################
a.wtd.var <- function(x, weights = NULL, normwt = FALSE, na.rm = TRUE)
{
    if (!length(weights)) {
        if (na.rm)
            x <- x[!is.na(x)]
        return(var(x))
    }
    if (na.rm) {
        s <- !is.na(x + weights)
        x <- x[s]
        weights <- weights[s]
    }
    if (normwt)
        weights <- weights * length(x)/sum(weights)
    xbar <- sum(weights * x)/sum(weights)
    sum(weights * ((x - xbar)^2))/(sum(weights) - 1)
}
#####################################

#####################################
# without Data.Frame or with
ismitdata=FALSE
if(is.numeric(x)){
data<-NULL
xname<-gsub('[(]{0}[A-Za-z0-9]*[$]', '' ,deparse(substitute(x)))
}
if(is.character(x)){
ismitdata=TRUE
if(ismitdata){
if(!is.null(data)){ if(!is.data.frame(data))  stop("(data) must be a data.frame!") }
}
xname<-x
if(!is.null(group)){
groupname<-group
group<-eval(parse(text=paste("data$",group, sep='')))
}
x<-eval(parse(text=paste("data$",x, sep='')))
}
if(is.character(w)){
w<-eval(parse(text=paste("data$",w, sep='')))
}
if(is.null(w)) w <- rep(1, length(x))


g=1
# Errors
if(!is.null(group))  {
if(!is.factor(group))   group<-as.factor(group)
if(nlevels(group )>6) stop('To many levels in group!')
g<-group
if(is.null(vnames)) vnames <- levels(g)
}

if(!is.numeric(x))     stop('x is not numeric!')
#####################################




#####################################
# Colors
if(is.null(tcol)  & wall==0)   tcol<-1
if(is.null(tcol)  & wall!=0)   tcol<-rgb(0.1,0.1,0.25)
if(is.null(bgcol) & wall==0)   bgcol<-1
if(is.null(bgcol) & wall!=0)   bgcol<-'#DBE0E8'
if(is.null(col) & wall==0 & is.null(group))   col<-'gray75'
if(is.null(col) & wall!=0 & is.null(group))   col<-rgb(0.6,0.6,0.75)


if(!is.null(group) & is.null(col)) {
if(nlevels(group)>1)  col <- epade:::a.getcol.ade(nlevels(group))
}


if(is.null(alpha) & wall==0) alpha<-1
if(is.null(alpha) & wall!=0 & is.null(group)) alpha<-1
if(is.null(alpha) & wall!=0 & nlevels(group)==2) alpha<-0.5
if(is.null(alpha) & wall!=0 & nlevels(group)>2 & nlevels(group)<=4 ) alpha<-0.33
if(is.null(alpha) & wall!=0 & nlevels(group)>4 ) alpha<-0.25


col<- epade:::a.alpha.ade(col, alpha)


if(is.null(lcol))  lcol<- tcol

col2<-epade:::a.alpha.ade(col, 1)
if(wall==0 & is.null(density)) density<-16
if(wall!=0 & is.null(density)) density<-NA
if(is.null(angle)){
if(nlevels(group)==2)  angle<- c(45, 135)
if(nlevels(group)==3)  angle<- c(45, 0, 135)
if(nlevels(group)==4)  angle<- c(45, 135, 0, 90)
if(nlevels(group)>4)  angle<- ((1:nlevels(group))*30)-30
}

#####################################

yrange<-range(x, na.rm=TRUE)



if(is.null(group))  xrange<-range(x, na.rm=TRUE)
if(!is.null(group)) xrange<-range(x[!is.na(g)], na.rm=TRUE)
if(!is.null(xlim))  xrange<- xlim
if(is.null(xlim) & (!is.numeric(breaks) | length(breaks)==1))   xlim<- xrange
if(is.null(xlim) & (is.numeric(breaks) & length(breaks)>1))     xlim<- range(breaks)
if(!is.null(ylim)){
if(length(ylim)==1) ylim <- c(0, ylim)
}


if(is.null(ylab)) ylab <-'Density'
if(is.null(xlab)) xlab <-xname




################################################################################
################################################################################
if(wall==0){
################################################################################
par(col.main=tcol)
par(col.axis=tcol)
par(col.lab =tcol)


################################################################################
if(is.null(group)){



# Erzeugen Hist
HH<-weighted.hist(x[which(!is.na(x))], w[which(!is.na(x))], breaks=breaks, plot=FALSE, freq=F)
HH$density <- HH$density/diff(HH$breaks)



# Limits
if(is.null(ylim) & !freq) ylim<-c(0, max(HH$density, na.rm=TRUE)+(max(HH$density, na.rm=TRUE)/10))
if(is.null(ylim) &  freq) ylim<-c(0, max(HH$counts, na.rm=TRUE) +(max(HH$counts,  na.rm=TRUE)/10))

# Plot
plot(0,0, col=rgb(0,0,0,0), main=main, xlab=xlab, ylab=ylab, ylim=ylim, xlim=xlim, axes=F)
if(bars & !freq)  rect(HH$breaks[-length(HH$breaks)] , 0, HH$breaks[-1], HH$density,  angle=angle, density=density, col=col, border=epade:::a.coladd.ade(col, -75), lwd=1)
if(bars &  freq)  rect(HH$breaks[-length(HH$breaks)] , 0, HH$breaks[-1], HH$counts,   angle=angle, density=density, col=col, border=epade:::a.coladd.ade(col, -75), lwd=1)
abline(v=v, h=h, col=lcol, lty=lty, lwd=lwd)

# Ticks
if(is.null(xticks))                      axis(1, labels=TRUE, tick=TRUE, col=bgcol, lwd.ticks=1, col.ticks=bgcol)
if(!is.null(xticks) & length(xticks)==1) axis(1, at=pretty(x, n=xticks), tick=TRUE, col=bgcol, lwd.ticks=1, col.ticks=bgcol)
if(!is.null(xticks) & length(xticks)>1 ) axis(1, at=xticks, tick=TRUE, col=bgcol, lwd.ticks=1, col.ticks=bgcol)
axis(2, labels=TRUE, tick=TRUE, col=bgcol, lwd.ticks=1, col.ticks=bgcol)


# norm Lines
xfit<-seq(min(x, na.rm=TRUE),max(x, na.rm=TRUE),length=length(HH$mids)*25)
yfit<-dnorm(xfit, mean=a.wtd.mean(x[which(!is.na(x))],w[which(!is.na(x))], na.rm=TRUE),sd=sqrt(a.wtd.var(x[which(!is.na(x))], w[which(!is.na(x))], na.rm=TRUE)))
yfit<- yfit
if(freq) yfit <-(yfit*sum(HH$counts))
if(norm) lines(xfit, yfit, col=epade:::a.alpha.ade(tcol, 1), lwd=lwd, lty=2)


# Density
w2<- w[which(!is.na(x))]/(sum(w[which(!is.na(x))]))
dd<-density(x[which(!is.na(x))], bw='nrd',na.rm =TRUE, weights =w2 )
dd$y<- (dd$y)
if(freq) dd$y <-(dd$y*sum(HH$counts))
if(kern & !bars) polygon(dd, col=col, lwd=lwd, border=epade:::a.coladd.ade(col, -75), angle=angle, density=density)
if(kern & bars ) lines(dd,   col=epade:::a.alpha.ade(tcol, 1), lwd=lwd)


box(col=bgcol)
}
################################################################################
################################################################################

if(!is.null(group)){


# Erzeugen Hist
HH<-weighted.hist(x[which(!is.na(x) & !is.na(g))], w[which(!is.na(x) & !is.na(g))], breaks=breaks, plot=FALSE, freq=F)
HH$density <- HH$density/diff(HH$breaks)

# Limits
if(is.null(ylim)){
for(i in 1:nlevels(g)){
HH2<-weighted.hist(x[which(g==levels(g)[i] & !is.na(x))], w[which(g==levels(g)[i] & !is.na(x))], breaks=HH$breaks, plot=FALSE, freq=F)
HH2$density <- HH2$density/diff(HH2$breaks)
if(!freq) ylim<-c(0, max(ylim, max(HH2$density, na.rm=T)+max(HH2$density, na.rm=T)/8, na.rm=T))
if(freq)  ylim<-c(0, max(ylim, max(HH2$counts,  na.rm=T)+max(HH2$counts,  na.rm=T)/8, na.rm=T))
}
}


# Plot
plot(0,0, col=rgb(0,0,0,0), main=main, xlab=xlab, ylab=ylab, ylim=ylim, xlim=xlim, axes=F)
abline(v=v, h=h, col=lcol, lty=lty, lwd=lwd)


# Ticks
if(is.null(xticks))                      axis(1, labels=TRUE, tick=TRUE, col=bgcol, lwd.ticks=1, col.ticks=bgcol)
if(!is.null(xticks) & length(xticks)==1) axis(1, at=pretty(x, n=xticks), tick=TRUE, col=bgcol, lwd.ticks=1, col.ticks=bgcol)
if(!is.null(xticks) & length(xticks)>1 ) axis(1, at=xticks, tick=TRUE, col=bgcol, lwd.ticks=1, col.ticks=bgcol)
axis(2, labels =TRUE, tick=TRUE, col=bgcol, lwd.ticks=1, col.ticks=bgcol)


#Legend
legend(legendon, legend=vnames, pt.cex=1.3 ,col = col2,  pch=15, box.col=bgcol,  border=bgcol, text.col=tcol, pt.bg=col2, text.width=max(strwidth(vnames,font = 2)))
legend(legendon, legend=vnames, pt.cex=1.3 ,col = bgcol, pch=0,  box.col=bgcol,  bg=epade:::a.alpha.ade(1,0), border=bgcol, text.col=tcol, text.width=max(strwidth(vnames,font = 2)))



# Gruppen inrun
#######################
for(i in 1:nlevels(g)){
xsub<- x[which(g==levels(g)[i] & !is.na(x))]
wsub<- w[which(g==levels(g)[i] & !is.na(x))]
b<-  max(diff(HH$breaks))


# Bars
HH2<-weighted.hist(xsub, wsub, breaks=HH$breaks, plot=FALSE, freq=F)
HH2$density <- HH2$density/diff(HH2$breaks)
if(bars & !freq) rect(HH2$breaks[-length(HH2$breaks)] , 0, HH2$breaks[-1], HH2$density,  angle=angle[i], density=density, col=col[i], border=epade:::a.alpha.ade(epade:::a.coladd.ade(col[i], -100),1), lwd=1)
if(bars & freq)  rect(HH2$breaks[-length(HH2$breaks)] , 0, HH2$breaks[-1], HH2$counts,   angle=angle[i], density=density, col=col[i], border=epade:::a.alpha.ade(epade:::a.coladd.ade(col[i], -100),1), lwd=1)



# norm Lines
xfit<-seq(min(xsub, na.rm=TRUE),max(xsub, na.rm=TRUE),length=100)
yfit<-dnorm(xfit, mean=a.wtd.mean(xsub,wsub, na.rm=TRUE),sd=sqrt(a.wtd.var(xsub,wsub, na.rm=TRUE)))
yfit<- yfit
if(freq) yfit <-(yfit*sum(HH2$counts))
if(norm) lines(xfit, yfit, col=epade:::a.alpha.ade(col[i], 1), lwd=lwd, lty=2)


# Density
w2<- wsub/sum(wsub)
dd<-density(xsub, bw='nrd',na.rm =TRUE, weights =w2)
dd$y<- dd$y
if(freq) dd$y <-(dd$y*sum(HH2$counts))
if(kern & !bars) polygon(dd, col=col[i], lwd=lwd, border=epade:::a.coladd.ade(col[i], -75), angle=angle[i], density=density)
if(kern & bars ) lines(dd,  col=epade:::a.alpha.ade(col[i], 1), lwd=lwd)

}
#######################


box(col=bgcol)
}
}
################################################################################
################################################################################




################################################################################
################################################################################
if(wall==1){
################################################################################
par(col.main=tcol)
par(col.axis=tcol)
par(col.lab =tcol)


################################################################################
if(is.null(group)){

# Erzeugen Hist
HH<-weighted.hist(x[which(!is.na(x))], w[which(!is.na(x))], breaks=breaks, plot=FALSE, freq=F)
HH$density <- HH$density/diff(HH$breaks)


# Limits
if(is.null(ylim) & !freq) ylim<-c(0, max(HH$density, na.rm=TRUE)+(max(HH$density, na.rm=TRUE)/10))
if(is.null(ylim) &  freq) ylim<-c(0, max(HH$counts, na.rm=TRUE) +(max(HH$counts,  na.rm=TRUE)/10))

# Plot
plot(0,0, col=rgb(0,0,0,0), main=main, xlab=xlab, ylab=ylab, ylim=ylim, xlim=xlim, axes=F)
polygon( c(par('usr')[c(1,1,2,2)]), par('usr')[c(3,4,4,3)], col=bgcol, border=FALSE)
if(is.null(xticks))                      a1<-axis(1, labels=TRUE, tick=TRUE, col=rgb(1,1,1), lwd.ticks=1, col.ticks=tcol)
if(!is.null(xticks) & length(xticks)==1) a1<-axis(1, at=pretty(x, n=xticks), tick=TRUE, col=rgb(1,1,1), lwd.ticks=1, col.ticks=tcol)
if(!is.null(xticks) & length(xticks)>1 ) a1<-axis(1, at=xticks, tick=TRUE, col=rgb(1,1,1), lwd.ticks=1, col.ticks=tcol)

a2<-axis(2, col=rgb(1,1,1), col.ticks=tcol, lwd.ticks=1)
abline(v=a1, h=a2, lty=1, col=rgb(1,1,1), lwd=1)
abline(v=v, h=h, col=tcol, lty=lty, lwd=lwd)

if(bars & !freq)  rect(HH$breaks[-length(HH$breaks)] , 0, HH$breaks[-1], HH$density,  angle=angle, density=density, col=col, border=epade:::a.coladd.ade(col, -75), lwd=1)
if(bars &  freq)  rect(HH$breaks[-length(HH$breaks)] , 0, HH$breaks[-1], HH$counts,   angle=angle, density=density, col=col, border=epade:::a.coladd.ade(col, -75), lwd=1)



# norm Lines

xfit<-seq(min(x, na.rm=TRUE),max(x, na.rm=TRUE),length=length(HH$mids)*25)
yfit<-dnorm(xfit, mean=a.wtd.mean(x[which(!is.na(x))],w[which(!is.na(x))], na.rm=TRUE),sd=sqrt(a.wtd.var(x[which(!is.na(x))], w[which(!is.na(x))], na.rm=TRUE)))
yfit<- yfit
if(freq) yfit <-(yfit*sum(HH$counts))
if(norm) lines(xfit, yfit, col=epade:::a.alpha.ade(tcol, 1), lwd=lwd, lty=2)


# Density
w2<- w[which(!is.na(x))]/(sum(w[which(!is.na(x))]))
dd<-density(x[which(!is.na(x))], bw='nrd',na.rm =TRUE, weights =w2 )
dd$y<- (dd$y)
if(freq) dd$y <-(dd$y*sum(HH$counts))
if(kern & !bars) polygon(dd, col=col, lwd=lwd, border=epade:::a.coladd.ade(col, -75), angle=angle, density=density)
if(kern & bars ) lines(dd,   col=epade:::a.alpha.ade(tcol, 1), lwd=lwd)


box(lwd=2, col=rgb(1,1,1))
}
################################################################################
################################################################################

if(!is.null(group)){


# Erzeugen Hist
HH<-weighted.hist(x[which(!is.na(x) & !is.na(g))], w[which(!is.na(x) & !is.na(g))], breaks=breaks, plot=FALSE, freq=F)
HH$density <- HH$density/diff(HH$breaks)

# Limits
if(is.null(ylim)){
for(i in 1:nlevels(g)){
HH2<-weighted.hist(x[which(g==levels(g)[i] & !is.na(x))], w[which(g==levels(g)[i] & !is.na(x))], breaks=HH$breaks, plot=FALSE, freq=F)
HH2$density <- HH2$density/diff(HH2$breaks)
if(!freq) ylim<-c(0, max(ylim, max(HH2$density, na.rm=T)+max(HH2$density, na.rm=T)/8, na.rm=T))
if(freq)  ylim<-c(0, max(ylim, max(HH2$counts,  na.rm=T)+max(HH2$counts,  na.rm=T)/8, na.rm=T))
}
}


# Plot
plot(0,0, col=rgb(0,0,0,0), main=main, xlab=xlab, ylab=ylab, ylim=ylim, xlim=xlim, axes=F)
polygon( c(par('usr')[c(1,1,2,2)]), par('usr')[c(3,4,4,3)], col=bgcol, border=FALSE)
if(is.null(xticks))                      a1<-axis(1, labels=TRUE, tick=TRUE, col=rgb(1,1,1), lwd.ticks=1, col.ticks=tcol)
if(!is.null(xticks) & length(xticks)==1) a1<-axis(1, at=pretty(x, n=xticks), tick=TRUE, col=rgb(1,1,1), lwd.ticks=1, col.ticks=tcol)
if(!is.null(xticks) & length(xticks)>1 ) a1<-axis(1, at=xticks, tick=TRUE, col=rgb(1,1,1), lwd.ticks=1, col.ticks=tcol)

a2<-axis(2, col=rgb(1,1,1), col.ticks=tcol, lwd.ticks=1)
abline(v=a1, h=a2, lty=1, col=rgb(1,1,1), lwd=1)
abline(v=v, h=h, col=tcol, lty=lty, lwd=lwd)
legend(legendon, legend=vnames, pt.cex=1.5 ,col = col2, pch=15, bg=bgcol, border=rgb(1,1,1), box.col=rgb(1,1,1), box.lwd=2, text.col=tcol, text.width=max(strwidth(vnames,font = 2)))
legend(legendon, legend=vnames, pt.cex=1.5 ,col = rgb(1,1,1), pch=0, bg=epade:::a.alpha.ade(1, 0), border=rgb(1,1,1), box.col=epade:::a.alpha.ade(1, 0), box.lwd=2, text.col=epade:::a.alpha.ade(1, 0), text.width=max(strwidth(vnames,font = 2)))





# Gruppen inrun
#######################
for(i in 1:nlevels(g)){
xsub<- x[which(g==levels(g)[i] & !is.na(x))]
wsub<- w[which(g==levels(g)[i] & !is.na(x))]
b<-  max(diff(HH$breaks))


# Bars
HH2<-weighted.hist(xsub, wsub, breaks=HH$breaks, plot=FALSE, freq=F)
HH2$density <- HH2$density/diff(HH2$breaks)
if(bars & !freq) rect(HH2$breaks[-length(HH2$breaks)] , 0, HH2$breaks[-1], HH2$density,  angle=angle[i], density=density, col=col[i], border=epade:::a.alpha.ade(epade:::a.coladd.ade(col[i], -100),1), lwd=1)
if(bars & freq)  rect(HH2$breaks[-length(HH2$breaks)] , 0, HH2$breaks[-1], HH2$counts,   angle=angle[i], density=density, col=col[i], border=epade:::a.alpha.ade(epade:::a.coladd.ade(col[i], -100),1), lwd=1)



# norm Lines
xfit<-seq(min(xsub, na.rm=TRUE),max(xsub, na.rm=TRUE),length=100)
yfit<-dnorm(xfit, mean=a.wtd.mean(xsub,wsub, na.rm=TRUE),sd=sqrt(a.wtd.var(xsub,wsub, na.rm=TRUE)))
yfit<- yfit
if(freq) yfit <-(yfit*sum(HH2$counts))
if(norm) lines(xfit, yfit, col=epade:::a.alpha.ade(col[i], 1), lwd=lwd, lty=2)


# Density
w2<- wsub/sum(wsub)
dd<-density(xsub, bw='nrd',na.rm =TRUE, weights =w2)
dd$y<- dd$y
if(freq) dd$y <-(dd$y*sum(HH2$counts))
if(kern & !bars) polygon(dd, col=col[i], lwd=lwd, border=epade:::a.coladd.ade(col[i], -75), angle=angle[i], density=density)
if(kern & bars ) lines(dd,  col=epade:::a.alpha.ade(col[i], 1), lwd=lwd)

}
#######################


box(lwd=2, col=rgb(1,1,1))
}
}
################################################################################
################################################################################



################################################################################
################################################################################
if(wall==2){
################################################################################
par(col.main=tcol)
par(col.axis=tcol)
par(col.lab =tcol)


################################################################################
if(is.null(group)){

# Erzeugen Hist
HH<-weighted.hist(x[which(!is.na(x))], w[which(!is.na(x))], breaks=breaks, plot=FALSE, freq=F)
HH$density <- HH$density/diff(HH$breaks)


# Limits
if(is.null(ylim) & !freq) ylim<-c(0, max(HH$density, na.rm=TRUE)+(max(HH$density, na.rm=TRUE)/10))
if(is.null(ylim) &  freq) ylim<-c(0, max(HH$counts, na.rm=TRUE) +(max(HH$counts,  na.rm=TRUE)/10))

# Plot
plot(0,0, col=rgb(0,0,0,0), main=main, xlab=xlab, ylab=ylab, ylim=ylim, xlim=xlim, axes=F)

# Ticks
if(is.null(xticks))                      a1<-axis(1, labels=TRUE, tick=TRUE, col=rgb(1,1,1), lwd.ticks=1, col.ticks=epade:::a.coladd.ade(bgcol, -75))
if(!is.null(xticks) & length(xticks)==1) a1<-axis(1, at=pretty(x, n=xticks), tick=TRUE, col=rgb(1,1,1), lwd.ticks=1, col.ticks=epade:::a.coladd.ade(bgcol, -75))
if(!is.null(xticks) & length(xticks)>1 ) a1<-axis(1, at=xticks, tick=TRUE, col=rgb(1,1,1), lwd.ticks=1, col.ticks=epade:::a.coladd.ade(bgcol, -75))
a2<-axis(2, col=rgb(1,1,1), col.ticks=epade:::a.coladd.ade(bgcol, -75), lwd.ticks=1)
abline(v=a1, h=a2, lty=1, col=bgcol, lwd=1)
abline(v=v, h=h, col=lcol, lty=lty, lwd=lwd)

if(bars & !freq)  rect(HH$breaks[-length(HH$breaks)] , 0, HH$breaks[-1], HH$density,  angle=angle, density=density, col=col, border=epade:::a.coladd.ade(col, -75), lwd=1)
if(bars &  freq)  rect(HH$breaks[-length(HH$breaks)] , 0, HH$breaks[-1], HH$counts,   angle=angle, density=density, col=col, border=epade:::a.coladd.ade(col, -75), lwd=1)

# norm Lines
xfit<-seq(min(x, na.rm=TRUE),max(x, na.rm=TRUE),length=length(HH$mids)*25)
yfit<-dnorm(xfit, mean=a.wtd.mean(x[which(!is.na(x))],w[which(!is.na(x))], na.rm=TRUE),sd=sqrt(a.wtd.var(x[which(!is.na(x))], w[which(!is.na(x))], na.rm=TRUE)))
yfit<- yfit
if(freq) yfit <-(yfit*sum(HH$counts))
if(norm) lines(xfit, yfit, col=epade:::a.alpha.ade(tcol, 1), lwd=lwd, lty=2)


# Density
w2<- w[which(!is.na(x))]/(sum(w[which(!is.na(x))]))
dd<-density(x[which(!is.na(x))], bw='nrd',na.rm =TRUE, weights =w2 )
dd$y<- (dd$y)
if(freq) dd$y <-(dd$y*sum(HH$counts))
if(kern & !bars) polygon(dd, col=col, lwd=lwd, border=epade:::a.coladd.ade(col, -75), angle=angle, density=density)
if(kern & bars ) lines(dd,   col=epade:::a.alpha.ade(tcol, 1), lwd=lwd)


box(col=epade:::a.coladd.ade(bgcol, -75))
}
################################################################################
################################################################################

if(!is.null(group)){


# Erzeugen Hist
HH<-weighted.hist(x[which(!is.na(x) & !is.na(g))], w[which(!is.na(x) & !is.na(g))], breaks=breaks, plot=FALSE, freq=F)
HH$density <- HH$density/diff(HH$breaks)

# Limits
if(is.null(ylim)){
for(i in 1:nlevels(g)){
HH2<-weighted.hist(x[which(g==levels(g)[i] & !is.na(x))], w[which(g==levels(g)[i] & !is.na(x))], breaks=HH$breaks, plot=FALSE, freq=F)
HH2$density <- HH2$density/diff(HH2$breaks)
if(!freq) ylim<-c(0, max(ylim, max(HH2$density, na.rm=T)+max(HH2$density, na.rm=T)/8, na.rm=T))
if(freq)  ylim<-c(0, max(ylim, max(HH2$counts,  na.rm=T)+max(HH2$counts,  na.rm=T)/8, na.rm=T))
}
}


# Plot
plot(0,0, col=rgb(0,0,0,0), main=main, xlab=xlab, ylab=ylab, ylim=ylim, xlim=xlim, axes=F)
if(is.null(xticks))                      a1<-axis(1, labels=TRUE, tick=TRUE, col=rgb(1,1,1), lwd.ticks=1, col.ticks=epade:::a.coladd.ade(bgcol, -75))
if(!is.null(xticks) & length(xticks)==1) a1<-axis(1, at=pretty(x, n=xticks), tick=TRUE, col=rgb(1,1,1), lwd.ticks=1, col.ticks=epade:::a.coladd.ade(bgcol, -75))
if(!is.null(xticks) & length(xticks)>1 ) a1<-axis(1, at=xticks, tick=TRUE, col=rgb(1,1,1), lwd.ticks=1, col.ticks=epade:::a.coladd.ade(bgcol, -75))
a2<-axis(2, col=rgb(1,1,1), col.ticks=epade:::a.coladd.ade(bgcol, -75), lwd.ticks=1)
abline(v=a1, h=a2, lty=1, col=bgcol, lwd=1)
abline(v=v, h=h, col=lcol, lty=lty, lwd=lwd)


#Legend
legend(legendon, legend=vnames, pt.cex=1.5 ,col = col2, pch=15, bg=rgb(1,1,1), box.col=epade:::a.coladd.ade(bgcol, -75), box.lwd=1, text.col=tcol, text.width=max(strwidth(vnames,font = 2)))
legend(legendon, legend=vnames, pt.cex=1.5 ,col = epade:::a.coladd.ade(bgcol, -75), pch=0,  bg=epade:::a.alpha.ade(1,0), box.col=epade:::a.alpha.ade(1,0), box.lwd=1, text.col=epade:::a.alpha.ade(1,0), text.width=max(strwidth(vnames,font = 2)))



# Gruppen inrun
#######################
for(i in 1:nlevels(g)){
xsub<- x[which(g==levels(g)[i] & !is.na(x))]
wsub<- w[which(g==levels(g)[i] & !is.na(x))]
b<-  max(diff(HH$breaks))


# Bars
HH2<-weighted.hist(xsub, wsub, breaks=HH$breaks, plot=FALSE, freq=F)
HH2$density <- HH2$density/diff(HH2$breaks)
if(bars & !freq) rect(HH2$breaks[-length(HH2$breaks)] , 0, HH2$breaks[-1], HH2$density,  angle=angle[i], density=density, col=col[i], border=epade:::a.alpha.ade(epade:::a.coladd.ade(col[i], -100),1), lwd=1)
if(bars & freq)  rect(HH2$breaks[-length(HH2$breaks)] , 0, HH2$breaks[-1], HH2$counts,   angle=angle[i], density=density, col=col[i], border=epade:::a.alpha.ade(epade:::a.coladd.ade(col[i], -100),1), lwd=1)



# norm Lines
xfit<-seq(min(xsub, na.rm=TRUE),max(xsub, na.rm=TRUE),length=100)
yfit<-dnorm(xfit, mean=a.wtd.mean(xsub,wsub, na.rm=TRUE),sd=sqrt(a.wtd.var(xsub,wsub, na.rm=TRUE)))
yfit<- yfit
if(freq) yfit <-(yfit*sum(HH2$counts))
if(norm) lines(xfit, yfit, col=epade:::a.alpha.ade(col[i], 1), lwd=lwd, lty=2)


# Density
w2<- wsub/sum(wsub)
dd<-density(xsub, bw='nrd',na.rm =TRUE, weights =w2)
dd$y<- dd$y
if(freq) dd$y <-(dd$y*sum(HH2$counts))
if(kern & !bars) polygon(dd, col=col[i], lwd=lwd, border=epade:::a.coladd.ade(col[i], -75), angle=angle[i], density=density)
if(kern & bars ) lines(dd,  col=epade:::a.alpha.ade(col[i], 1), lwd=lwd)

}
#######################


box(col=epade:::a.coladd.ade(bgcol, -75))
}
}
################################################################################
################################################################################



################################################################################
################################################################################
if(wall==3){
################################################################################
par(col.main=tcol)
par(col.axis=tcol)
par(col.lab =tcol)




################################################################################
if(is.null(group)){

# Erzeugen Hist
HH<-weighted.hist(x[which(!is.na(x))], w[which(!is.na(x))], breaks=breaks, plot=FALSE, freq=F)
HH$density <- HH$density/diff(HH$breaks)


# Limits
if(is.null(ylim) & !freq) ylim<-c(0, max(HH$density, na.rm=TRUE)+(max(HH$density, na.rm=TRUE)/10))
if(is.null(ylim) &  freq) ylim<-c(0, max(HH$counts, na.rm=TRUE) +(max(HH$counts,  na.rm=TRUE)/10))

# Plot
plot(0,0, col=rgb(0,0,0,0), main=main, xlab=xlab, ylab=ylab, ylim=ylim, xlim=xlim, axes=F)
polygon( c(par('usr')[c(1,1,2,2)]), par('usr')[c(3,4,4,3)], col=bgcol, border=FALSE)
if(is.null(xticks))                      a1<-axis(1, labels=TRUE, tick=TRUE, col=rgb(1,1,1), lwd.ticks=1, col.ticks=epade:::a.coladd.ade(bgcol, -50))
if(!is.null(xticks) & length(xticks)==1) a1<-axis(1, at=pretty(x, n=xticks), tick=TRUE, col=rgb(1,1,1), lwd.ticks=1, col.ticks=epade:::a.coladd.ade(bgcol, -50))
if(!is.null(xticks) & length(xticks)>1 ) a1<-axis(1, at=xticks, tick=TRUE, col=rgb(1,1,1), lwd.ticks=1, col.ticks=epade:::a.coladd.ade(bgcol, -50))
a2<-axis(2, col=rgb(1,1,1), col.ticks=epade:::a.coladd.ade(bgcol, -50), lwd.ticks=1)
abline(v=a1, h=a2, lty=1, col=epade:::a.coladd.ade(bgcol, -50), lwd=1)
abline(v=v, h=h, col=lcol, lty=lty, lwd=lwd)


if(bars & !freq)  rect(HH$breaks[-length(HH$breaks)] , 0, HH$breaks[-1], HH$density,  angle=angle, density=density, col=col, border=epade:::a.coladd.ade(col, -75), lwd=1)
if(bars &  freq)  rect(HH$breaks[-length(HH$breaks)] , 0, HH$breaks[-1], HH$counts,   angle=angle, density=density, col=col, border=epade:::a.coladd.ade(col, -75), lwd=1)

# norm Lines
b<-  max(diff(HH$breaks))
xfit<-seq(min(x, na.rm=TRUE),max(x, na.rm=TRUE),length=length(HH$mids)*25)
yfit<-dnorm(xfit, mean=a.wtd.mean(x[which(!is.na(x))],w[which(!is.na(x))], na.rm=TRUE),sd=sqrt(a.wtd.var(x[which(!is.na(x))], w[which(!is.na(x))], na.rm=TRUE)))
yfit<- yfit
if(freq) yfit <-(yfit*sum(HH$counts))
if(norm) lines(xfit, yfit, col=epade:::a.alpha.ade(tcol, 1), lwd=lwd, lty=2)


# Density
w2<- w[which(!is.na(x))]/(sum(w[which(!is.na(x))]))
dd<-density(x[which(!is.na(x))], bw='nrd',na.rm =TRUE, weights =w2 )
dd$y<- (dd$y)
if(freq) dd$y <-(dd$y*sum(HH$counts))
if(kern & !bars) polygon(dd, col=col, lwd=lwd, border=epade:::a.coladd.ade(col, -75), angle=angle, density=density)
if(kern & bars ) lines(dd,   col=epade:::a.alpha.ade(tcol, 1), lwd=lwd)


box(col=epade:::a.coladd.ade(bgcol, -50))
}
################################################################################
################################################################################



if(!is.null(group)){


# Erzeugen Hist
HH<-weighted.hist(x[which(!is.na(x) & !is.na(g))], w[which(!is.na(x) & !is.na(g))], breaks=breaks, plot=FALSE, freq=F)
HH$density <- HH$density/diff(HH$breaks)

# Limits
if(is.null(ylim)){
for(i in 1:nlevels(g)){
HH2<-weighted.hist(x[which(g==levels(g)[i] & !is.na(x))], w[which(g==levels(g)[i] & !is.na(x))], breaks=HH$breaks, plot=FALSE, freq=F)
HH2$density <- HH2$density/diff(HH2$breaks)
if(!freq) ylim<-c(0, max(ylim, max(HH2$density, na.rm=T)+max(HH2$density, na.rm=T)/8, na.rm=T))
if(freq)  ylim<-c(0, max(ylim, max(HH2$counts,  na.rm=T)+max(HH2$counts,  na.rm=T)/8, na.rm=T))
}
}



# Plot
plot(0,0, col=rgb(0,0,0,0), main=main, xlab=xlab, ylab=ylab, ylim=ylim, xlim=xlim, axes=F)
polygon( c(par('usr')[c(1,1,2,2)]), par('usr')[c(3,4,4,3)], col=bgcol, border=FALSE)
if(is.null(xticks))                      a1<-axis(1, labels=TRUE, tick=TRUE, col=rgb(1,1,1), lwd.ticks=1, col.ticks=epade:::a.coladd.ade(bgcol, -50))
if(!is.null(xticks) & length(xticks)==1) a1<-axis(1, at=pretty(x, n=xticks), tick=TRUE, col=rgb(1,1,1), lwd.ticks=1, col.ticks=epade:::a.coladd.ade(bgcol, -50))
if(!is.null(xticks) & length(xticks)>1 ) a1<-axis(1, at=xticks, tick=TRUE, col=rgb(1,1,1), lwd.ticks=1, col.ticks=epade:::a.coladd.ade(bgcol, -50))
a2<-axis(2, col=rgb(1,1,1), col.ticks=epade:::a.coladd.ade(bgcol, -50), lwd.ticks=1)
abline(v=a1, h=a2, lty=1, col=epade:::a.coladd.ade(bgcol, -50), lwd=1)
abline(v=v, h=h, col=lcol, lty=lty, lwd=lwd)
legend(legendon, legend=vnames, pt.cex=1.5 ,col = col2, pch=15, bg=rgb(1,1,1), box.col=epade:::a.coladd.ade(bgcol, -50), box.lwd=1, text.col=tcol, text.width=max(strwidth(vnames,font = 2)))
legend(legendon, legend=vnames, pt.cex=1.5 ,col = epade:::a.coladd.ade(bgcol, -50), pch=0 , bg=epade:::a.alpha.ade(1,0), box.col=epade:::a.alpha.ade(1,0), box.lwd=1, text.col=epade:::a.alpha.ade(1,0), text.width=max(strwidth(vnames,font = 2)))


# Gruppen inrun
#######################
for(i in 1:nlevels(g)){
xsub<- x[which(g==levels(g)[i] & !is.na(x))]
wsub<- w[which(g==levels(g)[i] & !is.na(x))]
b<-  max(diff(HH$breaks))


# Bars
HH2<-weighted.hist(xsub, wsub, breaks=HH$breaks, plot=FALSE, freq=F)
HH2$density <- HH2$density/diff(HH2$breaks)
if(bars & !freq) rect(HH2$breaks[-length(HH2$breaks)] , 0, HH2$breaks[-1], HH2$density,  angle=angle[i], density=density, col=col[i], border=epade:::a.alpha.ade(epade:::a.coladd.ade(col[i], -100),1), lwd=1)
if(bars & freq)  rect(HH2$breaks[-length(HH2$breaks)] , 0, HH2$breaks[-1], HH2$counts,   angle=angle[i], density=density, col=col[i], border=epade:::a.alpha.ade(epade:::a.coladd.ade(col[i], -100),1), lwd=1)



# norm Lines
xfit<-seq(min(xsub, na.rm=TRUE),max(xsub, na.rm=TRUE),length=100)
yfit<-dnorm(xfit, mean=a.wtd.mean(xsub,wsub, na.rm=TRUE),sd=sqrt(a.wtd.var(xsub,wsub, na.rm=TRUE)))
yfit<- yfit
if(freq) yfit <-(yfit*sum(HH2$counts))
if(norm) lines(xfit, yfit, col=epade:::a.alpha.ade(col[i], 1), lwd=lwd, lty=2)


# Density
w2<- wsub/sum(wsub)
dd<-density(xsub, bw='nrd',na.rm =TRUE, weights =w2)
dd$y<- dd$y
if(freq) dd$y <-(dd$y*sum(HH2$counts))
if(kern & !bars) polygon(dd, col=col[i], lwd=lwd, border=epade:::a.coladd.ade(col[i], -75), angle=angle[i], density=density)
if(kern & bars ) lines(dd,  col=epade:::a.alpha.ade(col[i], 1), lwd=lwd)

}
#######################


box(col=epade:::a.coladd.ade(bgcol, -50))
}
}
################################################################################
################################################################################



################################################################################
################################################################################
if(wall==4){
################################################################################
par(col.main=tcol)
par(col.axis=tcol)
par(col.lab =tcol)


################################################################################
if(is.null(group)){

# Erzeugen Hist
HH<-weighted.hist(x[which(!is.na(x))], w[which(!is.na(x))], breaks=breaks, plot=FALSE, freq=F)
HH$density <- HH$density/diff(HH$breaks)


# Limits
if(is.null(ylim) & !freq) ylim<-c(0, max(HH$density, na.rm=TRUE)+(max(HH$density, na.rm=TRUE)/10))
if(is.null(ylim) &  freq) ylim<-c(0, max(HH$counts, na.rm=TRUE) +(max(HH$counts,  na.rm=TRUE)/10))

# Plot
plot(0,0, col=rgb(0,0,0,0), main='', xlab='', ylab='', ylim=ylim, xlim=xlim, axes=F)
polygon( c(par('usr')[c(1,1,2,2)]), par('usr')[c(3,4,4,3)], col=bgcol, border=FALSE)
if(is.null(xticks))                      a1<-axis(1, labels=TRUE, tick=TRUE, col=rgb(1,1,1), lwd.ticks=1, col.ticks=tcol)
if(!is.null(xticks) & length(xticks)==1) a1<-axis(1, at=pretty(x, n=xticks), tick=TRUE, col=rgb(1,1,1), lwd.ticks=1, col.ticks=tcol)
if(!is.null(xticks) & length(xticks)>1 ) a1<-axis(1, at=xticks, tick=TRUE, col=rgb(1,1,1), lwd.ticks=1, col.ticks=tcol)
a2<-axis(2, col=rgb(1,1,1), col.ticks=tcol, lwd.ticks=1)
abline(v=a1, h=a2, lty=1, col=rgb(1,1,1), lwd=1)
if(bars & !freq)  rect(HH$breaks[-length(HH$breaks)] , 0, HH$breaks[-1], HH$density,  angle=angle, density=density, col=col, border=epade:::a.coladd.ade(col, -75), lwd=1)
if(bars &  freq)  rect(HH$breaks[-length(HH$breaks)] , 0, HH$breaks[-1], HH$counts,   angle=angle, density=density, col=col, border=epade:::a.coladd.ade(col, -75), lwd=1)
abline(v=v, h=h, col=lcol, lty=lty, lwd=lwd)
# Outer
par(xpd=TRUE)
polygon(epade:::a.glc(side=c(2,2,4,4), line=c(0,0,0,0)), epade:::a.glc(side=3, line=c(0, 2.75,  2.75, 0)), col=tcol, border=rgb(1,1,1))

if(is.expression(ylab))                           polygon( epade:::a.glc(side=2, line=c(3.5, 3.5, 2, 2)), epade:::a.glc(side=c(1, 3, 3, 1), line=0), col=bgcol, border=rgb(1,1,1))
if(is.character(ylab))  if(ylab!='' & ylab!=' ')  polygon( epade:::a.glc(side=2, line=c(3.5, 3.5, 2, 2)), epade:::a.glc(side=c(1, 3, 3, 1), line=0), col=bgcol, border=rgb(1,1,1))

if(is.expression(xlab))                           polygon( epade:::a.glc(side=c(2, 2, 4, 4), line=0),     epade:::a.glc(side=1, line=c(4, 2.5, 2.5, 4)), col=bgcol, border=rgb(1,1,1))
if(is.character(xlab))  if(xlab!='' & xlab!=' ')  polygon( epade:::a.glc(side=c(2, 2, 4, 4), line=0),     epade:::a.glc(side=1, line=c(4, 2.5, 2.5, 4)), col=bgcol, border=rgb(1,1,1))



text(epade:::a.glc(side=0), epade:::a.glc(side=3, line=1),    labels=main, cex = 1.25, font=2, col=rgb(1,1,1), adj=c(0.5,0))
text(epade:::a.glc(side=0), epade:::a.glc(side=1, line=3.5), labels=xlab, cex = 1.1,  font=2, col=tcol, adj=c(0.5,0))
text(epade:::a.glc(side=2, line=2.5), epade:::a.glc(side=5), labels=ylab, cex = 1.1,   font=2,  col=tcol, adj=c(0.5,0), srt=90)
par(xpd=FALSE)

# norm Lines
xfit<-seq(min(x, na.rm=TRUE),max(x, na.rm=TRUE),length=length(HH$mids)*25)
yfit<-dnorm(xfit, mean=a.wtd.mean(x[which(!is.na(x))],w[which(!is.na(x))], na.rm=TRUE),sd=sqrt(a.wtd.var(x[which(!is.na(x))], w[which(!is.na(x))], na.rm=TRUE)))
yfit<- yfit
if(freq) yfit <-(yfit*sum(HH$counts))
if(norm) lines(xfit, yfit, col=epade:::a.alpha.ade(tcol, 1), lwd=lwd, lty=2)


# Density
w2<- w[which(!is.na(x))]/(sum(w[which(!is.na(x))]))
dd<-density(x[which(!is.na(x))], bw='nrd',na.rm =TRUE, weights =w2 )
dd$y<- (dd$y)
if(freq) dd$y <-(dd$y*sum(HH$counts))
if(kern & !bars) polygon(dd, col=col, lwd=lwd, border=epade:::a.coladd.ade(col, -75), angle=angle, density=density)
if(kern & bars ) lines(dd,   col=epade:::a.alpha.ade(tcol, 1), lwd=lwd)


box(col=rgb(1,1,1))
}
################################################################################
################################################################################

if(!is.null(group)){


# Erzeugen Hist
HH<-weighted.hist(x[which(!is.na(x) & !is.na(g))], w[which(!is.na(x) & !is.na(g))], breaks=breaks, plot=FALSE, freq=F)
HH$density <- HH$density/diff(HH$breaks)

# Limits
if(is.null(ylim)){
for(i in 1:nlevels(g)){
HH2<-weighted.hist(x[which(g==levels(g)[i] & !is.na(x))], w[which(g==levels(g)[i] & !is.na(x))], breaks=HH$breaks, plot=FALSE, freq=F)
HH2$density <- HH2$density/diff(HH2$breaks)
if(!freq) ylim<-c(0, max(ylim, max(HH2$density, na.rm=T)+max(HH2$density, na.rm=T)/8, na.rm=T))
if(freq)  ylim<-c(0, max(ylim, max(HH2$counts,  na.rm=T)+max(HH2$counts,  na.rm=T)/8, na.rm=T))
}
}


# Plot
plot(0,0, col=rgb(0,0,0,0), main='', xlab='', ylab='', ylim=ylim, xlim=xlim, axes=F)
polygon( c(par('usr')[c(1,1,2,2)]), par('usr')[c(3,4,4,3)], col=bgcol, border=FALSE)
a2<-axis(2, col=rgb(1,1,1), col.ticks=tcol, lwd.ticks=1)
if(is.null(xticks))                      a1<-axis(1, labels=TRUE, tick=TRUE, col=rgb(1,1,1), lwd.ticks=1, col.ticks=tcol)
if(!is.null(xticks) & length(xticks)==1) a1<-axis(1, at=pretty(x, n=xticks), tick=TRUE, col=rgb(1,1,1), lwd.ticks=1, col.ticks=tcol)
if(!is.null(xticks) & length(xticks)>1 ) a1<-axis(1, at=xticks, tick=TRUE, col=rgb(1,1,1), lwd.ticks=1, col.ticks=tcol)

abline(v=a1, h=a2, lty=1, col=rgb(1,1,1), lwd=1)
abline(v=v, h=h, col=lcol, lty=lty, lwd=lwd)
legend(legendon, legend=vnames, pt.cex=1.5 ,col = col2, pch=15, bg=tcol, border=rgb(1,1,1), box.col=rgb(1,1,1), box.lwd=1, text.col=rgb(1,1,1), text.width=max(strwidth(vnames,font = 2)))
legend(legendon, legend=vnames, pt.cex=1.5 ,col = rgb(1,1,1), pch=0, bg=epade:::a.alpha.ade(1, 0), border=rgb(1,1,1), box.col=epade:::a.alpha.ade(1, 0), box.lwd=2, text.col=epade:::a.alpha.ade(1, 0), text.width=max(strwidth(vnames,font = 2)))
# Outer
par(xpd=TRUE)
polygon(epade:::a.glc(side=c(2,2,4,4), line=c(0,0,0,0)), epade:::a.glc(side=3, line=c(0, 2.75,  2.75, 0)), col=tcol, border=rgb(1,1,1))


if(is.expression(ylab))                          polygon( epade:::a.glc(side=2, line=c(3.5, 3.5, 2, 2)), epade:::a.glc(side=c(1, 3, 3, 1), line=0), col=bgcol, border=rgb(1,1,1))
if(is.character(ylab))  if(ylab!='' & ylab!=' ') polygon( epade:::a.glc(side=2, line=c(3.5, 3.5, 2, 2)), epade:::a.glc(side=c(1, 3, 3, 1), line=0), col=bgcol, border=rgb(1,1,1))

if(is.expression(xlab))                          polygon( epade:::a.glc(side=c(2, 2, 4, 4), line=0),     epade:::a.glc(side=1, line=c(4, 2.5, 2.5, 4)), col=bgcol, border=rgb(1,1,1))
if(is.character(xlab))  if(xlab!='' & xlab!=' ') polygon( epade:::a.glc(side=c(2, 2, 4, 4), line=0),     epade:::a.glc(side=1, line=c(4, 2.5, 2.5, 4)), col=bgcol, border=rgb(1,1,1))



text(epade:::a.glc(side=0), epade:::a.glc(side=3, line=1),    labels=main, cex = 1.25, font=2, col=rgb(1,1,1), adj=c(0.5,0))
text(epade:::a.glc(side=0), epade:::a.glc(side=1, line=3.5), labels=xlab, cex = 1.1,  font=2, col=tcol, adj=c(0.5,0))
text(epade:::a.glc(side=2, line=2.5), epade:::a.glc(side=5),  labels=ylab, cex = 1.1,  font=2,  col=tcol, adj=c(0.5,0), srt=90)
par(xpd=FALSE)


# Gruppen inrun
#######################
for(i in 1:nlevels(g)){
xsub<- x[which(g==levels(g)[i] & !is.na(x))]
wsub<- w[which(g==levels(g)[i] & !is.na(x))]
b<-  max(diff(HH$breaks))


# Bars
HH2<-weighted.hist(xsub, wsub, breaks=HH$breaks, plot=FALSE, freq=F)
HH2$density <- HH2$density/diff(HH2$breaks)
if(bars & !freq) rect(HH2$breaks[-length(HH2$breaks)] , 0, HH2$breaks[-1], HH2$density,  angle=angle[i], density=density, col=col[i], border=epade:::a.alpha.ade(epade:::a.coladd.ade(col[i], -100),1), lwd=1)
if(bars & freq)  rect(HH2$breaks[-length(HH2$breaks)] , 0, HH2$breaks[-1], HH2$counts,   angle=angle[i], density=density, col=col[i], border=epade:::a.alpha.ade(epade:::a.coladd.ade(col[i], -100),1), lwd=1)



# norm Lines
xfit<-seq(min(xsub, na.rm=TRUE),max(xsub, na.rm=TRUE),length=100)
yfit<-dnorm(xfit, mean=a.wtd.mean(xsub,wsub, na.rm=TRUE),sd=sqrt(a.wtd.var(xsub,wsub, na.rm=TRUE)))
yfit<- yfit
if(freq) yfit <-(yfit*sum(HH2$counts))
if(norm) lines(xfit, yfit, col=epade:::a.alpha.ade(col[i], 1), lwd=lwd, lty=2)


# Density
w2<- wsub/sum(wsub)
dd<-density(xsub, bw='nrd',na.rm =TRUE, weights =w2)
dd$y<- dd$y
if(freq) dd$y <-(dd$y*sum(HH2$counts))
if(kern & !bars) polygon(dd, col=col[i], lwd=lwd, border=epade:::a.coladd.ade(col[i], -75), angle=angle[i], density=density)
if(kern & bars ) lines(dd,  col=epade:::a.alpha.ade(col[i], 1), lwd=lwd)

}
#######################


box(col=rgb(1,1,1))
}
}
################################################################################
################################################################################





################################################################################
################################################################################
if(wall==5){
################################################################################
par(col.axis=tcol)
par(col.lab=tcol)
par(col.main=tcol)
newmai<-rep(0, 4)
oldmai<-par('mai')
if(oldmai[2]<1) newmai[2]<- 1 - oldmai[2]
if(oldmai[3]>0.75 & oldmai[3]<=0.82) newmai[3]<- 0.75-oldmai[3]
if(oldmai[4]>0.25 & oldmai[4]<=0.42) newmai[4]<- 0.25-oldmai[4]
par(mai=(oldmai+newmai))

par(font=2)

################################################################################
if(is.null(group)){

# Erzeugen Hist
HH<-weighted.hist(x[which(!is.na(x))], w[which(!is.na(x))], breaks=breaks, plot=FALSE, freq=F)
HH$density <- HH$density/diff(HH$breaks)


# Limits
if(is.null(ylim) & !freq) ylim<-c(0, max(HH$density, na.rm=TRUE)+(max(HH$density, na.rm=TRUE)/10))
if(is.null(ylim) &  freq) ylim<-c(0, max(HH$counts, na.rm=TRUE) +(max(HH$counts,  na.rm=TRUE)/10))

# Plot
plot(0,0, col=rgb(0,0,0,0), main='', xlab='', ylab='', ylim=ylim, xlim=xlim, axes=F)
if(is.null(xticks))                      a1<-axis(1, labels=TRUE, tick=TRUE, col=rgb(1,1,1), lwd.ticks=1, col.ticks=tcol)
if(!is.null(xticks) & length(xticks)==1) a1<-axis(1, at=pretty(x, n=xticks), tick=TRUE, col=rgb(1,1,1), lwd.ticks=1, col.ticks=tcol)
if(!is.null(xticks) & length(xticks)>1 ) a1<-axis(1, at=xticks, tick=TRUE, col=rgb(1,1,1), lwd.ticks=1, col.ticks=tcol)

a2<-axis(2, col=rgb(1,1,1), col.ticks=tcol, lwd.ticks=1)
abline(v=a1, h=a2, lty=1, col=rgb(1,1,1), lwd=1)
if(bars & !freq)  rect(HH$breaks[-length(HH$breaks)] , 0, HH$breaks[-1], HH$density,  angle=angle, density=density, col=col, border=epade:::a.coladd.ade(col, -75), lwd=1)
if(bars &  freq)  rect(HH$breaks[-length(HH$breaks)] , 0, HH$breaks[-1], HH$counts,   angle=angle, density=density, col=col, border=epade:::a.coladd.ade(col, -75), lwd=1)
abline(v=v, h=h, col=lcol, lty=lty, lwd=lwd)
# Outer
par(xpd=TRUE)
polygon(epade:::a.glc(side=2, line=c(4.25, 4.25, 0, 0)), epade:::a.glc(side=3, line=c(0.6, 3, 3, 0.6)), col=bgcol,        border=tcol)
polygon(epade:::a.glc(side=c(2,2,4,4), line=c(0,0,0,0)), epade:::a.glc(side=3, line=c(0.6, 3, 3, 0.6)), col=rgb(1,1,1,0), border=tcol)
polygon(epade:::a.glc(side=4, line=c(0, 0 ,0.6, 0.6)),   epade:::a.glc(side=3, line=c(0.6, 3, 3, 0.6)), col=bgcol,        border=tcol)
polygon(epade:::a.glc(side=2, line=c(4.25, 4.25 ,3.65, 3.65)),  epade:::a.glc(side=c(1,3,3,1), line=c(2.6, 0.6, 0.6, 2.6)), col=bgcol,  border=tcol)
polygon(epade:::a.glc(side=4, line=c(0, 0 ,0.6, 0.6)), epade:::a.glc(side=c(1, 3, 3, 1), line=0), col=bgcol, border=tcol)
polygon(epade:::a.glc(side=2, line=c(4.25, 4.25, 0, 0)), epade:::a.glc(side=1, line=c(2.6, 4.5, 4.5, 2.6)), col=bgcol, border=tcol)
polygon(epade:::a.glc(side=c(2, 2, 4, 4), line=0), epade:::a.glc(side=1, line=c(2.6, 4.5, 4.5, 2.6)), col=rgb(1,1,1,0), border=tcol)
polygon(epade:::a.glc(side=4, line=c(0, 0, 0.6, 0.6)), epade:::a.glc(side=1, line=c(2.6, 4.5, 4.5, 2.6)), col=bgcol, border=tcol)
text(epade:::a.glc(side=0), epade:::a.glc(side=3, line=1.5),  labels=main, cex = 1.25, font=2, col=tcol, adj=c(0.5,0))
text(epade:::a.glc(side=0), epade:::a.glc(side=1, line=3.75), labels=xlab, cex = 1.1,  font=2, col=tcol, adj=c(0.5,0))
text(epade:::a.glc(side=2, line=2.5), epade:::a.glc(side=5), labels=ylab, cex = 1.1,   font=2,  col=tcol, adj=c(0.5,0), srt=90)

par(xpd=FALSE)
box(col=tcol)
par(xpd=FALSE)

# norm Lines
xfit<-seq(min(x, na.rm=TRUE),max(x, na.rm=TRUE),length=length(HH$mids)*25)
yfit<-dnorm(xfit, mean=a.wtd.mean(x[which(!is.na(x))],w[which(!is.na(x))], na.rm=TRUE),sd=sqrt(a.wtd.var(x[which(!is.na(x))], w[which(!is.na(x))], na.rm=TRUE)))
yfit<- yfit
if(freq) yfit <-(yfit*sum(HH$counts))
if(norm) lines(xfit, yfit, col=epade:::a.alpha.ade(tcol, 1), lwd=lwd, lty=2)


# Density
w2<- w[which(!is.na(x))]/(sum(w[which(!is.na(x))]))
dd<-density(x[which(!is.na(x))], bw='nrd',na.rm =TRUE, weights =w2 )
dd$y<- (dd$y)
if(freq) dd$y <-(dd$y*sum(HH$counts))
if(kern & !bars) polygon(dd, col=col, lwd=lwd, border=epade:::a.coladd.ade(col, -75), angle=angle, density=density)
if(kern & bars ) lines(dd,   col=epade:::a.alpha.ade(tcol, 1), lwd=lwd)


}
################################################################################
################################################################################

if(!is.null(group)){


# Erzeugen Hist
HH<-weighted.hist(x[which(!is.na(x) & !is.na(g))], w[which(!is.na(x) & !is.na(g))], breaks=breaks, plot=FALSE, freq=F)
HH$density <- HH$density/diff(HH$breaks)


# Limits
if(is.null(ylim)){
for(i in 1:nlevels(g)){
HH2<-weighted.hist(x[which(g==levels(g)[i] & !is.na(x))], w[which(g==levels(g)[i] & !is.na(x))], breaks=HH$breaks, plot=FALSE, freq=F)
HH2$density <- HH2$density/diff(HH2$breaks)
if(!freq) ylim<-c(0, max(ylim, max(HH2$density, na.rm=T)+max(HH2$density, na.rm=T)/8, na.rm=T))
if(freq)  ylim<-c(0, max(ylim, max(HH2$counts,  na.rm=T)+max(HH2$counts,  na.rm=T)/8, na.rm=T))
}
}


# Plot
plot(0,0, col=rgb(0,0,0,0), main='', xlab='', ylab='', ylim=ylim, xlim=xlim, axes=F)
if(is.null(xticks))                      a1<-axis(1, labels=TRUE, tick=TRUE, col=rgb(1,1,1), lwd.ticks=1, col.ticks=tcol)
if(!is.null(xticks) & length(xticks)==1) a1<-axis(1, at=pretty(x, n=xticks), tick=TRUE, col=rgb(1,1,1), lwd.ticks=1, col.ticks=tcol)
if(!is.null(xticks) & length(xticks)>1 ) a1<-axis(1, at=xticks, tick=TRUE, col=rgb(1,1,1), lwd.ticks=1, col.ticks=tcol)
a2<-axis(2, col=rgb(1,1,1), col.ticks=tcol, lwd.ticks=1)
abline(v=v, h=h, col=lcol, lty=lty, lwd=lwd)
legend(legendon, legend=vnames, pt.cex=1.5 ,col = col2, pch=15, bg=rgb(1,1,1), border=TRUE,     box.col=tcol, box.lwd=1, text.col=tcol, text.width=max(strwidth(vnames,font = 2)))
legend(legendon, legend=vnames, pt.cex=1.5 ,col = tcol, pch=0, bg=epade:::a.alpha.ade(1, 0), border=TRUE, box.col=tcol, box.lwd=1, text.col=epade:::a.alpha.ade(1, 0), text.width=max(strwidth(vnames,font = 2)))

# Outer
par(xpd=TRUE)
polygon(epade:::a.glc(side=2, line=c(4.25, 4.25, 0, 0)), epade:::a.glc(side=3, line=c(0.6, 3, 3, 0.6)), col=bgcol,        border=tcol)
polygon(epade:::a.glc(side=c(2,2,4,4), line=c(0,0,0,0)), epade:::a.glc(side=3, line=c(0.6, 3, 3, 0.6)), col=rgb(1,1,1,0), border=tcol)
polygon(epade:::a.glc(side=4, line=c(0, 0 ,0.6, 0.6)),   epade:::a.glc(side=3, line=c(0.6, 3, 3, 0.6)), col=bgcol,        border=tcol)
polygon(epade:::a.glc(side=2, line=c(4.25, 4.25 ,3.65, 3.65)),  epade:::a.glc(side=c(1,3,3,1), line=c(2.6, 0.6, 0.6, 2.6)), col=bgcol,  border=tcol)
polygon(epade:::a.glc(side=4, line=c(0, 0 ,0.6, 0.6)), epade:::a.glc(side=c(1, 3, 3, 1), line=0), col=bgcol, border=tcol)
polygon(epade:::a.glc(side=2, line=c(4.25, 4.25, 0, 0)), epade:::a.glc(side=1, line=c(2.6, 4.5, 4.5, 2.6)), col=bgcol, border=tcol)
polygon(epade:::a.glc(side=c(2, 2, 4, 4), line=0), epade:::a.glc(side=1, line=c(2.6, 4.5, 4.5, 2.6)), col=rgb(1,1,1,0), border=tcol)
polygon(epade:::a.glc(side=4, line=c(0, 0, 0.6, 0.6)), epade:::a.glc(side=1, line=c(2.6, 4.5, 4.5, 2.6)), col=bgcol, border=tcol)
text(epade:::a.glc(side=0), epade:::a.glc(side=3, line=1.5),  labels=main, cex = 1.25, font=2, col=tcol, adj=c(0.5,0))
text(epade:::a.glc(side=0), epade:::a.glc(side=1, line=3.75), labels=xlab, cex = 1.1,  font=2, col=tcol, adj=c(0.5,0))
text(epade:::a.glc(side=2, line=2.5), epade:::a.glc(side=5), labels=ylab, cex = 1.1,   font=2,  col=tcol, adj=c(0.5,0), srt=90)

par(xpd=FALSE)
box(lwd=1, col=tcol)


# Gruppen inrun
#######################
for(i in 1:nlevels(g)){
xsub<- x[which(g==levels(g)[i] & !is.na(x))]
wsub<- w[which(g==levels(g)[i] & !is.na(x))]
b<-  max(diff(HH$breaks))


# Bars
HH2<-weighted.hist(xsub, wsub, breaks=HH$breaks, plot=FALSE, freq=F)
HH2$density <- HH2$density/diff(HH2$breaks)
if(bars & !freq) rect(HH2$breaks[-length(HH2$breaks)] , 0, HH2$breaks[-1], HH2$density,  angle=angle[i], density=density, col=col[i], border=epade:::a.alpha.ade(epade:::a.coladd.ade(col[i], -100),1), lwd=1)
if(bars & freq)  rect(HH2$breaks[-length(HH2$breaks)] , 0, HH2$breaks[-1], HH2$counts,   angle=angle[i], density=density, col=col[i], border=epade:::a.alpha.ade(epade:::a.coladd.ade(col[i], -100),1), lwd=1)



# norm Lines
xfit<-seq(min(xsub, na.rm=TRUE),max(xsub, na.rm=TRUE),length=100)
yfit<-dnorm(xfit, mean=a.wtd.mean(xsub,wsub, na.rm=TRUE),sd=sqrt(a.wtd.var(xsub,wsub, na.rm=TRUE)))
yfit<- yfit
if(freq) yfit <-(yfit*sum(HH2$counts))
if(norm) lines(xfit, yfit, col=epade:::a.alpha.ade(col[i], 1), lwd=lwd, lty=2)


# Density
w2<- wsub/sum(wsub)
dd<-density(xsub, bw='nrd',na.rm =TRUE, weights =w2)
dd$y<- dd$y
if(freq) dd$y <-(dd$y*sum(HH2$counts))
if(kern & !bars) polygon(dd, col=col[i], lwd=lwd, border=epade:::a.coladd.ade(col[i], -75), angle=angle[i], density=density)
if(kern & bars ) lines(dd,  col=epade:::a.alpha.ade(col[i], 1), lwd=lwd)

}
#######################



}
}
################################################################################
################################################################################




################################################################################
################################################################################
if(wall==6){
################################################################################
par(col.main=tcol)
par(col.axis=tcol)
par(col.lab =tcol)


################################################################################
if(is.null(group)){

# Erzeugen Hist
HH<-weighted.hist(x[which(!is.na(x))], w[which(!is.na(x))], breaks=breaks, plot=FALSE, freq=F)
HH$density <- HH$density/diff(HH$breaks)


# Limits
if(is.null(ylim) & !freq) ylim<-c(0, max(HH$density, na.rm=TRUE)+(max(HH$density, na.rm=TRUE)/10))
if(is.null(ylim) &  freq) ylim<-c(0, max(HH$counts, na.rm=TRUE) +(max(HH$counts,  na.rm=TRUE)/10))

# Plot
plot(0,0, col=rgb(0,0,0,0), main=main, xlab=xlab, ylab=ylab, ylim=ylim, xlim=xlim, axes=F)
polygon( c(par('usr')[c(1,1,2,2)]), par('usr')[c(3,4,4,3)], col=bgcol, border=FALSE)
if(is.null(xticks))                      a1<-axis(1, labels=TRUE, tick=TRUE, col=rgb(1,1,1), lwd.ticks=3, col.ticks=epade:::a.coladd.ade(bgcol, -35))
if(is.null(xticks))                      a1<-axis(1, labels=TRUE, tick=TRUE, col=rgb(1,1,1), lwd.ticks=1, col.ticks=rgb(1,1,1))
if(!is.null(xticks) & length(xticks)==1) a1<-axis(1, at=pretty(x, n=xticks), tick=TRUE, col=rgb(1,1,1), lwd.ticks=3, col.ticks=epade:::a.coladd.ade(bgcol, -35))
if(!is.null(xticks) & length(xticks)==1) a1<-axis(1, at=pretty(x, n=xticks), tick=TRUE, col=rgb(1,1,1), lwd.ticks=1, col.ticks=rgb(1,1,1))
if(!is.null(xticks) & length(xticks)>1 ) a1<-axis(1, at=xticks, tick=TRUE, col=rgb(1,1,1), lwd.ticks=3, col.ticks=epade:::a.coladd.ade(bgcol, -35))
if(!is.null(xticks) & length(xticks)>1 ) a1<-axis(1, at=xticks, tick=TRUE, col=rgb(1,1,1), lwd.ticks=1, col.ticks=rgb(1,1,1))
a2<-axis(2, col=rgb(1,1,1), col.ticks=epade:::a.coladd.ade(bgcol, -35), lwd.ticks=3)
a2<-axis(2, col=rgb(1,1,1), col.ticks=rgb(1,1,1), lwd.ticks=1)

abline(v=a1, h=a2, lty=1, col=epade:::a.coladd.ade(bgcol, -35), lwd=3)
abline(v=a1, h=a2, lty=1, col=rgb(1,1,1), lwd=1)
if(bars & !freq)  rect(HH$breaks[-length(HH$breaks)] , 0, HH$breaks[-1], HH$density,  angle=angle, density=density, col=col, border=epade:::a.coladd.ade(col, -75), lwd=1)
if(bars &  freq)  rect(HH$breaks[-length(HH$breaks)] , 0, HH$breaks[-1], HH$counts,   angle=angle, density=density, col=col, border=epade:::a.coladd.ade(col, -75), lwd=1)
abline(v=v, h=h, col=tcol, lty=lty, lwd=lwd)


# norm Lines
b<-  max(diff(HH$breaks))
xfit<-seq(min(x, na.rm=TRUE),max(x, na.rm=TRUE),length=length(HH$mids)*25)
yfit<-dnorm(xfit, mean=a.wtd.mean(x[which(!is.na(x))],w[which(!is.na(x))], na.rm=TRUE),sd=sqrt(a.wtd.var(x[which(!is.na(x))], w[which(!is.na(x))], na.rm=TRUE)))
yfit<- yfit
if(freq) yfit <-(yfit*sum(HH$counts))
if(norm) lines(xfit, yfit, col=epade:::a.alpha.ade(tcol, 1), lwd=lwd, lty=2)


# Density
w2<- w[which(!is.na(x))]/(sum(w[which(!is.na(x))]))
dd<-density(x[which(!is.na(x))], bw='nrd',na.rm =TRUE, weights =w2 )
dd$y<- (dd$y)
if(freq) dd$y <-(dd$y*sum(HH$counts))
if(kern & !bars) polygon(dd, col=col, lwd=lwd, border=epade:::a.coladd.ade(col, -75), angle=angle, density=density)
if(kern & bars ) lines(dd,   col=epade:::a.alpha.ade(tcol, 1), lwd=lwd)


box(lwd=3, col=rgb(1,1,1))
box(lwd=1, col=epade:::a.coladd.ade(bgcol, -35))
}
################################################################################
################################################################################

if(!is.null(group)){


# Erzeugen Hist
HH<-weighted.hist(x[which(!is.na(x) & !is.na(g))], w[which(!is.na(x) & !is.na(g))], breaks=breaks, plot=FALSE, freq=F)
HH$density <- HH$density/diff(HH$breaks)

# Limits
if(is.null(ylim)){
for(i in 1:nlevels(g)){
HH2<-weighted.hist(x[which(g==levels(g)[i] & !is.na(x))], w[which(g==levels(g)[i] & !is.na(x))], breaks=HH$breaks, plot=FALSE, freq=F)
HH2$density <- HH2$density/diff(HH2$breaks)
if(!freq) ylim<-c(0, max(ylim, max(HH2$density, na.rm=T)+max(HH2$density, na.rm=T)/8, na.rm=T))
if(freq)  ylim<-c(0, max(ylim, max(HH2$counts,  na.rm=T)+max(HH2$counts,  na.rm=T)/8, na.rm=T))
}
}


# Plot
plot(0,0, col=rgb(0,0,0,0), main=main, xlab=xlab, ylab=ylab, ylim=ylim, xlim=xlim, axes=F)
polygon( c(par('usr')[c(1,1,2,2)]), par('usr')[c(3,4,4,3)], col=bgcol, border=FALSE)
if(is.null(xticks))                      a1<-axis(1, labels=TRUE, tick=TRUE, col=rgb(1,1,1), lwd.ticks=3, col.ticks=epade:::a.coladd.ade(bgcol, -35))
if(is.null(xticks))                      a1<-axis(1, labels=TRUE, tick=TRUE, col=rgb(1,1,1), lwd.ticks=1, col.ticks=rgb(1,1,1))
if(!is.null(xticks) & length(xticks)==1) a1<-axis(1, at=pretty(x, n=xticks), tick=TRUE, col=rgb(1,1,1), lwd.ticks=3, col.ticks=epade:::a.coladd.ade(bgcol, -35))
if(!is.null(xticks) & length(xticks)==1) a1<-axis(1, at=pretty(x, n=xticks), tick=TRUE, col=rgb(1,1,1), lwd.ticks=1, col.ticks=rgb(1,1,1))
if(!is.null(xticks) & length(xticks)>1 ) a1<-axis(1, at=xticks, tick=TRUE, col=rgb(1,1,1), lwd.ticks=3, col.ticks=epade:::a.coladd.ade(bgcol, -35))
if(!is.null(xticks) & length(xticks)>1 ) a1<-axis(1, at=xticks, tick=TRUE, col=rgb(1,1,1), lwd.ticks=1, col.ticks=rgb(1,1,1))
a2<-axis(2, col=rgb(1,1,1), col.ticks=epade:::a.coladd.ade(bgcol, -35), lwd.ticks=3)
a2<-axis(2, col=rgb(1,1,1), col.ticks=rgb(1,1,1), lwd.ticks=1)

abline(v=a1, h=a2, lty=1, col=epade:::a.coladd.ade(bgcol, -35), lwd=3)
abline(v=a1, h=a2, lty=1, col=rgb(1,1,1), lwd=1)

abline(v=v, h=h, col=tcol, lty=lty, lwd=lwd)
legend(legendon, legend=vnames, pt.cex=1.5 ,col = col2, pch=15, bg=bgcol, border=epade:::a.coladd.ade(bgcol, -35), box.col=rgb(1,1,1), box.lwd=3, text.col=tcol, text.width=max(strwidth(vnames,font = 2)))
legend(legendon, legend=vnames, pt.cex=1.5 ,col = rgb(1,1,1), pch=0, bg=epade:::a.alpha.ade(1, 0), border=epade:::a.coladd.ade(bgcol, -35), box.col=epade:::a.coladd.ade(bgcol, -35), box.lwd=1, text.col=epade:::a.alpha.ade(1, 0), text.width=max(strwidth(vnames,font = 2)))


# Gruppen inrun
#######################
for(i in 1:nlevels(g)){
xsub<- x[which(g==levels(g)[i] & !is.na(x))]
wsub<- w[which(g==levels(g)[i] & !is.na(x))]
b<-  max(diff(HH$breaks))


# Bars
HH2<-weighted.hist(xsub, wsub, breaks=HH$breaks, plot=FALSE, freq=F)
HH2$density <- HH2$density/diff(HH2$breaks)

if(bars & !freq) rect(HH2$breaks[-length(HH2$breaks)] , 0, HH2$breaks[-1], HH2$density,  angle=angle[i], density=density, col=col[i], border=epade:::a.alpha.ade(epade:::a.coladd.ade(col[i], -100),1), lwd=1)
if(bars & freq)  rect(HH2$breaks[-length(HH2$breaks)] , 0, HH2$breaks[-1], HH2$counts,   angle=angle[i], density=density, col=col[i], border=epade:::a.alpha.ade(epade:::a.coladd.ade(col[i], -100),1), lwd=1)



# norm Lines
xfit<-seq(min(xsub, na.rm=TRUE),max(xsub, na.rm=TRUE),length=100)
yfit<-dnorm(xfit, mean=a.wtd.mean(xsub,wsub, na.rm=TRUE),sd=sqrt(a.wtd.var(xsub,wsub, na.rm=TRUE)))
yfit<- yfit
if(freq) yfit <-(yfit*sum(HH2$counts))
if(norm) lines(xfit, yfit, col=epade:::a.alpha.ade(col[i], 1), lwd=lwd, lty=2)


# Density
w2<- wsub/sum(wsub)
dd<-density(xsub, bw='nrd',na.rm =TRUE, weights =w2)
dd$y<- dd$y
if(freq) dd$y <-(dd$y*sum(HH2$counts))
if(kern & !bars) polygon(dd, col=col[i], lwd=lwd, border=epade:::a.coladd.ade(col[i], -75), angle=angle[i], density=density)
if(kern & bars ) lines(dd,  col=epade:::a.alpha.ade(col[i], 1), lwd=lwd)

}
#######################


box(lwd=3, col=rgb(1,1,1))
box(lwd=1, col=epade:::a.coladd.ade(bgcol, -35))
}
}
################################################################################
################################################################################



}
