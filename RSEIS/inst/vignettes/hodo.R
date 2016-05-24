### R code from vignette source 'hodo.Rnw'

###################################################
### code chunk number 1: hodo.Rnw:66-67
###################################################
library(RSEIS)


###################################################
### code chunk number 2: hodo.Rnw:74-78
###################################################
data(GH)

swig(GH, SHOWONLY=TRUE)



###################################################
### code chunk number 3: hodo.Rnw:84-85
###################################################
print(GH$pickfile$STAS$name)


###################################################
### code chunk number 4: hodo.Rnw:91-101
###################################################

thesta = "CE1"
iwv  = which(GH$STNS==thesta & GH$COMPS=="V")
iwn  = which(GH$STNS==thesta & GH$COMPS=="N")
iwe  = which(GH$STNS==thesta & GH$COMPS=="E")


data = cbind(GH$JSTR[[iwv]], GH$JSTR[[iwn]], GH$JSTR[[iwe]])




###################################################
### code chunk number 5: hodo.Rnw:109-118
###################################################
ipphase = which(GH$pickfile$STAS$name==thesta & GH$pickfile$STAS$phase=="P" )
isphase = which(GH$pickfile$STAS$name==thesta & GH$pickfile$STAS$phase=="S" )


lat=GH$pickfile$STAS$lat[ipphase]
lon = GH$pickfile$STAS$lon[ipphase]
DAZ = rdistaz(lat, lon,   GH$pickfile$LOC$lat,GH$pickfile$LOC$lon )
 rbaz = grotseis(DAZ$baz, flip=FALSE)



###################################################
### code chunk number 6: hodo.Rnw:124-134
###################################################

plot( c(GH$pickfile$STAS$lon, GH$pickfile$LOC$lon) , c(GH$pickfile$STAS$lat, GH$pickfile$LOC$lat), type='n', xlab="LON", ylab="LAT")
points(GH$pickfile$STAS$lon, GH$pickfile$STAS$lat, pch=6)
points(GH$pickfile$LOC$lon, GH$pickfile$LOC$lat, pch=8)
text(GH$pickfile$STAS$lon, GH$pickfile$STAS$lat,  GH$pickfile$STAS$name, pos=3)

##  plot( c(lon, GH$pickfile$LOC$lon) , c(GH$pickfile$STAS$lat, GH$pickfile$LOC$lat))
##  text(lon, lat, labels="station")




###################################################
### code chunk number 7: hodo.Rnw:144-164
###################################################

x = DAZ$dist*sin(DAZ$baz*pi/180)
y = DAZ$dist*cos(DAZ$baz*pi/180)

plot(c(0,1.3*x), c(0,1.3*y), type='n', asp=1, xlab="E-W, km", ylab="N-S, km")
points(c(0,x), c(0,y), pch=c(3,6))
text(x,y, labels="station", pos=1)
text(x,y, labels=GH$pickfile$STAS$name[ipphase], pos=2)


text(0,0, labels="source", pos=3)

vecs  = rbind(c(0,0,1), c(0,1,0))
bvec  = vecs %*%  rbaz

bvec = .1*DAZ$dist*bvec
arrows(x,y, x+bvec[,2], y+bvec[,3], col=c("red", "blue"))

 



###################################################
### code chunk number 8: hodo.Rnw:170-181
###################################################

## data = cbind(GH$JSTR[[iwv]], GH$JSTR[[iwn]], GH$JSTR[[iwe]])

  vnelabs=c("Vertical", "North", "East")
  
rotlabs=c("Vertical", "Radial(away)", "Transvers(right)")

 xt=seq(from=0, by=GH$dt[iwv], length=length(GH$JSTR[[iwv]]))

 PLOT.MATN(data, tim=xt, dt=GH$dt[iwv], notes=vnelabs)



###################################################
### code chunk number 9: hodo.Rnw:189-195
###################################################

## data = cbind(GH$JSTR[[iwv]], GH$JSTR[[iwn]], GH$JSTR[[iwe]])
btemp  = data  %*%  rbaz

 PLOT.MATN(btemp, tim=xt, dt=GH$dt[iwv], notes=rotlabs)



###################################################
### code chunk number 10: hodo.Rnw:205-234
###################################################


i1 = match(GH$STNS[iwv],  GH$pickfile$STAS$name)


reft = list(jd=GH$info$jd[iwv], hr=GH$info$hr[iwv], mi=GH$info$mi[iwv], sec=GH$info$sec[iwv] )


ptim = list(jd=GH$pickfile$LOC$jd, hr=GH$pickfile$LOC$hr, mi=GH$pickfile$LOC$mi, sec=GH$pickfile$STAS$sec[ipphase] )
stim = list(jd=GH$pickfile$LOC$jd, hr=GH$pickfile$LOC$hr, mi=GH$pickfile$LOC$mi, sec=GH$pickfile$STAS$sec[isphase] )

t1 = secdifL( reft, ptim)
t2 = secdifL( reft, stim)

 PLOT.MATN(btemp, WIN=c(5,8) , tim=xt, dt=GH$dt[iwv], notes=rotlabs)
abline(v=t1, col='red', lty=2)
abline(v=t2, col='blue', lty=2)

mtext(side=3, at=t1, line=.1, text="Pwave", col='red')
mtext(side=3, at=t2, line=.1, text="Swave", col='blue')

pwin = c(t1-.02, t1+.09)
swin = c(t2-.01, t2+.1)

abline(v=pwin, col='red', lty=2)
abline(v=swin, col='blue', lty=2)





###################################################
### code chunk number 11: hodo.Rnw:240-247
###################################################
rbow=rainbow(140)[1:100]

atemp = btemp[xt>pwin[1]&xt<pwin[2]  ,]
##  PLOT.MATN(atemp,  tim=xt[xt>pwin[1]&xt<pwin[2]], dt=GH$dt[iwv], notes=rotlabs)

 hodogram(atemp, dt=GH$dt[iwv]  ,labs=rotlabs,  STAMP=thesta,  COL=rbow )



###################################################
### code chunk number 12: hodo.Rnw:252-258
###################################################

atemp = btemp[xt>swin[1]&xt<swin[2]  ,]
###  PLOT.MATN(atemp,  tim=xt[xt>swin[1]&xt<swin[2]], dt=GH$dt[iwv], notes=rotlabs)

 hodogram(atemp, dt=GH$dt[iwv]  ,labs=rotlabs,  STAMP=thesta,  COL=rbow )



