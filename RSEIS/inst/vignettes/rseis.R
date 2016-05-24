### R code from vignette source 'rseis.Rnw'

###################################################
### code chunk number 1: rseis.Rnw:113-119
###################################################

library(RSEIS)
library(RPMG)

options(continue = " ")



###################################################
### code chunk number 2: rseis.Rnw:128-131
###################################################
data(KH)
names(KH)



###################################################
### code chunk number 3: rseis.Rnw:138-147
###################################################

pdark=c("tomato3","royalblue","forestgreen","blueviolet","tan3",
"lightseagreen","deeppink","cyan3",
"bisque3","magenta1","lightsalmon3","darkcyan", "gold4" ,      
"darkorange1", "goldenrod4" , "seagreen3" ,  "lightpink3")

palette(pdark )




###################################################
### code chunk number 4: rseis.Rnw:150-158
###################################################


STDLAB = c("DONE",  "zoom in", "zoom out", "refresh", "restore",
 "XTR", "SPEC", "SGRAM" ,"WLET", "FILT",  "Pinfo")


swig(KH, SHOWONLY=FALSE, WIN = c(90, 200) , STDLAB=STDLAB)



###################################################
### code chunk number 5: rseis.Rnw:260-265
###################################################

data(GH)
numstas = length(GH$STNS)




###################################################
### code chunk number 6: rseis.Rnw:276-287
###################################################


verts = which(GH$COMPS == "V")

STDLAB = c("DONE", "QUIT",  "NEXT","PREV",  "zoom in", "zoom out", "refresh", "restore", "SavePF", 
"PickWin", "XTR", "SPEC", "SGRAM" ,"WLET", "FILT", "Pinfo", "WINFO", "PTS", "YPIX", "WPIX")


swig(GH,  sel=verts[c(1,2,3,4,6) ] ,STDLAB =STDLAB, WIN=c(4,15),  SHOWONLY=TRUE)




###################################################
### code chunk number 7: rseis.Rnw:294-299
###################################################

vertord = getvertsorder(GH$pickfile, GH)

swig(GH,  sel=vertord$sel[1:5] ,      STDLAB =STDLAB,  WIN=c(4, 15), SHOWONLY=FALSE)



###################################################
### code chunk number 8: rseis.Rnw:315-316
###################################################
names(GH$pickfile)


###################################################
### code chunk number 9: rseis.Rnw:320-321
###################################################
names(GH$pickfile$STAS)


###################################################
### code chunk number 10: rseis.Rnw:330-331
###################################################
data.frame(cbind(name=GH$pickfile$STAS$name, comp=GH$pickfile$STAS$comp, phase=GH$pickfile$STAS$phase, time=GH$pickfile$STAS$sec, lat=GH$pickfile$STAS$lat, lon=GH$pickfile$STAS$lon))


###################################################
### code chunk number 11: rseis.Rnw:334-335
###################################################
names(GH$pickfile$LOC)


###################################################
### code chunk number 12: rseis.Rnw:342-347
###################################################

apx = uwpfile2ypx(GH$pickfile)
swig(GH,  sel=vertord$sel[1:5] , APIX=apx, STDLAB =STDLAB,  WIN=c(4, 15), SHOWONLY=FALSE, velfile=VELMOD1D, )




###################################################
### code chunk number 13: rseis.Rnw:354-356
###################################################
PICK.DOC('WLET')



###################################################
### code chunk number 14: rseis.Rnw:360-362
###################################################
PICK.DOC()



###################################################
### code chunk number 15: rseis.Rnw:372-385
###################################################
data(sunspots)

AA = attributes(sunspots)
starttime=list(yr=AA$tsp[1], jd=1,mo=1,dom=1,hr=0,mi=0,sec=0)
ES = prep1wig(wig=sunspots, dt=1/12, sta="STA", comp="CMP", units="UNITS", starttime=starttime    )

EH=prepSEIS(ES)

STDLAB = c("DONE",  "zoom in", "zoom out", "refresh", "restore",
 "XTR", "SPEC", "SGRAM" ,"WLET", "FILT",  "Pinfo")

xx =  swig( EH, STDLAB = STDLAB,  SHOWONLY=TRUE)



###################################################
### code chunk number 16: rseis.Rnw:389-419
###################################################
a = list(y=EH$JSTR[[1]], dt=EH$dt[1])
 Mspec = mtapspec(a$y, a$dt, klen =1024 , MTP = list(kind = 1,
        nwin = 5, npi = 3, inorm = 0))
    f = Mspec$freq
    amp = Mspec$spec[1:length(f)]
    ma = amp
displ = ma
f1 = 0.01
f2 = 1/(2*EH$dt[1])

    flag = f >= f1 & f <= f2
plxy = "xy"

 plot(range(f[flag]),range(displ[flag]),type='n',log=plxy,axes=FALSE, xlab="Hz", ylab="Spec")
      lines(f[flag], displ[flag], col=1, lty=1)       
      axis(2, las=2)
      axis(1)
      box()
  ##     L = locator()
 ## rDUMPLOC(L)
  
L=list()
L$x=c(0.0939105803154482,0.183351679178275,0.350333785192083)
L$y=c(31166.8951116052,7536.36156151748,2370.29538151056)

abline(v=L$x, lty=2, col=grey(0.8) )
text(L$x, rep(max(range(displ[flag])), length(L$x)),  labels=round((1/L$x)), xpd=TRUE, srt=45, adj=c(0,0) )





###################################################
### code chunk number 17: rseis.Rnw:433-436
###################################################
data(OH)
xx =  swig( OH, sel=which(OH$COMPS == "V"), STDLAB = STDLAB,  SHOWONLY=TRUE)



###################################################
### code chunk number 18: rseis.Rnw:439-472
###################################################
a = list(y=OH$JSTR[[1]], dt=OH$dt[1])
 Mspec = mtapspec(a$y, a$dt, klen =1024 , MTP = list(kind = 1,
        nwin = 5, npi = 3, inorm = 0))
    f = Mspec$freq
    amp = Mspec$spec[1:length(f)]
    ma = amp
displ = ma
f1 = 0.01
f2 = 1/(2*EH$dt[1])

    flag = f >= f1 & f <= f2
plxy = "xy"

 plot(range(f[flag]),range(displ[flag]),type='n',log=plxy,axes=FALSE, xlab="Hz", ylab="Spec")
      lines(f[flag], displ[flag], col=1, lty=1)       
      axis(2, las=2)
      axis(1)
      box()
      
##   L = locator()
##   rDUMPLOC(L)

u = par("usr")
L=list()
L$x=c(0.242745864538716,0.423063797271721,0.447439169221329,0.529320129815409,0.100457805181183)
L$y=c(37.4392793655756,19.0231348719557,15.6564640841027,15.0862607475168,40.6986171075992)

abline(v=L$x, lty=2, col=grey(0.8) )
text(L$x, rep(max(range(displ[flag])), length(L$x)),  labels=round(10000*(1/L$x)), xpd=TRUE, srt=45, adj=c(0,0) )






###################################################
### code chunk number 19: rseis.Rnw:481-499
###################################################
 data(KH)
     dt = KH$dt[1]


     y =  KH$JSTR[[1]]
     y = y[1:50000]
y = y-mean(y)
     x =  seq(from=0, by=dt, length=length(y))




     fl=rep(1/100, 5)
     fh=1/c(4,3,2,1, .5)

     FILT.spread(x, y, dt, fl = fl, fh = fh, sfact = 1, WIN = NULL, PLOT = TRUE, TIT = NULL, TAPER = 0.05, POSTTAPER = 0.1)




