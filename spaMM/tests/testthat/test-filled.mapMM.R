cat("\ntest filled.mapMM:")
# filled.mapMM

data(blackcap)
bfit <- corrHLfit(migStatus ~ means+ Matern(1|longitude+latitude),data=blackcap,
                  HLmethod="ML",
                  ranFix=list(lambda=0.5537,phi=1.376e-05,rho=0.0544740,nu=0.6286311))

## showing add.map
if (require(maps)) { ## required for add.map=TRUE 
  filled.mapMM(bfit,add.map=TRUE,plot.axes={axis(1);axis(2)},
             plot.title=title(main="Inferred migration propensity of blackcaps",
                              xlab="longitude",ylab="latitude"))
}
## filled.mapMM takes a bit longer
filled.mapMM(bfit,nlevels=30,plot.axes={axis(1);axis(2)},
             plot.title=title(main="Inferred migration propensity of blackcaps",
                              xlab="longitude",ylab="latitude"))

data(Loaloa)  
lfit <- corrHLfit(cbind(npos,ntot-npos)~elev1+elev2+elev3+elev4+maxNDVI1+seNDVI
                  +Matern(1|longitude+latitude),HLmethod="HL(0,1)",data=Loaloa,
                  family=binomial(),ranFix=list(nu=0.5,rho=2.255197,lambda=1.075))   

## longer computation requiring interpolation of 197 points 
if (require(maps)) { ## required for add.map=TRUE 
  filled.mapMM(lfit,add.map=TRUE,plot.axes={axis(1);axis(2)},
             decorations=quote(points(pred[,coordinates],pch=15,cex=0.3)),
             plot.title=title(main="Inferred prevalence, North Cameroon",
                              xlab="longitude",ylab="latitude"))
}
# test of syntax, no expect_ yet

## local maximum issue
# obj <- function(v) {logLik(corrHLfit(formula=migStatus ~ 1 + Matern(1|latitude+longitude),
#                                      HLmethod='ML',data=blackcap,ranFix=as.list(v)))}
# mygrid <- expand.grid(rho=seq(0.2,0.5,0.025),nu=seq(16))
# apply(mygrid,1L,obj) -> blu
# spaMM.filled.contour(x=seq(0.2,0.5,0.025),y=seq(16),z=matrix(blu,ncol=16))

