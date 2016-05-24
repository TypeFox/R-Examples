### R code from vignette source 'PTM_MarkerFinder.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: PTM_MarkerFinder.Rnw:39-42
###################################################
library(protViz)
data(HexNAc)
str(HexNAc[[1]], nchar.max=30)


###################################################
### code chunk number 2: PTM_MarkerFinder.Rnw:62-64
###################################################
HexNAc_MarkerIons <- c(126.05495, 138.05495, 144.06552, 168.06552, 
    186.07608, 204.08665)


###################################################
### code chunk number 3: PTM_MarkerFinder.Rnw:70-82
###################################################
ptm.0<-cbind(AA="-",
	mono=0.0, avg=0.0, desc="unmodified", unimodAccID=NA)

ptm.1<-cbind(AA='N', 
	mono=317.122300, avg=NA, desc="HexNAc",
        unimodAccID=2)
             
ptm.2<-cbind(AA='M',
	mono=147.035400, avg=NA, desc="Oxidation",
        unimodAccID=1)
     
m<-as.data.frame(rbind(ptm.0, ptm.1, ptm.2))


###################################################
### code chunk number 4: PTM_MarkerFinder.Rnw:88-95
###################################################
s <- PTM_MarkerFinder(data=HexNAc, modification=m$mono, 
	modificationName=m$desc,
        minMarkerIntensityRatio=3,
        itol_ppm=20,
        mZmarkerIons=HexNAc_MarkerIons) 

s


###################################################
### code chunk number 5: PTM_MarkerFinder.Rnw:99-100 (eval = FALSE)
###################################################
## demo(PTM_MarkerFinder) 


###################################################
### code chunk number 6: PTM_MarkerFinder.Rnw:105-126
###################################################
op<-par(mfrow=c(2,2), mar=c(4,4,4,1)); 
dump <- lapply(split(s,s$query), 
    function(x){ plot(x$mZ, x$markerIonIntensity, 
        type='h',
        col='lightblue',
        cex=2,
        ylab='intensity', xlab='m/z',
        xlim=range(c(HexNAc_MarkerIons,  
            max(HexNAc_MarkerIons) 
                + 0.1 * (max(HexNAc_MarkerIons) - min(HexNAc_MarkerIons)), 
            min(HexNAc_MarkerIons) 
                - 0.1 * (max(HexNAc_MarkerIons) - min(HexNAc_MarkerIons)))),
            ylim=range(s$markerIonIntensity),
            log='y',
            main=paste("scan=", unique(x$scans),
                "/query=", unique(x$query), sep='')); 
            text(x$mZ, x$markerIonIntensity, 
                round(x$mZ,2),col='red',cex=0.7)
        }
    )
par(op)


###################################################
### code chunk number 7: PTM_MarkerFinder.Rnw:131-139
###################################################
d<-list(); d[[1]]<-HexNAc[[3]]; d[[2]]<-HexNAc[[4]];
d[[3]]<-HexNAc[[5]]
ss<-PTM_MarkerFinder(data=d, modification=m$mono, 
	modificationName=m$desc,
        minMarkerIntensityRatio=3,
        itol_ppm=20,
        mZmarkerIons=HexNAc_MarkerIons) 



###################################################
### code chunk number 8: PTM_MarkerFinder.Rnw:144-147
###################################################
w<-reshape(s[,c(1,7,3,4)], direction='wide', 
    timevar="markerIonMZ", idvar=c('scans','query'))
w


###################################################
### code chunk number 9: PTM_MarkerFinder.Rnw:151-154
###################################################
write.table(w, file="HexNAc_PTM_markerFinder.csv",
    sep=',', row.names=FALSE,col.names=TRUE, quote=FALSE)



