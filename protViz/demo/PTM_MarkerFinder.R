#R

# 2012-11-29
# Christian Panse <cp@fgcz.ethz.ch>
# Paolo Nanni

data(HexNAc)
#-- setting variable modification
#        ->HexNAc (N)
#                ->N 317.122300
#        ->Oxidation (M)
#                ->M 147.035400
#-- setting fixed modification
#        ->Carbamidomethyl (C)
#                ->C 160.030649


ptm.0<-cbind(AA="-",
    mono=0.0, avg=0.0, desc="unmodified", unimodAccID=NA)

ptm.1<-cbind(AA='N',
    mono=317.122300, avg=NA, desc="HexNAc",
        unimodAccID=2)

ptm.2<-cbind(AA='M',
    mono=147.035400, avg=NA, desc="Oxidation",
        unimodAccID=1)

m<-as.data.frame(rbind(ptm.0, ptm.1, ptm.2))

HexNAc_MarkerIons<-c(126.05495, 138.05495, 144.06552, 168.06552, 186.07608, 204.08665)

s <- PTM_MarkerFinder(data=HexNAc, modification=m$mono, 
    modificationName=m$desc,
    minMarkerIntensityRatio=3,
    itol_ppm=20,
    mZmarkerIons=HexNAc_MarkerIons) 


#w<-reshape(s[,c(1,7,3,4)], direction='wide', 
#    timevar="markerIonMZ", idvar=c('scans','query'))
#
#write.table(w, file="F175813_PTM_markerFinder.csv",
#    sep=',', row.names=FALSE,col.names=TRUE, quote=FALSE)

    # some statistics
    boxplot(markerIonPpmError~as.factor(scans),
        data=s,
        main='Overview plot: boxplot of marker ion errors (ppm) for each pPTM spectrum',
        xlab='scan number', ylab='ppm error',
        ylim=c(-5,5))
       
    abline(h=0.0,col='grey')

    boxplot(markerIonPpmError~s$markerIonMZ,
        data=s,
        main='Overview plot: boxplot of marker ion errors (ppm) for each pPTM spectrum',
        xlab='m/z marker ions', ylab='ppm error')
       
    abline(h=0.0,col='grey')

    boxplot(s$markerIonIntensity ~ s$markerIonMZ,
        log='y',
        main='Summary plot: boxplot of marker ion intensities from all pPTM spectra',
        xlab='markerIon m/z', 
        ylab='log10 based marker ion intensity')



op<-par(mfrow=c(2,2), mar=c(4,4,4,4)); 
dump <- lapply(split(s,s$query), 
    function(x){ plot(x$mZ, x$markerIonIntensity, 
        type='h',
        col='lightblue',
        cex=2,
        ylab='intensity', xlab='m/z',
        xlim=range(c(HexNAc_MarkerIons,  
            max(HexNAc_MarkerIons) + 0.1*(max(HexNAc_MarkerIons)-min(HexNAc_MarkerIons)), 
            min(HexNAc_MarkerIons)-0.1*(max(HexNAc_MarkerIons)-min(HexNAc_MarkerIons)))),
            ylim=range(s$markerIonIntensity),
            log='y',
            main=paste("scan=",unique(x$scans),"/query=", unique(x$query), sep='')); 
            text(x$mZ, x$markerIonIntensity,round(x$mZ,1),col='red',cex=0.7)
        }
    )
par(op)

op<-par(mfrow=c(2,2), mar=c(4,4,4,4)); 
dump <- lapply(split(s,s$query), 
    function(x){ pie(x$markerIonIntensity, 
        cex=0.7,
        round(x$markerIonMZ,2), 
        col='white',
        main=paste("scan=",unique(x$scans),"/query=", unique(x$query), sep=''))})
par(op)
