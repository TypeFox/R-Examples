#R

fetuin<-c('MK', 'SFVLLFCLAQLWGCHSIPLDPVAGYK',
'EPACDDPDTEQAALAAVDYINK',
'HLPR', 'GYK', 'HTLNQIDSVK', 'VWPR',
'RPTGEVYDIEIDTLETTCHVLDPTPLANCSVR',
'QQTQHAVEGDCDIHVLK', 'QDGQFSVLFTK',
'CDSSPDSAEDVR', 'K', 'LCPDCPLLAPLNDSR',
'VVHAVEVALATFNAESNGSYLQLVEISR',
'AQFVPLPVSVSVEFAVAATDCIAK',
'EVVDPTK', 'CNLLAEK', 'QYGFCK',
'GSVIQK', 'ALGGEDVR',
'VTCTLFQTQPVIPQPQPDGAEAEAPSAVPDAAGPTPSAAGPPVASVVVGPSVVAVPLPLHR',
'AHYDLR', 'HTFSGVASVESSSGEAFHVGK',
'TPIVGQPSIPGGPVR', 'LCPGR', 'IR', 'YFK', 'I')



plot(pim<-parentIonMass(fetuin), ylab="parent ion mass", main="FETUA_BOVIN")
hist(pim, main="Histogram of parent ion mass of FETUA_BOVIN", xlab="m/Z")

fi<-fragmentIon(fetuin)

peptides<-c('HTLNQIDSVK', 'ALGGEDVR', 'TPIVGQPSIPGGPVR')

pim<-parentIonMass(peptides)
fi<-fragmentIon(peptides)

par(mfrow=c(3,1)); 
for (i in 1:length(peptides)){
    plot(0,0,
        xlab='m/Z',
        ylab='',
        xlim=range(c(fi[i][[1]]$b,fi[i][[1]]$y)),
        ylim=c(0,1),
        type='n',
        axes=FALSE,
        sub=paste("fragment ions of", peptides[i], "/ parent ion mass =", pim[i], "Da"));
    box()
    axis(1,fi[i][[1]]$b,round(fi[i][[1]]$b,2))
    pepSeq<-strsplit(peptides[i],"")
    axis(3,fi[i][[1]]$b,pepSeq[[1]])

    abline(v=fi[i][[1]]$b, col='red',lwd=2) 
    abline(v=fi[i][[1]]$c, col='orange') 
    abline(v=fi[i][[1]]$y, col='blue',lwd=2)
    abline(v=fi[i][[1]]$z, col='cyan')
}

spec<-list(scans=1138,
    title="178: (rt=22.3807) [20080816_23_fetuin_160.RAW]",
    rtinseconds=1342.8402,
    charge=2,
    mZ=c(195.139940, 221.211970, 239.251780, 290.221750, 
    316.300770, 333.300050, 352.258420, 448.384360, 466.348830, 
    496.207570, 509.565910, 538.458310, 547.253380, 556.173940, 
    560.358050, 569.122080, 594.435500, 689.536940, 707.624790, 
    803.509240, 804.528220, 822.528020, 891.631250, 909.544400, 
    916.631600, 973.702160, 990.594520, 999.430580, 1008.583600, 
    1017.692500, 1027.605900),
    intensity=c(931.8, 322.5, 5045, 733.9, 588.8, 9186, 604.6,
    1593, 531.8, 520.4, 976.4, 410.5, 2756, 2279, 5819, 2.679e+05,
    1267, 1542, 979.2, 9577, 3283, 9441, 1520, 1310, 1.8e+04,
    587.5, 2685, 671.7, 3734, 8266, 3309)
)


par(mfrow=c(2,3));
m1<-psm('HTLNQIDSVK', spec,plot=TRUE)
hist(m1$mZ.Da.error,20)
abline(v=c(-5,5),col='grey')
hist(m1$mZ.ppm.error,20)
abline(v=c(-1000,1000),col='grey')

m2<-psm('ENNTHLLVSK', spec,plot=TRUE)
hist(m2$mZ.Da.error,20)
abline(v=c(-5,5),col='grey')
hist(m2$mZ.ppm.error,20)
abline(v=c(-1000,1000),col='grey')



### PEAKPLOT

        peptides<-c('HTLNQIDSVK', 'ALGGEDVR', 'TPIVGQPSIPGGPVR')

        pim<-parentIonMass(peptides)
        fi<-fragmentIon(peptides)

        op<-par(mfrow=c(3,1)); 
        for (i in 1:length(peptides)){
            plot(0,0,
                xlab='m/Z',
                ylab='',
                xlim=range(c(fi[i][[1]]$b,fi[i][[1]]$y)),
                ylim=c(0,1),
                type='n',
                axes=FALSE,
                sub=paste(peptides[i], "/", pim[i], "Da"));
            box()
            axis(1,fi[i][[1]]$b,round(fi[i][[1]]$b,1),las=2)
            axis(1,fi[i][[1]]$y,round(fi[i][[1]]$y,1),las=2)

            pepSeq<-strsplit(peptides[i],"")
            axis(3,fi[[i]]$b, paste("b", row.names(fi[[i]]),sep=''),las=2)
            axis(3,fi[[i]]$y, paste("y", row.names(fi[[i]]),sep=''),las=2)

            text(fi[i][[1]]$b, rep(0.3, nchar(peptides[i])), pepSeq[[1]],pos=3,cex=4, lwd=4, col="#aaaaaaaa")

            abline(v=fi[[i]]$b, col='red') 
            abline(v=fi[[i]]$y, col='blue',lwd=2)
        }
        par(op)


ptm.0<-cbind(AA="-",
    mono=0.0, avg=0.0, desc="unmodified", unimodAccID=NA)

ptm.616<-cbind(AA='S',
    mono=-27.010899, avg=NA, desc="Substituition",
    unimodAccID=616)

ptm.651<-cbind(AA='N',
    mono=27.010899, avg=NA, desc="Substituition",
    unimodAccID=651)


m<-as.data.frame(rbind(ptm.0, ptm.616, ptm.651))

genMod(c('TAFDEAIAELDTLNEESYK','TAFDEAIAELDTLSEESYK'), m$AA)

fi<-fragmentIon(c('TAFDEAIAELDTLSEESYK', 
    'TAFDEAIAELDTLNEESYK', 'TAFDEAIAELDTLSEESYK', 
    'TAFDEAIAELDTLNEESYK'),
        modified=c('0000000000000200000', 
        '0000000000000100000', '0000000000000000000', 
        '0000000000000000000'),
    modification=m$mono)

data(msms)
op<-par(mfrow=c(2,1), mar=c(4,4,4,4),cex=0.5)
peakplot('TAFDEAIAELDTLSEESYK',msms[[1]], fi=fi[[1]])
peakplot('TAFDEAIAELDTLNEESYK',msms[[1]], fi=fi[[4]])
peakplot('TAFDEAIAELDTLNEESYK',msms[[2]], fi=fi[[2]])
peakplot('TAFDEAIAELDTLSEESYK',msms[[2]], fi=fi[[3]])

