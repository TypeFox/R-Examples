#R

# $HeadURL: http://fgcz-svn.unizh.ch/repos/fgcz/testing/proteomics/R/protViz/R/PTM_MarkerFinder.R $
# $Id: PTM_MarkerFinder.R 6318 2014-04-02 13:02:51Z cpanse $
# $Date: 2014-04-02 15:02:51 +0200 (Wed, 02 Apr 2014) $


.PTM_MarkerFinder_writeMGF_Header <- function(mgfFilename){
    FILE <- file(mgfFilename, "a")
    writeLines(as.character(paste("# protViz packageVersion ", packageVersion('protViz'))), FILE)
    writeLines(paste("# ",date(),sep=''), FILE)
    close(FILE)
}

.PTM_MarkerFinder_writeMGF <- function(spec, mgfFilename, pattern=-1){
    FILE <- file(mgfFilename, "a")
    if (length(pattern)>0){
        writeLines(as.character(paste("# PATTERN ", pattern, sep=' ')), FILE)
    }
    writeLines("BEGIN IONS", FILE)
    writeLines(paste("TITLE=",as.character(spec$title),sep=''), FILE)
    writeLines(paste("PEPMASS=",as.character(spec$pepmass),sep=''), FILE)
    writeLines(paste("CHARGE=",as.character(spec$charge),'+',sep=''), FILE)
    writeLines(paste("SCANS=",as.character(spec$scans),sep=''), FILE)
    writeLines(paste("RTINSECONDS=",as.character(spec$rtinseconds),sep=''), FILE)
    writeLines(as.character(paste(spec$mZ, spec$intensity, sep=' ')), FILE)
    writeLines("END IONS", FILE)                              
    writeLines("", FILE)                              
    close(FILE)
}

.PTM_MarkerFinder_ <- function(peptideSequence, modification, modificationName){

    result <- as.character()

    n <- nchar(peptideSequence)
    seq <- strsplit(peptideSequence,'')[[1]]
    mod <- as.numeric((strsplit(modification,'')[[1]]))

    for (i in 1:n){
        result<-paste(result, seq[i], sep='')
        if (mod[i+1] > 0){
            result<-paste(result,'[', modificationName[mod[i+1]+1], ']', sep='')
        }
    }
    return(result)
}

.PTM_MarkerFinder_no_peakplot<-function(data, 
    modification,
    modificationName,
    mZmarkerIons, 
    minNumberIons=2, 
    itol_ppm=10, 
    minMarkerIntensityRatio=5){

    lapply(data, function(x){
            idx<-findNN(mZmarkerIons, x$mZ)
             ppm.error <- 1e-06 * itol_ppm * x$mZ[idx]
             b<-(abs(mZmarkerIons-x$mZ[idx]) < ppm.error)

             sum.mZmarkerIons.intensity <- sum(x$intensity[idx[b]])
             sum.intensity <- sum(x$intensity) - sum.mZmarkerIons.intensity
             percent.mZmarkerIons <- round(100 * (sum.mZmarkerIons.intensity / (sum.mZmarkerIons.intensity + sum.intensity)),1)

             if ((length(x$mZ[idx[b]]) >= minNumberIons) & percent.mZmarkerIons > minMarkerIntensityRatio ){
            split.screen(c(2,1))
            split.screen(c(1,4), 2)
            screen(1)
            plot(x$mZ, x$intensity, type='h', xlab='m/z',ylab='Intensity', xlim=c(0, max(x$mZ)), sub=paste(x$title))
            points(x$mZ[idx[b]], x$intensity[idx[b]], pch=22, col='#6CC417AA', bg='#6CC417AA', cex=1)
            legend("topright", paste(c('query', 'pepmass', 'charge'), c(x$id, x$pepmass, x$charge)), cex=0.75)

################################################################################
            screen(6)
            pie(x$intensity[idx[b]], 
                paste(round(100*x$intensity[idx[b]]/sum(x$intensity[idx[b]]),0),"% ", 
                    round(x$mZ[idx[b]], 2), sep=""), 
                main=list("relative marker ions %", cex=0.8),
                cex=0.75)

################################################################################
            screen(5)
            pie(c(sum.intensity, sum.mZmarkerIons.intensity), 
                c('Intensity', paste(percent.mZmarkerIons, '%',sep='')),
                col=c('white', '#6CC417AA'),
                main=list("% marker ion intensities", cex=0.8),
                cex=0.75)

################################################################################
            screen(3)
            #plot(mZmarkerIons, (mZmarkerIons-data[[i]]$mZ[idx]),
            plot(mZmarkerIons[b], 1e+06 * (mZmarkerIons[b] - x$mZ[idx[b]]) / x$mZ[idx[b]],
                axes=TRUE, type='p', xlab='m/z', ylab='ppm error',log='')# ylim=c(-max(ppm.error),max(ppm.error)))
            axis(3, mZmarkerIons[b], round(mZmarkerIons[b],1),las=2)
            abline(h=0.0, col='grey')

################################################################################
####### P A T T E R N P L O T ##################################################
            screen(4)
            plot(x$mZ[idx[b]],
                x$intensity[idx[b]],
                col='#6CC417AA', 
                lwd=2, 
                type='h',
                xlab='m/z',
                ylab='Intensity',
                axes=TRUE)
            axis(3, x$mZ[idx[b]],round(x$mZ[idx[b]],1), las=2)

            # readline("Press <Enter> to continue")
            close.screen(all.screens = TRUE)
        }
    close.screen(all.screens = TRUE)
    }
    )
}

# mZmarker which ion you want to see
PTM_MarkerFinder <- function(data, 
    modification, 
    modificationName, 
    mZmarkerIons, 
    minNumberIons=2, 
    itol_ppm=10, 
    minMarkerIntensityRatio=5,
    mgfFilename=-1, 
    PEAKPLOT=TRUE){

    if(!PEAKPLOT){
        .PTM_MarkerFinder_no_peakplot(data, 
            modification,
            modificationName,
            mZmarkerIons, 
            minNumberIons, 
            itol_ppm,
            minMarkerIntensityRatio)
        return (TRUE)
    }


    query.idx<-1:length(data)
    query.to.scan<-as.integer(as.character(lapply(data, function(x){
            if (length(x$scans) == 1){
                return(x$scans)
            }else{return (x$scans[1])}
            })))

    scan.to.query<-rep(-1, max(query.to.scan, na.rm=TRUE))
    scan.to.query[query.to.scan[query.idx]]<-query.idx


    if (mgfFilename != -1){
        unlink(mgfFilename)
        .PTM_MarkerFinder_writeMGF_Header(mgfFilename)
    }

    rr<-numeric()

    for (i in 2:length(data)){

        idx<-findNN(mZmarkerIons, data[[i]]$mZ)

# ppm.itol.cutoff
        ppm.error <- 1e-06 * itol_ppm * data[[i]]$mZ[idx]
        b<-(abs(mZmarkerIons-data[[i]]$mZ[idx]) < ppm.error)

        sum.mZmarkerIons.intensity <- sum(data[[i]]$intensity[idx[b]])
        sum.intensity <- sum(data[[i]]$intensity) - sum.mZmarkerIons.intensity
        percent.mZmarkerIons <- round(100 * (sum.mZmarkerIons.intensity / (sum.mZmarkerIons.intensity + sum.intensity)),1)

        if ((length(data[[i]]$mZ[idx[b]]) >= minNumberIons) 
            & percent.mZmarkerIons > minMarkerIntensityRatio ){

           if (mgfFilename != -1){

            .PTM_MarkerFinder_writeMGF(data[[i]], mgfFilename, 
                pattern=paste(data[[i]]$mZ[idx[b]], data[[i]]$intensity[idx[b]]))
            }
           r <- cbind(scans=data[[i]]$scans, 
                mZ=data[[i]]$mZ[idx[b]],
                markerIonMZ=mZmarkerIons[b],
                markerIonIntensity=data[[i]]$intensity[idx[b]], 
                markerIonMzError=mZmarkerIons[b]-data[[i]]$mZ[idx[b]],
                markerIonPpmError=1e+06 * (mZmarkerIons[b]-data[[i]]$mZ[idx[b]])/data[[i]]$mZ[idx[b]],
                query=i,
                pepmass=data[[i]]$pepmass
                )

           rr<-rbind(rr,r)

            split.screen(c(3,1))
            split.screen(c(1,5), 2)
            screen(1)

    ####### P E A K P L O T ########################################################
            if (!is.na(data[[i]]$mascotScore) ){
                fi<-fragmentIon(sequence=data[[i]]$peptideSequence, 
                    FUN=defaultIon,
                    modified=substr(data[[i]]$modification, 2, nchar(data[[i]]$modification)-1),
                    modification=modification)

                fi.by<-as.data.frame(cbind(b=fi[[1]]$b, y=fi[[1]]$y))

                peakplot(data[[i]]$peptideSequence, 
                    spec=data[[i]], fi=fi.by, ion.axes=FALSE,  
                    main=paste("scantype: HCD / peptide sequence: ", 
                        .PTM_MarkerFinder_(data[[i]]$peptideSequence, data[[i]]$modification, modificationName), sep=''),
                    xlim=c(0,max(data[[i]]$mZ)))
            }else{
                par(cex=0.75)
                plot(data[[i]]$mZ, data[[i]]$intensity, type='h',
                    xlab='m/z',ylab='Intensity',
                    main=paste("scantype: HCD"),
                    xlim=c(0, max(data[[i]]$mZ)),
                    sub=paste(i, data[[i]]$title))
            }

            points(data[[i]]$mZ[idx[b]],
                data[[i]]$intensity[idx[b]], 
                pch=22,
                col='#6CC417AA', 
                bg='#6CC417AA', 
                cex=1)

            legend("topright", paste(c('protein', 'm/z', 'charge', 'ionScore', 'query#'), 
                c(data[[i]]$proteinInformation, 
                data[[i]]$pepmass, data[[i]]$charge, data[[i]]$mascotScore,i)), cex=0.75)

                #c(data[[i]]$proteinInformation, 
        

################################################################################
####### P E A K P L O T ########################################################
            screen(3)
            
            j<-1;
            PLOTFLAG<-FALSE
            while (j<5){ 
                #k<-i+j
                k <- scan.to.query[query.to.scan[i] + j]

                if (k < 0){break;}

            if (!is.na(data[[k]]$mascotScore) & abs(data[[k]]$pepmass - data[[i]]$pepmass) < 1){


           if (mgfFilename != -1){
                .PTM_MarkerFinder_writeMGF(data[[k]], mgfFilename)
           }

            fi<-fragmentIon(sequence=data[[k]]$peptideSequence, 
                FUN=defaultIon,
                modified=substr(data[[k]]$modification, 2, nchar(data[[k]]$modification)-1),
                modification=modification)

            fi.cyz<-as.data.frame(cbind(y=fi[[1]]$y, c=fi[[1]]$c, z=fi[[1]]$z))

            p<-peakplot(data[[k]]$peptideSequence, spec=data[[k]], 
                main=paste("scantype: ETD / peptide sequence: ",
                        .PTM_MarkerFinder_(data[[k]]$peptideSequence, data[[k]]$modification, modificationName), sep=''),
                fi=fi.cyz,
                itol=0.6,
                    xlim=c(0, max(data[[i]]$mZ)),
                ion.axes=FALSE)
            PLOTFLAG<-TRUE

            legend("topleft", paste(c('protein', 'm/z', 'charge', 'ionScore', 'query#'), 
                c(data[[k]]$proteinInformation, 
                data[[k]]$pepmass, data[[k]]$charge, data[[k]]$mascotScore, k)), cex=0.75)
            break;
            }
                j <- j+1
            }

            if (PLOTFLAG == FALSE){
                plot(0,0,type='n', xlab='',ylab='', axes=FALSE)
                text(0,0,"no peptide assignment", cex=3, col="lightgrey")
            }
################################################################################
            screen(8)
            plot(0,0,axes=FALSE,type='n', xlab='',ylab='')
################################################################################
            screen(7)
            pie(data[[i]]$intensity[idx[b]], 
                paste(round(100*data[[i]]$intensity[idx[b]]/sum(data[[i]]$intensity[idx[b]]),0),"% ", 
                    round(data[[i]]$mZ[idx[b]], 2), sep=""), 
                main=list("relative marker ions %", cex=0.8),
                cex=0.75)

################################################################################
            screen(6)
            pie(c(sum.intensity, sum.mZmarkerIons.intensity), 
                c('Intensity', paste(percent.mZmarkerIons, '%',sep='')),
                col=c('white', '#6CC417AA'),
                main=list("% marker ion intensities", cex=0.8),
                cex=0.75)

################################################################################
            screen(4)
            #plot(mZmarkerIons, (mZmarkerIons-data[[i]]$mZ[idx]),
            plot(mZmarkerIons[b], 1e+06 * (mZmarkerIons[b]-data[[i]]$mZ[idx[b]])/data[[i]]$mZ[idx[b]],
                axes=TRUE,type='p', xlab='m/z',ylab='ppm error',log='')# ylim=c(-max(ppm.error),max(ppm.error)))
            axis(3, mZmarkerIons[b], round(mZmarkerIons[b],1),las=2)
            abline(h=0.0, col='grey')

################################################################################
####### P A T T E R N P L O T ##################################################
            screen(5)
            plot(data[[i]]$mZ[idx[b]],
                data[[i]]$intensity[idx[b]],
                col='#6CC417AA', 
                lwd=2, 
                type='h',
                xlab='m/z',
                ylab='Intensity',
                axes=TRUE)
            axis(3,data[[i]]$mZ[idx[b]],round(data[[i]]$mZ[idx[b]],1), las=2)

            # readline("Press <Enter> to continue")
            close.screen(all.screens = TRUE)
        }
    }
    close.screen(all.screens = TRUE)

    return(as.data.frame(rr))
}


# code for screening multiple modifications
# by Paolo <paolo.nanni@fgcz.uzh.ch> and Christian <cp@fgcz.ethz.ch>, March 2013
#
# processed using R version 2.15.3 (2013-03-01)
# protViz packageVersion  0.1.41
# Linux fgcz-c-071 2.6.32-5-amd64 #1 SMP Mon Jan 16 16:22:28 UTC 2012 x86_64 GNU/Linux
# Intel(R) Xeon(R) CPU  X5650  @ 2.67GHz

.PTM_MarkerFinder_summaryPlot<-function(s, config){
    plot(0,0,xlab='', ylab='', axes=FALSE, type='n')
    legend(-0.95,+0.5, paste(c(names(config), 'R', 'protViz.version', 'create date'), 
        c(config, R.version.string, as.character(packageVersion('protViz')), date()), 
        sep=": "), cex=2.0)

    if (nrow(s) > 1){
        op<-par(mfrow=c(4,5), mar=c(4,4,3,1)); 
        dump <- lapply(split(s,s$query), 

    function(x){ plot(x$mZ, x$markerIonIntensity, 
        type='h',
        col='lightblue',
        cex=2,
        ylab='intensity', xlab='m/z',
        ylim=range(s$markerIonIntensity),
        log='y',
        main=paste("scan=",unique(x$scans),"/query=", unique(x$query), sep='')); 
        text(x$mZ, x$markerIonIntensity,round(x$mZ,1),col='red',cex=0.7)})

par(op)

op<-par(mfrow=c(4,5), mar=c(1,1,3,1)); 
dump <- lapply(split(s,s$query), 
    function(x){ pie(x$markerIonIntensity, 
        cex=0.7,
        round(x$markerIonMZ,2), 
        col='white',
        main=paste("scan=",unique(x$scans),"/query=", unique(x$query), sep=''))})
par(op)

op<-par(mfrow=c(1,1), mar=c(5,5,5,5)); 

boxplot(s$markerIonIntensity ~ s$markerIonMZ,
    log='y',
    main='summary marker ion intensity distribution',
    xlab='markerIon m/z', ylab='log10 based marker ion intensity')
box()

boxplot(markerIonPpmError~as.factor(scans),
    data=s, 
    main='marker ions PPM error summary',
    xlab='scan number', ylab='ppm error')
abline(h=0.0,col='grey')

par(mfrow=c(1,1));
barplot(tapply(s$markerIonIntensity, s$scans, FUN=sum),log='y',
    main='sum marker ion intensities summary', 
    xlab='scan numbers',
    ylab='sum marker ion intensities [log10]'
)
}
}


PTM_MarkerFinder_util<-function(dataFileName,   
    mZmarkerIons=c(126.0550, 138.0550, 144.0655, 168.0655, 186.0761, 204.0866), 
    minMarkerIntensityRatio=2, 
    minNumberIons=2, 
    itol_ppm=20, 
    write_csv=FALSE){

    config<-list(
        dataFileName=dataFileName,
        MarkerIon=mZmarkerIons,
        minMarkerIntensityRatio=minMarkerIntensityRatio,
        minNumberIons=minNumberIons,
        itol_ppm=itol_ppm)

    load(paste(config$dataFileName, "RData", sep='.'))

    pdf(paste(config$dataFileName, "pdf", sep='.'), 19, 12)
    s<-PTM_MarkerFinder(data=get(config$dataFileName), 
        modification=get(paste(config$dataFileName, "modification.mass", sep=".")),
        modificationName=get(paste(config$dataFileName, "modification.description", sep = '.')),
        minMarkerIntensityRatio=config$minMarkerIntensityRatio,
        minNumberIons=config$minNumberIons,
        itol_ppm=config$itol_ppm,
        mgfFilename=paste(config$dataFileName,".mgf", sep=''),
        mZmarkerIons=config$MarkerIon)
    dev.off()

# overwrite the it again.
    pdf(paste(config$dataFileName, "pdf", sep='.'), 19, 12)
    .PTM_MarkerFinder_summaryPlot(s, config)

    s<-PTM_MarkerFinder(data=get(config$dataFileName), 
        modification=get(paste(config$dataFileName, "modification.mass", sep=".")),
        modificationName=get(paste(config$dataFileName, "modification.description", sep = '.')),
        minMarkerIntensityRatio=config$minMarkerIntensityRatio,
        minNumberIons=config$minNumberIons,
        itol_ppm=config$itol_ppm,
        mZmarkerIons=config$MarkerIon)

    dev.off()

    if ( nrow(s) > 1 & write_csv==TRUE ){
        w<-reshape(s[, c(1,7,8,3,4)], direction='wide', timevar="markerIonMZ", idvar=c('scans','query','pepmass'))

        write.table(w, file=paste(config$dataFileName, "csv", sep='.'),
            sep=',', row.names=FALSE, col.names=TRUE, quote=FALSE)
    }

#    return(s)
}

# ToDo adapt the parameter mZmarkerIon in all function 

