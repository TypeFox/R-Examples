#R

# $HeadURL: http://fgcz-svn.unizh.ch/repos/fgcz/testing/proteomics/R/protViz/R/peakplot.R $
# $Id: peakplot.R 6683 2014-09-18 06:52:41Z cpanse $
# $Date: 2014-09-18 08:52:41 +0200 (Thu, 18 Sep 2014) $

# TODO: make a class peakplot
# peakplot.print
# peakplot.plot

.peakplot.putlabel<-function(MASS, INTENSITY, LABEL, 
    l.col="green", 
    delta=0, 
    yMin, 
    maxIntensity=max(INTENSITY)) {

	noise<-seq(0, 2 * delta, length=length(MASS))

	segments(MASS, 
		INTENSITY+0.03*maxIntensity,
		MASS, 
		1.1*maxIntensity,
        lty=2,
		#yMin+noise,
		col=l.col, 
		pch=5, 
		cex=0.5,
		lwd=0.5)

	segments(MASS, INTENSITY, MASS, rep(0, length(MASS)), lwd=1.5, col=l.col)


	text(MASS, yMin+noise, round(MASS,2), 
		cex=0.50,
		pos=1,
		#offset=0.0,
		col="red",
		srt=90)
}

.peakplot.label<-function(spec, match, yMax){

    bin.n<-10
    lab.n<-3
    bin<-seq(min(spec$mZ), max(spec$mZ), length=(bin.n+1))
    bin.min<-bin[seq(1,length(bin)-1)];
    bin.max<-bin[seq(2,length(bin))];
    bin.range<-1:lab.n

    # filtering the ions
    LABEL.abc<-(abs(match$mZ.Da.error) < 0.6) & (regexpr("[abc].*", match$label) > 0)
    LABEL.xyz<-(abs(match$mZ.Da.error) < 0.6) & (regexpr("[xyz].*", match$label) > 0)

    points(spec$mZ[match$idx[LABEL.abc]], spec$intensity[match$idx[LABEL.abc]], col="black", cex=1.0, pch=22)
    points(spec$mZ[match$idx[LABEL.xyz]], spec$intensity[match$idx[LABEL.xyz]], col="blue", cex=1.0, pch=22)

    for (i in c(1:length(bin)-1)) {
        bin.label <- bin.min[i] <= spec$mZ[match$idx] & spec$mZ[match$idx] < bin.max[i]

        bin.label.abc <- bin.label & LABEL.abc 
        bin.label.xyz <- bin.label & LABEL.xyz 

        bin.label.filter <- bin.label.abc | bin.label.xyz

        if (length(bin.label.filter)>0){ 

        sum.intensity <- sum(spec$intensity[match$idx[bin.label.filter]], na.rm = TRUE)

        if (sum.intensity > 0){
            d<-(yMax - max(spec$intensity[match$idx[bin.label.filter]], na.rm=TRUE)) / 2

            if (length(bin.label.abc) > 0){
                if (sum(spec$intensity[match$idx[bin.label.abc]], na.rm=TRUE) > 0)
            .peakplot.putlabel(MASS=spec$mZ[match$idx[bin.label.abc]], 
                INTENSITY=spec$intensity[match$idx[bin.label.abc]], 
                LABEL=match$label[bin.label.abc],
                l.col="black",
                yMin=1.1 * d + max(spec$intensity[match$idx[bin.label.abc | bin.label.xyz]]), 
                delta=d, 
                maxIntensity=max(spec$intensity))
            }

            if (length(bin.label.xyz) > 0){
                if ( sum(spec$intensity[match$idx[bin.label.xyz]], na.rm=TRUE) > 0)
            .peakplot.putlabel(MASS=spec$mZ[match$idx[bin.label.xyz]], 
                INTENSITY=spec$intensity[match$idx[bin.label.xyz]], 
                LABEL=match$label[bin.label.xyz],
                l.col="blue",
                yMin=d + max(spec$intensity[match$idx[bin.label.xyz | bin.label.xyz]]), 
                delta=d, 
                maxIntensity=max(spec$intensity))
            }
        }
        }
    }
}

.peakplot.pie <- function(spec, match){ 

    LABEL.abc<-abs(match$mZ.Da.error < 0.6) & (regexpr("[abc].*", match$label) > 0)
    LABEL.xyz<-abs(match$mZ.Da.error < 0.6) & (regexpr("[xyz].*", match$label) > 0)

    i.abc<-spec$intensity[match$idx[LABEL.abc]]
    i.xyz<-spec$intensity[match$idx[LABEL.xyz]]

    l.abc<-match$label[LABEL.abc]
    l.xyz<-match$label[LABEL.xyz]

    i.rest<-sum(spec$intensity)-sum(i.abc)-sum(i.xyz)

    pie(c(i.abc,i.xyz,i.rest), c(l.abc, l.xyz, "rest"), col=c(rep("blue",length(i.abc)), rep("grey",length(i.abc)), "white"))
}                      

peakplot <- function(peptideSequence,
    spec, 
    FUN=defaultIon, 
    fi=fragmentIon(peptideSequence, FUN=FUN)[[1]],
    main=NULL,
    sub=paste(peptideSequence, spec$title, sep=" / "),
    xlim=range(spec$mZ, na.rm=TRUE),
    ylim=range(spec$intensity, na.rm=TRUE),
    itol=0.6,
    pattern.abc="[abc].*",
    pattern.xyz="[xyz].*",
    ion.axes=TRUE){ 

    n<-nchar(peptideSequence)

    m<-psm(peptideSequence, spec, FUN, fi=fi, plot=FALSE)

    max.intensity<-max(spec$intensity, na.rm=TRUE)
    yMax <- 1.0 * max.intensity

    plot(spec$mZ, spec$intensity,
        xlab='m/z',
        ylab='Intensity',
        type='h',
        main=main,
        xlim=xlim,
        ylim=c(0, 1.2 * yMax),
        sub=sub,
        axes='F'
    ) 

    LABEL.abc<-(abs(m$mZ.Da.error) < itol) & (regexpr(pattern.abc, m$label) > 0)
    LABEL.xyz<-(abs(m$mZ.Da.error) < itol) & (regexpr(pattern.xyz, m$label) > 0)

    if (ion.axes){
        if (length(m$idx[LABEL.abc]) > 0){
            axis(1, spec$mZ[m$idx[LABEL.abc]], m$label[LABEL.abc],las=2)
        }
        axis(2)
        if (length(m$idx[LABEL.xyz]) > 0){
            axis(3, spec$mZ[m$idx[LABEL.xyz]], m$label[LABEL.xyz], col.axis='blue', las=2)
        }
    }else{
        axis(1)
        axis(2)
        a.at <- spec$mZ[m$idx[LABEL.abc | LABEL.xyz]]
        a.label <- m$label[LABEL.abc | LABEL.xyz]

        if (length(a.at)>0) {
            axis(3,a.at, a.label, col.axis='black', las=2)
        } else {
            print ("WARNING")
            print (a.at)
            print (a.label)
        }
        box()
    }
    axis(4,seq(0,yMax,length=6), seq(0,100,length=6))

    .peakplot.label(spec=spec, match=m, yMax=yMax)

    return(m)
}                      
