#R

# $HeadURL: http://fgcz-svn.unizh.ch/repos/fgcz/testing/proteomics/R/protViz/R/de_novo.R $
# $Id: de_novo.R 6178 2014-02-27 09:33:30Z cpanse $
# $Date: 2014-02-27 10:33:30 +0100 (Thu, 27 Feb 2014) $

.de_novo_rev <- function (x, ax=0){

    #return (-(x - ax)+ax)
    return (-x + (2 * ax))

}


.de_novo<-function(mZ, intensity, deltamass, label, error=0.1, precursor=NA, topN=20, sub=''){

    cm<-rainbow(length(deltamass), alpha=0.75)

    y.max<-max(intensity)
    precursor2<-precursor/2

    plot(mZ, intensity, type='h', 
        xlab='mZ', 
        ylab='intensity', 
        log='', 
        xlim=c(0, 1*max(mZ)),
        ylim=c(0, 2 * y.max),
        sub=sub,
        col='black',
        main='protViz de-novo debug plot')


    #c<-combn(1:length(deltamass), 2)

    n<-length(mZ)
    filter<-rep(FALSE, n)
    filter[order(intensity, decreasing=TRUE)[1:topN]] <- TRUE

    abline(v=precursor/2, col='red')

    mZ.mirror <- .de_novo_rev(mZ, precursor2)
    nn<-findNN(mZ.mirror, mZ)
    filter[abs(mZ - mZ.mirror[nn]) < error] <- TRUE;

    mZ<-mZ[filter]
    intensity<-intensity[filter]
    mZ.mirror <- .de_novo_rev(mZ, precursor2)


    segments(mZ.mirror, 2*y.max, mZ.mirror, 2*y.max-intensity)

    abline(v=mZ, lty=2, col='#00AA0099')

    #abline(v=mZ[mZ - .de_novo_rev(mZ, precursor2)[nn] > error], lty=2, col='#AA000099')

    d<-as.matrix(dist(mZ, diag = FALSE, upper = FALSE))
    idx <- order(d) 

    for (i in 1:length(deltamass)){

        if (deltamass[i] > 0){

            r<-which(abs(d - deltamass[i]) < error,  arr.ind = TRUE)

            if (length(r) > 0){

                x <- r[,1]
                y <- r[,2]

                points(mZ[c(x,y)], intensity[c(x, y)], col='#05052511', bg='#05050511', pch=22)

                text(round(c(mZ[r[,2]],mZ[r[,1]])), 
                    round(c(intensity[r[,2]], intensity[r[,1]])), 
                    round(c(mZ[r[,2]],mZ[r[,1]]),2), 
                    cex=0.75, srt=90, col='#FF0000AA', pos=4)

                segments(x0=mZ[x], 
                    y0=0.5*(intensity[x] + intensity[y]),
                    x1=mZ[y], 
                    y1=0.5*(intensity[x] + intensity[y]),
                    col="#11111111",
                    lwd=15
                    )

                text(0.5 * (mZ[x] + mZ[y]), 
                     0.5 * (intensity[x] + intensity[y]), 
                    label[i], 
                    col=cm[i], 
                    cex=1)
            }
        }
    }              
}

de_novo<-function(data){

    AA.name<-c('A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P','Q','R','S','T','U', 'V', 'W', 'X', 'Y', 'Z')
    AA<-c(71.037110, 114.534940, 160.030649, 115.026940, 129.042590, 147.068410, 57.021460, 137.058910, 113.084060, 0.0, 128.094960, 113.084060, 131.040480, 114.042930, 0.000000, 97.052760, 128.058580, 156.101110, 87.032030, 101.047680, 150.953630, 99.068410, 186.079310, 111.000000, 163.063330, 128.550590)
    .de_novo(data$mZ, data$intensity, AA, AA.name, error=0.1, 
        precursor=(data$pepmass[1] * data$charge)-1.008,
        topN=30, sub=data$title)

}
