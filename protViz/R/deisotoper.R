#R

# $HeadURL: http://fgcz-svn.unizh.ch/repos/fgcz/testing/proteomics/R/protViz/R/deisotoper.R $
# $Id: deisotoper.R 6178 2014-02-27 09:33:30Z cpanse $
# $Date: 2014-02-27 10:33:30 +0100 (Thu, 27 Feb 2014) $

deisotoper <- function(data,
    Z=1:4,
    isotopPatternDF=averagine,
    massError=0.005,
    plot=FALSE){

    colormap=rainbow(length(Z), alpha=0.5)

    val <- lapply(data, function(x){
        if ( length(x$mZ) > 1 ){
            # TODO check if sorted 
            mZ.idx <- order(x$mZ)
            x$mZ <- x$mZ[mZ.idx]
            x$intensity <- x$intensity[mZ.idx]

            out <- .Call("deisotoper_main", x$mZ, x$intensity, Z, averagine, massError, PACKAGE="protViz")
            out$group <- lapply(out$group, function(x){ x[x==-1] <- NA; return(x)})

            if(plot){
                op<-par(mfrow=c(1,1), mar=c(4,4,4,4))
                plot(x$mZ, x$intensity, col='grey',type ='h', log='', main=x$title)

                for (i in 1:length(Z)){
                    mapply(function(xx, ss){
                        if (ss > 0.80 & length(xx)>2){
                        points(x$mZ[xx+1], x$intensity[xx+1], type='h', col=colormap[i], lwd=max(Z)-i)

                        iso.mean<-mean(x$mZ[xx+1])

                        text(x$mZ[min(xx+1)], x$intensity[min(xx+1)], round(x$mZ[min(xx+1)],2), srt=0, cex=0.50, pos=3, col=colormap[i]) 
                        text(iso.mean, x$intensity[min(xx+1)], round(iso.mean,2), srt=0, cex=0.50, pos=4, col=colormap[i]) 
                        }
                    }, out$result[[i]], out$score[[i]])
                }
                .deisotoperUtilPlot(x, out, Z, colormap)
            }
        }
        return(out)
    })
    return(val)
}

.deisotoperUtilPlot <- function(x, out, Z, colormap){
    op<-par(mfrow=c(2, 2), mar=c(3,3,3,1))

    lapply(out$group, function(g){
        for (i in 1:length(Z)){
            if (! is.na(g[i])){
                idx<-(out$result[[i]][[g[i] + 1]]) + 1

                my.title <- paste("mass=", x$mZ[min(idx)],sep='')
                plot(xx<-x$mZ[min(idx):max(idx)], yy<-x$intensity[min(idx):max(idx)], 
                    type='h', 
                    axes=FALSE,
                    col='grey', 
                    xlim=c(x$mZ[min(idx)]-0.5, x$mZ[max(idx)]+0.5),
                    ylim=c(0, max(x$intensity[min(idx):max(idx)])),
                    xlab="mZ", 
                    main=paste("mass=", x$mZ[min(idx)],sep=''),
                    ylab="intensity", 
                )

                text(xx,yy,xx,pos=4,cex=0.75)
                axis(1, x$mZ[min(idx):max(idx)], round(x$mZ[min(idx):max(idx)], 2))
                axis(2)
                break;
            }
        }

        score <- rep(NA, length(Z))
        score1 <- rep(NA, length(Z))
        cscore1 <- rep(NA, length(Z))
        cscore <- rep(NA, length(Z))

        for (i in 1:length(Z)){
            if (! is.na(g[i])){
                idx<-(out$result[[i]][[g[i]+1]])
                cscore[i]<-round(out$score[[i]][[g[i]+1]],2)
                cscore1[i]<-round(out$score1[[i]][[g[i]+1]],2)

                averagine.mZ <- x$mZ[idx+1]
                averagine.intensity <- averagine[, findNN_(min(x$mZ[idx]) * Z[i], as.double(colnames(averagine)))]

                L2.intensity <- sqrt(sum((x$intensity[idx+1])^2))
                score[i] <- round((x$intensity[idx+1] / L2.intensity) %*% averagine.intensity[1:length(averagine.mZ)], 2)

                idx2<-(idx+1)[2:length(idx)]

                L2.intensity2 <- sqrt(sum((x$intensity[idx2])^2))
                score1[i] <- round((x$intensity[idx2] / L2.intensity) %*% averagine.intensity[1:(length(averagine.mZ)-1)], 2)

                points(x$mZ[idx+1], x$intensity[idx+1], 
                    type='h', 
                    col='#AAAAAAAA', lwd=3)

                points(a.x<-averagine.mZ, a.y<-max(x$intensity[idx+1])*averagine.intensity[1:length(averagine.mZ)], col=colormap[i], pch=25, lwd=3)
                text(a.x, a.y, length(idx), pos=4, col=colormap[i])
            }
        }

        legend("topright", paste("c",Z,'=(',cscore,", ",cscore1,")", sep=''), pch=25, col=colormap, title='C++score', cex=1.0)

        box()
    })
}
