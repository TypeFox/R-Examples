plotNScores <-
function(sgRSEA.out){
        stdTmat = sgRSEA.out$stdTmat
        stdNullTmat = sgRSEA.out$stdNullTmat
        par(mfrow=c(2,1))
        hist(stdNullTmat[,1], breaks='FD', freq=F, xlim=range(stdTmat[,1]), xlab='null dist',
        main=NULL)
        hist(stdTmat[,1], breaks='FD', freq=F, xlim=range(stdTmat[,1]), xlab='observed',
        main=NULL  )
        title(main='Normalized Gene Scores', line=-1.5, outer=TRUE)
        }
