plotScores <-
function(sgRSEA.out, m){
        obsTmat = sgRSEA.out$Tmat
        nullTmat = sgRSEA.out$NullTmat
        mvec = obsTmat[,2]
        cat('m values:')
        print(sort(unique(mvec)))
        if (all(unique(mvec)!=m)) stop( paste('There is no m=',m,'.',sep='') )
        null.T = nullTmat[(nullTmat[,2]==m), 1]
        obs.T = obsTmat[ (obsTmat[,2]==m), 1 ]
        par(mfrow=c(2,1))
        hist(null.T, breaks='FD', freq=F, xlim=range(obs.T), main=NULL   )
        hist(obs.T, breaks='FD', freq=F, xlim=range(obs.T), main=NULL   )
        title( main=paste('Gene Scores, m=',m,sep=''), line=-1.5, outer=TRUE )
        }
