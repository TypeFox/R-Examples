PermPval <-
function(obs.T, null.T){
        pooledT = cbind( c(obs.T, null.T), c(rep(0, length(obs.T)), rep(1, length(null.T)) )
        )
        pTs = pooledT[ order(pooledT[,1], decreasing=T), ]
        pTs1 = cbind( pTs, cumsum(pTs[,2]) )
        pTs.obs = pTs1[ (pTs[,2]==0),]

        pvec1 = (pTs.obs[,3]+1)/(1+length(null.T) )
        pvec2 = ( length(null.T) - pTs.obs[,3] + 1) / (1+ length(null.T) )

        pmat = cbind( NScore=pTs.obs[,1], p.value.pos=pvec1, p.value.neg=pvec2)
        return(pmat)
        }
