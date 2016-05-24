Tmean.Null.Tss2 <-
function(gnames, null.Z,  multiplier=30){

    multi.sim=NULL
    for (k in 1:multiplier){
        resampled.Z = sample(null.Z, size=length(gnames), replace=FALSE)
        re.geneZ = cbind(gnames, resampled.Z)
        sim.Txx =  T.mean(geneZ=re.geneZ)
        multi.sim = rbind(multi.sim, sim.Txx)
        }
    return(multi.sim)
    }
