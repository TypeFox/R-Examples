stdofT <-
function( tnmvec, meanarr, sdarr ){
        mvalue = tnmvec[2]
        mum = meanarr[ which(meanarr[,2]==mvalue),1]
        sigmam = sdarr[ which(sdarr[,2]==mvalue),1]
        stdtn = (tnmvec[1] - mum)/sigmam
        return( c(stdtn, mvalue) )
        }
