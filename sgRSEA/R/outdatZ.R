outdatZ <-
function(dat, Zvec, Txx){
    Tvec <- Txx[,1]
    mvec <- Txx[,2]
    mvec1 = inverse.rle( list(lengths=mvec, values=mvec))
    Tvec1 = inverse.rle( list(lengths=mvec, values=Tvec))

    datZ <- dat
    datZ$Z = Zvec
    datZ$T = Tvec1
    datZ$m = mvec1
    return(datZ)
    }
