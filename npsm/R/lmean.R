lmean <-
function(z,p) {

        ind<-p*length(z)
        i0<-floor(ind)

        if( i0 < 1 ) {
                Lp<-min(z)
        } else {
                z<-sort(z,partial=1:(i0+1))
                i1<-ind-i0
                Lp<-(sum(z[1:i0])+i1*z[i0+1])/ind
        }

        Lp

}
