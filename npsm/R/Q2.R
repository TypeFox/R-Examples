Q2 <-
function(z) {

        L05<-lmean(z,0.05)
        L50<-lmean(z,0.5)

        U05<-umean(z,0.05)
        U50<-umean(z,0.5)

        (U05-L05)/(U50-L50)

}
