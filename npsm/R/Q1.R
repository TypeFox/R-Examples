Q1 <-
function(z) {

        U05<-umean(z,0.05)
        M50<-mean(z,trim=0.25)
        L05<-lmean(z,0.05)

        (U05-M50)/(M50-L05)

}
