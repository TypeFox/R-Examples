STTtoShapeParameter <-
function(STT) {

        if(STT> 0.4) stop("only works for STT <= 0.4")
        f<-function(w,STT) abs(w-qgamma(1-STT,shape=w))
        out<-optimize(f,lower=0,upper=2,STT=STT)$minimum
        return(out)

}
