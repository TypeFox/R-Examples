varTrans <- function(asr.object){
    if(asr.object$sigma2 == 1){
       vars <- diag(aiFun(asr.object))
    } else{
       Rcomp <- which(asr.object$gammas == 1.00)
       AI <- aiFun(asr.object)
       comps <- asr.object$gammas * asr.object$sigma2
       vars <- c(((asr.object$gammas[-Rcomp]^2) * diag(AI)[Rcomp] + comps[Rcomp]^2 * diag(AI)[-Rcomp] + 2*asr.object$gammas[-Rcomp]*comps[Rcomp]*AI[Rcomp, -Rcomp]), diag(AI)[Rcomp])
      }
vars
}

