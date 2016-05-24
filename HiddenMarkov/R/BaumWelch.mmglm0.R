BaumWelch.mmglm0 <- function (object, control = bwcontrol(), ...){
    object <- as.dthmm(object)
    object <- BaumWelch.dthmm(object, control)
    return(as.mmglm0(object))
}

