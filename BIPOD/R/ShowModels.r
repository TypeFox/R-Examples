ShowModels <- function(Model){
    if(missing(Model)){stop("Model must be either 'OU', 'CIR', 'FHN' or 'FHN5'.")}
    if(!(Model %in% c("OU","FHN","FHN5","CIR"))){
        stop("Model must be either 'OU', 'CIR', 'FHN' or 'FHN5'.")
    }
    if(Model=="OU"){
        tmp <-
            list(Model =
                 rbind(paste("dX_t = ( - drift1 * X_t - drift2 * Y_t + drift5 ) dt + diff1 dB1_t"),
                 paste("dY_t = ( - drift3 * X_t - drift4 * Y_t + drift6 ) dt + diff2 dB2_t")),Ndrift=6,Ndiff=2)
    }
    if(Model=="CIR"){
        tmp <-
            list(Model =
    rbind(paste("dX_t = ( - drift1 * X_t - drift2 * Y_t + drift5 ) dt + diff1 * sqrt( X_t ) dB1_t"),
          paste("dY_t = ( - drift3 * X_t - drift4 * Y_t + drift6 ) dt + diff2 * sqrt( Y_t ) dB2_t")),
                 Ndrift=6,Ndiff=2)
    }
    if(Model=="FHN"){
        tmp <-
            list(Model = rbind(paste("dX_t = ( drift1 * ( X_t - X_t^3 - Y_t ) + drift2 ) dt + diff1 dB1_t"),
                  paste("dY_t = -( drift3 * X_t - Y_t + drift4 ) dt + diff2 dB2_t")),Ndrift=4,Ndiff=2)
    }
    if(Model=="FHN5"){
        tmp <-
            list(Model = rbind(paste("dX_t = ( - drift1 * X_t^3 + drift2 * ( X_t - Y_t ) + drift3 ) dt + diff1 dB1_t"),
                 paste("dY_t =  drift4 * X_t - Y_t + drift5 ) dt + diff2 dB2_t")),Ndrift=5,Ndiff=2)
    }
    return(tmp)
}
