######### Check the parameters are valid ####################################

skewhypCheckPars <- function(param){

    param <- as.numeric(param)

    if (length(param) !=4){
        case <- "error"
        errMessage <- "param vector must contain 4 values"
    }else{
        mu <- param[1]
        delta <- param[2]
        beta <- param[3]
        nu <- param[4]

        case <- "normal"
        errMessage <- ""

        if( delta <= 0){
            case <- "error"
            errMessage <- "Delta must be greater than zero"
        }
        if( nu < 0 ){
            case <- "error"
            errMessage <- "Nu must be greater than zero"
        }
    }
    result <- list(case = case, errMessage = errMessage)
    return(result)

}
