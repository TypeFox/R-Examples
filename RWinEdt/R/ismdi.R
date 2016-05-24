ismdi <- function(){
    return(mdi = as.logical(.C("ismodemdi", as.integer(0), PACKAGE = "RWinEdt")[[1]]))
}
