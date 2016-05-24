"lkhlmomco" <- function(f,para) {
    if(! are.par.valid(para)) return()
    Ku <- lrzlmomco(1-f,para)
    return(1 - Ku)
}
