"bfrlmomco" <- function(f,para) {
    if(! are.par.valid(para)) return()
    Bu <- lrzlmomco(f,para)/f
    return(Bu)
}
