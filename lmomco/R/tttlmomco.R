"tttlmomco" <- function(f,para) {
    if(! are.par.valid(para)) return()
    Qu <- par2qua(f,para,paracheck=FALSE)
    Ru <- rrmlmomco(f,para)
    return(Qu - f*Ru)
}
