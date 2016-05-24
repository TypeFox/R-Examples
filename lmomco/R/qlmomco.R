"qlmomco" <-
function(f,para) {
    if(! are.par.valid(para)) return()
    x <- par2qua(f,para,paracheck=FALSE)
    return(x)
}
