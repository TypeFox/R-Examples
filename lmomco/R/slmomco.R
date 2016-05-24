"slmomco" <-
function(x,para) {
    if(! are.par.valid(para)) return()
    f <- par2cdf(x,para,paracheck=FALSE)
    return(1 - f)
}
