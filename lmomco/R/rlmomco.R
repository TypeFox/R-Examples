"rlmomco" <-
function(n,para) {
    if(! are.par.valid(para)) return()
    f <- runif(n)
    x <- par2qua(f,para,paracheck=FALSE)
    return(x)
}
