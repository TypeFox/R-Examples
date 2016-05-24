"riglmomco" <- function(f,para) {
    if(! are.par.valid(para)) return()
    Qu <- par2qua(f,para,paracheck=FALSE)
    Theta1 <- sapply(f, function(u) { return(rreslife.lmoms(u,para, nmom=1)$lambdas[1]) })
    Gu <- 1 - (Theta1/Qu)
    return(Gu)
}
