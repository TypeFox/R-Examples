"pdfwei" <-
function(x,para) {
    if(! are.parwei.valid(para)) return()
    ZETA <- para$para[1]
    B    <- para$para[2]
    D    <- para$para[3]

    K  <- 1/D
    A  <- B/D
    XI <- ZETA - B 
    gev.para <- list(type = 'gev', para = c(XI,A,K))
    return(pdfgev(-x,gev.para))
}
