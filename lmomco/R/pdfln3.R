"pdfln3" <-
function(x,para) {
    if(! are.parln3.valid(para)) return()
    ZETA <-     para$para[1]
    U    <- exp(para$para[2])
    A    <-     para$para[3]

    gnopara <- vec2par(c(ZETA + U, U*A, -A), type="gno")
    return(pdfgno(x,gnopara))
}
