spodt.tree <-
function(object)     #ancienne fonction arbre()
{
    nf <- nb.feuilles(object@racine)
    cplx <- complexite(object@racine)

    L <- c(0, 2 * nf)
    H <- c(-6 * cplx - 2, 0)

    u <- 2
    v <- -1
    mrg <- 1/30
    
    R <- (abs(v*cplx))/(u*nf)

    l <- 2*u/3-mrg
    h <- 5*v
    ht <- 6*v
    
    hf <- -6 * (cplx - 1)

    plot.new()
    par(cex=0.65, oma=c(0,0,0,0), mar=c(0,0,0,0), mai=c(0,0,0,0), fig=c(0,1,0,1))
    plot.window(xlim=L, ylim=H)

    #rect(nf-l, h, nf+l, 0)
    text(nf, h/10, labels=object@racine@n)
    text(nf, 3*h/10, labels=round(object@racine@m, digits=3))
    text(nf, 5*h/10, labels=round(object@racine@v, digits=3))
    #segments(nf-l, 3*h/5, nf+l, 3*h/5)
    text(nf, 7*h/10, labels=round(object@racine@R2, digits=3))
    if (class(object@racine) == "sp.spodt")
    {
        if (object@racine@coeff[1] == 0)
        {
            t <- paste("y<=",round(object@racine@coeff[2],digits=1))
        }
        else if (object@racine@coeff[1] == 1)
        {
            if (object@racine@coeff[2] > 0)
            {
                t <- paste("y<=","x+",round(object@racine@coeff[2],digits=1))
            }
            if (object@racine@coeff[2] < 0)
            {
                t <- paste("y<=","x",round(object@racine@coeff[2],digits=1))
            }
            if (object@racine@coeff[2] == 0)
            {
                t <- paste("y<=","x")
            }
        }
        else
        {
            if (object@racine@coeff[2] > 0)
            {
                t <- paste("y<=",round(object@racine@coeff[1],digits=1),"x+",round(object@racine@coeff[2],digits=1))
            }
            if (object@racine@coeff[2] < 0)
            {
                t <- paste("y<=",round(object@racine@coeff[1],digits=1),"x",round(object@racine@coeff[2],digits=1))
            }
            if (object@racine@coeff[2] == 0)
            {
                t <- paste("y<=",round(object@racine@coeff[1],digits=1),"x")
            }
        }
    }
    if (class(object@racine) == "vql.spodt")
    {
        t <- paste(object@racine@vrbl,"=",object@racine@mod)
    }
    if (class(object@racine) == "vqt.spodt")
    {
        t <- paste(object@racine@vrbl,"<=",round(object@racine@seuil, 2))
    }
    text(nf, 9*h/10, labels=t)
    
    nfg <- nb.feuilles(object@racine@fg)
    nfd <- nb.feuilles(object@racine@fd)
    
    segments(nf, h, nfg, ht)
    segments(nf, h, 2*nfg+nfd, ht)
    
    tracer.noeud(object@racine@fg, nfg, ht,   0,   u, v, l , h, ht, hf)
    tracer.noeud(object@racine@fd, nfd, ht, 2*nfg, u, v, l , h, ht, hf)
}
