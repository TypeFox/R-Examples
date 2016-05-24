tracer.noeud <-
function(noeud, x, y, d, u, v, l, h, ht, hf)
{
    x <- d + x

    if (class(noeud) != "f.spodt")
    {
        #rect(x - l, y + h, x + l, y)
        text(x, y + h/10, labels=noeud@n)
        text(x, y + 3*h/10, labels=round(noeud@m, digits=2))
        text(x, y + 5*h/10, labels=round(noeud@v, digits=2))
        #segments(x - l, y + 3*h/5 -1/30, x + l, y + 3*h/5 -1/30)
        text(x, y + 7*h/10, labels=round(noeud@R2, digits=2))
        
        nfg <- nb.feuilles(noeud@fg)
        nfd <- nb.feuilles(noeud@fd)    

        segments(x, y+h, d+nfg, y+ht)
        segments(x, y+h, d+2*nfg+nfd, y+ht)
        if (class(noeud) == "sp.spodt")
        {
            if (noeud@coeff[1] == 0)
            {
                t <- paste("y=",round(noeud@coeff[2],digits=1))
            }
            else
            {
                if (noeud@coeff[2] > 0)
                {
                    t <- paste("y<=",round(noeud@coeff[1],digits=1),"x+",round(noeud@coeff[2],digits=1))
                }
                if (noeud@coeff[2] < 0)
                {
                    t <- paste("y<=",round(noeud@coeff[1],digits=1),"x",round(noeud@coeff[2],digits=1))
                }
                if (noeud@coeff[2] == 0)
                {
                    t <- paste("y<=",round(noeud@coeff[1],digits=1),"x")
                }
            }
        }
        if (class(noeud) == "vql.spodt")
        {
            t <- paste(noeud@vrbl,"=",noeud@mod)
        }
        if (class(noeud) == "vqt.spodt")
        {
            t <- paste(noeud@vrbl,"<=",round(noeud@seuil, 2))
        }
        text(x, y + 9*h/10, labels=t)
        tracer.noeud(noeud@fg, nfg, y+ht,   d,     u, v, l , h, ht, hf)
        tracer.noeud(noeud@fd, nfd, y+ht, d+2*nfg, u, v, l , h, ht, hf)
    }
    else
    {
        h <- h - v
        segments(x, y, x, hf)
        #rect(x-l/2, hf+h, x+l/2, hf)
        text(x, hf + h/8, noeud@id)
        #segments(x - l/2, hf + h/5 -1/30, x + l/2, hf + h/5 -1/30)
        text(x, hf + 3*h/8, labels=noeud@n)
        text(x, hf + 5*h/8, labels=round(noeud@m, digits=1))
        text(x, hf + 7*h/8, labels=round(noeud@v, digits=2))
    } 
}
