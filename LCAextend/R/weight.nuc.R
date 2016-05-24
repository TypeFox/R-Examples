weight.nuc <-
function(connect,spouse.connect,children.connect,status,probs,fyc,p.ybarF.c,ww,w,res.upward)
{
    p.yF.c <- res.upward$p.yF.c
    sum.child <- res.upward$sum.child

    p.spouse <- p.post.found(spouse.connect,status,probs,fyc)
    if(length(children.connect)>0) for(child in children.connect)
    {
        for(c.connect in 1:(length(probs$p)+1)) for(c.spouse in 1:(length(probs$p)+1))
        {
            p.child <- p.post.child(child,c.connect,c.spouse,status,probs,fyc)
            p.child.complet <- matrix(0,nrow=2,ncol=length(probs$p)+1)
            for(s in 1:2) p.child.complet[s,] <- p.yF.c[child,]*p.child[s,]

            bro.child <- setdiff(children.connect,child)
            prod.bro <- 1
            if(length(bro.child)>0) prod.bro <- prod(sum.child[bro.child,c.connect,c.spouse])
            for(s1 in 1:2) for(s2 in 1:2) for(si in 1:2)
            ww[child,s1,s2,si,c.connect,c.spouse,] <- p.ybarF.c[connect,s1,c.connect]*p.spouse[s2,c.spouse]*p.child.complet[si,]*prod.bro
                
        }
        w[child,,] <- apply(ww[child,,,,,,],c(3,6),sum)
    }
    w[connect,,] <- apply(ww[children.connect[1],,,,,,],c(1,4),sum)
    w[spouse.connect,,] <- apply(ww[children.connect[1],,,,,,],c(2,5),sum)

    res <- list("ww"=ww,"w"=w)
    res
}

