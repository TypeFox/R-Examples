upward.connect <-
function(connect,spouse.connect,children.connect,status,probs,p.yF.c,fyc,sum.child)
{
    p.spouse <- p.post.found(spouse.connect,status,probs,fyc)
    uu <- matrix(0,nrow=length(probs$p)+1,ncol=length(probs$p)+1)
    for(c.connect in 1:(length(probs$p)+1)) for(c.spouse in 1:(length(probs$p)+1))
    {
        u <- 1
        if(length(children.connect)>0) for(child in children.connect)
        {
            p.child <- p.post.child(child,c.connect,c.spouse,status,probs,fyc)
            p.child.complet <- matrix(0,nrow=2,ncol=length(probs$p)+1)
            for(s in 1:2) p.child.complet[s,] <- p.yF.c[child,]*p.child[s,]
            sum.child[child,c.connect,c.spouse] <- sum(p.child.complet)
            u <- u*sum.child[child,c.connect,c.spouse]
        }
        uu[c.connect,c.spouse] <- u*sum(p.spouse[,c.spouse])
    }
    p.yF.c[connect,] <- apply(uu,1,sum)
    res <- list("sum.child"=sum.child,"p.yF.c"=p.yF.c)
    res
}

