downward.connect <-
function(connect,parent1,parent2,bro.connect,status,probs,fyc,p.ybarF.c,res.upward)
{
    sum.child <- res.upward$sum.child
    p.parent1 <- p.ybarF.c[parent1,,]
    p.parent2 <- p.post.found(parent2,status,probs,fyc)

    u <- array(0,dim=c(rep(2,3),rep(length(probs$p)+1,3)))
    for(c.parent1 in 1:(length(probs$p)+1)) for(c.parent2 in 1:(length(probs$p)+1))
    {
        p.child <- p.post.child(connect,c.parent1,c.parent2,status,probs,fyc)
        prod.bro <- 1
        if(length(bro.connect)>0) prod.bro <- prod(sum.child[bro.connect,c.parent1,c.parent2])
        for(s1 in 1:2) for(s2 in 1:2) for(si in 1:2)
        u[s1,s2,si,c.parent1,c.parent2,] <- prod.bro*p.parent1[s1,c.parent1]*p.parent2[s2,c.parent2]*p.child[si,]
    }
    p.ybarF.c[connect,,] <- apply(u,c(3,6),sum)
    p.ybarF.c
}

