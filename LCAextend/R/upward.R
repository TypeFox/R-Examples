upward <-
function(id,dad,mom,status,probs,fyc,peel)
{
    p.yF.c <- matrix(1,nrow=length(id),ncol=length(probs$p)+1)
    sum.child <- array(0,c(length(id),length(probs$p)+1,length(probs$p)+1))
    for(generat in 1:peel$generation)
    {
        connects <- peel$peel.connect[generat,]
        connects <- connects[connects>0]
        for(connect in connects)
        {
            spouse.connect <- peel$couple[peel$couple[,1]==connect,2]
            children.connect <- union(id[dad==connect],id[mom==connect])
            res.upward.connect <- upward.connect(connect,spouse.connect,children.connect,status,probs,p.yF.c,fyc,sum.child)
            p.yF.c <- res.upward.connect$p.yF.c
            sum.child <- res.upward.connect$sum.child
        }
    }
    res <- list("p.yF.c"=p.yF.c,"sum.child"=sum.child)
    res
}

