downward <-
function(id,dad,mom,status,probs,fyc,peel,res.upward)
{
    ww <- array(0,dim=c(length(id),rep(2,3),rep(length(probs$p)+1,3)))
    w <- array(0,dim=c(length(id),2,length(probs$p)+1))
    p.ybarF.c <- array(1,dim=c(length(id),2,length(probs$p)+1))

    generat <- peel$generation
    connect <- peel$peel.connect[generat,]
    connect <- connect[connect>0]
    p.ybarF.c[connect,,] <- p.post.found(connect,status,probs,fyc)
    spouse.connect <- peel$couple[peel$couple[,1]==connect,2]
    children.connect <- union(id[dad==connect],id[mom==connect])
    res.nuc <- weight.nuc(connect,spouse.connect,children.connect,status,probs,fyc,p.ybarF.c,ww,w,res.upward)
    ww <- res.nuc$ww
    w <- res.nuc$w

    while(generat>1)
    {
        generat <- generat-1
        connects <- peel$peel.connect[generat,]
        connects <- connects[connects>0]
        for(connect in connects)
        {
            parent1.connect <- intersect(peel$peel.connect[generat+1,],c(dad[id==connect],mom[id==connect]))
            parent2.connect <- setdiff(c(dad[id==connect],mom[id==connect]),parent1.connect)
            bro.connect <- union(id[dad==parent1.connect],id[mom==parent1.connect])
            bro.connect <- setdiff(bro.connect,connect)

            p.ybarF.c <- downward.connect(connect,parent1.connect,parent2.connect,bro.connect,status,probs,fyc,p.ybarF.c,res.upward)

            spouse.connect <- peel$couple[peel$couple[,1]==connect,2]
            children.connect <- union(id[dad==connect],id[mom==connect])
            res.nuc <- weight.nuc(connect,spouse.connect,children.connect,status,probs,fyc,p.ybarF.c,ww,w,res.upward)
            ww <- res.nuc$ww
            w <- res.nuc$w
        }
    }
    res <- list("ww"=ww,"w"=w)
    res
}

