
simultoccur.l<-function(data,occur,ln=FALSE,...) {
    n.obs<-dim(data)[2]
    if(ln) {
        v<-0
    } else {
        v<-1
    }
    for(j in 1:n.obs) {
        z<-data[,j]
        p<-occur[,j]
        lc<-prod.mubz(p,z,...)
        vs<-v.star(length(z),z,...)
        if(ln) {
            v<-v+log(lc)-vs
        } else {
            v<-v*lc*exp(-vs)
        }
    }
    v
}

maxlik.simultoccur<-function(...) {
    maxlik.f(f=simultoccur.l,...)
}

dens.grid.simultoccur<-function(...) {
    dens.grid(f=simultoccur.l,...)
}

to.classes<-function(x) {
    f<-factor(x)
    levels(f)<-1:max(as.numeric(levels(f)))
    as.numeric(f)
}

maxblocks<-function(y,n.blocks=50) {
    t.max<-dim(y)[2]
    n<-dim(y)[1]
    block.n<-1+floor( ((1:t.max)-1) * n.blocks/t.max)
    ym<-matrix(nrow=n,ncol=n.blocks)
    yi<-matrix(nrow=n,ncol=n.blocks)
    for(i in 1:n.blocks) {
        ym[,i]<-apply(y[,block.n==i],1,max)/sum(block.n==i)
        yi[,i]<-to.classes(apply(y[,block.n==i],1,which.max))
    }
    list(normalized.max=ym,classes.max=yi)
}

