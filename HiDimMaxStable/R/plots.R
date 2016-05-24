sfOuter<-function(x,y,f,...) {
    params<-matrix(c(rep(x,length(y)),rep(y,rep(length(x),length(y)))),ncol=2)
    results<-snowfall::sfApply(params,1,function(p,...) f(p[1],p[2],...),...)
    matrix(results,ncol=length(y))
}

dens.grid<-function(f,data,params,seqx,seqy,
                    ln=TRUE,parallel=FALSE,...) {
    check.parallel(parallel)
    if(parallel) {
        dens<-sfOuter(seqx,seqy,
                      Vectorize(function(x,y,data,ln,params,...) {
                          p<-params
                          p[is.na(p)]<-c(x,y)
                          f(data,ln=ln,params=p,...)
                      }, vectorize.args=c("x","y")),
                      data=data,ln=ln,params=params,...)
    } else {
        dens<-outer(seqx,seqy,Vectorize(function(x,y) {
            p<-params
            p[is.na(p)]<-c(x,y)
            f(data,ln=ln,params=p,...)
        }))
    }
    dg<-list(params=params,seqx=seqx,seqy=seqy,dens=dens)
    class(dg)<-'densgrid'
    dg
}

persp.densgrid<-function(x,...) {
    persp(x$seqx,x$seqy,x$dens,...)
}

plot3d.densgrid<-function(dg,...) {
    perspcolor(dg$seqx,dg$seqy,dg$dens,
               xlab="param_1",ylab="param_2",zlab="Log-density",...)
}

maxgrid<-function(dg) {
    p<-dg$params
    l<-length(dg$seqx)
    ii<-which.max(dg$dens)
    i<-((ii-1) %% l)+1
    j<-((ii-1) %/% l)+1
    pmax<-c(dg$seqx[i],dg$seqy[j])
    p[is.na(p)]<-pmax
    list(p=p,xy=pmax,value=dg$dens[i,j])
}

perspcolor<-function(x,y,z,n.colors=50,tight=TRUE,...) {
    if (requireNamespace("rgl", quietly = TRUE)) {
    
        nrz <- nrow(z)
        ncz <- ncol(z)

        if(tight) {
            color <- rep(terrain.colors(n.colors),n.colors:1)
            n.colors<-sum(n.colors:1)
        } else {
            color <- terrain.colors(n.colors)
        }
    
        zfacet <- z[-1, -1] + z[-1, -ncz] + z[-nrz, -1] + z[-nrz, -ncz]
    
        facetcol <- cut(zfacet, n.colors,labels=F)
        col <- rbind(1, cbind(matrix(color[facetcol],
                                     nrz-1, ncz-1), 1))
        rgl::persp3d(x, y, z, theta=130, phi=25, expand=1, col=col,
                     ticktype="detailed",
                     axes=TRUE,smooth=FALSE,lit=FALSE,...)
    } else {
        stop("Package rgl must be installed.")
    }
}

