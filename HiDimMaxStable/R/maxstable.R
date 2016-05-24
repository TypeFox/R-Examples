
prod.mubz<-function(partition,z,...) {
    n.blocks<-max(partition)
    prod(Vectorize(function(j) mubz(b=partition==j,z=z,...))(1:n.blocks))
}

sumprod.mubz<-function(z,...) {
    partitions<-setparts(length(z))
    sum(apply(partitions,2,function(part)
              prod.mubz(part,z,...)))
}

maxstable.l<-function(z,ln=FALSE,...) {
    if(ln) {
        log(sumprod.mubz(z,...))-v.star(length(z),z,...)
    } else {
        sumprod.mubz(z,...)*exp(-v.star(length(z),z,...))
    }
}

maxstable.l.clusters<-function(data,clusters=rep(1,dim(data)[1]),
                               ln=FALSE,spatial=NULL,...) {
    n.obs<-dim(data)[2]
    if(ln) {
        v<-0
    } else {
        v<-1
    }
    for(j in 1:n.obs) {
        z<-data[,j]
        lc<-Vectorize(function(i) {
            if(sum(clusters==i)<=1) {
                ifelse(ln,0,1)
            } else {
                spatial.selected<-NULL
                if(!is.null(spatial)) {
                    spatial.selected<-spatial
                    spatial.selected$sites<-spatial.selected$sites[clusters==i,]
                }
                maxstable.l(z[clusters==i],spatial=spatial.selected,ln=ln,...)
            }
        })(1:max(clusters))
        if(ln) {
            v<-v+sum(lc)
        } else {
            v<-v*prod(lc)
        }
    }
    v
}

build.clusters.spatial<-function(xy,max.size=5,plot=FALSE) {
    n.site<-dim(xy)[1]

    # utilise kmeans pour regrouper les points en clusters. On
    # commence avec un nombre de clusters n.site/max.size, et on augmente ce
    # nombre tant que la taille d'un des clusters dépasse max.size
    
    n.clust<-max(1,floor(n.site/max.size))
    cl <- kmeans(xy, n.clust, nstart=10)
    while (max(cl$size) > max.size) {
        n.clust<-n.clust+1
        cl <- kmeans(xy, n.clust, nstart=10)
    }

    # représentation graphique

    if(plot) {
        colnames(xy) <- c("x", "y")
        # points avec couleur en fonction du cluster
        plot(xy, col = cl$cluster)
        # centres des clusters
        points(cl$centers, col = 1:n.clust, pch = 8, cex=1)
    }
    
    cl$cluster
}

dens.grid.maxstable<-function(...) {
    dens.grid(f=maxstable.l.clusters,...)
}

maxlik.maxstable<-function(...) {
    maxlik.f(f=maxstable.l.clusters,...)
}

forme.csv<-function(...) {
    paste(c(...),sep=",",collapse=",")
}

simul.maxlik.maxstable<-function(i=0,
                                 params=c(0.5,1.5),n.obs=20,
                                 n.site=50,xy.size=2,
                                 results.csv=NULL,
                                 model="whitmat",
                                 model.f=switch(model,
                                     "whitmat"=spatialWhittleMatern,
                                     "cauchy"=spatialCauchy,
                                     "powexp"=spatialPowerexp,
                                     "brown"=spatialPower
                                     ),
                                 category=switch(model,
                                     "whitmat"="normal",
                                     "cauchy"="normal",
                                     "powexp"="normal",
                                     "brown"="lnormal"
                                     ),
                                 iter.max=200,reltol=0.000001,
                                 xy=NULL,
                                 ...) {
    if(is.null(results.csv)) {
        results.csv<-paste("MCMLcomp-",model,"-",n.site,"-",n.obs,"-",xy.size,"-",params[1],"-",params[2],".csv",sep="")
    }

    if(is.null(xy)) {
        xy<-matrix(runif(2 * n.site, 0, xy.size), ncol = 2)
    } else {
        n.site<-dim(xy)[1]
    }
    if(model=="brown"){
        data<-t(SpatialExtremes::rmaxstab(n.obs, xy, model,
                         range = params[1], smooth = params[2]))

        f<-SpatialExtremes::fitmaxstab(t(data),xy,model,
                      start=list(range=params[1],smooth=params[2]),
                      control=list(maxit=iter.max,reltol=reltol),...)
    } else {
        data<-t(SpatialExtremes::rmaxstab(n.obs, xy, model,
                         nugget = 0, range = params[1], smooth = params[2]))

        f<-SpatialExtremes::fitmaxstab(t(data),xy,model,
                      start=list(range=params[1],smooth=params[2]),nugget=0,
                      control=list(maxit=iter.max,reltol=reltol),...)
    }
    f.res<-c(f$param["range"],f$param["smooth"],f$convergence,f$counts["function"])
 
    cl<-build.clusters.spatial(xy)
    ml<-maxlik.maxstable(data,c(NA,NA),start=params,iterlim=iter.max,reltol=reltol,
                          category=category,
                          spatial=list(sites=xy,family=model.f),
                          clusters=cl)

    if(!is.null(results.csv)) {
        if(!file.exists(results.csv)) {
            cat(forme.csv("model","n.site","xy.size","n.obs","range","smooth","est.range","est.smooth","convergence","iterations","est.range(pairwise)","est.smooth(pairwise)","convergence","iterations"),
                "\n",sep="",file=results.csv)
        }
        cat(forme.csv(model,n.site,xy.size,n.obs,params,ml$estimate,ml$message,ml$iterations,f.res),
            "\n",sep="",file=results.csv,append=TRUE)
    }

}

simul.maxlik.maxstable.runs<-function(n=100,parallel=TRUE,...) {
    if(parallel) {
        snowfall::sfSapply(1:n,function(i)
            simul.maxlik.maxstable(i,...))
    } else {
        for(i in 1:n) {
            simul.maxlik.maxstable(...)
        }
    }
    TRUE
}
