
excess.l<-function(data,ln=FALSE,...) {
    n.obs<-dim(data)[2]
    m<-dim(data)[1]

    if(ln) {
        v<-0
    } else {
        v<-1
    }

    rb<-v.star(m=m,...)

    for(i in 1:n.obs) {
        z<-data[,i]
        ra<-mubz(b=z>1,z=z,...)

        if(ln) {
            v<-v+log(ra)-log(rb)
        } else {
            v<-v*ra/rb
        }
    }

    v
}

maxlik.excess<-function(...) {
    maxlik.f(f=excess.l,...)
}

# get the increasing ranks from the data vector x 
vector.ranks<-function(x) {
    r<-numeric(length(x))
    r[sort(x,index.return=T)$ix]<-1:length(x)
    r
}

# for each row from the matrix y, get the normalized ranks from the
# data row.
matrix.ranks<-function(y) {
    n.obs<-dim(y)[2]
    m<-dim(y)[1]

    for(i in 1:m) {
        y[i,]<-(vector.ranks(y[i,])-1)/n.obs
    }
    y
}

excess.censor<-function(z,t=10) {
    # Transform empirical distributions to unit Pareto
    z<-1/(1-matrix.ranks(z))
    #
    z<-z/t
    # Extracts rows with at least one coordinate over 1
    vanish<-apply(z<=1,2,all)
    z<-z[,!vanish,drop=FALSE]
    # Censor coordinates <=1 to 1
    z[z<1]<-1
    z
}

dens.grid.excess<-function(...) {
    dens.grid(f=excess.l,...)
}

# random deviates for Clustered max-stable distribution
rCMS<-function(copulas,margins,classes,params,n=100) {
    u<-matrix(NA,ncol=n,nrow=length(classes))
    for(i in 1:max(classes)) {
        d<-sum(classes==i)
        copi<-onacopulaL(copulas[[i]],list(params[i*2-1],1:d))
        u[classes==i,]<-margins[[i]]@q(t(
              rCopula(n,copi)
              ),params[2*i])
    }
    t( t(u)*1/runif(n) )
}

# Needs MASS for mvrnorm
rSchlatherExcess<-function(n=500,spatial,params) {
    cov.matrix<-2*pi*spatial.cor.matrix(c(params,1),spatial)
    n.site<-dim(spatial$sites)[1]
    -1/log(runif(n)) * MASS::mvrnorm(n,rep(0,n.site),cov.matrix)
}
