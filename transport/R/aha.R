# Computes optimal transport or wasserstein distance between n*m matrices a and b
# where sum(a)==sum(b) must be satisfied. [NO THIS IS FIXED] The algorithm treats matrix a as a
# measure on [0,m-1]x[0,n-1] which is constant in each pixel [i,i+1)x[j,j+1),
# while matrix b is interpreted as a discrete measure. In particular,
# aha(a,a,wasser=T) != 0 if a!=0.
#
# nscales: determines the number of scales to be used
# scmult: determines the multiplier (w.r.t. the number of support points) to get from one scale to the next
#         in the case of nscales>1
# maxit/factr: control the quality of approximation, see 'optim'.
# wasser: if true, compute the wasserstein distance instead of optimal transport.
# wasser.spt: number of support points used to compute the discretisation of b.
aha <- function(a,b,nscales=1,scmult=2,factr=1e+05,maxit=10000,wasser=FALSE,wasser.spt=NA,approx=FALSE,...) {
    n <- dim(a)[1] # y
    m <- dim(a)[2] # x

    #seed <- floor(runif(1,0,10000))
    #print(seed)
    #set.seed(seed)

    if (all(c("x","y","mass") %in% names(b))) {
        stopifnot(all(b$x>=0) && all(b$y>=0) && all(b$x<=m) && all(b$y<=n))
        target_mass <- b$mass
        x <- b$x
        y <- b$y
    } else {
        stopifnot(dim(a)==dim(b))
        target_mass <- as.vector(b)
        x <- as.vector(mapply(function(k) rep(k,n),1:m))-0.5
        y <- rev(rep(1:n,m))-0.5
    }

    x <- x + runif(length(x),-1e-5,1e-5)
    y <- y + runif(length(y),-1e-5,1e-5)

    # permutate target measure points
    # for faster power diagram generation
    perm <- sample(1:length(x),prob=target_mass+1e-7)
    x <- x[perm]
    y <- y[perm]
    target_mass <- target_mass[perm]

    pixel_density <- as.vector(a)
    rect <- c(0,0,m,n)

    if (!is.na(wasser.spt)) {
        wasser <- TRUE
    }

    if (is.na(wasser.spt) || wasser.spt>length(x)) {
        wasser.spt <- length(x) 
    }

    pixel_density <- pixel_density/sum(pixel_density)
    target_mass <- target_mass/sum(target_mass)
    
    # uses lloyd's algorithm to determine n points (x,y,m) such that
    # sum(m)=sum(m0) and which are close to (x0,y0,m0) w.r.t. wasserstein distance
    # p will repesent the mapping between the points of (x,y) and (x0,y0), i.e.
    # p[i]=j means (x0[i],y0[i]) is mapped to (x[j],y[j])
    decompose <- function(n,x0,y0,m0) {
        n0 <- length(x0)
        if (n==n0) {
            return(list(x=x0,y=y0,m=m0,p=0:(n-1)))
        }

        # sample initial cluster coordinates from (x0,y0,m0)
        if (length(m0[m0>0])>=n0) {
            i <- sample.int(n0,size=n,prob=m0)
        } else {
            i <- sample.int(n0,size=n,prob=m0+1e-7)
        }

        return(.C("decompose_c",as.integer(n),x=as.double(x0[i]),y=as.double(y0[i]),
                       m=as.double(m0[i]),as.integer(n0),as.double(x0),as.double(y0),
                       as.double(m0),p=integer(n0),as.double(0.01),PACKAGE="transport"))
    }

    # recursively decompose the target measure and apply bfgs to get
    # optimal weight vectors
    multiscale <- function(n0,x0,y0,m0,depth) {
        f <- function (w) {
            #plot(power_diagram(x0,y0,w,rect))
            -.C("aha_phi", as.integer(n0), as.double(x0), as.double(y0), as.double(w), 
               as.double(pixel_density), as.double(m0), as.integer(!approx), res=double(1),PACKAGE="transport")$res
        }
        g <- function (w) {
            #plot(power_diagram(x0,y0,w,rect))
            -.C("aha_dphi", as.integer(n0), as.double(x0), as.double(y0), as.double(w),
               as.double(pixel_density), as.double(m0), as.integer(!approx), res=double(n0),PACKAGE="transport")$res
        }

        if (depth<nscales && n0>scmult) {
            v <- decompose(floor(n0/scmult),x0,y0,m0)
            w <- multiscale(floor(n0/scmult),v$x,v$y,v$m,depth+1)
            res <- optim(w[v$p+1],f,g,method="L-BFGS-B",control=list(maxit=maxit,factr=factr,...))
            #if (res$convergence!=0) { print(res$message) }
            return (res$par)
        } else {
            res <- optim(rep(0,n0),f,g,method="L-BFGS-B",control=list(maxit=maxit,factr=factr,...))
            #if (res$convergence!=0) { print(res$message) }
            return (res$par)
        }
    }

    .C("aha_init",as.integer(n),as.integer(m),as.double(rect),PACKAGE="transport")
    if (wasser) {
        v <- decompose(wasser.spt,x,y,target_mass)
        w <- multiscale(wasser.spt,v$x,v$y,v$m,1)
        #plot(power_diagram(v$x,v$y,w,c(0,0,m,n)))

        res <- .C("aha_wasserstein",as.integer(wasser.spt),as.double(v$x),as.double(v$y),as.double(w),
                  as.double(pixel_density), res=double(1),PACKAGE="transport")

        # error bound
        error <- sqrt(target_mass%*%((x-v$x[v$p+1])^2+(y-v$y[v$p+1])^2))

        .C("aha_free",PACKAGE="transport")

        return(data.frame(wasser.dist=res$res,error.bound=error))
    } else {
        w <- multiscale(length(x),x,y,target_mass,1)
        #plot(power_diagram(x,y,w,rect=rect))

        tmemsize <- .C("aha_compute_transport", as.integer(length(x)), as.double(x), as.double(y),
                       as.double(w), as.double(as.vector(a)), res=integer(1),PACKAGE="transport")$res
        res <- .C("aha_get_transport", as.integer(tmemsize), from=double(tmemsize), 
                  to=double(tmemsize), mass=double(tmemsize),PACKAGE="transport")[2:4]
        .C("aha_free",PACKAGE="transport")

        tp <- data.frame(from=1+res$from,to=perm[1+res$to],mass=res$mass)
        if (!("mass" %in% names(b))) {
            tp <- tp[tp$from!=tp$to,]
        }
        return(tp)
    }
}

transport_apply <- function(a,tplan) {
    n <- dim(a)[1]
    m <- dim(a)[2]
    av <- as.vector(a)
    for (i in seq(1,dim(tplan)[1])) {
        av[tplan$to[i]] <- av[tplan$to[i]]+tplan$mass[i]
        av[tplan$from[i]] <- av[tplan$from[i]]-tplan$mass[i]
    }
    return(matrix(av,n,m))
}

transport_error <- function(a,b,tplan) {
    if (all(c("x","y","mass") %in% names(b))) {
        return(sum(abs(aggregate(tplan$mass,by=list(tplan$to),sum)[2]-b$mass)))
    } else { 
        return(sum(abs(transport_apply(a,tplan)-b)))
    }
}

# Computes the power diagram of weigted points (x,y,w) in R^2,
# intersected with the rectangle 'rect', which defaults to
# rect=c(min(x),min(y),max(x),max(y))
power_diagram <- function (xi,eta,w,rect=NA) {

    stopifnot(length(xi)==length(eta),length(eta)==length(w))

    if (!identical(cbind(xi,eta),unique(cbind(xi,eta)))) {
        stop("input data must consist of distinct points")
    }

    if (is.na(rect[1])) {
        rect <- c(min(xi),max(xi),min(eta),max(eta))
    }

    if (length(rect)!=4) {
        stop("rectangle format is c(xmin,xmax,ymin,ymax)")
    }

    n <- length(xi)

    # get cells
    cell_sizes <- .C("compute_power_diagram", res = integer(n), as.integer(n), 
                     as.double(xi), as.double(eta), as.double(w), as.double(rect[c(1,3,2,4)]), 
                     PACKAGE="transport")$res
    memory <- sum(cell_sizes)
    res <- .C("get_power_diagram", as.integer(memory), x = double(memory), y = double(memory),
              PACKAGE="transport")

    # format cells
    cells <- as.list(rep(NA,n))
    j <- 1
    for (i in 1:n) {
        m <- cell_sizes[i]
        if (m>2) {
            cells[[i]] <- cbind(x=res$x[seq(j,j+m-1)],y=res$y[seq(j,j+m-1)])
            j <- j + m
        } else {
            cells[[i]] <- NA
        }
    }

    pd <- list(sites=data.frame(xi=xi,eta=eta,w=w), cells=cells, rect=rect)
    class(pd) <- c("power_diagram")
    return(pd)
}

plot.power_diagram <- function(x, weights=FALSE, ...) {
    stopifnot(class(x) == "power_diagram")
    pd <- x
    segmentize <- function(pg) {
      if (is.na(pg[1])) {
        return(matrix(0,0,4))
      } else {
        n <- dim(pg)[1]
        res <- cbind(pg[,1], pg[,2], pg[c(2:n,1),1], pg[c(2:n,1),2])
        return(res)
      }
    }

    mcircle <- function(x,y,r) {
      nx <- length(x)
      phi <- seq(0,2*pi,length.out=200)
      xer <- function(phi,r) {r*cos(phi)}
      yer <- function(phi,r) {r*sin(phi)}
      cx <- outer(phi,r,xer)
      cy <- outer(phi,r,yer)
      xx <- matrix(x,200,nx,byrow=TRUE) + cx
      yy <- matrix(y,200,nx,byrow=TRUE) + cy
      matplot(xx,yy,col=grey(0.5), type="l", lty=1, add=TRUE)
    }

	temp <- lapply(pd$cells, segmentize)
    temp2 <- do.call(rbind, temp)
    rect <- pd$rect
    plot(c(rect[1],rect[2]),c(rect[3],rect[4]),type="p",asp=1,axes=FALSE,xaxs="i",yaxs="i",xlab="",ylab="",pch="",...)
    sites <- pd$sites[!is.na(pd$cells),]
    points(sites[,1], sites[,2], pch=20)
    if (any(is.na(pd$cells))) { 
        hidden_sites <- pd$sites[is.na(pd$cells),]
        points(hidden_sites[,1], hidden_sites[,2], pch=20, col=grey(0.7))
    }
    segments(temp2[,1],temp2[,2],temp2[,3],temp2[,4],col="blue")
    if (weights) {
        mcircle(sites[,1], sites[,2], sqrt(sapply(sites[,3],function(w) max(0,w))))
    }
}
