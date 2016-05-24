cshell <- function (x, centers, iter.max = 100, verbose = FALSE,
                    dist = "euclidean", method = "cshell",
                    m=2, radius= NULL) 
{
    x <- as.matrix(x)  
    xrows <- dim(x)[1]
    xcols <- dim(x)[2]
    xold <- x
    perm <- sample(xrows)
    x <- x[perm, ]
    ## initial values are given
    if (is.matrix(centers)) 
        ncenters <- dim(centers)[1]
    else {
        ## take centers random vectors as initial values
        ncenters <- centers
        centers <- x[rank(runif(xrows))[1:ncenters], ]+0.001
    }

    ## initialize radius
    if (missing(radius))
        radius <- rep(0.2,ncenters)
    else
        radius <- as.double(radius)

    dist <- pmatch(dist, c("euclidean", "manhattan"))
    if (is.na(dist)) 
        stop("invalid distance")
    if (dist == -1) 
        stop("ambiguous distance")
  
    method <- pmatch(method, c("cshell"))
    if (is.na(method)) 
        stop("invalid clustering method")
    if (method == -1) 
        stop("ambiguous clustering method")
    
    initcenters <- centers
    ## dist <- matrix(0, xrows, ncenters)
    ## necessary for empty clusters
    pos <- as.factor(1 : ncenters)
    rownames(centers) <- pos
    iter <- integer(1)

    flag <- integer(1)

  
    retval <- .C("cshell",
                 xrows = as.integer(xrows),
                 xcols = as.integer(xcols), 
                 x = as.double(x),
                 ncenters = as.integer(ncenters), 
                 centers = as.double(centers), 
                 iter.max = as.integer(iter.max),
                 iter = as.integer(iter), 
                 verbose = as.integer(verbose),
                 dist = as.integer(dist-1), 
                 U = double(xrows*ncenters),
                 UANT = double(xrows*ncenters),
                 m = as.double(m),
                 ermin = double(1),
                 radius = as.double(radius),
                 flag = as.integer(flag),
                 PACKAGE = "e1071")

    centers <- matrix(retval$centers, ncol = xcols, dimnames = dimnames(initcenters))
  
  
    radius <- as.double(retval$radius)
    U <- retval$U
    U <- matrix(U, ncol=ncenters)
    UANT <- retval$UANT
    UANT <- matrix(UANT, ncol=ncenters)

    iter <- retval$iter
    flag <- as.integer(retval$flag)
  
    ## Optimization part
    while (((flag == 1) || (flag==4)) && (iter<=iter.max)) {
    
        flag <- 3
    
        system <- function (spar=c(centers,radius), x, U, m, i) {
            k <- dim(x)[1]
            d <- dim(x)[2]
            nparam<-length(spar)
            
            v<-spar[1:(nparam-1)]
            r<-spar[nparam]
            
            ##distance matrix x_k - v_i
            distmat <- t(t(x)-v)
            
            ##norm from x_k - v_i
            normdist <- distmat[,1]^2
            for (j in 2:d)
                normdist<-normdist+distmat[,j]^2
            normdist <- sqrt(normdist)
            
            ##equation 5
            op <- sum( (U[,i]^m) * (normdist-r) )^2
            ##equation 4
            equationmatrix <- ((U[,i]^m) * (1-r/normdist))*distmat
            ## <FIXME KH 2005-01-14>
            ## This had just apply(), but optim() really needs a scalar
            ## fn. 
            ## What do we really want here?
            op<- op+sum(apply(equationmatrix, 2, sum)^2)
            ## </FIXME>
      
        }
    
        for (i in 1:ncenters) {
            spar <- c(centers[i,],radius[i])
            npar <- length(spar)
      
            optimres <- optim(spar, system, method="CG", x=x, U=U, m=m, i=i)
            centers[i,] <- optimres$par[1:(npar-1)]
            radius[i] <- optimres$par[npar]
        }
    
    
        retval <- .C("cshell",
                     xrows = as.integer(xrows),
                     xcols = as.integer(xcols), 
                     x = as.double(x),
                     ncenters = as.integer(ncenters), 
                     centers = as.double(centers), 
                     iter.max = as.integer(iter.max),
                     iter = as.integer(iter-1), 
                     verbose = as.integer(verbose),
                     dist = as.integer(dist-1), 
                     U = as.double(U),
                     UANT = as.double(UANT),
                     m = as.double(m),
                     ermin = double(1),
                     radius = as.double(radius),
                     flag = as.integer(flag),
                     PACKAGE = "e1071")
    
        flag<-retval$flag
        if (retval$flag!=2)
            flag<-1
    
    
        centers <- matrix(retval$centers, ncol = xcols,
                          dimnames = dimnames(initcenters))
    
        radius <- as.double(retval$radius)
        U <- retval$U
        U <- matrix(U, ncol=ncenters)
        UANT <- retval$UANT
        UANT <- matrix(UANT, ncol=ncenters)
    
        iter <- retval$iter
    }
  
    centers <- matrix(retval$centers, ncol = xcols,
                      dimnames = list(pos, colnames(initcenters)))
    
    U <- matrix(retval$U, ncol = ncenters,
                dimnames = list(rownames(x), 1 : ncenters))
    U <- U[order(perm),]  
    clusterU <- apply(U, 1, which.max)
    
    clustersize <- as.integer(table(clusterU))
    radius <- as.double(retval$radius)
    
    retval <- list(centers = centers, radius=radius,
                   size = clustersize, cluster = clusterU,
                   iter = retval$iter - 1, membership=U,
                   withinerror = retval$ermin,
                   call = match.call())
    
    class(retval) <- c("cshell", "fclust")
    return(retval)
}
  

#predict.cshell <- function( clobj, x){
  
#  xrows<-dim(x)[1]
#  xcols<-dim(x)[2]
#  ncenters <- clobj$ncenters
#  cluster <- integer(xrows)
#  clustersize <- integer(ncenters)
#  f <- clobj$m
#  radius <- clobj$radius

#  if(dim(clobj$centers)[2] != xcols){
#    stop("Number of variables in cluster object and x are not the same!")
#  }

  
#  retval <- .C("cshell_assign",
#               xrows = as.integer(xrows),
#               xcols = as.integer(xcols),
#               x = as.double(x),
#               ncenters = as.integer(ncenters),
#               centers = as.double(clobj$centers),
#               dist = as.integer(clobj$dist-1),
#               U = double(xrows*ncenters),
#               f = as.double(f),
#               radius = as.double(radius))

  

#  U <- retval$U
#  U <- matrix(U, ncol=ncenters)

#  clusterU <- apply(U,1,which.max)
#  clustersize <- as.integer(table(clusterU))
     

#  clobj$iter <- NULL
#  clobj$cluster <- clusterU
#  clobj$size <- retval$clustersize
#  clobj$membership <- U
  
#  return(clobj)
#}
