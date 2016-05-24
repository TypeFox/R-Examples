trace.variog <-
function(coords, L2norm, bin=FALSE, max.dist, uvec="default", breaks="default", nugget.tolerance){

  # Argument validation
  if(is.null(coords)) stop("coords is not an optional parameter")
  if(ncol(coords)!=2) stop("coords must be an n x 2 matrix")
  if(!isSymmetric(L2norm)) stop("L2norm must be a symmetric matrix")
  if(sum(diag(L2norm))!=0) stop("each element of the diagonal of L2norm must zero")

  # Euclidian distance among sites
  Eu.d <- as.matrix(dist(coords, method="euclidean"))

  # The variogram values are the point to point differences between all 
  # the curves in L2norm. The variogram values are obtained by 
  # removing the upper part of matrix L2norm, excluding the diagonal and
  # sequantially passing the result by rows to an array
  vtemp <- array(as.dist(L2norm))

  # The distance values are obtained similarly to the variogram values
  utemp <- array(as.dist(Eu.d))

  # The maxim distance value is defined temporarily to be used 
  # in .define.bins and as a definite value max.dist in cloud
  if(missing(max.dist)){
    umax <- max(utemp)
  }else{
    umax <- max(utemp[utemp < max.dist])
  }

  # |bin|<---lag--->|bin|<---lag--->|bin|<---lag--->|bin|<---lag--->|bin|

  if(bin){
    # If there is no value for nugget.tolerance or 
    # if it is smaller than a given value
    # then nugget.tolerance is defined and 
    # it is indicated that it must be removed from the first bin value
    if(missing(nugget.tolerance) || nugget.tolerance < 1e-11){
      nugget.tolerance <- 1e-12
      nt.ind <- FALSE
    }else{
      if(mode(nugget.tolerance) != "numeric") stop("nugget.tolerance must be numeric")
      nt.ind <- TRUE
    }
    # If the minimum distance is minor than the nugget.tolerance then it is
    # indicated that the first bin should NOT be removed
    min.dist <- min(utemp)
    if(min.dist < nugget.tolerance){
      nt.ind <- TRUE
    }
    # The bins are defined
    dbins <- .define.bins(max.dist=umax, nugget.tolerance=nugget.tolerance, uvec=uvec, breaks=breaks)
    # The max.dist value is defined
    if(missing(max.dist)){
      max.dist <- max(dbins$bins.lim)
    }
    # The first bin is adjusted
    if(dbins$bins.lim[1] < 1e-16){
      dbins$bins.lim[1] <- 0
    }
    if(!nt.ind){
      dbins$uvec <- dbins$uvec[-1]
      dbins$bins.lim <- dbins$bins.lim[-1]
    }

    u <- rep(NA, length(dbins$bins.lim)-1)
    v <- rep(NA, length(dbins$bins.lim)-1)
    for(i in seq(1,length(dbins$bins.lim)-1)){
      # vector defining points belong which lag
      # ToDo: specify what happens if a point is exactly on the edge of the lag
      lagpoints <- dbins$bins.lim[i]<utemp & utemp<dbins$bins.lim[i+1]
      # If there is any point inside the lag then it is used in the variogram
      if(any(lagpoints)){
        u[i] <- mean(utemp[lagpoints])
        v[i] <- mean(vtemp[lagpoints])
      }
    }
    output.type <- "bin"
    emp.trace.vari <- list(bins.lim=dbins$bins.lim, nugget.tolerance=nugget.tolerance)
#    emp.trace.vari <- list(bins.lim=dbins$bins.lim, nugget.tolerance=nugget.tolerance, uvec=dbins$uvec)
  }else{
    u <- utemp
    v <- vtemp
    max.dist <- umax
    output.type <- "cloud"
    emp.trace.vari <- list()
  }

#  plot(u, v)

  # ToDo: what to do when bins are too small, to manage with nugget.tolerance

  # The trace.variogram object to be returned is created
  # ToDo: When the variogram type is cloud, 
  #       objects 'Eu.d' and 'L2norm' are returned, but they are augmented 
  #       versions of 'u' and 'v' objects.
  #       The Eu.d object is then used to compute the semivariogram range.
  #       Eu.d and L2norm can be obtained from u and v by doing a cycle
  #       iterating over all values
  # [a, b, c]
  # and creating a matrix with the structure 
  # 0 a b
  # a 0 c
  # b c 0
  emp.trace.vari <- c(emp.trace.vari, list(u=u, v=v, output.type=output.type, max.dist=max.dist, Eu.d=Eu.d, L2norm=L2norm))

  # 'variogram' class is assigned to the returned object. 
  # Then it can be used as any variomodel in geoR
  class(emp.trace.vari) <- "variogram"

  return(emp.trace.vari)
}

".define.bins" <-
  function(max.dist, uvec = "default", breaks = "default", nugget.tolerance)
  {
    if(all(breaks ==  "default")){
      if (all(uvec == "default")) uvec <- 13
      if (mode(uvec) == "numeric"){
        if(length(uvec) == 1){
          bins.lim <- seq(0, max.dist, l = uvec+1)
          bins.lim <- c(0, nugget.tolerance, bins.lim[bins.lim >  nugget.tolerance])
          uvec <- 0.5 * (bins.lim[-1] + bins.lim[-length(bins.lim)])
        }
        else{
          uvec <- c(0, uvec)
          nvec <- length(uvec)
          d <- 0.5 * diff(uvec[2:nvec])
          bins.lim <- c(0, (uvec[2:(nvec - 1)] + d), (d[nvec - 2] + uvec[nvec]))
          bins.lim <- c(0, nugget.tolerance, bins.lim[bins.lim >  nugget.tolerance])
        }
      }
      else stop("argument uvec can only take a numeric vector")
    }
    else{
      if(mode(breaks) != "numeric") stop("argument breaks can only take a numeric vector")
      else bins.lim <- breaks
      bins.lim <- c(0, nugget.tolerance, bins.lim[bins.lim >  nugget.tolerance])
      uvec <- 0.5 * (bins.lim[-1] + bins.lim[-length(bins.lim)])
    }
    return(list(uvec = uvec, bins.lim = bins.lim))
  }
