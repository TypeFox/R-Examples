Rc <- function(x, ...){
  UseMethod("Rc")
}  


base.Rc <-function(data,  closest.coords,   type="curve"){
  
  data <- as.matrix(data)
  if (missing(closest.coords)) {
            stop("closest.coords needs to be provided.")
      }
  if (closest.coords[1]=="none"){
       stop("The provided vector of fitted values is empty.")
  }   
  
  n <- dim(data)[1]
  d <- dim(data)[2]

  # if weights argument copy in preamble here.
  
  all.dist <- vecdist(data, closest.coords)
  
  if (type=="curve"){
    pc.dist <- rep(0, n)
    pc <- princomp(data)
        # princomp does not support weights; so no observation weights possible here.
        # would need to implement by hand using cov.wt
    t.max <- max(apply(data, 2, function(dat) {
        diff(range(dat))
    }))
    t.line <- seq(-t.max, t.max, length = 1001)
    pc.line <- matrix(0, length(t.line), d)
    for (t in 1:length(t.line)) {
        pc.line[t, ] <- pc$center + pc$loadings[, 1] * t.line[t]
    }
    
    for (i in 1:n) {
        pc.dist[i] <- mindist(pc.line, data[i, ])$mindist
    }
    Ac <- mean(abs(all.dist))/mean(abs(pc.dist))
    # Ac<- weighted.mean( abs(all.dist), w=weights)/weighted.mean(abs(pc.dist), w=weights)
  } else if (type=="points"){
     m.dist <- rep(0, n)
     m <- colMeans(as.data.frame(data))
     m.dist <- sqrt(d)* distancevector(data, m)   
     Ac <- mean(abs(all.dist))/mean(abs(m.dist))
     # Ac<-weighted.mean(abs(all.dist, w=weights)/weighted.mean(abs(m.dist), w=weights)
  }
    Rc <- 1 - Ac
    return(Rc)
}
  


Rc.ms <- function(x,...){
  base.Rc(x$data, x$cluster.center[x$closest.label,], type="points")
}  

  



Rc.lpc <- function(x,...){
  object<- x
  if (class(object)=="lpc"){
     data <-object$data
     closest.coords <- lpc.spline(object, project=TRUE)$closest.coords
     # weights <- object$weights
  } else if (class(object)=="lpc.spline"){
     # weights <-object$lpcobject$weights
     if (object$Rc=="none"){
        data <- object$lpcobject$data
        if (object$closest.coords[1]=="none"){
                 closest.coords<-lpc.spline(object$lpcobject, project=TRUE)$closest.coords
        } else {
           closest.coords<-object$closest.coords
        }           
    } else {
        return(object$Rc)
    }
  }   else {
      stop("invalid object class.")
  }
   
  R <- base.Rc(data, closest.coords, type="curve")
  #  R <- Rc(data, closest.coords, weights, type="curve")
  return(R)
  }





Rc.lpc.spline <- Rc.lpc
