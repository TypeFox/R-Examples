infoTheoreticGCM <- function(g, dist=NULL, coeff="lin", infofunct="sphere", lambda=1000, custCoeff=NULL, alpha=0.5, prec=53, flag.alpha=FALSE){
  if (prec > 53)
    require("Rmpfr")

  allowed.functionals <- c("sphere","pathlength","vertcent","degree")

  # check if g is a graphNEL object
  if(class(g)[1]!="graphNEL")
    stop("'g' must be a 'graphNEL' object")
  if(is.null(dist))
    dist <- distanceMatrix(g)

  # assign the coefficients
  ci <- .infoTheoreticCoeff(coeff, custCoeff, max(dist))

  # check the information functional
  if (!infofunct %in% allowed.functionals) {
    stop(paste("unknown information functional '",infofunct,"'; must be one of ", paste(allowed.functionals,collapse=", "), sep=""))
  }

  # actually calculate the functional
  nan.message <- "check your parameters"
  if (infofunct == allowed.functionals[1]) {          # "sphere"
    fvi <- .functionalJSphere(g, dist=dist, ci=ci, alpha=alpha, flag=flag.alpha)
  } else if (infofunct == allowed.functionals[2]) {   # "pathlength"
    fvi <- .functionalPathlength(g, dist=dist, ci=ci)
  } else if (infofunct == allowed.functionals[3]) {   # "vertcent"
    fvi <- .functionalLocalProperty(g, dist=dist, ci=ci)
  } else if (infofunct == allowed.functionals[4]) {   # "degree"
    fvi <- .functionalDegreeDegree(g, dist=dist, ci=ci, alpha=alpha, prec=prec)
    nan.message <- "please try higher values for 'prec'"
  }

  .infoTheoretic(fvi, lambda, nan.message)
}  

.functionalJSphere <- function(g, dist, ci, alpha, flag){
  #calculate Spheres
  nam <- nodes(g)
  Sj <- lapply(nam,function(n){
    table(dist[n,],exclude=0)
  })
  names(Sj) <- nam
  
  fvi <- sapply(Sj,function(s,pc=ci){
     sum(s*ci[1:length(s)])
  })
  names(fvi) <- nam
  if(flag==FALSE){
    return (fvi)
  }else{
    return (0.5^fvi)
  }
}


.functionalPathlength <- function(g, dist,ci){
  ig <- .G2IG(g)
  vs <- igraph::V(ig)
  lvs <- length(vs$name)
  fvi <- rep(0,lvs)
  nam <- nodes(g)
  #determine number of all possible shortest path
  for(n in 1:lvs){
    #f <- vs[vs$name==n] 
    asp <- igraph::get.all.shortest.paths(ig,from=n)$res
    lvi <- table(sapply(asp,length)-1,exclude=0)
    fvi[n] <- sum(lvi*ci[1:length(lvi)])
  }
  names(fvi) <- nam
  return (fvi)
}

.functionalLocalProperty <- function(g, dist,ci){

  ig <- .G2IG(g)
  vs <- igraph::V(ig)
  lvs <- length(vs$name)
  fvi <- rep(0,lvs)
  nam <- nodes(g)
  #determine number of all possible shortest path
  for(n in 1:lvs){
    #f <- vs[vs$name==n] 
    asp <- igraph::get.all.shortest.paths(ig,from=n)$res #changes between igraph0 and igraph >> $res
    lvi <- table(sapply(asp,length)-1,exclude=0)
    tmp.sum <- sapply(1:max(as.numeric(names(lvi))),function(lpl){
      sum(1:lpl)
    })
    betai <- 1/tmp.sum*lvi
    fvi[n] <- sum(betai*ci[1:length(lvi)])
  }
  names(fvi) <- nam
  return (fvi)
}

.functionalDegreeDegree <- function(g, dist, ci, alpha, prec){
  m <- adjacencyMatrix(g)
  deg <- graph::degree(g)
  size <- numNodes(g)

  expts <- .C("quacn_degdeg_exponents",
    as.integer(m), as.integer(deg), as.integer(size),
    as.double(ci), as.integer(length(ci)),
    double(size))[[6]]

  fvi <- alpha ^ if (prec > 53) {
    Rmpfr::mpfr(expts, prec)}
  else {
    expts
  }
  names(fvi) <- nodes(g)

  fvi
}
