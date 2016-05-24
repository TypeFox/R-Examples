.potentialKernelTypes <- c("EDKernel", "GKernel", "EKernel", "TriangleKernel", "VarGKernel")

.potential_depths <- function(ddalpha, objects, class = 0){
  
  if (ddalpha$ignoreself) stop("ignoreself not supported")  # todo ignore only if objects were NULL, @ what class do the points belong to? when seperately scaled
  
  # count all potential on the untransformed data
  if (class == 0) { 
    a = ddalpha$kernel.bandwidth
    data <- NULL
    cardinalities <- c()
    for (i in 1:ddalpha$numPatterns){
      data <- rbind(data, ddalpha$patterns[[i]]$points)
      cardinalities <- c(cardinalities, ddalpha$patterns[[i]]$cardinality)
    }
    if (is.null(objects))
      objects = data
  }
  # count potential w.r.t the given class
  else {
    a = ddalpha$kernel.bandwidth[class]
    data <- ddalpha$patterns[[class]]$transformer( ddalpha$patterns[[class]]$points )
    cardinalities <- c(ddalpha$patterns[[class]]$cardinality)      
    objects <- ddalpha$patterns[[class]]$transformer(objects)
    
    #   if (is.null(objects)){
    #     for (i in (1:ddalpha$numPatterns)){
    #        objects <- rbind(objects, ddalpha$patterns[[i]]$points)
    #     }
    #      objects <- ddalpha$patterns[[class]]$transformer(objects)
    #   } else {
    #     ignoreself = F
    #    }
  }
  
  kernelType = ddalpha$kernel
  if (is.character(kernelType))
    kernelType =  switch (kernelType, EDKernel = 1, GKernel = 2, EKernel = 3, TriangleKernel = 4, 1)
  
  points <- as.vector(t(data))
  numPoints <- sum(cardinalities)
  dimension <- ncol(data)
  points2 <- as.vector(t(objects))
  numPoints2 <- nrow(objects)
  classes <- length(cardinalities)
  
  depth <- .C("PotentialDepthsCount", 
              as.double(points), 
              as.integer(numPoints), 
              as.integer(dimension), 
              as.integer(classes),
              as.integer(cardinalities), 
              as.double(points2), 
              as.integer(numPoints2), 
              as.integer(kernelType), 
              as.double(a), 
              as.integer(ddalpha$ignoreself),
              depth=double(numPoints2*classes))$depth
  
  d = as.matrix(depth)
  dim(d)<-c(numPoints2,classes)
  
  return(d)
}

.potential_depths_wrt <- function(ddalpha, objects_){
  
  d <- NULL
  
  # w.r.t. each class
  for (cls in 1:ddalpha$numPatterns){
    
    data <- ddalpha$patterns[[cls]]$transformer( ddalpha$patterns[[cls]]$points )
    cardinalities <- c(ddalpha$patterns[[cls]]$cardinality)
    
    objects = objects_
    if (is.null(objects_))  # count for all data
      for (i in 1:ddalpha$numPatterns){
        objects <- rbind(objects, ddalpha$patterns[[i]]$points)
      }
    
    objects <- ddalpha$patterns[[cls]]$transformer(objects)
    
    kernelType = ddalpha$kernel
    if (is.character(kernelType))
      kernelType =  switch (kernelType, EDKernel = 1, GKernel = 2, EKernel = 3, TriangleKernel = 4, 1)
    
    points <- as.vector(t(data))
    numPoints <- sum(cardinalities)
    dimension <- ncol(data)
    points2 <- as.vector(t(objects))
    numPoints2 <- nrow(objects)
    classes <- length(cardinalities)
    
    depth <- .C("PotentialDepthsCount", 
                as.double(points), 
                as.integer(numPoints), 
                as.integer(dimension), 
                as.integer(classes),
                as.integer(cardinalities), 
                as.double(points2), 
                as.integer(numPoints2), 
                as.integer(kernelType), 
                as.double(ddalpha$kernel.bandwidth[cls]),
                as.integer(ddalpha$ignoreself),
                depth=double(numPoints2*classes))$depth
    
    d = cbind(d, depth)
  }
  
  return(d)
}

depth.potential <- function(x, data, pretransform = "1Mom", kernel = "GKernel", kernel.bandwidth = NULL, mah.parMcd = 0.75){
  if (!is.numeric(x)){
    stop("Argument \"x\" should be numeric")
  }
  if (!is.matrix(x)) {
    if(is.vector(x))
      x <- matrix(x, nrow=1)
    if(is.data.frame(x))
      x <- data.matrix(x)
  }
  if (!(is.matrix(data) && is.numeric(data)
        || is.data.frame(data) && prod(sapply(data, is.numeric))) 
      || ncol(data) < 2){
    stop("Argument \"data\" should be a numeric matrix of at least 2-dimensional data")
  }
  if(is.data.frame(data))
    data <- data.matrix(data)
  if (ncol(x) != ncol(data)){
    stop("Dimensions of the arguments \"x\" and \"data\" should coincide")
  }
  if (ncol(data) + 1 > nrow(data)){ #?
    stop("To few data points")
  }
  
  if(!is.null(pretransform)){
    if (pretransform == "1Mom" || pretransform == "NMom")
      mm <- mah.moment(data)
    else if (pretransform == "1MCD" || pretransform == "NMCD")
      mm <- mah.mcd(data, mah.parMcd)              
    
    transformer <- MahMomentTransformer(mm$mu, mm$b)
    data <- transformer(data)
    x <- transformer(x)
  }
  
  kernelType = kernel
  if (is.character(kernelType))
    kernelType =  switch (kernelType, EDKernel = 1, GKernel = 2, EKernel = 3, TriangleKernel = 4, 1)
  
  if (is.null(kernel.bandwidth)) {    # use the rule of thumb
    kernel.bandwidth = nrow(data) ^ (-2/(ncol(data)+4))
  }
  else{
    if (length(kernel.bandwidth) != 1 || is.na(kernel.bandwidth) || kernel.bandwidth == 0)
      stop("Argument \"kernel.bandwidth\" has invaid format.")
  }
  
  points <- as.vector(t(data))
  numPoints <- nrow(data)
  dimension <- ncol(data)
  points2 <- as.vector(t(x))
  numPoints2 <- nrow(x)
  cardinalities = numPoints
  classes <- 1
  
  ignoreself = F
  
  depth <- .C("PotentialDepthsCount", 
              as.double(points), 
              as.integer(numPoints), 
              as.integer(dimension), 
              as.integer(classes),
              as.integer(cardinalities), 
              as.double(points2), 
              as.integer(numPoints2), 
              as.integer(kernelType), 
              as.double(kernel.bandwidth),
              as.integer(ignoreself),
              depth=double(numPoints2*classes))$depth
  
  return(depth)
}




depth.space.potential <- function(data, cardinalities, pretransform = "NMom", kernel = "GKernel", kernel.bandwidth = NULL, mah.parMcd = 0.75){
  if (!(is.matrix(data) && is.numeric(data)
        || is.data.frame(data) && prod(sapply(data, is.numeric))) 
      || ncol(data) < 2){
    stop("Argument \"data\" should be a numeric matrix of at least 2-dimensional data")
  }
  if(is.data.frame(data))
    data <- data.matrix(data)
  if (!is.vector(cardinalities, mode = "numeric") 
      || is.na(min(cardinalities)) 
      || sum(.is.wholenumber(cardinalities)) != length(cardinalities) 
      || min(cardinalities) <= 0 
      || sum(cardinalities) != nrow(data)){
    stop("Argument \"cardinalities\" should be a vector of cardinalities of the classes in \"data\" ")
  }
  if (sum(cardinalities < ncol(data) + 1) != 0){
    stop("Not in all classes sufficiently enough objetcs")
  }
  
  d <- NULL
  
  needtransform = F
  
  if (pretransform == "1Mom" || pretransform == "1MCD"){
    if (pretransform == "1Mom")
      mm <- mah.moment(data)
    else   # "1MCD"
      mm <- mah.mcd(data, mah.parMcd)              
    
    transformer <- MahMomentTransformer(mm$mu, mm$b)
    data <- transformer(data)
    
    if (is.null(kernel.bandwidth)) {    # use the rule of thumb
      kernel.bandwidth = nrow(data) ^ (-2/(ncol(data)+4))
    }
    else{
      if (length(kernel.bandwidth) != 1 || is.na(kernel.bandwidth) || kernel.bandwidth == 0)
        stop("Argument \"kernel.bandwidth\" has invaid length, Zero or NA elements.")
    }
    
  } else if (pretransform == "NMom" || pretransform == "NMCD"){
    needtransform = T
    
    if (is.null(kernel.bandwidth)) {    # use the rule of thumb
      #separately calculated later
    }
    else{
      if (!is.numeric(kernel.bandwidth)
          ||!(is.vector(kernel.bandwidth) || is.list(kernel.bandwidth))){
        stop("Argument \"kernel.bandwidth\" has invaid format.")
      }
      if (length(kernel.bandwidth) == 1)
        kernel.bandwidth = rep(kernel.bandwidth, length(cardinalities))
      if (sum(!is.na(kernel.bandwidth)) != length(cardinalities) || sum(kernel.bandwidth != 0) != length(cardinalities)){
        stop("Argument \"kernel.bandwidth\" has invaid length, Zero or NA elements.")
      } 
    }
  }
  
  if(needtransform)
    # w.r.t. each class
    for (cls in 1:length(cardinalities)){
      pattern <- data[(1 + sum(cardinalities[0:(cls - 1)])):sum(cardinalities[1:cls]),]
      if(is.null(kernel.bandwidth)) band <- NULL
      else band <- kernel.bandwidth[cls]
      
      depth <- depth.potential(data, pattern, pretransform, kernel, band, mah.parMcd)
      d <- cbind(d, depth)
    }
  else # not w.r.t.
  {
    points <- as.vector(t(data))
    numPoints <- sum(cardinalities)
    dimension <- ncol(data)
    points2 <- as.vector(t(data))
    numPoints2 <- nrow(data)
    classes <- length(cardinalities)
    
    kernelType = kernel
    if (is.character(kernelType))
      kernelType =  switch (kernelType, EDKernel = 1, GKernel = 2, EKernel = 3, TriangleKernel = 4, 1)
    
    ignoreself = F
    
    depth <- .C("PotentialDepthsCount", 
                as.double(points), 
                as.integer(numPoints), 
                as.integer(dimension), 
                as.integer(classes),
                as.integer(cardinalities), 
                as.double(points2), 
                as.integer(numPoints2), 
                as.integer(kernelType), 
                as.double(kernel.bandwidth), 
                as.integer(ignoreself),
                depth=double(numPoints2*classes))$depth
    
    d = as.matrix(depth)
    dim(d)<-c(numPoints2,classes)
  }
  
  return(d)
}
