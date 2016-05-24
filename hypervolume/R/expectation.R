randsphere <- function(m,n,r)
{
  X = matrix(data=rnorm(m*n),nrow=m,ncol=n)
  s2 = rowSums(X * X)
  
  X = X * matrix(rep(r*(pgamma(s2/2,n/2)^(1/n))/sqrt(s2),n),ncol=n)
  
  return(X)
}

nball_volume <- function(n, r)
{
  return(pi^(n/2) * r^n / gamma(n/2+1))
}

expectation_maximal <- function(input, ...)
{
  if(class(input)=="Hypervolume")
  {
    return(input)
  }
  else
  {
    return(hypervolume(input, name=sprintf("Maximal expectation for %s", deparse(substitute(input))[1]), ...))
  }
}

expectation_box <- function(input, npoints=NULL, userandom=FALSE)
{
  if (class(input)=="Hypervolume")
  {
    if (userandom==TRUE)
    {
      data <- input@RandomUniformPointsThresholded
    }
    else 
    {
      data <- input@Data
    }
  }
  else
  {
    data <- input
    
    
  }
  
  minv <- apply(data,2,min)
  maxv <- apply(data,2,max)
  
  volume = prod(maxv - minv)
  
  if (is.null(npoints))
  {
    if (class(input)=="Hypervolume")
    {
      npoints = input@PointDensity * volume
    }
    else 
    {
      npoints <- 10 * 10^ncol(input)
    }
    cat(sprintf('Choosing npoints=%.0f (use a larger value for more accuracy.)\n',npoints))
  }
  
  density = npoints / volume
  
  result <- matrix(NA, nrow=npoints, ncol=length(minv), dimnames=list(NULL, dimnames(data)[[2]]))
  
  for (i in 1:length(minv)) # across each dimension
  {
    result[,i] <- runif(npoints, minv[i], maxv[i])
  }
  
  hv_box <- new("Hypervolume",
                Name=sprintf("Box expectation for %s", ifelse(class(input)=="Hypervolume", input@Name, deparse(substitute(data))[1])),
                Data = as.matrix(data),
                RandomUniformPointsThresholded=result, 
                Dimensionality=ncol(data), 
                Volume=volume,
                PointDensity=density, 
                DisjunctFactor=NaN,
                Bandwidth= NaN, 
                RepsPerPoint=floor(npoints / nrow(data)), 
                QuantileThresholdDesired=0, 
                QuantileThresholdObtained=0,
                ProbabilityDensityAtRandomUniformPoints = rep(1, npoints)
  )  
  
  return(hv_box)
}


expectation_ball <- function(input, npoints=NULL, userandom=FALSE)
{
  if (class(input)=="Hypervolume")
  {
    if (userandom==TRUE)
    {
      data <- input@RandomUniformPointsThresholded
    }
    else 
    {
      data <- input@Data
    }
  }
  else
  {
    data <- input
    
    
  }
  
  center <- apply(data,2,mean,na.rm=T)
  data_centered <- sweep(data, 2, center, "-")
  radii = sqrt(rowSums(data_centered^2))
  maxradius = max(radii)
  
  volume = nball_volume(n=ncol(data),r=maxradius)
  
  if (is.null(npoints))
  {
    if (class(input)=="Hypervolume")
    {
      npoints = input@PointDensity * volume
    }
    else 
    {
      npoints <- 10 * 10^ncol(data)
    }
    cat(sprintf('Choosing npoints=%.0f (use a larger value for more accuracy.)\n',npoints))
  }
  
  density = npoints / volume
  
  result <- randsphere(m=npoints, n=ncol(data), r=maxradius)
  result <- sweep(result, 2, center,"+")
  dimnames(result) <- list(NULL, dimnames(data)[[2]])
  
  hv_ball <- new("Hypervolume",
                 Name=sprintf("Ball expectation for %s", ifelse(class(input)=="Hypervolume", input@Name, deparse(substitute(data))[1])),
                 Data = as.matrix(data),
                 RandomUniformPointsThresholded=result, 
                 Dimensionality=ncol(data), 
                 Volume=volume,
                 PointDensity=density, 
                 DisjunctFactor=NaN,
                 Bandwidth= NaN, 
                 RepsPerPoint=floor(npoints / nrow(data)), 
                 QuantileThresholdDesired=0, 
                 QuantileThresholdObtained=0,
                 ProbabilityDensityAtRandomUniformPoints = rep(1, npoints)
  )  
  
  return(hv_ball)
}



expectation_convex <- function(input, npoints_inhull=NULL, npoints_onhull=NULL, check_memory=TRUE, userandom=FALSE)
{  
  if (class(input)=="Hypervolume")
  {
    if (userandom==TRUE)
    {
      data <- input@RandomUniformPointsThresholded
    }
    else 
    {
      data <- input@Data
    }
  }
  else
  {
    data <- input
    
  }  
  
  if (is.null(npoints_onhull))
  {
    npoints_onhull <- min(nrow(data),floor(10^sqrt(ncol(data))))
    cat(sprintf('Choosing npoints_onhull=%.0f (use a larger value for more accuracy.)\n',npoints_onhull))
  }
  
  numconvexhullpoints <- min(nrow(data),npoints_onhull)
  
  data_reduced <- data[sample(nrow(data), numconvexhullpoints, prob=rowSums(scale(data, center=TRUE, scale=FALSE)^2)),]
  
  if (ncol(data) > 5)
  {
    warning(sprintf("Algorithm may be very slow on high dimensional data (n>5: here, n=%d)",ncol(data)))
  }
  
  # FIND THE CONVEX HULL of the reduced data  
  convexhull <- geometry::convhulln(data_reduced,options="FA")
  hull_matrix <- convexhull$hull #extract coordinates of vertices
  hull_volume <- convexhull$vol # extract volume
  
  cat(sprintf("Convex hull calculated with %.0f simplices.\nCalculation of inequality constraints will require allocation of %.0f double-precision numbers.\n", nrow(hull_matrix), nrow(hull_matrix)^2))
  
  if (check_memory)
  {
    stop('Set check_memory=F to continue.\n')
  }	
  
  # make two proposals for npoints_inhull
  proposal_npoints_inhull_hv <- Inf
  proposal_npoints_inhull_nohv <- Inf
  
  if (class(input)=="Hypervolume" && is.null(npoints_inhull))
  {
    proposal_npoints_inhull_hv = input@PointDensity * hull_volume
    cat(sprintf('Density match proposal: npoints_inhull=%.0f\n',proposal_npoints_inhull_hv))
  }  

  proposal_npoints_inhull_nohv = floor(1000*10^sqrt(ncol(data)))
  cat(sprintf('Auto match proposal: npoints_inhull=%.0f\n',proposal_npoints_inhull_nohv))    

  if (is.null(npoints_inhull))
  {
    npoints_inhull <- min(proposal_npoints_inhull_hv, proposal_npoints_inhull_nohv)
    cat(sprintf('Choosing npoints_inhull=%.0f (use a larger value for more accuracy.)\n',npoints_inhull)) 
  }

  # REJECTION SAMPLING FROM BOUNDING BOX TO FILL IN CONVEX POLYTOPE APPROACH  
  ntrials <- 1000 # number of candidate points to propose at a time (algorithm always works regardless of this choice)
  done = FALSE
  inpoints <- NULL
  niter <- 0
  
  # generate a bounding box, and test if each point is in the convex hull or not
  while(!done)
  {
    # try a set of test points
    testpoints <- expectation_box(data_reduced, npoints=ntrials)
    
    chullin <- (inhull_compiled(testpts= testpoints@RandomUniformPointsThresholded, calpts=data_reduced, hull=hull_matrix)==1)
    # figure out which are 'in'
    inpoints_temp <- testpoints@RandomUniformPointsThresholded[chullin,]
    
    niter <- niter + 1
    inpoints <- rbind(inpoints, inpoints_temp)
    
    done = nrow(inpoints) >= npoints_inhull
    
    cat(sprintf('Rejection sampling: iteration %.0f - %.0f / %.0f points accepted\n', niter, nrow(inpoints), npoints_inhull))
  }
  
  dimnames(inpoints) <- list(NULL,dimnames(data)[[2]])
  
  # MAKE a new hypervolume with the convex hull shape and volume
  hv_chull <- new("Hypervolume",
                  Name=sprintf("Convex expectation for %s", ifelse(class(input)=="Hypervolume", input@Name, deparse(substitute(data))[1])),
                  Data=as.matrix(data),
                  RandomUniformPointsThresholded=inpoints, 
                  Dimensionality=ncol(inpoints), 
                  Volume=hull_volume,
                  DisjunctFactor=NaN,
                  PointDensity = nrow(inpoints) / hull_volume,
                  Bandwidth= NaN, 
                  RepsPerPoint=NaN, 
                  QuantileThresholdDesired=0, 
                  QuantileThresholdObtained=0,
                  ProbabilityDensityAtRandomUniformPoints = rep(1, nrow(inpoints))
  )
  
  return(hv_chull)
}