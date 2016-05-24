constantDensitySampling <- function(x, polygon.id='pID', ...) {
  
  # retain proj4 information
  p4s <- proj4string(x)
  
  ## TODO: enable parallel processing
  ## must load packages in each cluster
  #   cl <- makeCluster(getOption("cl.cores", 2))
  #   res <- parLapply(cl=cl, X=slot(x, 'polygons'), fun=sample.by.poly, n.pts.per.ac=n.pts.per.ac, min.samples=min.samples, p4s=p)
  #   stopCluster(cl)
  
  # sample and return a list, one element / valid polygon
  res <- lapply(slot(x, 'polygons'), sample.by.poly, p4s=p4s, ...)
  
  # this happens when there aren't enough sample points based on min.samples
  # check for NULL in this list-- cases where it was too difficult to place a point
  null.items <- which(sapply(res, is.null))
  if(length(null.items) > 0) {
    message('some polygons were too small for the requested number of samples')
    res <- res[-null.items]
  }
    
  # when rbind-ing the result there is an error related to duplicate rownames
  # reset with running counter
  total <- 0
  for(i in seq_along(res)) {
    
    # number of points
    n <- length(res[[i]])
    
    # new row names seq
    i.seq <- seq(from=(total + 1), to=(total + n))
    
    # reset rownames
    dimnames(res[[i]]@coords)[[1]] <- as.character(i.seq)
    
    # increment counter
    total <- total + n
  }
  
  # convert into a single SP object
  if(length(res) > 0)
    res <- do.call('rbind', res)
  
  # upgrade to SPDF with polygon ID
  pID.df <- sp::over(res, x)[, polygon.id, drop=FALSE]
  res <- SpatialPointsDataFrame(res, data=pID.df)

  return(res)
}


# sample by polygon, must be from a projected CRS
# p: Polygon object
# n: number of points per acre (results will be close)
# min.samples: minimum requested samples / polygon
# iterations: number of sampling "tries"
# p4s: proj4string assigned to SpatialPoints object
sample.by.poly <- function(p, n.pts.per.ac=1, min.samples=5, sampling.type='hexagonal', iterations=10, p4s=NULL) {
  # convert _projected_ units to acres
  ac.i <- p@area * 2.47e-4
  # determine number of points based on requested density
  n.samples <- round(ac.i * n.pts.per.ac)
  
  # polygon must be at least large enough to support requested number of samples
  if(n.samples >= min.samples) {
    # trap errors caused by bad geometry
    s.i <- try(spsample(p, n=n.samples, type=sampling.type, iter=iterations), silent=TRUE)
    
    # errors? return NULL
    if(class(s.i) == 'try-error')
      return(NULL)
    
    # assign original proj4 string
    if(!is.null(p4s) & !is.null(s.i))
      proj4string(s.i) <- p4s
    
    return(s.i)
  }
  
  # not enough samples, return NULL
  else
    return(NULL)
}
