
#########################################################
# extract urban coordinates from land-cover data
#########################################################

coordinate.list <- function(raster, cell.class, compare=""){
  
  if(compare=="g"){compare <- function(v){return(which(v > cell.class))}}
  else{ if(compare=="s"){compare <- function(v){return(which(v < cell.class))}}
        else{compare <- function(v){return(which(v %in% cell.class))}}
  }
  
  
  all <- NULL
  coordinates <- NULL
  
  bs <- blockSize(raster)
  
  print("get coordinates:")
  # create progress bar
  pb <- txtProgressBar(min = 0, max = bs$n, style = 3)
    
  for (i in 1:bs$n) {
      v <- getValues(raster, row=bs$row[i], nrows=bs$nrows[i] )
      v <- compare(v)
      
      if(length(v)>0){
        coordinates <- xyFromCell(object=raster, cell= (raster@ncols * (bs$row[i]-1)) +v)
      }
      all <- rbind(all, coordinates)
      coordinates <- NULL
      setTxtProgressBar(pb, i)
    }
  
  close(pb)
  
  return(all)
}

#########################################################
# cluster urban coordinates
#########################################################


ccaR.order <- function(m) { m <- m[order(-m[,2]),]; m <- m[order(m[,1]),]; return(m); }

ccaRm <- function(data="", d=500, res.x=NULL, res.y=NULL, cell.class, compare)
{
  
  # If raster get coordinate list
  if(class(data)=="RasterLayer"){
    res.x <- xres(data)
    res.y <- yres(data)
    data <- coordinate.list(data, cell.class, compare)
  }
  
  # preprocessing
  n1 <- length(data[,2]);
  data <- ccaR.order(data);
  data <- data[order(data[,2]),];
  print("Sorting... Done");
  data[,1] <- data[,1]*pi/180;
  data[,2] <- data[,2]*pi/180;
   
  # calculate exit condition/distance
  max_abs <- max(abs(data[,2]));
  min_dist <- acos(sin(max_abs)*sin(max_abs)+cos(max_abs)*cos(max_abs)*cos(res.y*pi/180))*6371000;
  min_dist_h <- acos(sin(0)*sin(res.x*pi/180)+cos(0)*cos(res.x*pi/180)*cos(0))*6371000;
 
  # create place holder for cluster id
  data <- cbind(data, rep(0,length(data[,1])))
  data <- as.matrix(data)
  data <- as.numeric(data[,1:3]);
 
  step_w <- floor(1+(d/min_dist));
  step_h <- ceiling(d/min_dist_h);
  print("Start Clustering...");
  out <- .C("ccaRevT", m=as.numeric(data), n=as.integer(n1), l=as.numeric(d),
              step_w=as.integer(step_w), step_h=as.integer(step_h),
              res_x=as.numeric(res.x*pi/180), res_y=as.numeric(res.y*pi/180),
              w=as.integer(array(0,n1)));
  rm(data)
  print("Clustering... Done");
  mat1 <- matrix(out$m, ncol=3, byrow=FALSE);
  m1 <- array(0, max(mat1[,3]));
  m1 <- .C("ccaSumT", m=as.numeric(mat1), m3=as.integer(mat1[,3]), mm=as.numeric(m1), n=as.integer(n1));
  print("Summary... Done");
   
  mat1[,1] <- mat1[,1]*180/pi;
  mat1[,2] <- mat1[,2]*180/pi;
  mat1 <- as.data.frame(mat1)
  colnames(mat1) <- c("long","lat","cluster_id")
  return(list(cluster=mat1,size=m1$mm));
}


ccaRd <- function(data="", d=500, cell.class, compare)
{
  
  # If raster get coordinate list
  if(class(data)=="RasterLayer"){
     data <- coordinate.list(data, cell.class, compare)
  }
  
  # preprocessing
  n1 <- length(data[,2]);
  #data <- ccaR.order(data);
  data <- data[order(data[,1]),];
  print("Sorting... Done");
  
  # create place holder for cluster id
  data <- cbind(data, rep(0,length(data[,1])))
  data <- as.matrix(data)
  data <- as.numeric(data[,1:3]);
  
  print("Start Clustering...");
  out <- .C("ccaRev", m=as.numeric(data), n=as.integer(n1), l=as.numeric(d),w=as.integer(array(0,n1)));
  rm(data)
  print("Clustering... Done");
  mat1 <- matrix(out$m, ncol=3, byrow=FALSE);
  m1 <- array(0, max(mat1[,3]));
  m1 <- .C("ccaSum", m=as.numeric(mat1), m3=as.integer(mat1[,3]), mm=as.numeric(m1), n=as.integer(n1));
  print("Summary... Done");
  
  mat1 <- as.data.frame(mat1)
  colnames(mat1) <- c("long","lat","cluster_id")
  return(list(cluster=mat1,size=m1$mm));
}

#########################################################
# population assignment to clusters 
#########################################################

assign.data <- function(cluster, points, dist=1000){

  # conversion from meter to degree
  wgs84 <- 6371000
  scale <- 360 / (2*pi*wgs84)
  y_dist <- dist * scale
  
  points <- cbind(points, cluster_id=rep(0, length(points[, 1])))
  
  pop_cluster <- NULL
  
  for(j in 1:length(points[, 1])){
    #print(j)
    x <- points$long[j]
    y <- points$lat[j]
    x_dist <- dist * scale * (1 / cos( y*pi / 180 ))
    temp <- which((x-x_dist) < cluster[, 1] &
                    cluster[, 1] < (x+x_dist) &
                    (y-y_dist) < cluster[, 2] &
                    cluster[, 2] < (y+y_dist))
    if(length(temp)>0){
      pdist <- pointDistance(c(x, y), matrix(c(cluster[temp, 1],cluster[temp, 2]),
                                             ncol=2, byrow=FALSE), longlat=TRUE)
      if(min(pdist)< dist){
        index <- which(pdist[]==min(pdist))
        cluster_id <- cluster[temp[index], 3]
        points$cluster_id[j] <- cluster_id
        }
    }
  }
  return(points)
}

#########################################################
# buffer
#########################################################

osc.buffer <- function(input, width)
{	
  if(class(input)=="RasterLayer"){
  m <- matrix(input[], nrow=input@nrows, byrow=TRUE)
  returnraster <- T}else{
    m <- input
    returnraster <- F
  }
	m1 <- .C("ccaBuffED", m=as.integer(m), nr=as.integer(dim(m)[1]), nc=as.integer(dim(m)[2]), sz=as.integer(width))
	m1 <- m1$m
#	m1[which(m1<0)] <- -1
	m1 <- matrix(m1, nrow=dim(m)[1])
	if(returnraster){
	  input[] <- raster(m1)[]
	  return(input)
	  }else{
	    return(m1)
	  }
}
