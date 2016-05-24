distmatrix <- function(date, type="mindist", tolerance=0.1, useGW=T) {

  # check input
  if (!inherits(date, "Date")) {
    stop("date is not of type Date")
  }
  
  if (date < as.Date("1946-1-1") | date > as.Date("2015-6-30")) {
    stop("Specified date is out of range")
  }
  
  if (!(type %in% c("mindist", "capdist", "centdist"))) {
  	stop("Wrong type argument. Possible values: mindist, capdist, centdist")
  }
  
  if (tolerance<0) {
  	stop("Tolerance must be >=0")
  }
  
  # where to look for the dataset
  path <- paste(system.file(package = "cshapes"), "shp/cshapes.shp", sep="/")
  
  # load the dataset
  cshp.full <- readShapePoly(path, proj4string=CRS("+proj=longlat +ellps=WGS84"), IDvar="FEATUREID")
    
  # select all relevant polygons
  if (useGW) {
    cshp.full <- cshp.full[cshp.full$GWCODE>=0,]
    startdate <- as.Date(paste(cshp.full$GWSYEAR, cshp.full$GWSMONTH, cshp.full$GWSDAY, sep="-"))
    enddate <- as.Date(paste(cshp.full$GWEYEAR, cshp.full$GWEMONTH, cshp.full$GWEDAY, sep="-"))
  } else {
    cshp.full <- cshp.full[cshp.full$COWCODE>=0,]
    startdate <- as.Date(paste(cshp.full$COWSYEAR, cshp.full$COWSMONTH, cshp.full$COWSDAY, sep="-"))
    enddate <- as.Date(paste(cshp.full$COWEYEAR, cshp.full$COWEMONTH, cshp.full$COWEDAY, sep="-"))
  }	
  cshp.full <- cshp.full[!is.na(startdate) & !is.na(enddate),]
  startdate <- startdate[!is.na(startdate)]
  enddate <- enddate[!is.na(enddate)]
  cshp.part <- cshp.full[startdate <= date & enddate >= date,]
  
  # compute pairwise distances
  if (useGW) {
    cshp.part <- cshp.part[order(cshp.part$GWCODE),]
    ccodes <- cshp.part$GWCODE
  } else {
    cshp.simple <- cshp.part[order(cshp.part$COWCODE),]
    ccodes <- cshp.part$COWCODE
  }
    
  resultmatrix <- matrix(0, nrow=length(ccodes), ncol=length(ccodes))
  colnames(resultmatrix) <- ccodes
  rownames(resultmatrix) <- ccodes
  
  if (type=="mindist") {
    
    # simplify the polygons
    cshp.simple <- thinnedSpatialPoly(cshp.part, tolerance, minarea=0, avoidGEOS=T)
    
    for (c1 in 1:(length(ccodes)-1)) {
      for (c2 in (c1+1):length(ccodes)) {
        
        # compute distance
        dist <- cshp.mindist(cshp.simple[c1,], cshp.simple[c2,])
        resultmatrix[c1,c2] <- dist
        resultmatrix[c2,c1] <- dist 
      }
    }
  } else {
    
    for (c1 in 1:(length(ccodes)-1)) {
      for (c2 in (c1+1):length(ccodes)) {
  
        if (type=="capdist") {
          dist <- cshp.capdist(cshp.part[c1,], cshp.part[c2,])
        }
        if (type=="centdist") {
          dist <- cshp.centdist(cshp.part[c1,], cshp.part[c2,])
        }
        
        resultmatrix[c1,c2] <- dist
        resultmatrix[c2,c1] <- dist
      }
    }    
  }
  resultmatrix
}

cshp.mindist <- function(polygon1, polygon2) {
  
  # create matrices containing all points of the polygons
  p1 <- ldply(polygon1@polygons[[1]]@Polygons, function(y) {y@coords})
  p2 <- ldply(polygon2@polygons[[1]]@Polygons, function(y) {y@coords})
  
  # use spDists function to compute distances between all pairs of points
  min(spDists(as.matrix(p1), as.matrix(p2), longlat=T))
}

cshp.centdist <- function(polygon1, polygon2) {
  
  # use spDists function to compute distances between centroids
  spDistsN1(coordinates(polygon1), coordinates(polygon2), longlat=T)
}

cshp.capdist <- function(polygon1, polygon2) {
  
  # use spDists function to compute distances between centroids
  spDistsN1(matrix(c(polygon1$CAPLONG, polygon1$CAPLAT), ncol=2), matrix(c(polygon2$CAPLONG, polygon2$CAPLAT), ncol=2), longlat=T)
}
