fuse <- function(map, regions, name) {
  
  regions <- regions
  testmap <- map
  surroundingNames <-(attributes(map)$surrounding)
  surroundingNames2 <-unlist(attributes(map)$surrounding)
  ind <- c()
  indexrem <- NULL
  j <- 0
  regisurrounding <- c()
  for(i in 1:length(surroundingNames)) {
    
    if(length(surroundingNames[[i]])>0) {
      j <- j+1
      if(any(regions==surroundingNames2[j])) {
        if(any(regions==names(map)[i])) {
          
        } else {
          ind <- c(ind,j)
          regisurrounding <- c(regisurrounding,names(map)[i])
        }
      }
    }
  }

  polys <- list()
  id <- c()
  for(i in 1:length(regions)) {
    check <- which(names(map)==regions[i])
    for(l in 1:length(check)) {
      id <- c(id,paste("p",1+length(id),sep=""))
      if(any(indexrem==check[l])) {
        
      }else{
        indexrem <- c(indexrem,check[l])
      }
    }
    if(length(which(names(map)==regions[i]))==0) {
      warning(paste(regions[i]," not contained in map"),sep="")
    }
    
  }
  
  for(i in 1:length(indexrem)) {
    
    polys <- c(polys,sp::Polygons(list(sp::Polygon(rbind(map[[indexrem[i]]][dim(map[[indexrem[i]]])[1],],map[[indexrem[i]]]))),id[i]))
  }
  
  testmap <- map[-indexrem]
  ntemp <- length(testmap)
  #combine seleted regions to one single polygon, add this region to map again
  spatpol <- sp::SpatialPolygons(polys)
  combine <- (maptools::unionSpatialPolygons(spatpol,rep(1,length(slot(spatpol, "polygons")))))
  regionadd <- slot(combine,"polygons")
  regionadd <- (regionadd[[1]]@Polygons)
  for(i in 1:length(regionadd)) {
    testmap[[ntemp+i]] <- regionadd[[i]]@coords
    names(testmap)[ntemp+i] <- name
  }
  
  
  
  minima <- sapply(map, function(x) {
    apply(x, 2, min)
  })
  maxima <- sapply(map, function(x) {
    apply(x, 2, max)
  })
  minimum <- apply(minima, 1, min)
  maximum <- apply(maxima, 1, max)
  x.range <- maximum[1] - minimum[1]
  y.range <- maximum[2] - minimum[2]
  height2width <- round(y.range/x.range, digits = 2)
  surrounding <- replicate(n = length(testmap), expr = character())
  if(length(regisurrounding)>0) {
    for(i in 1:length(regisurrounding)) {
      ind2 <- which(names(testmap)==regisurrounding[i])
      surrounding[[ind2]] <- surroundingNames2[ind[i]]
    }
  }
  regi <- unique(names(testmap))
  
  attributes(testmap) <- list(names=names(testmap),height2width = height2width, class="bnd",
                          surrounding = surrounding, regions = regi)
  
  return(testmap)
  
}

