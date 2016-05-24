
findOverlap = function(areas1, areas2, debug.level = 1, ...) {
  params = getRtopParams(...)
  olim = params$olim
  partialOverlap = params$partialOverlap
  ndim = length(areas1@polygons)
  t0 = proc.time()[[3]]
  if (!missing(areas2))  {
    mdim = length(areas2@polygons)
    sym = FALSE
  } else {
    mdim = ndim
    areas2 = areas1
    sym = TRUE
  }
  overlap = matrix(0,nrow = ndim,ncol = mdim)
  iareas = lapply(areas1@polygons, FUN = function(x) x@area)
  jareas = lapply(areas2@polygons, FUN = function(x) x@area)
  if (partialOverlap && requireNamespace("rgeos")) {
    for (ia in 1:(ndim-sym)) {
      ifi = ifelse(sym,ia+1,1)
      for (ib in ifi:mdim) {
        tmp = rgeos::gIntersection(areas1[ia,],areas2[ib,])
        if (!is.null(tmp) && inherits(tmp, "SpatialPolygons")) {
          oarea = sum(sapply(tmp@polygons, FUN = function(x) x@area))
          if (oarea > min(iareas[[ia]], jareas[[ib]]) * olim) {
            overlap[ia,ib] = oarea
            if (sym) overlap[ib,ia] = overlap[ia,ib]
          }
        }
      }
    }
    if (sym) diag(overlap) = iareas
  } else {
    if (partialOverlap) stop("partialOverlap of areas requested, but rgeos is not installed")
    ndim = length(areas1@polygons)
    t0 = proc.time()[[3]]
    ptdim = 25
    pts1 = SpatialPoints(coordinates(areas1))
    if (!missing(areas2))  {
      mdim = length(areas2@polygons)
      pts2 = SpatialPoints(coordinates(areas2))
      sym = FALSE
    } else {
      mdim = ndim
      areas2 = areas1
      pts2 = pts1
      sym = TRUE
    }
  
    plist2 = list()
    for (ib in 1:mdim) {
      poly2 = areas2@polygons[[ib]]
      plist2[[ib]] = SpatialPolygons(list(poly2)) 
    }
    nnover = 0
    t1 = proc.time()[[3]]
    if (debug.level > 1) print(t1-t0)
    for (ia in 1:(ndim-sym)) {
      t1 = proc.time()[[3]]
      poly1 = areas1@polygons[[ia]]
      SP1 = SpatialPolygons(list(poly1))
      a1 = SP1@polygons[[1]]@Polygons[[1]]@area
#    pls_SP1 = slot(SP1, "polygons")
#    lb_SP1 <- lapply(pls_SP1, function(x) as.double(bbox(x)))

      pt1 = pts1[ia,]
      ifi = ifelse(sym,ia+1,1)
      t2 = proc.time()[[3]]
      for (ib in ifi:mdim) {
        SP2 = plist2[[ib]]
        a2 = areas2@polygons[[ib]]@Polygons[[1]]@area
   #   pls_SP2 = slot(SP2, "polygons")
#      lb_SP2 <- lapply(pls_SP2, function(x) as.double(bbox(x)))
        if (max(unlist(commonArea(SP1,SP2))) > 0.1) {
          if (a2 < a1) {
            pt2 = pts2[ib,]
            nover = over(pt2, SP1)
            if (!is.na(nover)) overlap[ia,ib] = min(a1,a2)
            nnover = nnover + 1
          } else {
            nover = over(pt1, SP2)
            if (!is.na(nover)) overlap[ia,ib] = min(a1,a2)
            nnover = nnover + 1
          }
          if (sym) overlap[ib,ia] = overlap[ia,ib]
        }
      }
      t3 = proc.time()[[3]]
      if (debug.level > 1) print(paste(ia,round(sqrt(a1),2),round(t2-t1,3), round(t3-t2,3), round(t3-t0,3)))
    }
    if (debug.level > 1) print(paste("nnover",nnover)) 
    if (sym) diag(overlap) = a1
  }
  overlap
}





findVarioOverlap = function(vario) {
  overlap = function(a1,a2,dist) {
    ad = sqrt(a1)/2
    ad[2] = sqrt(a2)/2
    if (ad[1]+ad[2] > dist) {
      Srl = list()
      for (i in 1:2) {
        pt1 = c(0,ifelse(i==1,0,dist))
        x1 = pt1[1]-ad[i]
        x2 = pt1[1]+ad[i]
        y1 = pt1[2]-ad[i]
        y2 = pt1[2]+ad[i]
        boun = data.frame(x=c(x1,x2,x2,x1,x1),y=c(y1,y1,y2,y2,y1))
        Srl[[i]] = Polygon(SpatialPoints(boun))
      }
      cArea = commonArea(Srl[[1]],Srl[[2]])
     } else cArea = 0
     cArea[[1]]*a1
   }
   mapply(FUN = overlap,vario$a1,vario$a2,vario$dist)
}



bbArea = function(bb) {
  xd = bb[[3]]-bb[[1]]
  yd = bb[[4]]-bb[[2]]
  abs(xd) * abs(yd)
}
commonArea = function(objecti,objectj) {
  bi = bbox(objecti)
  bj = bbox(objectj)
  iarea = bbArea(bi)
  jarea = bbArea(bj)
  sdim = sqrt((iarea+jarea)/2)
  bl = list()
  for (i in 1:2) bl[[i]] =  max(bi[[i]],bj[[i]])
  for (i in 3:4) bl[[i]] =  min(bi[[i]],bj[[i]])
  if (bl[[3]] >= bl[[1]] & bl[[4]] >= bl[[2]]) {
    larea = bbArea(bl)
  } else {
    larea = 0
  }
  ilarea = larea/iarea
  jlarea = larea/jarea
  return(list(ilarea,jlarea))
}

#bOver = findOverlap(observations)
