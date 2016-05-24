`earth.poly` <-
function(lats) {
    #tests if lats are a matrix or spatial points
    if (class(lats)=='SpatialPoints') lats<-coordinates(lats)
    lats<-unique(lats)
    rownames.lats<-rownames(lats)
    nr <- nrow(lats)
    #bearing matrix == bm == pairwise bearing matrix
    bm <- matrix(, nr, nr)
    for (i in 1:nr) {
        for (j in 1:nr) {
            if (i == j)
                bm[i, j] <- -1
            else bm[i, j] <- earth.bear(long1 = lats[i, 1], lat1 = lats[i, 2], long2 = lats[j, 1], lat2 = lats[j, 2])
        }
    }
    #quadrant matrix == qm == which quadrants the bearings to other points fall in to
    qm <- matrix(, nr, 4)
    for (i in 1:nr) {
        qm[i, 1] <- length(bm[1, ][bm[i, ] >= 0]) - length(bm[i,][bm[i, ] >= 90])
        qm[i, 2] <- length(bm[1, ][bm[i, ] >= 90]) - length(bm[i,][bm[i, ] >= 180])
        qm[i, 3] <- length(bm[1, ][bm[i, ] >= 180]) - length(bm[i,][bm[i, ] >= 270])
        qm[i, 4] <- length(bm[1, ][bm[i, ] >= 270])
    }
    id <- numeric(nr)
    #identifying which points have bearings in greater than 1 quadrant
    #those that have bearings in less than 2 quadrants must be vertices (ie <180 degrees)
    for (i in 1:nr) {
        for (j in 1:4) {
            if (length(qm[i, ][qm[i, ] > 0]) <= 2)
                id[i] <- i
            else {
                if (length(qm[i, ][qm[i, ] > 0]) == 3) {
                  if (length(qm[i, 2:3][qm[i, 2:3] > 0]) == 2) {
                    if ((max(bm[i, ]) - min(bm[i, ][bm[i, ] >=
                      0])) < 180)
                      id[i] <- i
                  }
                  else {
                    if ((360 - min(bm[i, ][bm[i, ] > 180]) +
                      max(bm[i, ][bm[i, ] < 180])) < 180)
                      id[i] <- i
                  }
                }
            }
        }
    }
    #mi == points on the perimeter
    mi <- id[id > 0]
    #the starting/primary point that will be used for the calculations
    prpt <- mi[1]
    #all points on the perimeter except the starting point
    snpts <- mi[-1]
    #the bearing for all the points from the primary point
    pt <- bm[prpt, snpts]
    #to get coordinates between -180 and 0 to between 180 and 360
    if ((max(pt) - min(pt)) > 180) {
        pt[pt < 180] <- pt[pt < 180] + 360
    }
    #ranked bearings, to give the points in the proper order so that the sequence of points is in order
    rpt <- rank(pt)
    ind <- numeric(length(pt))
    for (i in 1:length(pt)) ind[rpt[i]] <- snpts[i]
    ta <- numeric(1)
    for (i in 1:(length(pt) - 1)) {
        k <- i + 1
        ta <- ta + earth.tri(long1 = lats[prpt, 1], lat1 = lats[prpt, 2], long2 = lats[ind[i], 1], lat2 = lats[ind[i], 2], long3 = lats[ind[k], 1], lat3 = lats[ind[k], 2])
    }
    verts <- c(prpt, ind)
    return(list(area = ta, vertices = verts))
}
