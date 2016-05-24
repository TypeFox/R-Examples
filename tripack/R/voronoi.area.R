voronoi.area <- function(voronoi.obj)
{
  ## Compute the area of each Voronoi polygon.
  ## If the area of a polygon cannot be computed, NA is returned.
  ##
  ## TODO: currently, the list of Voronoi vertices (vs) of each site
  ## is found, but then discarded.  They could be reused for other
  ## calls?
  
  nsites <- length(voronoi.obj$tri$x)
  areas <- double(nsites)
  for (i in 1:nsites) {
    vs <- voronoi.findvertices(i, voronoi.obj)
    if (length(vs) > 0) {
      areas[i] <- voronoi.polyarea( voronoi.obj$x[vs], voronoi.obj$y[vs])
    } else {
      areas[i] <- NA
    }
  }
  areas
}


voronoi.findvertices <- function(site, vor) {
  ## Helper function.
  ## Return the ordered list of Voronoi vertices for site number SITE
  ## in the Voronoi tesselation.
  
  p <- cbind(vor$p1, vor$p2, vor$p3)
  a <- which(p == site, arr.ind=TRUE)
  vertices <- a[,1]                     #list of the vertice indexes.
  triples <- p[a[,1],]
  triples
  ## Now remove the entries that are not site.

  ## Need to take transpose, as `which' runs down by column, rather
  ## than by row, and we want to keep rows together.
  triples <- t(triples)
  pairs <- triples[ which (triples!= site)]
  m <- matrix(pairs, ncol=2, byrow=TRUE)

  ## Now go through the list of sites and order the vertices.  We
  ## build up the list of vertices in the vector `orderedvs'.  This
  ## vector is truncated to the exact size at the end of the function.


  ## To order the vertices of the Voronoi polygon associated with a
  ## site, we first find all vertices that are associated with a site.
  ## These will come in threes, from the array `triples'.  We then
  ## remove the site number itself from the triples to come up with a
  ## list of pairs.  e.g. trying to find the vertices for site 6:
  
  ## sites     v number
  ## 3 9 6     6
  ## 6 4 3     2
  ## 9 6 7     3
  ## 6 7 4     9
  ##
  ## remove the `6':
  ## sites   v number
  ## 3 9     6
  ## 4 3     2
  ## 9 7     3
  ## 7 4     9
  
  ## and then starting with site 3, we find each subsequent site.
  ## i.e. 3 then 9 (output v 6), then 7 (output v 3), then 4 (output v
  ## 9) then 3 (output v 2).  We are now back to the starting site so
  ## the ordered list of vertices is 6, 3, 9, 2.
  
  orderedvs <- integer(30); vnum <- 1
  orderedvs[vnum] <- vertices[1]; vnum <- 1 + vnum
  firstv <- m[1,1];   nextv <- m[1,2]; m[1,] <- -1; #blank 1st row out.
  looking <- TRUE
  while (looking) {
    ##cat(paste("looking for ", nextv, "\n"))
    t <- which(m == nextv, arr.ind=TRUE)
    if (length(t) == 0) {               #could check length(t) != 1
      ## cannot compute area...
      vnum <- 1; looking <- FALSE
    } else {
      t.row <- t[1,1]
      t.col <- t[1,2]
      orderedvs[vnum] <- vertices[t.row]; vnum <- 1 + vnum
      othercol <- (3 - t.col)            #switch 1 to 2 and vice-versa.
      nextv <- m[ t.row, othercol]
      m[t.row,] <- -1                    #blank this row out.
      if (nextv == firstv) looking <- FALSE
    }
  }

  orderedvs[1:vnum-1]                   #truncate vector to exact length.
}

  
voronoi.polyarea <- function (x, y)
{
  ## Return the area of the polygon given by the points (x[i], y[i]).
  ## Absolute value taken in case coordinates are clockwise.
  ## Taken from the Octave implementation.
  ## Helper function.
  r <- length(x)
  p <- matrix(c(x, y), ncol=2, nrow=r)
  p2 <- matrix( c(y[2:r], y[1],  -x[2:r], -x[1]), ncol=2, nrow=r)
  a <- abs(sum (p * p2 ) / 2)
}
