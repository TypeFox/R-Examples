#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# this function was taken from Simon Wood's MGCV package
# the main difference between  this and and the equivalent function in 
# BayesX is that one point neighbours
# Polygons to neighbour function
polys2nb <- function(polys) 
{
  ## polys is a list of polygons. i.e. 
  ## polys[[i]] is 2 column matrix defining 
  ## polygons for ith area (NA separated). 
  ## The function returns list of neightbours 
  ## for each area.
  ## Bounding box speed up from a comment in spdep package help.
  ## WARNING: neighbours defined by sharing 
  ## vertices. So one having vertices on another's line-segment 
  ## is not detected! 
  ## MS what are the implications? MS
  ##-----------------------------------------------------------------------------
  # start  polys2nb
  n.poly <- length(polys) ## total numer of polygons
  ## work through list of list of polygons, computing bounding boxes
  a.ind <- p.ind <- lo1 <- hi1 <- lo2 <- hi2 <- rep(0,n.poly)
  k <- 0
  for (i in 1:n.poly) {
    ## bounding box limits...
    polys[[i]] <- polys[[i]][!is.na(rowSums(polys[[i]])),] ## strip NA's
    lo1[i] <- min(polys[[i]][,1])
    lo2[i] <- min(polys[[i]][,2])
    hi1[i] <- max(polys[[i]][,1])
    hi2[i] <- max(polys[[i]][,2])
    ## strip out duplicates
    polys[[i]] <- uniquecombs(polys[[i]])
  }
  ##---------------------------------------------------------------------------- 
  ## now work through finding neighbours....
  nb <- list() ## nb[[k]] is vector indexing neighbours of k
  for (i in 1:length(polys)) nb[[i]] <- rep(0,0) ## setting to zeros MS
  
  for (k in 1:n.poly) 
  { ## work through poly list looking for neighbours
    ol1 <- (lo1[k] <= hi1 & lo1[k] >= lo1)|(hi1[k] <= hi1 & hi1[k] >= lo1)|
      (lo1 <= hi1[k] & lo1 >= lo1[k])|(hi1 <= hi1[k] & hi1 >= lo1[k])
    ol2 <- (lo2[k] <= hi2 & lo2[k] >= lo2)|(hi2[k] <= hi2 & hi2[k] >= lo2)|
      (lo2 <= hi2[k] & lo2 >= lo2[k])|(hi2 <= hi2[k] & hi2 >= lo2[k])
    ol <- ol1&ol2
    ol[k] <- FALSE
    ind <- (1:n.poly)[ol] ## index of potential neighbours of poly k
    ## co-ordinates of polygon k...
    cok <- polys[[k]]
    if (length(ind)>0) for (j in 1:length(ind)) {
      co <- rbind(polys[[ind[j]]],cok) 
      cou <- uniquecombs(co)
      n.shared <- nrow(co) - nrow(cou)
      ## if there are common vertices add area from which j comes
      ## to vector of neighbour indices 
      if (n.shared>0) nb[[k]] <- c(nb[[k]],ind[j]) 
    }
  }
  for (i in 1:length(polys)) nb[[i]] <- unique(nb[[i]])
  names(nb) <- names(polys)
  list(nb=nb,xlim=c(min(lo1),max(hi1)),ylim=c(min(lo2),max(hi2)))
}
##------------------------------------------------------------------------------
## polys2nb ends here
##------------------------------------------------------------------------------
##------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
## neibourhood  to precision
##------------------------------------------------------------------------------
## This function it needed for checking especially as at the moment do not check the 
## neighbour information
## FDB: This functions needs information about x: the factor containing the areas, then we can have levels(k)
nb2prec <- function(neighbour,x,area=NULL)
{
  if (!is(x, "factor")) stop("x must be a factor")  
  k <- area
  if (is.null(k))
    k <- factor(levels(x),levels=levels(x)) # default knots = all regions in the data
  else{
    if (class(area)=="character") k <- as.factor(k)
    if (!(class(k)=="character"||class(k)=="factor")) stop("area must be a factor or a chacacter vector")
  }
  if (length(levels(x))>length(levels(k))) stop("MRF basis dimension set too high")
  if (sum(!levels(x)%in%levels(k))) stop("data contain regions that are not contained in the area specification") 
  x <- factor(x,levels=levels(k))
  nfv <- nlevels(x)
  a.name <- names(neighbour$nb)
  if (all.equal(sort(a.name), sort(levels(k))) != TRUE) 
    stop("mismatch between neighbour/polys supplied area names and data area names")
  np <- nfv
  G <- matrix(0, np, np)
  rownames(G) <- colnames(G) <- levels(k) 
  for (i in 1:np) 
  {
    ind <- neighbour$nb[[i]]
    lind <- length(ind)
    G[a.name[i], a.name[i]] <- lind
    if (lind > 0) 
      for (j in 1:lind) G[a.name[i], a.name[ind[j]]] <- -1
  }
  if (sum(G != t(G)) > 0) 
    stop("Something wrong with auto- penalty construction")
  if (all.equal(sort(a.name),a.name)!=TRUE) 
  { ## re-order penalty to match X
    G <- G[levels(x),]
    G <- G[,levels(x)]
  }
  G
}## end of function nb2prec
##------------------------------------------------------------------------------
##------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------
#This function was adapted from Simon Wood - mgcv package!!

draw.polys <-function( polys, 
                       object = NULL, 
                       scheme = NULL,
                   swapcolors = FALSE,
                        n.col = 100,...)
  {
## to get the range of all polygons    
  for (i in 1:length(polys)) {
    yr <- range(polys[[i]][, 2], na.rm = TRUE)
    xr <- range(polys[[i]][, 1], na.rm = TRUE)
    if (i == 1) {
      ylim <- yr
      xlim <- xr
    }
    else {
      if (yr[1] < ylim[1]) 
        ylim[1] <- yr[1]
      if (yr[2] > ylim[2]) 
        ylim[2] <- yr[2]
      if (xr[1] < xlim[1]) 
        xlim[1] <- xr[1]
      if (xr[2] > xlim[2]) 
        xlim[2] <- xr[2]
    }
  }
  ## of no object just plot the polygons
  mar <- par("mar")
  oldpar <- par(mar = c(2, mar[2], 2, 1))
  if (is.null(object)) {
    plot(0, 0, ylim = ylim, xlim = xlim, xaxt = "n", yaxt = "n", 
         type = "n", bty = "n", ylab = "", xlab = "",...)
    for (i in 1:length(polys)) {
      polygon(polys[[i]], col = NA)
    }
  }
  else 
  { # if object is defined  we must two alternatives
    ## i) it is MRF object
    ## ii) a list which defines the values and the areas    
    if(class(object)=="MRF")
      {
      y.y <- object$beta 
     #  x.x <- object$x
       }
      else
       {
         if (!is.vector(object))  stop("object class should be MRF or a vector with names matching the areas in the polys")
         else
         { 
           y.y <- object
          # x.x <- object[[2]]  
         }
            
        }

    npolys <- names(polys)
   nobject <- names(y.y)
    if (is.null(nobject)) stop("the object do not have names")
    else (!is.null(nobject) && !is.null(npolys)) 
      {
      if (!all(sort(nobject)%in% sort(npolys))) 
        stop("object names and polys names must match")
      }
    
    y.y <- y.y[npolys]
      #fv1 <- tapply(y, x, mean)
      #fv <- object$beta  
    xmin <- xlim[1]
    xlim[1] <- xlim[1] - 0.1 * (xlim[2] - xlim[1])
    n.col <- n.col
    if (is.null(scheme)||scheme=="gray")
      newscheme <- gray(0:n.col/n.col)
    else if (scheme == "heat"){
      newscheme <- heat.colors(n.col + 1)
    }
    else if (scheme == "rainbow")
      newscheme <- rainbow(n.col+1)
    else if(scheme == "terrain")
      newscheme <- terrain.colors(n.col+1)
    else if(scheme == "topo")
      newscheme <- topo.colors(n.col+1)
    else if(scheme=="cm")
      newscheme <- cm.colors(n.col+1)
    else {scheme=scheme
          ramp <- colorRamp(c(scheme, "white"))
          newscheme <-  rgb(ramp(seq(0, 1, length = n.col)), maxColorValue = 255)
    }
    
    if(swapcolors==TRUE){
      if((scheme=="heat")||(scheme=="rainbow")||(scheme=="terrain")||(scheme=="topo")||(scheme=="cm")) 
        newscheme=rev(newscheme)
      else stop("swapcolors just work for few options. Please, see help file.")
    }
    zlim <- range(pretty(y.y))
    for (i in 1:length(polys)) polys[[i]][, 2] <- zlim[1] + (zlim[2] - 
                                                           zlim[1]) * (polys[[i]][, 2] - ylim[1])/(ylim[2] - ylim[1])
    ylim <- zlim
    plot(0, 0, ylim = ylim, xlim = xlim, type = "n", xaxt = "n", 
         bty = "n", xlab = "", ylab = "",...)
    for (i in 1:length(polys)) {
      coli <- round((y.y[i] - zlim[1])/(zlim[2] - zlim[1]) * 
                      n.col) + 1
      polygon(polys[[i]], col = newscheme[coli])
    }
    xmin <- min(c(axTicks(1), xlim[1]))
    dx <- (xlim[2] - xlim[1]) * 0.05
    x0 <- xmin - 2 * dx
    x1 <- xmin + dx
    dy <- (ylim[2] - ylim[1])/n.col
    poly <- matrix(c(x0, x0, x1, x1, ylim[1], ylim[1] + 
                       dy, ylim[1] + dy, ylim[1]), 4, 2)
    for (i in 1:n.col) {
      polygon(poly, col = newscheme[i], border = NA)
      poly[, 2] <- poly[, 2] + dy
    }
    poly <- matrix(c(x0, x0, x1, x1, ylim[1], ylim[2], ylim[2], 
                     ylim[1]), 4, 2)
    polygon(poly, border = "black")
  }
  par(oldpar)
}
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# a function to get a neighbour from a S4 shape file and transfer it a format we can use  
nb2nb <-function(neighbour.nb)
{
  if (!is(neighbour.nb,"nb")) stop("the object is not a nb class")
  newList<-list()
  LL<-length(neighbour.nb)
  naMes <- attr(neighbour.nb, "region.id")
  for (i in 1:LL)
  {
    newList[[as.character(naMes[i])]] <- neighbour.nb[[i]]    
  }
  list(nb=newList)
}


#------------------------------------------------------------------------------
#a function to get polygons from a S4 shape file and transfer it a format we can use  
polys2polys<- function(object,neighbour.nb){#neighbour.nb is given because we need the names according the region.id
  newlist<-list()
  naMes <- attr(neighbour.nb, "region.id")
  for(i in 1:length(object)){
    newlist[[as.character(naMes[i])]] = object[[i]]@Polygons[[1]]@coords
  }
  #as.character(naMes[i]) is writen because we want the polys with the right names
  #then in drawmap the areas will match with the fitted values
  list(newlist)
}
#-------------------------------------------------------------------------------