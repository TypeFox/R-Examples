################################################################################
#### optimiser
################################################################################
#'Function for optimising a landscape genetic analysis based on resistance layers
#'
#'@param landscape resistance layer as a raster object (if not projected dimensions are assumed to be in meters!!) 
#'@param nlocations the number of locations
#'@param mindist minimal distance in map units (meter if not specified)
#'@param fixed n fixed locations, provided as a data.frame with dimenstions nx2
#'@param method least cost path algorithm provided by the gdistance package. Options are "leastcost", "SPDistance" and "commute". see function costdistances.  
#'@param NN number of next neighboring cells to be used, default is 8. see function costdistances.  
#'@param iter number of iterations that should be used to find an optimised design. Try initially the default and if it runs fast (depends on the type of costdistance an d dimenstions of landscape), increase to higher values.
#'@param retry number of retries if optimisation is not possible in a certain iteration (due to mindist and/or fixed locations). Defaults to 10, which should be sufficient in most cases. 
#'@param mask a raster object that masks positions which are not available. Areas which are not to be used for locations are coded as NA (not available), all other values are treated as suitable for placing locations.
#'@param plot if true, some plots are presented that show the histogramm, ecdf of the best (and the worst scenario).
#'@return a list object with two slots. The first slot is called opt and contains iter optimisation values (sd.detour and sd.cost) in a iter x 2 dimenstional data.frame. The second slot is called scenario and contains the corrosponding 1:iter scenarios, given as coordinates in a data.frame of dimensions nlocations x 2. Both slots are ordered in decreasing order of sd.detour values. So the best scenario is at 1 and the worst is at position iter. 
#'@description opt.landgen places iter times nlocations in a given landscape. For each scenario the pairwise costdistances and Euclidean distances are calculated and the standard deviation of detour (see Gruber et al. in prep) is calculated. This metric evaluates the scenerio in their ability to find a significant effect for the given resistance layer on population structure (based on the causal modelling approach of Cushman et al.). To allow for more realistic designs previously locations, a minimal distance between locations and a mask that indicates "forbidden" areas can be specified.  

opt.landgen <- function(landscape,  nlocations, mindist=0, fixed = NULL, method, NN=8, iter =100, retry=10, mask=NULL, plot=TRUE)
{
  scenario <- list()
  opt <- data.frame(sd.cost=NA, sd.detour=NA)#, sd.detour2=NA)
  if (!is.null(fixed))
  {
    specified   <- nrow(fixed)
    
  } else specified = 0
  
  nadd <- nlocations - specified
  
  if (is.na(crs(landscape))) crs(landscape) <-"+proj=merc +units=m" #already proojected data ????
    
  
  for (it in 1:iter)
  {
    r1 <-landscape
    values(r1) <- NA
    
    
    
    #mask #NA values are cutted, everything else is left
    if (!is.null(mask))
    {
      r1 <- mask
      values(r1)<- ifelse(is.na(values(r1)),1,NA)
    }
    retryc <- retry
    ##Define x and y locations
    # random no mindist specified
    if (mindist==0) {
      r1inv <- r1
      values(r1inv) <- ifelse(is.na(values(r1inv)),1,NA)
      rp<-coordinates(as(r1inv,"SpatialPointsDataFrame") )
      choosep <- sample(1:nrow(rp),nadd)
      xs <- rp[choosep,1]
      ys <- rp[choosep,2]
    }
    
    # random with minddist>0
    if (mindist>0)
    {
      
      if (!is.null(fixed))
      {
        rd <- rasterize(fixed, r1,1)
        rd <- buffer(rd, mindist)
        r1 <- sum(r1,rd, na.rm=T)
        values(r1) <- ifelse(values(r1)>0,1,NA)
      }
      left <- sum(is.na(values(r1)))
      if (left==0) {cat("There is no area left to place the required number of locations after placing the fixed locations. Reduce mindist or the amount of fixed locations.\n"); stop()}
      xc <- NA
      yc <- NA
      i <- 1
      rback <- r1
      while (i<=nadd)
      {
        #choose point from left over areas
        r1inv <- r1
        values(r1inv) <- ifelse(is.na(values(r1inv)),1,NA)
        rp<-coordinates(as(r1inv,"SpatialPointsDataFrame") )
        choosep <- sample(1:nrow(rp),1)
        xs <- rp[choosep,1]
        ys <- rp[choosep,2]
        xc[i] <- xs
        yc[i] <- ys
        i <- i+1
        #add new point to mask    
        rd <- rasterize(cbind(xs,ys), r1,1)
        rd <- buffer(rd, mindist)
        r1 <- sum(r1,rd, na.rm=T)
        values(r1) <- ifelse(values(r1)>0,1,NA)
        left <- sum(is.na(values(r1)))
        if (left==0 & i<=nadd) 
        {
          retryc <- retryc - 1
          cat(paste("No area left after ",i,"points at iteration",it,". I go back and try again.", retryc, "retries left...\n"))
          i <- 1
          r1 <- rback
          if (retryc<1) {cat("Could not find any good combination, reduce mindist or increase retry or reduce number of fixed locations.\n");return(list(r1,xc,yc))}
        }
      }
      xs <- xc
      ys <- yc
    }
    
    locs <- cbind(xs,ys)
    rownames(locs) <- LETTERS[1:nadd]
    
    #put fixed and locs together
    locs <- rbind(fixed, locs)
    
    #create a costdistance matrix
    cost.mat <- costdistances(landscape, locs, method, NN)
    eucl.mat <- as.matrix(dist(locs))
    
    sd.cost <- sd(as.dist(cost.mat))
    sd.detour = sd(resid(lm(as.dist(cost.mat)~as.dist( eucl.mat))))
    #sd.detour2 <- sd(as.dist(cost.mat)-as.dist(eucl.mat))
    
    opt[it,] <-c(sd.cost, sd.detour)#, sd.detour2)
    scenario[[it]] <- locs
  } #end of iter loop
  if (plot)
  {
  op <- par(mfrow=c(2,2), mai=c(0.5,0.5,0.2,0.2))
  opt.val <- which.max(opt$sd.detour)
  locs.opt <- scenario[[opt.val]]
  locs <- locs.opt
  cost.mat <- costdistances(landscape, locs, method, NN)
  eucl.mat <- as.matrix(dist(locs))
  detour = resid(lm(as.dist(cost.mat)~as.dist( eucl.mat)))
  
  hist(detour, main="Distrubtion of detour", col="darkgreen")
  
  plot(ecdf(opt$sd.detour), main="Cummulative Density of sd.detour")
  abline(v=opt[opt.val,"sd.detour"] ,col="blue")
  
  #best versus worst case...
  opt.val <- which.max(opt$sd.detour)
  locs.opt <- scenario[[opt.val]]
  locs <- locs.opt
  plot(landscape, main=paste("best:",round(opt[opt.val,"sd.detour"],2) ))
  points(locs[,1], locs[,2], pch=16, cex=1.2, col=c(rep("blue", specified), rep("black", nadd)))
  text(locs[,1],locs[,2], row.names(locs), cex=1)
  
  opt.val <- which.min(opt$sd.detour)
  locs.opt <-  scenario[[opt.val]]
  locs <- locs.opt
  plot(landscape, main=paste("worst:",round(opt[opt.val,"sd.detour"],2) ))
  points(locs[,1], locs[,2], pch=16, cex=1.2,  col=c(rep("blue", specified), rep("black", nadd)))
  text(locs[,1],locs[,2], row.names(locs), cex=1)
  par(op)
  }
  ord <- order(opt$sd.detour, decreasing = TRUE)
  list(opt=opt[ord,], scenario = scenario[ord])
}

