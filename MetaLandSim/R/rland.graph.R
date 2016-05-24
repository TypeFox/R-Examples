rland.graph <-
function(mapsize, dist_m, areaM, areaSD, Npatch, disp, plotG)
  {
  if(areaM^2<3*areaSD^2){
  SDmax <- round(sqrt(areaM^2/3),3)
  Mmin <- round(sqrt(3*areaSD^2),3)
  stop(call. = FALSE,cat(paste("For a mean patch size of ",areaM," hectares the areaSD parameter must be \nlower than ",SDmax,".","\n",
  "\nAlternatively, for a standard deviation of ",areaSD, " the areaM parameter \nmust be higher than ", Mmin,".\n",
  "\nOtherwise patches with negative areas might be produced.\n",sep="")))
  }
    reg_pts <- rSSI(r=dist_m, n=Npatch, win = owin(c(0,mapsize),c(0, mapsize)),
                    giveup = 1000, x.init=NULL)
    reg_pts1 <- as.data.frame(reg_pts)
    if(areaSD > 0){
      areas_Ha <- matrix(nrow=Npatch,ncol=1)
	  b <- (sqrt(areaSD^2*12)+(areaM*2))/2
	  a <- (areaM*2) - b
      x1 <- runif(Npatch,a,b)
	  areas_Ha[,1] <- x1
    }
	if(areaSD == 0){
      areas_Ha <- matrix(nrow=Npatch,ncol=1)
      x1 <- rep(areaM, Npatch)
	  areas_Ha[,1] <- x1
	}
	radius <- sqrt((areas_Ha*10000)/pi)	
	dframe <- cbind(reg_pts1, areas_Ha, radius)
    grouping <- hclust(dist(dframe[, 1:2], method = "euclidean"), "single")
    clusters <- cutree(grouping, h=disp) 
    rg2 <- cbind(dframe, clusters)
    col1 <- rainbow(max(rg2[, 5]))
    col2 <- as.data.frame(col1)
    col2[, 2] <- seq(1:max(rg2[, 5]))
    col3 <- merge(rg2, col2, by.x = "clusters", by.y = "V2")
    col5 <- nndist (reg_pts1)
    col4 <- data.frame(col3$x, col3$y, col3$areas_Ha, col3$radius, col3$clusters, col3$col1)
    col4 <- cbind(col4,col5)
    ids <- 1:nrow(col4)
    col4 <- cbind (col4,ids)
    names(col4)[names(col4)=="col3.x"] <- "x"
    names(col4)[names(col4)=="col3.y"] <- "y"
    names(col4)[names(col4)=="col3.areas_Ha"] <- "areas"
    names(col4)[names(col4)=="col3.radius"] <- "radius"
    names(col4)[names(col4)=="col3.clusters"] <- "cluster"
    names(col4)[names(col4)=="col3.col1"] <- "colour"
    names(col4)[names(col4)=="col5"] <- "nneighbour"
    names(col4)[names(col4)=="ids"] <- "ID"
    if (plotG==TRUE)
      {
        plot(col4[,1], col4[,2], xlim=c(0,mapsize), ylim=c(0,mapsize),
             pch=19, xlab="X", ylab="Y", col=col4[,6])
        symbols(dframe[, 1], dframe[, 2], circles = dframe[, 4], col="deepskyblue4", add = TRUE, inches = FALSE)
        distxy <- pairdist (dframe[,1], dframe[,2])
        d <- distxy < disp
        ind <- which(d, arr.ind = TRUE, useNames = FALSE)
        x0 <- as.vector(reg_pts1[as.vector(ind[,1]),1])
        x1 <- as.vector(reg_pts1[as.vector(ind[,2]),1])
        y0 <- as.vector(reg_pts1[as.vector(ind[,1]),2])
        y1 <- as.vector(reg_pts1[as.vector(ind[,2]),2])
        segments(x0, y0, x1, y1,col = "goldenrod3")
      }
    minimum.distance <- dist_m
    mean.area2 <- mean(col4$areas)
    SD.area2 <- sd(col4$areas)
    number.patches2 <- nrow(col4)
    dispersal <- disp
    rland.out <- list(mapsize=mapsize, minimum.distance=dist_m, mean.area=mean.area2,
                      SD.area=SD.area2, number.patches=number.patches2,dispersal=disp,
                      nodes.characteristics=col4)
    class(rland.out) <- "landscape"
    return(rland.out)
  }
