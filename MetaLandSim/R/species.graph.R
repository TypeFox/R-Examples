species.graph <-
  function(rl, method="percentage", parm, nsew="none", plotG=TRUE)
  {
    if (class(rl)!="landscape") 
    {
      stop(paste(rl, " should be an object of class class 'landscape'.", sep=""), call. = FALSE)
    }
    d1 <- rl$nodes.characteristics
    if(any(names(d1)=="age")){
	d1<-d1
    }else{
	  nrow1 <- nrow(d1)
      ages1 <- rep(1,nrow1)
      d1 <- cbind(d1,ages1)
      names(d1)[9] <- "age"
	      }
    mapsize <- rl$mapsize
    sp1 <- SpatialPoints(d1[, 1:2])
    dist_pair <- matrix.graph(rl,mat="centr_distance")
    dp1 <- as.data.frame(dist_pair)
    age <- rl$nodes.characteristics$age
    d2 <- cbind(d1,dp1)
    names(d2)[10:ncol(d2)] <- as.character(d2$ID)
    if(method == "click")
    {
      select <- select.spatial(data=sp1, digitize = FALSE, pch = "+", rownames = FALSE)
    }
    if(method == "percentage")
    {
      vector_sel <- c(as.numeric(rownames(d2)))
      size_vecsel <- length(vector_sel)
      if(nsew == "none")
      {
        nr <- round((size_vecsel*parm)/100)
        select <- sample(vector_sel, nr, replace = FALSE, prob = NULL)
      }
      if(nsew == "N")
      {
        d5 <- rev(sort(d2[, 2]))
        ny <- round((size_vecsel*parm)/100)
        py <- d5[1:ny]
        mN <- match (d2[, 2], py)
        dfN <- cbind (d2, mN)
        select <- as.numeric(rownames(na.omit(dfN[])))
      }
      if(nsew == "S")
      {
        d5 <- sort(d2[, 2])      
        ny <- round((size_vecsel*parm)/100) 
        py <- d5[1:ny] 
        mN <- match (d2[, 2], py)
        dfN <- cbind (d2, mN) 
        select <- as.numeric(rownames(na.omit(dfN[]))) 
      }
      if(nsew == "E")
      {
        d5 <- rev(sort(d2[, 1])) 
        nx <- round((size_vecsel*parm)/100) 
        px <- d5[1:nx] 
        mN <- match (d2[, 1], px)
        dfN <- cbind (d2, mN) 
        select <- as.numeric(rownames(na.omit(dfN[]))) 
      }
      if(nsew == "W")
      {
        d5 <- sort(d2[, 1])      
        nx <- round((size_vecsel*parm)/100) 
        px <- d5[1:nx] 
        mN <- match (d2[, 1], px)
        dfN <- cbind (d2, mN) 
        select <- as.numeric(rownames(na.omit(dfN[]))) 
      }
    }
    if(method == "number")
    {
      vector_sel <- c(as.numeric(rownames(d2))) 
      size_vecsel <- length(vector_sel)
      if(nsew == "none")
      {
        select <- sample(vector_sel, parm, replace = FALSE, prob = NULL)
      }
      if(nsew == "N")
      {
        d5 <- rev(sort(d2[, 2])) 
        ny <- parm 
        py <- d5[1:ny] 
        mN <- match (d2[, 2], py)
        dfN <- cbind (d2, mN) 
        select <- as.numeric(rownames(na.omit(dfN[]))) 
      }
      if(nsew == "S")
      {
        d5 <- sort(d2[, 2]) 
        ny <- parm 
        py <- d5[1:ny] 
        mN <- match (d2[, 2], py)
        dfN <- cbind (d2, mN) 
        select <- as.numeric(rownames(na.omit(dfN[]))) 
      }
      if(nsew == "E")
      {
        d5 <- rev(sort(d2[, 1])) 
        nx <- parm 
        px <- d5[1:nx] 
        mN <- match (d2[, 1], px)
        dfN <- cbind (d2, mN) 
        select <- as.numeric(rownames(na.omit(dfN[]))) 
      }
      if(nsew == "W")
      {
        d5 <- sort(d2[, 1]) 
        nx <- parm 
        px <- d5[1:nx] 
        mN <- match (d2[, 1], px)
        dfN <- cbind (d2, mN) 
        select <- as.numeric(rownames(na.omit(dfN[]))) 
      }
    }
    vec1 <- rep(0, nrow(dist_pair))
    df0 <- as.data.frame(vec1)
    df0[select,] <- 1
    d3 <- cbind (d2, df0)
    names(d3)[names(d3) == "vec1"] <- "species"
    d4 <- d3 
    if(plotG == TRUE)
    {
      colour <- ifelse(d4$species == 0, "red", "forestgreen")
      circ <- symbols(x=d4[, 1], y =d4[, 2], circles=d4[, 4],
                      inches=FALSE, add=FALSE, fg=NULL, bg=colour,
                      xlab="X", ylab="Y", main=NULL)
    }
    mapsize <- as.numeric(rl[[1]])
    minimum.distance <- as.numeric(rl[[2]])
    mean.area <- as.numeric(rl[[3]])
    SD.area <- as.numeric(rl[[4]])
    number.patches <- as.numeric(rl[[5]])
    dispersal <- as.numeric(rl[[6]])
    neigh <- d4[,10:(nrow(d4)+9)]
    nodes <- d4[,c("x", "y", "areas", "radius", "cluster", "colour", "nneighbour",
                   "ID","age", "species")]
    if(any(names(rl$nodes.characteristics)=="age")){
       nodes <- nodes
	} else nodes <- nodes[!(names(nodes) %in% c("age"))]
    species.out <- list(mapsize=mapsize, minimum.distance=minimum.distance, 
                        mean.area=mean.area, SD.area=SD.area, number.patches=number.patches,
                        dispersal=dispersal, distance.to.neighbours=neigh,
                        nodes.characteristics=nodes)
    class(species.out) <- "metapopulation"
    return(species.out)
  }