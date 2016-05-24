plot_graph <-
function (rl, species, links)
  {
	if (class(rl)!="landscape" & class(rl)!="metapopulation") 
  {
  stop(paste(rl, " should be either, an object of class class 'landscape' or 'metapopulation'.", sep=""), call. = FALSE)
  }
    
	if(class(rl)=="metapopulation" & species==FALSE)
  {
	stop("With an object of class 'metapopulation' the 'species' parameter has to be TRUE.\n")
  }
	
	if(class(rl)=="metapopulation" & links==TRUE)
  {
	stop("With an object of class 'metapopulation' the 'links' parameter has to be FALSE.\n")
  }
	
	if(class(rl)=="landscape" & species==TRUE)
  {
	stop("With an object of class 'landscape' the 'species' parameter has to be FALSE.\n")
  }
	
    dframe <- rl$nodes.characteristics
    mapsize <- rl$mapsize
    disp <- rl$dispersal
    distx <- pairdist(dframe[,1:2])
    diag(distx) <- NA
    dist_min <- min(distx,na.rm=TRUE)
    if(dist_min > disp) links <- "FALSE"
    if(species == FALSE)
      {
        if(nrow(rl$nodes.characteristics) > 1)  rl <- cluster.id(rl)
        if(nrow(rl$nodes.characteristics) == 1) rl <- rl
        plot(dframe[,1], dframe[,2], xlim=c(min(dframe[,1]),min(dframe[,1])+mapsize), ylim=c(min(dframe[,2]),min(dframe[,2])+mapsize),
             pch=19, xlab="X", ylab="Y", col=dframe[,6])
        symbols(dframe[, 1], dframe[, 2], circles = dframe[, 4], col="deepskyblue4", add = TRUE, inches = FALSE)
      }
    if(species == TRUE)
      {
        cores <- vector(length=nrow(dframe))
        for(i in 1:(nrow(dframe)))
          {
            cores[i] <- ifelse(dframe$species[i] == 0,"grey","green")
          }
        plot(dframe[,1], dframe[,2], xlim=c(min(dframe[,1]),min(dframe[,1])+mapsize), ylim=c(min(dframe[,2]),min(dframe[,2])+mapsize),
             pch=20, xlab="X", ylab="Y", col=cores)
        symbols(dframe[, 1], dframe[, 2], circles = dframe[, 4], col="black", add = TRUE, inches = FALSE)
      }
    if(nrow(rl$nodes.characteristics) > 1)
      { 
        if(links==TRUE)
          {
            df_edges <- edge.graph(cluster.id(rl))
            x0 <- df_edges$XA
            x1 <- df_edges$XB
            y0 <- df_edges$YA
            y1 <- df_edges$YB
            segments(x0, y0, x1, y1, col="goldenrod3")
          }
      }
    if(species == TRUE)
      {
        if(links == TRUE)
          {
            segments(x0, y0, x1, y1, col="red")
          }
      }
    if(species == FALSE)
      {
        if(links == TRUE)
          {
            segments(x0, y0, x1, y1, col="goldenrod3")
          }
      }
  }
