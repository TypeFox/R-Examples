simulate_graph <-
  function (rl, rlist, simulate.start, method, parm,
            nsew="none", succ="none", param_df,kern, conn, 
            colnz, ext,beta1,b, c1, c2, z, R)
  {
    if(class(rl)=="landscape" && simulate.start==FALSE) stop("When the object 'rl' is a landscape (class 'landscape')\n the parameter 'simulate.start' must be set to TRUE!")
    if(class(rl)=="metapopulation" && simulate.start==TRUE) stop("When the object 'rl' is an occupied landscape (class 'metapopulation')\n the parameter 'simulate.start' must be set to FALSE!")
    
	span <- length(rlist)
    max_age <- span
    rlist <- patch_age(rlist)
    rl <- attr_age_unique_landscape(rland=rl,span_age=rlist,position=1)
    metpop.list <- as.list(rep("", span))
    turnover.list <- as.list(rep("", span))
    if(simulate.start==TRUE)
    {
      sp_0 <- species.graph(rl,method=method,parm=parm,nsew=nsew,plotG=FALSE)
    }
    if (simulate.start==FALSE)
    {
	  n1 <- rl$nodes.characteristics
	  n1 <- cbind(n1[,1:8],n1[,10],n1[,9])
	  names(n1)[names(n1)=="n1[, 9]"] <- "species"
	  names(n1)[names(n1)=="n1[, 10]"] <- "age"
      
	  mapsize1 <- rl$mapsize
      minimum.distance1 <- rl$minimum.distance
      mean.area1 <- rl$mean.area
      SD.area1 <- rl$SD.area
      number.patches1 <- rl$number.patches
      dispersal1 <- rl$dispersal
      neigh1 <- rl$distance.to.neighbours
	  
      sp_0 <- list(mapsize=mapsize1, minimum.distance=minimum.distance1, 
                        mean.area=mean.area1, SD.area=SD.area1, number.patches=number.patches1,
                        dispersal=dispersal1, distance.to.neighbours=neigh1,
                        nodes.characteristics=n1)
      class(sp_0) <- "metapopulation"
	    
    }
    sp_1 <- sp_0$nodes.characteristics
    metpop.list[[1]] <- sp_1
    turnover.list[[1]] <- 0
    for(i in 2:span)
    {
      prec.sp <- metpop.list[[i-1]]
      mapsize <- as.numeric(rl[[1]])
      minimum.distance <- as.numeric(rl[[2]])
      mean.area <- mean(prec.sp$areas)
      SD.area <- sd(prec.sp$areas)
      number.patches <- nrow(prec.sp)
      dispersal <- as.numeric(rl[[6]])
      neigh <- as.data.frame(pairdist(prec.sp[, 1:2]))
      names(neigh) <- prec.sp$ID
      rownames(neigh) <- prec.sp$ID
      prec.sp$nneighbour <- nndist(prec.sp[, 1:2])
      prec.sp_1 <- list(mapsize=mapsize, minimum.distance=minimum.distance,
                        mean.area=mean.area, SD.area=SD.area,
                        number.patches=number.patches,dispersal=dispersal,
                        distance.to.neighbours=neigh,nodes.characteristics=prec.sp)
      class(prec.sp_1) <- "metapopulation"
      out_0 <- spom(prec.sp_1, kern, conn, colnz, ext, param_df, beta1,
                    b, c1, c2, z, R, succ, max_age)
      turnover.list[[i]] <- ((out_0$turnover*100)/nrow(out_0$nodes.characteristics))
      out_1 <- out_0$nodes.characteristics[, -c(10,12)]
      names(out_1)[names(out_1)=="species2"] <- "species"
      lands_i <- rlist[[i]]
      out_2 <- merge_order(lands_i, out_1, by.x = "ID", by.y = "ID",sort=FALSE,keep_order=TRUE,all.x=TRUE,all.y=TRUE)
      out_3 <- out_2[, c(1:9,18)]
      out_3 <- na.omit(out_3)
      if(any(is.na(out_2[, 10:18]))==TRUE)
      {
        out_4 <- out_2[is.na(out_2$species),]
        out_4 <- out_4[,c(1:9,18)]
        out_4[,10] <- rep(0,nrow(out_4))
        out_3 <- rbind(out_3,out_4)
      }
      out_4 <- data.frame(out_3$x.x, out_3$y.x, out_3$areas.x, out_3$radius.x,
                          out_3$cluster.x, out_3$colour.x, out_3$nneighbour.x,
                          out_3$ID, out_3$age, out_3$species)
      names(out_4)[names(out_4)=="out_3.x.x"] <- "x"
      names(out_4)[names(out_4)=="out_3.y.x"] <- "y"
      names(out_4)[names(out_4)=="out_3.areas.x"] <- "areas"
      names(out_4)[names(out_4)=="out_3.radius.x"] <- "radius"
      names(out_4)[names(out_4)=="out_3.cluster.x"] <- "cluster"
      names(out_4)[names(out_4)=="out_3.colour.x"] <- "colour"
      names(out_4)[names(out_4)=="out_3.nneighbour.x"] <- "nneighbour"
      names(out_4)[names(out_4)=="out_3.ID"] <- "ID"
      names(out_4)[names(out_4)=="out_3.species"] <- "species"
      names(out_4)[names(out_4)=="out_3.age"] <- "age"
      metpop.list[[i]] <- out_4
    }
    turnover <-as.numeric(turnover.list)
    output <- list(metpop.list,turnover=turnover)
    return(output)
}