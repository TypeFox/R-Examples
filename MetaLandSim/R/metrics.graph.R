metrics.graph <-
function(rl, metric)
  {
	if (class(rl)!="landscape") 
  {
  stop(paste(rl, " should be an object of class class 'landscape'.", sep=""), call. = FALSE)
  }
  
  result <- c()
  
    if("NC" %in% metric)
      {
        result <- c(result,NC = components.graph(rl))
      }
    if("LNK" %in% metric)
      {
        result <- c(result, LNK = nrow(edge.graph(rl)))
      }
    if("SLC" %in% metric)
      {
        df0 <- rl$nodes.characteristics
        NC <- components.graph(rl)
        area_sum <- rep(NA, NC)
        for(i in 1:NC)
          {
            component <- df0[df0$cluster==i, ]
            area_sum[i] <- sum(component$areas)
          }
        result <- c(result, SLC = max(area_sum))
      }
    if("MSC" %in% metric)
      {
        df0 <- rl$nodes.characteristics
        NC <- components.graph(rl)
        area_sum <- rep(NA, NC)
        for(i in 1:NC)
          {
            component <- df0[df0$cluster==i, ]
            area_sum[i] <- sum(component$areas)
          }
        result <- c(result, MSC = mean(area_sum))
      }
    if("HI" %in% metric)
      {
        m0 <- min_distance(rl)
		m1 <- m0[upper.tri(m0)]
        m2 <- 1/m1
        result <- c(result, HI = sum(m2[!is.infinite(m2)]))
      }
    if("NH" %in% metric)
      {
        df0 <- rl$nodes.characteristics
        Ncomps <- length(unique(df0$cluster))
        H_vector <- rep(NA, Ncomps)
        for(b in 1:Ncomps)
          {
            component1 <- df0[df0$cluster==b, ]
            nnodesComp1 <- nrow(component1)
            r1 <- list(mapsize=rl$mapsize, minimum.distance=rl$minimum.distance,
                       mean.area=rl$mean.area, SD.area=rl$SD.area, number.patches=nnodesComp1,
                       dispersal=rl$dispersal, nodes.characteristics=component1)
			class(r1) <- "landscape"
			m_r1 <- min_distance(r1)
			m1 <- m_r1[upper.tri(m_r1)]
            m2 <- 1/m1
            H <- sum(m2[!is.infinite(m2)])
            H_vector[b] <- H
          }
        H_chain_vector <- rep(NA, Ncomps)
        H_planar_vector <- rep(NA, Ncomps)
        for(r in 1:Ncomps)
          {
            component <- df0[df0$cluster==r, ]
            nnodesComp <- nrow(component)
            H_chain_vector_2 <- rep(NA, nnodesComp-1)
            H_chain_vector_2[1] <- (nnodesComp-1)
            H_chain_vector_2[nnodesComp-1] <- 1/(nnodesComp-1)
            for(f in 2:nnodesComp-2)
              {
                H_chain_vector_2[f] <- (nnodesComp-f)/f
              }
            H_chain_c <- sum(H_chain_vector_2)
            H_chain_vector[r] <- H_chain_c
            H_planar <- ((nnodesComp*(nnodesComp+5))/4)-3
            H_planar_vector[r] <- H_planar
          }
        H <- sum(H_vector)
        H_chain <- sum(H_chain_vector)
        H_planar <- sum(H_planar_vector)
        result <- c(result, NH = (H-H_chain)/(H_planar-H_chain))
      }
    if("ORD" %in% metric)
      {
        a <- cluster.graph(rl)
        result <- c(result, ORD = max(a[, 2]))
      }
    if("GD" %in% metric)
      {
        result <- c(result, GD = max(min_distance(rl), na.rm = TRUE))
      }
    if("CCP" %in% metric)
      {
        df0 <- rl$nodes.characteristics
        NC <- components.graph(rl)
        Ac <- sum(df0$areas)
        r0 <- rep(NA, NC)
        for(i in 1:NC)
          {
            df1 <- df0[df0$cluster==i, ]
            ci <- sum(df1$areas)
            r0[i] <- (ci/Ac)^2
          }
        result <- c(result, CCP = sum (r0))
      }
    if("LCP" %in% metric)
      {
        NC <- components.graph(rl)
        df0 <- rl$nodes.characteristics
        r0_vec <- rep(NA, NC)
        AL <- ((rl$mapsize)^2)/10000
        for(i in 1:NC)
          {
            df1 <- df0[df0$cluster==i, ]
            ci <- sum(df1$areas)
            r0_vec[i] <- (ci/AL)^2
          }
        result <- c(result, LCP = sum(r0_vec))
      }
    if("CPL" %in% metric)
      {
        m0 <- min_distance(rl)
		m1 <- m0[upper.tri(m0)]
        L <- length(m1)
        result <- c(result, CPL = sum(m1)/L)
      }
    if("ECS" %in% metric)
      {
        NC <- components.graph(rl)
        df0 <- rl$nodes.characteristics
        ai <- rep(NA, NC)
        for(i in 1:NC)
          {
            df1 <- df0[df0$cluster==i, ]
            ai[i] <- sum((df1$areas)^2)
          }
        a_num <- sum(ai)
        a <- sum(df0$areas)
        result <- c(result, ECS = a_num/a)
      }
    if("AWF" %in% metric)
      {
        disp <- rl$dispersal
        d0 <- rl$nodes.characteristics
        nnodes <-rl$number.patches
        distP <- pairdist (d0[,1:2])
        ext <- log(0.05)/disp
        pij <- as.data.frame(exp(ext*(distP)))
        names(pij) <- d0$ID
        rownames(pij) <- d0$ID
        paths <- as.data.frame(which(min_distance(rl)!=0, arr.ind = TRUE, useNames = FALSE))
        for(i in 1: nrow(paths)) paths[i, 3] <- d0$areas[d0$ID %in% paths[i, 1]]
        for(i in 1: nrow(paths)) paths[i, 4] <- d0$areas[d0$ID %in% paths[i, 2]]
        for(i in 1: nrow(paths))paths[i, 5] <- pij[paths[i,1],paths[i,2]]
        result <- c(result, AWF = 2*(sum(paths[, 3]*paths[, 4]*paths[, 5])))
      }
    if("IIC" %in% metric)
      {
        dist_tp <- min_distance(rl)
        df0 <- rl$nodes.characteristics
        paths <- as.data.frame(which(min_distance(rl)!=0, arr.ind = TRUE, useNames = FALSE))
        paths <- rbind(as.data.frame(which(dist_tp > 0, arr.ind = TRUE, useNames = FALSE)),
                       cbind(min(paths):max(paths), min(paths):max(paths)))
        names(paths)[names(paths) == "V1"] <- "nodo_A"
        names(paths)[names(paths) == "V2"] <- "nodo_B"
        for(i in 1:nrow(paths)) paths[i, 3] <- dist_tp[paths[i, 1],paths[i, 2]]
        names(paths)[names(paths) == "V3"] <- "dist_topo"
        for(f in 1:nrow(paths))paths[f, 4] <- df0$areas[df0$ID %in% paths[f, 1]]
        names(paths)[names(paths) == "V4"] <- "area_A"
        for(x in 1:nrow(paths))paths[x, 5] <- df0$areas[df0$ID %in% paths[x, 2]]
        names(paths)[names(paths) == "V5"] <- "area_B"
        paths[, 6] <- (paths[, 4]*paths[, 5])/(1+(2*paths[, 3]))
        Al2 <- (((rl$mapsize)^2)/10000)^2
        result <- c(result, IIC = sum(paths[, 6])/Al2)
      }
    if("PC" %in% metric)
      {
        d0 <- rl$nodes.characteristics
        disp <- rl$dispersal
        m2 <- as.matrix(d0[, 1:2])
        m3 <- pairdist(m2)
        m3[m3>disp] <- 0
        m3[m3 == 0] <- NA
        dmin <- min_distance(rl)
        connect <- as.data.frame(which(dmin!=0, arr.ind = TRUE, useNames = FALSE))
        paths <- rbind(connect, cbind(min(connect):max(connect), min(connect):max(connect)))
        nr_conn <- nrow(paths)
        dists <- as.data.frame(pairdist(m2))
        colnames(dists) <- d0$ID
        rownames(dists) <- d0$ID
        for(r in 1:nr_conn) paths[r, 3] <- dists[paths[r, 1],paths[r, 2]]
        for(r in 1:nr_conn)paths[r, 4] <- dmin[paths[r, 1],paths[r, 2]]
        colnames(paths)[1] <- "node A"
        colnames(paths)[2] <- "node B"
        colnames(paths)[3] <- "distance"
        colnames(paths)[4] <- "top_distance"
        m4 <- allShortestPaths(m3)
        m5 <- rep(NA, nr_conn)
        for(i in 1:nr_conn)
          {
            if (paths$top_distance[i]>1)
              {
                a <- extractPath(m4, paths[i, 1], paths[i, 2])
                b <- length (a)-1
                b2 <- length(a)
                c1 <- as.data.frame(matrix(NA, b, 4))
                a2 <- rep(NA, (2*(b2-2))+2)
                a2[1] <- a[1]
                a2[length(a2)] <- a[length(a)]
                a3 <- rep(a[2:(length(a)-1)] ,each=2)
                for(s in 1:length(a3)) a2[s+1] <- a3[s]
                newpairs <- as.data.frame(matrix(a2, ncol=2, byrow=TRUE))
                ext <- log(0.05)/disp
                for(g in 1: nrow(newpairs))newpairs[g, 3] <- dists[newpairs[g, 1],newpairs[g, 2]]
                for(h in 1: nrow(newpairs))newpairs[h, 4] <- exp(ext*(newpairs[h, 3]))
                colnames(newpairs)[1] <- "node_A"
                colnames(newpairs)[2] <- "node_B"
                colnames(newpairs)[3] <- "distance"
                colnames(newpairs)[4] <- "probability"
                d <- prod(newpairs[, 4])
                m5[i] <- d
              }
            if(paths$top_distance[i] == 1)
              {
                ext <- log(0.05)/disp
                m5[i] <- exp(paths$distance[i]*ext)
              }
            if(paths[i, 1] == paths[i, 2])
              {
                m5[i] <- 1
              }
          }
        paths2 <- cbind(paths, m5)
        colnames(paths2)[5] <- "pij*"
        for(i in 1:nrow(paths2))paths2[i, 6] <- d0$areas[d0$ID %in% paths2[i, 1]]
        for(i in 1:nrow(paths2))paths2[i, 7] <- d0$areas[d0$ID %in% paths2[i, 2]]
        colnames(paths2)[6] <- "areaA"
        colnames(paths2)[7] <- "areaB"
        paths2[, 8] <- paths2[, 5]*paths2[, 6]*paths2[, 7]
        Al2 <- ((rl$mapsize)^2)/10000
        result <- c(result, PC = (sum(paths2[, 8]))/Al2^2)
      }
    return(round(result,3))
  }

  #Nada?