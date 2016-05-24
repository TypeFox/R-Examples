span.graph <-
  function (rl, span=100, par1="none", par2=NULL, par3=NULL, par4=NULL, par5=NULL)
  {
    if (class(rl)!="landscape")
    {
      stop(paste(rl, " should be an object of class class 'landscape'.", sep=""), call. = FALSE)
    }
    rll_0 <- rl 
    mapsize2 <- rl$mapsize
    dist_m2 <- rl$minimum.distance
    areaM2 <- rl$mean.area
    areaSD2 <- rl$SD.area
    Npatch2 <- rl$number.patches
    disp2 <- rl$dispersal
    rl <- rl$nodes.characteristics
    ID2 <- rl$ID
    rland.list <- as.list(rep(NA, span))
    rland.list[[1]] <- rl
    if(par1 == "hab")
    { 
      vect1 <- c()
      for(i in 2:span)
      { 
        n_select <- (nrow(na.omit(rland.list[[i-1]]))-round((nrow(na.omit(rland.list[[i-1]]))*par2)/100))
        if(n_select == nrow(rland.list[[i-1]]))vect1 <- c(vect1,i)
        else vect1 <- c(vect1,0)
        g1 <- na.omit(rland.list[[i-1]][sample(nrow(na.omit(rland.list[[i-1]])), 
                                               n_select, replace=FALSE),])
        g2 <- g1[sort.list(as.numeric(rownames(g1))),] 
        rownames(g2) <- 1:nrow(g2)
        rland.list[[i]] <- g2
        if (nrow(g2) == 0)
        {
          warning ("Unable to generate all the landscapes.\nSome landscapes would have no habitat patches.") 
          break
        }
      }
      if (sum(vect1)!=0)
      {
        
        vect1 <- vect1[vect1!=0]   
        index1 <- min(vect1)-1
        message("After time step ",index1," all landscapes had the same number of patches.","\n")
      }
    } 
    if(par1 == "dincr")
    {
      idist <- par2 
      md2 <- seq(from=(dist_m2+idist), by=idist, length.out=(span-1))
      md3 <- c(dist_m2, md2)    
      md <- vector(length=span) 
      for (i in 2:span)
      {
        md <- md3[i] 
        df1 <- rland.list[[i-1]] 
        ndist <- nndist(df1[, 1:2]) 
        df2 <- cbind(df1, ndist) 
        vec_s <- df2$ndist < md
        df3 <- cbind(df2, vec_s) 
        df4 <- df3[df3$vec_s == "FALSE",]
        df5 <- df4[, -c(7,8)]
        names(df5)[names(df5) == "ndist"] <- "nneighbour" 
        df6 <- df5[, -8] 
        df6_1 <- nndist(df6[, 1:2])
        df6_2 <- cbind(df6, df6_1)
        df6_3 <- df6_2[, -7]
        names(df6_3)[names(df6_3) == "df6_1"] <- "nneighbour" 
        match1 <- rl[, c(1, 8)]
        df7 <- merge_order(df6_3, match1, by.x="x", by.y="x",
                           sort=FALSE, keep_order=TRUE)
        rland.list[[i]] <- df7
      }
    }
    if(par1 == "darea")
    { 
      for(b in 2:span)
      {
        rl0 <- rland.list[[b-1]] 
        areaM <-(rl0[, "areas"]) 
        a_dim <- (areaM*par2)/100
        areaM2 <- areaM-a_dim
        radius <- sqrt((areaM2*10000)/pi)
        rl0$areas <- areaM2
        rl0$radius <- radius 
        rl1 <- subset(rl0, rl0[, 3] > par3)
        if(nrow(rl1)==0)
        {
          message("After time step",b+1, "the landscape has no patches!", "\n")
          break		  
        }
        rland.list[[b]] <- rl1
      }
    }
    if(par1 == "none")
    {
      rland.list <- rep(list(rl), span)
    }
    if(par1=="stoc")
    { 
      addpoints1 <- function (rl,nr,IDmax)
      {
        if(nrow(rl$nodes.characteristics)!=0){   
          if (nr!=0)
          {
            mapsize3 <- rl$mapsize
            dist_m3 <- rl$minimum.distance
            areaM3 <- rl$mean.area
            areaSD3 <- rl$SD.area
            Npatch3 <- rl$number.patches
            disp3 <- rl$dispersal
            rl <- rl$nodes.characteristics
            ID2 <- rl$ID
            rl2 <- rl[,1:2]
            wind <- owin(xrange=c(0,mapsize3), yrange=c(0,mapsize3))
            suppressWarnings(pts_0 <- as.ppp(rl2, W = wind,fatal=TRUE))
            if(npoints(pts_0)+nr >1)  pts_1 <- rSSI(r=dist_m3, n = npoints(pts_0)+nr, win = wind,giveup = 1000, x.init=pts_0)
            if(npoints(pts_0)+nr==1)	pts_1 <- rSSI(r=dist_m3, n = npoints(pts_0)+nr, win = wind, giveup = 1000)		
            df_pts0 <- as.data.frame(coords(pts_1))
            df_pts1 <- as.data.frame(coords(pts_0))
            df_pts2 <- df_pts0[!duplicated(rbind(df_pts1, df_pts0))[nrow(df_pts1) + 
                                                                      1:nrow(df_pts0)],]
            nrow_0 <- nrow(rl)
            na_lines <- as.data.frame(matrix(NA,nrow=nr,ncol=ncol(rl)))
            colnames(na_lines) <- colnames(rl)
            if(nrow(na_lines)>1)rownames(na_lines) <- max(as.numeric(rownames(rl)))+1:nrow(na_lines)
            if(nrow(na_lines)==1)rownames(na_lines)<-1
            rl <- rbind(rl,na_lines)
            rl[(nrow_0+1):nrow(rl),1:2] <- df_pts2
            b1 <- (sqrt(areaSD3^2*12)+(areaM3*2))/2
            a <- (areaM3*2) - b1
            areas0 <- runif(nr,a,b1)
            radius0 <- sqrt((areas0*10000)/pi)
            rl[(nrow_0+1):nrow(rl),"areas"] <- areas0
            rl[(nrow_0+1):nrow(rl),"radius"] <- radius0
            if(nrow(rl)>1)new_ID <- (IDmax+1):(IDmax+nr)
            if(nrow(rl)==1)new_ID <- IDmax+1
            rl[,8] <- as.character (c(ID2,new_ID))
            if(nrow(rl)>1){
              grouping <- hclust(dist(rl[,1:2],method = "euclidean"),"single")
              clusters <- as.data.frame(cutree(grouping, h=disp2))[,1]
            }
            if(nrow(rl)==1)clusters <- 1
            rl[,"cluster"] <- clusters
            col1 <- rainbow(max(rl[,5]))
            col2 <- as.data.frame(col1)
            col2[,2] <- seq(1:max(rl[,5]))
            col3 <- merge_order(rl, col2, by.x = "cluster", 
                                by.y = "V2",sort=FALSE,keep_order=TRUE)[,9]
            rl[,"colour"] <- col3
            if(nrow(rl)>1)col4 <- nndist (rl[,1:2])
            if(nrow(rl)==1) col4 <- 0
            rl[,"nneighbour"] <- col4
            rland.out <- list(mapsize=mapsize3, minimum.distance=dist_m3,
                              mean.area=areaM3, SD.area=areaSD3,
                              number.patches=nrow(rl),dispersal=disp3,
                              nodes.characteristics=rl)
            class(rland.out) <- "landscape"
            if(nrow(rl)>1)rland.out <- cluster.id(rland.out)
          }
          if (nr==0)
          {
            rland.out <- rl
          }
        }
        else if(nrow(rl$nodes.characteristics)==0){
          if(nr==0){
          mapsize3 <- rl$mapsize
          dist_m3 <- rl$minimum.distance
          areaM3 <- rl$mean.area
          areaSD3 <- rl$SD.area
          disp3 <- rl$dispersal
          
          rland.out <-  list(mapsize=mapsize3, minimum.distance=dist_m3,
                               mean.area=areaM3, SD.area=areaSD3,
                               number.patches=nr,dispersal=disp3,
                               nodes.characteristics=rl$nodes.characteristics)
            
          class(rland.out) <- "landscape"
          }
          if(nr!=0 && nr!=1){
            mapsize3 <- rl$mapsize
            dist_m3 <- rl$minimum.distance
            areaM3 <- rl$mean.area
            areaSD3 <- rl$SD.area
            disp3 <- rl$dispersal
            
            rland.out <- rland.graph(mapsize=mapsize3, dist_m=dist_m3, areaM=areaM3, 
                         areaSD=areaSD3, Npatch=nr,
                         disp=disp3, plotG=FALSE)
            
            rland.out$nodes.characteristics$ID <- (IDmax+1):(IDmax+nrow(rland.out$nodes.characteristics))
            
          }
          if(nr==1){
            mapsize3 <- rl$mapsize
            dist_m3 <- rl$minimum.distance
            areaM3 <- rl$mean.area
            areaSD3 <- rl$SD.area
            disp3 <- rl$dispersal
            rl0 <- rl$nodes.characteristics
            
            wind <- owin(xrange=c(0,mapsize3), yrange=c(0,mapsize3))
            pts_1 <- rSSI(r=dist_m3, n = 1, win = wind, giveup = 1000)		
            df_pts0 <- as.data.frame(coords(pts_1))
            rl0[1,1] <- df_pts0[1,1]
            rl0[1,2] <- df_pts0[1,2]
            
            b1 <- (sqrt(areaSD3^2*12)+(areaM3*2))/2
            a <- (areaM3*2) - b1
            areas0 <- runif(nr,a,b1)
            rl0[1,3] <- areas0
            
            radius0 <- sqrt((areas0*10000)/pi)
            rl0[1,4] <- radius0 
            
            rl0[1,5] <- 1
            
            colour0 <- rainbow(nr)
            rl0[1,6] <- colour0
            rl0[1,7] <- 0
            rl0[1,8] <- IDmax+1
            
            rland.out <-  list(mapsize=mapsize3, minimum.distance=dist_m3,
                               mean.area=areaM3, SD.area=areaSD3,
                               number.patches=nr,dispersal=disp3,
                               nodes.characteristics=rl0)
            
            class(rland.out) <- "landscape"
          }
        }
        
        if(nrow(rland.out$nodes.characteristics)!=0) rownames(rland.out$nodes.characteristics) <- 1:nrow(rland.out$nodes.characteristics)
        
    rland.out$mean.area <- rl$mean.area
    rland.out$SD.area <- rl$SD.area
		
        return(rland.out)
      }
      
      removepoints1 <- function (rl,nr)
      {
        if(nrow(rl$nodes.characteristics)==nr) {
          
          rl$nodes.characteristics <- rl$nodes.characteristics[-(1:nrow(rl$nodes.characteristics)),]
          
          return(rl)
        }
        if(nrow(rl$nodes.characteristics)==0) {
          
          rll_1 <- rl$nodes.characteristics
        }  
        else{
          mapsize4 <- rl$mapsize
          dist_m4 <- rl$minimum.distance
          areaM4 <- rl$mean.area
          areaSD4 <- rl$SD.area
          Npatch4 <- rl$number.patches
          disp4 <- rl$dispersal
          rl_01 <- rl$nodes.characteristics
          ID2 <- rl_01$ID
          nr_select <- nrow(rl_01)-nr
          rl_1 <- rl_01[sample(1:nrow(rl_01), nr_select,replace=FALSE),]
          rl_2 <- rl_1[sort.list(as.numeric(rownames(rl_1))),]
          rl_3 <- list(mapsize=mapsize4, minimum.distance=dist_m4,
                       mean.area=areaM4, SD.area=areaSD4, number.patches=nrow(rl_2),
                       dispersal=disp4,nodes.characteristics=rl_2)
          class(rl_3) <- "landscape"
          rl_4 <- cluster.id(rl_3)
          rownames(rl_4$nodes.characteristics) <- 1:nrow(rl_4$nodes.characteristics)
          rl_4$SD.area <- areaSD4
          rl_4$mean.area <- areaM4
          return(rl_4)
        }
      }
      
      for(i in 2:span){
        tot <- nrow(rland.list[[i-1]])
        rland.list3 <- rland.list[1:length(unlist(lapply(rland.list,nrow)))]
        IDmax <- max(do.call(rbind, rland.list3)[,8])
        create <- tot*(par2/100)
        destroy <- tot*(par3/100)
        p2 <- rpois(1, lambda=create)
        p3 <- rpois(1, lambda=destroy)
        
        rl_m1 <- rland.list[[i-1]]
        
        rl_m1_0 <- list(mapsize=mapsize2, minimum.distance=dist_m2,
                        mean.area=areaM2,SD.area=areaSD2, number.patches=nrow(rl_m1),
                        dispersal=disp2,nodes.characteristics=rl_m1)
        class(rl_m1_0) <- "landscape" 
        if(p3>nrow(rl_m1_0$nodes.characteristics)) p3 <- nrow(rl_m1_0$nodes.characteristics)
        rll_r <- removepoints1(rl_m1_0,nr=p3)
        suppressWarnings(rll_a <- addpoints1(rll_r,nr=p2, IDmax))
        if(nrow(rll_a$nodes.characteristics)!=nrow(rll_r$nodes.characteristics)+p2){
          stop("During simulation one landscape has more patches than it can
               hold!\nConsider making 'dist_m=0' or increasing 'mapsize' in the input landscape.")
        }
        rll_1 <- rll_a$nodes.characteristics
        
        if(nrow(rll_1)==1) rll_1$nneighbour <- NA
        
        rland.list[[i]] <- rll_1
        }
    }
    if (par1=="ncsd")
    {
      removepoints2 <- function (rl,nr,parameter)
      {
        mapsize2 <- rl$mapsize
        dist_m2 <- rl$minimum.distance
        areaM2 <- rl$mean.area
        areaSD2 <- rl$SD.area
        Npatch2 <- rl$number.patches
        disp2 <- rl$dispersal
        rl <- rl$nodes.characteristics
        ID2 <- rl$ID
        ymax <- mapsize2*(1-parameter/100)
        df1_2 <- rl[rl$y < ymax,]
        if(nr > nrow(df1_2))
        {
          ID_delete <- as.numeric(as.vector(df1_2$ID))
          rl1 <- subset(rl,! rl[,8] %in%  ID_delete)
          rl3 <- list(mapsize=mapsize2, minimum.distance=dist_m2, 
                      mean.area=mean(rl1$areas),SD.area=sd(rl1$areas), number.patches=nrow(rl1),
                      dispersal=disp2,nodes.characteristics=rl1)
          class(rl3) <- "landscape"    
          rl2 <- cluster.id(rl3)
        }
        else
        {
          deleted <- sample(rownames(df1_2),nr)
          `%ni%` = Negate(`%in%`)
          df1_3 <- rl[rownames(rl) %ni% deleted,]
          rl1 <- list(mapsize=mapsize2, minimum.distance=dist_m2, 
                      mean.area=mean(df1_3$areas),SD.area=sd(df1_3$areas), number.patches=nrow(df1_3),
                      dispersal=disp2,nodes.characteristics=df1_3)
          class(rl1) <- "landscape"
          rl2 <- cluster.id(rl1)
        }
        rownames(rl2$nodes.characteristics) <- 1:nrow(rl2$nodes.characteristics)
        rl2$mean.area <- areaM2
        rl2$SD.area <- areaSD2
        return(rl2)
      }
      addpoints2 <- function (rl,nr, parameter)
      {
        if (nr!=0)
        {
          mapsize2 <- rl$mapsize
          dist_m2 <- rl$minimum.distance
          areaM2 <- rl$mean.area
          areaSD2 <- rl$SD.area
          Npatch2 <- rl$number.patches
          disp2 <- rl$dispersal
          rl <- rl$nodes.characteristics
          ID2 <- rl$ID
          ymin <- mapsize2*(parameter/100)
          rl2 <- rl[,1:2]
          wind <- owin(xrange=c(0,mapsize2), yrange=c(ymin,mapsize2))
          suppressWarnings(pts_0 <- as.ppp(rl2, W = wind,fatal=TRUE))
          pts_1 <- rSSI(r=dist_m2, n = npoints(pts_0)+nr, win = wind, 
                        giveup = 1000, x.init=pts_0)
          df_pts0 <- as.data.frame(coords(pts_1))
          df_pts1 <- as.data.frame(coords(pts_0))
          df_pts2 <- df_pts0[!duplicated(rbind(df_pts1, df_pts0))[nrow(df_pts1) + 
                                                                    1:nrow(df_pts0)],]
          nrow_0 <- nrow(rl)
          na_lines <- as.data.frame(matrix(NA,nrow=nr,ncol=ncol(rl)))
          colnames(na_lines) <- colnames(rl)
          rownames(na_lines) <- max(as.numeric(rownames(rl)))+1:nrow(na_lines)
          rl <- rbind(rl,na_lines)
          rl[(nrow_0+1):nrow(rl),1:2] <- df_pts2
		  b1 <- (sqrt(areaSD2^2*12)+(areaM2*2))/2
          a <- (areaM2*2) - b1
          areas0 <- runif(nr,a,b1)
		  radius0 <- sqrt((areas0*10000)/pi)
          rl[(nrow_0+1):nrow(rl),"areas"] <- areas0
          rl[(nrow_0+1):nrow(rl),"radius"] <- radius0
          new_ID <- (max(ID2)+1):(max(ID2)+nr)
          rl[,8] <- as.character (c(ID2,new_ID))
          grouping <- hclust(dist(rl[,1:2],method = "euclidean"),"single")
          clusters <- as.data.frame(cutree(grouping, h=disp2))[,1]
          rl[,"cluster"] <- clusters
          col1 <- rainbow(max(rl[,5]))
          col2 <- as.data.frame(col1)
          col2[,2] <- seq(1:max(rl[,5]))
          col3 <- merge_order(rl, col2, by.x = "cluster", 
                              by.y = "V2",sort=FALSE,keep_order=TRUE)[,9]
          rl[,"colour"] <- col3
          col4 <- nndist (rl[,1:2])
          rl[,"nneighbour"] <- col4
          rland.out <- list(mapsize=mapsize2, minimum.distance=dist_m2, 
                            mean.area=mean(rl$areas), SD.area=sd(rl$areas), number.patches=nrow(rl),
                            dispersal=disp2,nodes.characteristics=rl)
          class(rland.out) <- "landscape"
          rland.out$nodes.characteristics$ID <- as.numeric(rland.out$nodes.characteristics$ID)
          rland.out <- cluster.id(rland.out)
        }
        if (nr==0)
        {
          rland.out <- rl
        }
        rownames(rland.out$nodes.characteristics) <- 1:nrow(rland.out$nodes.characteristics) 
        rland.out$mean.area <- rl$mean.area
        rland.out$SD.area <- rl$SD.area
        return(rland.out)
      }
      for(i in 2:span)
      {
        rl0 <- rland.list[[i-1]]
        n_TOT <- nrow(rl0)
        p2 <- n_TOT*(par5/100)
        p3 <- n_TOT*(par4/100)
        n_D <- rpois(1, lambda=p2)
        n_C <- rpois(1, lambda=p3)
        rland0_0 <- list(mapsize=mapsize2, minimum.distance=dist_m2, 
                         mean.area=areaM2,SD.area=areaSD2, number.patches=nrow(rl0),
                         dispersal=disp2,nodes.characteristics=rl0)
        class(rland0_0) <- "landscape"
        rl0_d <- removepoints2(rland0_0, n_D, parameter=par3)
        suppressWarnings(rland0_2 <- addpoints2(rl0_d, n_C, parameter=par2))
        rland0_3 <- cluster.id(rland0_2)
        out <- rland0_3$nodes.characteristics
        rland.list[[i]] <- out
      }
    }
    if (par1=="aggr")
    {
      f5 <- function(df1,parameter)
      {
        vector1 <- sin(df1$x*pi/parameter+df1$y*pi/parameter)+1
        return(vector1)
      }
      f6 <- function(a,...) round(a*(1-(par2/100)))
      vector_kept <- rep(0,span)
      vector_kept[1] <- nrow(rl)
      vector_kept <- Reduce(f6,vector_kept,accumulate=TRUE)
      for(i in 2:span)
      {
        rl0 <- rland.list[[i-1]]
        n_K <- vector_kept[i]
        rl0[,9] <- f5(rl0,par3)
        kept_order <- rl0[,c("ID","V9")]
        kept_order2 <- kept_order[order(-kept_order["V9"]),]
        kept_patches <- as.numeric(as.vector(kept_order2[1:n_K,1]))
        rl0_2 <- rl0[1:8]
        rl0_3 <- subset(rl0_2, rl0_2[,8] %in%  kept_patches)
        rland.1 <- list(mapsize=mapsize2, minimum.distance=dist_m2, 
                        mean.area=mean(rl0_3$areas),SD.area=sd(rl0_3$areas), number.patches=nrow(rl0_3),
                        dispersal=disp2,nodes.characteristics=rl0_3)
        class(rland.1) <- "landscape"
        rland.2 <- cluster.id(rland.1)
        rl0_4 <- rland.2$nodes.characteristics
        rownames(rl0_4) <- 1:nrow(rl0_4)
        rland.list[[i]] <- rl0_4
      }
      unlist(lapply(rland.list,nrow))
    }
    if (par1 == "stoc2")
    {
      addpoints3 <- function (rl,nr)
      {
        if (nr!=0)
        {
          mapsize2 <- rl$mapsize
          dist_m2 <- rl$minimum.distance
          areaM2 <- rl$mean.area
          areaSD2 <- rl$SD.area
          Npatch2 <- rl$number.patches
          disp2 <- rl$dispersal
          rl <- rl$nodes.characteristics
          ID2 <- rl$ID
          rl2 <- rl[,1:2]
          wind <- owin(xrange=c(0,mapsize2), yrange=c(0,mapsize2))
          suppressWarnings(pts_0 <- as.ppp(rl2, W = wind,fatal=TRUE))
          pts_1 <- rSSI(r=dist_m2, n = npoints(pts_0)+nr, win = wind, 
                        giveup = 1000, x.init=pts_0)
          df_pts0 <- as.data.frame(coords(pts_1))
          df_pts1 <- as.data.frame(coords(pts_0))
          df_pts2 <- df_pts0[!duplicated(rbind(df_pts1, df_pts0))[nrow(df_pts1) + 
                                                                    1:nrow(df_pts0)],]
          nrow_0 <- nrow(rl)
          na_lines <- as.data.frame(matrix(NA,nrow=nr,ncol=ncol(rl)))
          colnames(na_lines) <- colnames(rl)
          rownames(na_lines) <- max(as.numeric(rownames(rl)))+1:nrow(na_lines)
          rl <- rbind(rl,na_lines)
          rl[(nrow_0+1):nrow(rl),1:2] <- df_pts2
		  b1 <- (sqrt(areaSD2^2*12)+(areaM2*2))/2
          a <- (areaM2*2) - b1
          areas0 <- runif(nr,a,b1)
		  radius0 <- sqrt((areas0*10000)/pi)
          rl[(nrow_0+1):nrow(rl),"areas"] <- areas0
          rl[(nrow_0+1):nrow(rl),"radius"] <- radius0
          new_ID <- (max(ID2)+1):(max(ID2)+nr)
          rl[,8] <- as.character (c(ID2,new_ID))
          grouping <- hclust(dist(rl[,1:2],method = "euclidean"),"single")
          clusters <- as.data.frame(cutree(grouping, h=disp2))[,1]
          rl[,"cluster"] <- clusters
          col1 <- rainbow(max(rl[,5]))
          col2 <- as.data.frame(col1)
          col2[,2] <- seq(1:max(rl[,5]))
          col3 <- merge_order(rl, col2, by.x = "cluster", 
                              by.y = "V2",sort=FALSE,keep_order=TRUE)[,9]
          rl[,"colour"] <- col3
          col4 <- nndist (rl[,1:2])
          rl[,"nneighbour"] <- col4
          rland.out <- list(mapsize=mapsize2, minimum.distance=dist_m2, 
                            mean.area=mean(rl$areas), SD.area=sd(rl$areas),
                            number.patches=nrow(rl),dispersal=disp2,
                            nodes.characteristics=rl)
          class(rland.out) <- "landscape"
          rland.out <- cluster.id(rland.out)
        }
        if (nr==0)
        {
          rland.out <- rl
        }
        rownames(rland.out$nodes.characteristics) <- 1:nrow(rland.out$nodes.characteristics) 
        return(rland.out)
      }
      removepoints3 <- function (rl,nr)
      {
        mapsize2 <- rl$mapsize
        dist_m2 <- rl$minimum.distance
        areaM2 <- rl$mean.area
        areaSD2 <- rl$SD.area
        Npatch2 <- rl$number.patches
        disp2 <- rl$dispersal
        rl_0 <- rl$nodes.characteristics
        ID2 <- rl_0$ID
        nr_select <- nrow(rl_0)-nr
        rl_1 <- rl_0[sample(1:nrow(rl_0), nr_select,replace=FALSE),]
        rl_2 <- rl_1[sort.list(as.numeric(rownames(rl_1))),]
        rl_3 <- list(mapsize=mapsize2, minimum.distance=dist_m2, 
                     mean.area=mean(rl_2$areas), SD.area=sd(rl_2$areas), number.patches=nrow(rl_2),
                     dispersal=disp2,nodes.characteristics=rl_2)
        class(rl_3) <- "landscape"
        rl_4 <- cluster.id(rl_3)
        rownames(rl_4$nodes.characteristics) <- 1:nrow(rl_4$nodes.characteristics)
        return(rl_4)
      }
      for(m in 2:span)
      {
        tot <- nrow(rland.list[[m-1]])
        p2 <- round((tot*par2)/100)
        rl_m1 <- rland.list[[m-1]]
        rl_m1_0 <- list(mapsize=mapsize2, minimum.distance=dist_m2, 
                        mean.area=mean(rl_m1$areas),SD.area=sd(rl_m1$areas), number.patches=nrow(rl_m1),
                        dispersal=disp2,nodes.characteristics=rl_m1)
        class(rl_m1_0) <- "landscape"    
        rll_r <- removepoints3(rl_m1_0,nr=p2)
        rll_a <- addpoints3(rll_r,nr=p2)
        rll_1 <- rll_a$nodes.characteristics
        rland.list[[m]] <- rll_1
      }
    }
    rland.list2 <- rland.list[1:length(unlist(lapply(rland.list,nrow)))]
    if (length(rland.list2) < span)
    {
      cat("Only ",length(rland.list2)," landscape(s) generated.","\n")
    }
    
    for(i in 1:length(rland.list2)){
      df01 <- rland.list2[[i]]
      if(nrow(df01)!=0){
      xy <- df01[,1:2]
      d1 <- pairdist(xy)
      d1[d1 == 0] <- NA
      nneigh <- rep(NA,nrow(d1))
      for(j in 1:length(nneigh)){
        r1 <- d1[j,]
        val1 <- suppressWarnings(min(r1, na.rm=TRUE))
        if(is.finite(val1)) nneigh[j]<-min(r1, na.rm=TRUE)
        if(!is.finite(val1)) nneigh[j]<-0
      }
      }
      if(nrow(df01)!=0)df01[,7] <- nneigh
      rland.list2[[i]] <- df01
    }
    
return(rland.list2)

}