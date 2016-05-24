pred.strk <- function (temp,
                       zcol=1,
                       newdata, # only STFDF
                       pred.id= 'tempPred',
                       zero.tol=0,
                       dynamic.cov=c(1,2),
                       static.cov=c(1,2),
                       reg.coef=list(tmean= c(-0.126504415,0.4051734447,0.4943247727,0.0001837527,-0.0189207588), # Intercept, temp_geo,modis,dem,twi 
                                     tmin = c(-0.9825601517,0.5672140021,0.3344561638, 0.0003119777,-0.0243629638),
                                     tmax = c(1.7873573081,0.350228076, 0.5569091092, 0.0002571338,-0.0012988123) ) [[ 'tmean']] ,
                       vgm.model=list(tmean=vgmST("sumMetric",
                                                  space=vgm( 14.13, "Sph", 5903, 1.933),
                                                  time =vgm(0, "Sph",  0.1, 0),
                                                  joint=vgm(9.06, "Sph", 2054, 0.474),
                                                  stAni=497.9),
                                      tmin = vgmST("sumMetric",
                                                  space=vgm( 22.682, "Sph", 5725, 3.695),
                                                  time =vgm(0, "Sph",  0.1, 0),
                                                  joint=vgm(9.457, "Sph",1888, 1.67),
                                                  stAni=485),
                                     tmax = vgmST("sumMetric",
                                                  space=vgm( 8.31, "Sph", 4930, 2.872),
                                                  time =vgm(0, "Sph",  0.1, 0),
                                                  joint=vgm(11.175, "Sph", 2117, 1.75),
                                                  stAni=527) ) [['tmean']] ,
                       tiling= FALSE,
                       ntiles=64,
                       parallel.processing = FALSE,
                       cpus=3,
                       sp.nmax=18,
                       time.nmax=2,
                       fast= FALSE,
                       computeVar=FALSE,
                       do.cv= FALSE,
                       only.cv = FALSE,
                       out.remove = FALSE,
                       threshold.res=15,
                       progress=TRUE){
 
  temp <- rm.dupl(temp, zcol,zero.tol)
  gg <- newdata
  time <- gg@time
  # dynamic predictors  
  if(!is.null(dynamic.cov)) {
   dyn<- as.data.frame ( gg@data[,dynamic.cov] )
  }else{ dyn<-NA } 
  if(only.cv){ do.cv=TRUE}
  # static predictors
  if(!is.null(static.cov)  ) {
    sta <- lapply(static.cov, function(i) rep(gg@sp@data[,i],length(time)) )
    sta <- as.data.frame( do.call(cbind,sta) )
    names(sta)<- names(gg@sp@data[,static.cov])
  }else{ sta<- NA }
  
  df= cbind(dyn,sta)
  df <- df[,colSums(is.na(df))<nrow(df)]
  
  if(is.null(df)){
     gg$tlm<-0 
     temp$tres <- temp@data[,zcol]
     } else {
       
       gg$tlm= reg.coef[1] + as.matrix(df )  %*%  reg.coef[-1] #regression model-trend
       gg$tlm<- as.vector(gg$tlm)
       nrowsp <- length(temp@sp)
       gg@sp=as(gg@sp,'SpatialPixelsDataFrame')
       ov <- sapply(1:length(time),function(i) over(temp@sp,as(gg[,i,'tlm'],'SpatialPixelsDataFrame' ) ) )
       ov <- do.call('cbind',ov)
       ov <- as.vector(ov)
       
       t1 <- which(as.character(index(time[1])) == as.character( index(temp@time)) )
       t2 <- which(as.character(index(time[ length(time) ])) == as.character(index(temp@time)) )
       
       
         temp <- temp[,t1:t2, drop=F]
       
     
       
       temp$tlm <- ov
       temp$tres <- temp@data[,zcol]- temp$tlm #residuals
       
       # count NAs per stations
       numNA <- apply(matrix(temp@data[,'tres'],
                             nrow=nrowsp,byrow=F), MARGIN=1,
                      FUN=function(x) sum(is.na(x)))
       
       # Remove stations out of covariates
       rem <- numNA != length(time)
       
         temp <-  temp[rem,drop=F]
      
       
       
#        row.names(temp@sp) <- 1: nrow(temp@sp)
       
       
     } # end of regression 
  if( sp.nmax > nrow(temp@sp) ) {sp.nmax <- nrow(temp@sp) }
  
  i_1 <- (1:length(time)) - ceiling(time.nmax/2)
  i_1[i_1<1]=1        
  ip1 <- i_1 + floor(time.nmax/2)
  ip1[ip1>length(time)] <- length(time)  
  
  
  if(parallel.processing) {
    sfInit ( parallel = parallel.processing , cpus =cpus)
    sfLibrary(package="gstat", character.only=TRUE)
    sfExport("vgm.model" )
    sfExport( "computeVar" )
    sfExport( "i_1" )
    sfExport( "ip1" )
    sfExport( "time" )
    sfExport( "temp" )
  }
  
#   

  gg@sp=as(gg@sp,'SpatialPointsDataFrame')
  row.names(gg@sp) = 1:nrow(gg@sp)
  
  remove <- rep(FALSE,length(temp@sp))
  robs=NA
  idsOUT= NA
  
  if(do.cv){
    
    N_POINTS <- length(temp@sp@coords[,1])
    
    if(parallel.processing) {
    sfExport( "temp" )
    sfExport( "N_POINTS" )
    sfExport( "sp.nmax" )   }
    
    
    cv = as.list ( rep(NA, N_POINTS)) 
    if (progress)   pb <- txtProgressBar(style = 3,char= sprintf("cv ") , max=N_POINTS)
    
    cv =   (if (parallel.processing) sfLapply else lapply )(1:N_POINTS, function(i) 
    {
      st= temp@sp
      st$dist=spDists(temp@sp,temp@sp[i,])
      
      tmp_st<-st[ order(st$'dist') , ]
      
      local_t= row.names(tmp_st[2:sp.nmax,] ) # remove target station
      
      xxx = as.list ( rep(NA, length(time) ) )
      for( ii in 1:length(time) ) {
        xxx[[ii]]=krigeST(as.formula("tres~1"),
                          data=as(temp[local_t, i_1[ii]:ip1[ii],'tres',drop=F], "STSDF"),
                          newdata=STF(as(temp@sp[i,],"SpatialPoints"),
                                      temp@time[ii],  
                                      temp@endTime[ii]),     
                          modelList=vgm.model,
                          computeVar=computeVar)@data[,1]
        
      } # end of  for
      
      ret=unlist(xxx) 
      if (progress)  setTxtProgressBar(pb, i )
      ret } ) # end of (if (parallel.processing) sfLapply else lapply )
    if (progress)  close(pb)
    cv <- do.call(rbind,cv)
    cv <- as.vector(cv)
    cv.temp <- temp
    cv.temp$pred.cv <- cv + cv.temp$tlm
    cv.temp$resid.cv <- cv.temp$pred.cv  - cv.temp@data[,zcol]
    
    temp<- cv.temp
    
    bb <- cv.temp
    bb <- as.data.frame(bb)
    bb$sp.ID <-as.character(bb$sp.ID )
    bb$abs.res <-abs(bb$resid.cv )
    bb <- bb[order(bb$abs.res,  decreasing = TRUE),]
    idsOUT <- unique(bb[which( bb$abs.res > threshold.res),]$sp.ID) 
#     idsOUT <- as.numeric (idsOUT) 
    #     stplot( temp[idsOUT,,c('tempc','pred.cv')] , mode='tp')

    if(out.remove==TRUE & length(idsOUT)!=0){
      if (progress)   pb <- txtProgressBar(style = 3,char= sprintf("remove OUT. "), max=length(temp@sp@coords[,1]) )
      while(length(idsOUT)>0){ 
        
        rm.station = which(row.names(temp@sp)==idsOUT[1])
        remove [rm.station ] <- TRUE
        
        rm.station = which(row.names(cv.temp@sp)==idsOUT[1])
        cv.temp <- cv.temp[ - rm.station , drop=F]
        
        N_POINTS <- length(cv.temp@sp@coords[,1])
        
        if(parallel.processing) {
          sfExport( "cv.temp" )
          sfExport( "N_POINTS" )
          sfExport( "sp.nmax" )  }
        
        
        cv = as.list ( rep(NA, N_POINTS)) 
        
        cv =   (if (parallel.processing) sfLapply else lapply )(1:N_POINTS, function(i) 
        {
          st= cv.temp@sp
          st$dist=spDists(cv.temp@sp,cv.temp@sp[i,] )
          
          tmp_st<-st[ order(st$'dist') , ]
          if( sp.nmax > nrow(tmp_st) ) {sp.nmax <- nrow(tmp_st) }
          local_t= row.names(tmp_st[2:sp.nmax,] ) # remove target station
          
          xxx = as.list ( rep(NA, length(time) ) )
          for( ii in 1:length(time) ) {
            xxx[[ii]]=krigeST(as.formula("tres~1"),
                              data=as(temp[local_t, i_1[ii]:ip1[ii],'tres',drop=F], "STSDF"),
                              newdata=STF(as(cv.temp@sp[i,],"SpatialPoints"),
                                          cv.temp@time[ii],  
                                          cv.temp@endTime[ii]),     
                              modelList=vgm.model,
                              computeVar=computeVar)@data[,1]
            
          } # end of  for

          if (progress)  setTxtProgressBar(pb, i )
          unlist(xxx)  } ) # end of (if (parallel.processing) sfLapply else lapply )
        
        cv <- do.call(rbind,cv)
        cv <- as.vector(cv)
        cv.temp$pred.cv <- cv + cv.temp$tlm
        cv.temp$resid.cv <- cv.temp$pred.cv  - cv.temp@data[,zcol]
        
        bb <- cv.temp
        bb <- as.data.frame(bb)
        bb$sp.ID <-as.character(bb$sp.ID )
        bb$abs.res <-abs(bb$resid.cv )
        bb <- bb[order(bb$abs.res,  decreasing = TRUE),]
        idsOUT <- unique(bb[which( bb$abs.res > threshold.res),]$sp.ID) 
#         idsOUT <- as.numeric (idsOUT) 
     
      } # while
      if (progress)  close(pb)
      remove <- temp@sp[remove,]
      robs <- temp[remove,drop=F]
      temp <- cv.temp
    }# remove
    
    
    if(only.cv){
      sfStop()
      return(list(pred=NA,cv =cv.temp, out= idsOUT, remst =remove, remobs=robs) )
    }
    
  }else{ 
    cv.temp =NA
    idsOUT=NA}
  
  if (!tiling) { #prediction
    
    temp.local<-temp[,,'tres',drop=F]
    
    if(parallel.processing) {
      sfExport( "temp.local" )
      sfExport( "gg" )    } 
    
    if (progress)   pb <- txtProgressBar(style = 3,char= sprintf("pred krigeST ") , max=length(time) )
    xxx<- (if (parallel.processing & length(time)>1) sfLapply else lapply )(1:length(time), function(i) {

      obs=temp.local[,i_1[i]:ip1[i],'tres', drop=F]

      pred.var=1
      if(computeVar==TRUE){ pred.var=2}
      
        
      ret= krigeST(as.formula("tres~1"),
              data=as(obs, "STSDF"), # [,1:6] # in this case I must to limit for the fist few days
              newdata=STF(as(gg@sp,"SpatialPoints"),
                          temp.local@time[i],    # [3]    
                          #srb@data,  # [3,1]
                          temp.local@endTime[i]),    # [3]    
              modelList=vgm.model,
              computeVar=computeVar)@data[,pred.var]
      if (progress)  setTxtProgressBar(pb, i )
      ret } )
    if (progress)  close(pb)

    res=do.call(cbind,xxx)  
    res= as.vector(res)
    gg$temp.pred<- res + gg$tlm
  
    stfdf <-gg


  }else{  # do tiling
    dimnames(gg@sp@coords)[[2]] <- c('x','y')
    xy=as.data.frame(gg@sp)
    xy= xy[row.names(gg@sp),]
    nxy= floor(sqrt(ntiles) )
    xy$xg=as.character( cut(xy$x,nxy,labels=paste("x",1:nxy,sep="") ) )
    xy$yg=as.character( cut(xy$y,nxy,labels=paste("y",1:nxy,sep="") )  )
    #  xy=xy[order(-xy$yg, xy$xg),]
    xy$g=as.factor(paste(xy$xg,xy$yg,sep="") )
    xy$xg=NULL
    xy$yg=NULL
    coordinates(xy) = ~ x+y
    xy$index=1:nrow(xy)
    xy@proj4string <- gg@sp@proj4string
    g_list <- split(xy, xy$g)
    
    # g_list= lapply(g_list,function(i) g[i,])
    
    # mAKE CHUNKS OF INITIAL DATA
    Mpoint=data.frame(x=mean(gg@sp@coords[,1]),y=mean(gg@sp@coords[,2]) )
    coordinates(Mpoint)=~x+y
    Mpoint@proj4string <- gg@sp@proj4string
    st=temp@sp
    
    
#     # Find nearest 500 stations
#     #
#     st$dist=spDists(temp@sp,Mpoint)
#     tmp_st<-st[ order(st$'dist') ,]
#     local1= row.names(tmp_st[1:(sp.nmax+1),] )
    
    
    # Middle point for each chunk
    Mpts= (if (parallel.processing) sfLapply else lapply )(g_list,function(i) {
      Mpoint=data.frame(x=mean(i@coords[,1]),y=mean(i@coords[,2]) )
      coordinates(Mpoint)=~x+y
      Mpoint@proj4string <- temp@sp@proj4string
      Mpoint   })
    
    # Local obs.
#     temp.local <-temp[local1, ,'tres']
    temp.local <-temp[ , ,'tres',drop=F]
    st=temp@sp
    
    if(parallel.processing) {
      sfExport( "temp.local" )
      sfExport( "st" )
      sfExport( "Mpts" )
      sfExport( "g_list" )
      sfExport( "sp.nmax" ) }

    if (progress)   pb <- txtProgressBar(style = 3,char= sprintf("pred krigeST ") , max=length(g_list))

    res =   (if (parallel.processing) sfLapply else lapply )(1:length(g_list), function(i) 
    {
      pred.var = 1
      if(computeVar==TRUE){ pred.var=2}
      
      st$dist=spDists(temp.local@sp,Mpts[[i]])
      
      tmp_st<-st[ order(st$'dist') ,]
      
      local_t= row.names(tmp_st[1:sp.nmax,] )
      
      xxx = as.list ( rep(NA, length(time) ) )
      for( ii in 1:length(time) ) {
        xxx[[ii]]=krigeST(as.formula("tres~1"),
                          data=as(temp[local_t, i_1[ii]:ip1[ii],'tres',drop=F], "STSDF"), 
                          newdata=STF(as(g_list[[i]],"SpatialPoints"),
                                      temp.local@time[ii],  
                                      temp.local@endTime[ii]),     
                          modelList=vgm.model,
                          computeVar=computeVar)@data[,pred.var]
      } # end of  for
      
      ret=do.call(cbind,xxx)
      if (progress)  setTxtProgressBar(pb, i )
      ret } ) # end of (if (parallel.processing) sfLapply else lapply )
    if (progress)  close(pb)
    res=do.call(rbind,res)
    cc=do.call(rbind,g_list)
    g=cc
    
      stfdf= gg[as.character(g$index),drop=F]    
      res <- as.numeric(res)
      stfdf$temp.pred <- res +stfdf$tlm
    
    if(!fast){
      stfdf1= stfdf
      
      dimnames(gg@sp@coords)[[2]] <- c('x','y')
      xy=as.data.frame(gg@sp)
      xy= xy[row.names(gg@sp),]
      nxy= floor(sqrt(ntiles+ ntiles) )
      xy$xg=as.character( cut(xy$x,nxy,labels=paste("x",1:nxy,sep="") ) )
      xy$yg=as.character( cut(xy$y,nxy,labels=paste("y",1:nxy,sep="") )  )
      #  xy=xy[order(-xy$yg, xy$xg),]
      xy$g=as.factor(paste(xy$xg,xy$yg,sep="") )
      xy$xg=NULL
      xy$yg=NULL
      coordinates(xy) = ~ x+y
      xy@proj4string <- gg@sp@proj4string
      xy$index=1:nrow(xy)
      g_list <- split(xy, xy$g)
      
      # g_list= lapply(g_list,function(i) g[i,])
      
      # mAKE CHUNKS OF INITIAL DATA
      Mpoint=data.frame(x=mean(gg@sp@coords[,1]),y=mean(gg@sp@coords[,2]) )
      coordinates(Mpoint)=~x+y
      Mpoint@proj4string <- gg@sp@proj4string
      st=temp@sp

      
      
      # Middle point for each chunk
      Mpts= (if (parallel.processing) sfLapply else lapply )(g_list,function(i) {
        Mpoint=data.frame(x=mean(i@coords[,1]),y=mean(i@coords[,2]) )
        coordinates(Mpoint)=~x+y
        Mpoint@proj4string <- temp@sp@proj4string
        Mpoint   })
      
      # Local obs.
      #     temp.local <-temp[local1, ,'tres']
      temp.local <-temp[ , ,'tres', drop = FALSE]
      st=temp.local@sp
      if(parallel.processing) {
        sfExport( "temp.local" )
        sfExport( "st" )
        sfExport( "Mpts" )
        sfExport( "g_list" )
        sfExport( "sp.nmax" )
        sfExport( "computeVar" ) }
      if (progress)   pb <- txtProgressBar(style = 3,char= sprintf("pred krigeST "), max=length(g_list) )
      
      res =   (if (parallel.processing) sfLapply else lapply )(1:length(g_list), function(i) 
      {
        pred.var = 1
        if(computeVar==TRUE){ pred.var=2}
        
        st$dist=spDists(temp@sp,Mpts[[i]] ) ###kkk
        
        tmp_st<-st[ order(st$'dist') ,]
        
        local_t= row.names(tmp_st[1:sp.nmax,] )
        
        xxx = as.list ( rep(NA, length(time) ) )
        for( ii in 1:length(time) ) {
          xxx[[ii]]=krigeST(as.formula("tres~1"),
                            data=as(temp[local_t, i_1[ii]:ip1[ii],'tres',drop=F], "STSDF"), 
                            newdata=STF(as(g_list[[i]],"SpatialPoints"),
                                        temp.local@time[ii],  
                                        temp.local@endTime[ii]),     
                            modelList=vgm.model,
                            computeVar=computeVar)@data[,pred.var]
        
        } # end of  for
        
        ret = do.call(cbind,xxx)  
        if (progress)  setTxtProgressBar(pb, i )
        ret } ) # end of (if (parallel.processing) sfLapply else lapply )
      if (progress)  close(pb)
      res=do.call(rbind,res)
      cc=do.call(rbind,g_list)
      g=cc
      
      stfdf= gg[as.character(g$index), drop=F]    
      res <- as.numeric(res)
      stfdf$temp.pred <- res +stfdf$tlm
      
      stfdf<- stfdf[ row.names(stfdf1@sp)  , drop=F]
      
      stfdf$temp.pred <- 0.5*stfdf$temp.pred + 0.5*stfdf1$temp.pred
      
      
    } # end of fast
    
  } # and of tiling else do tiling
  
  stfdf@sp <- as(stfdf@sp, class(newdata@sp) )
   
  stfdf <- STFDF(stfdf@sp,stfdf@time, data.frame(stfdf$temp.pred), endTime=stfdf@endTime)
  
  names(stfdf@data) <- pred.id 
  
if (parallel.processing){
  sfStop() 
}
  

  return(list(pred=stfdf,cv =cv.temp, out= idsOUT, remst =remove, remobs=robs) )
}
  
    

  
  
  

