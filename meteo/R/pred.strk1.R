pred.strk1 <- function (obs,
                        stations,
                        newdata,
                        zero.tol=0,
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
                        computeVar=FALSE,
                        out.remove=FALSE,
                        threshold.res=15,
                        returnList=FALSE){
  

  
  temp <- meteo2STFDF(obs,stations)
  temp <- rm.dupl(temp, 1,zero.tol)
  pr=as.data.frame ( temp@data[,2:length(temp@data)] )  #   covariates
  temp$tlm= reg.coef[1] + as.matrix(pr)  %*%  reg.coef[-1]
  temp$tlm<- as.vector(temp$tlm)
  temp$tres <- temp@data[,1]- temp$tlm
  
  nrowsp <- length(temp@sp)
  timeT = temp@time
  # count NAs per stations
  numNA <- apply(matrix(temp@data[,'tres'],
                        nrow=nrowsp,byrow=F), MARGIN=1,
                 FUN=function(x) sum(is.na(x)))
  
  # Remove stations out of covariates
  rem <- numNA != length(timeT)
  
  temp <-  temp[rem,drop=F]
  
  
  
  #out.remove cv
  remove <- rep(FALSE,length(temp@sp))
  robs=NA
  idsOUT= NA
  cv.temp =NA
  
  if(out.remove){
    
    
    N_POINTS <- length(temp@sp@coords[,1])
    
    cv = as.list ( rep(NA, N_POINTS)) 
        
    cv = lapply (1:N_POINTS, function(i) 
    {
      st= temp@sp
      
      #   attr(vgm.model,"temporal unit") = "days"
      xxx = as.list ( rep(NA, length(time) ) )
      for( ii in 1:length(time) ) {
        xxx[[ii]]=krigeST(as.formula("tres~1"),
                          data=temp[-i, ,'tres',drop=F], 
                          newdata=STF(as(temp@sp[i,],"SpatialPoints"),
                                      temp@time[ii],  
                                      temp@endTime[ii]),     
                          modelList=vgm.model,
                          computeVar=computeVar)@data[,1]
        
      } # end of  for
      
      ret=unlist(xxx) 
      ret } ) # end of (if (parallel.processing) sfLapply else lapply )
    
    cv <- do.call(rbind,cv)
    cv <- as.vector(cv)
    cv.temp <- temp
    cv.temp$pred.cv <- cv + cv.temp$tlm
    cv.temp$resid.cv <- cv.temp$pred.cv  - cv.temp@data[,1]
    
    temp<- cv.temp
    
    bb <- cv.temp
    bb <- as.data.frame(bb)
    bb$sp.ID <-as.character(bb$sp.ID )
    bb$abs.res <-abs(bb$resid.cv )
    bb <- bb[order(bb$abs.res,  decreasing = TRUE),]
    idsOUT <- unique(bb[which( bb$abs.res > threshold.res),]$sp.ID) 
    #     idsOUT <- as.numeric (idsOUT) 
    #     stplot( temp[idsOUT,,c('tempc','pred.cv')] , mode='tp')
    
    if(length(idsOUT)!=0){
      
      while(length(idsOUT)>0){ 
        
        rm.station = which(row.names(temp@sp)==idsOUT[1])
        remove [rm.station ] <- TRUE
        
        rm.station = which(row.names(cv.temp@sp)==idsOUT[1])
        cv.temp <- cv.temp[ - rm.station , drop=F]
        
        N_POINTS <- length(cv.temp@sp@coords[,1])
        
        cv = as.list ( rep(NA, N_POINTS)) 
        
        cv =  lapply(1:N_POINTS, function(i) 
        {
          st= cv.temp@sp
                    
          xxx = as.list ( rep(NA, length(time) ) )
          for( ii in 1:length(time) ) {
            xxx[[ii]]=krigeST(as.formula("tres~1"),
                              data=temp[-i, ,'tres', drop=F],
                              newdata=STF(as(cv.temp@sp[i,],"SpatialPoints"),
                                          cv.temp@time[ii],  
                                          cv.temp@endTime[ii]),     
                              modelList=vgm.model,
                              computeVar=computeVar)@data[,1]
            
          } # end of  for
          
          
          unlist(xxx)  } ) # end of (if (parallel.processing) sfLapply else lapply )
        
        cv <- do.call(rbind,cv)
        cv <- as.vector(cv)
        cv.temp$pred.cv <- cv + cv.temp$tlm
        cv.temp$resid.cv <- cv.temp$pred.cv  - cv.temp@data[,1]
        
        bb <- cv.temp
        bb <- as.data.frame(bb)
        bb$sp.ID <-as.character(bb$sp.ID )
        bb$abs.res <-abs(bb$resid.cv )
        bb <- bb[order(bb$abs.res,  decreasing = TRUE),]
        idsOUT <- unique(bb[which( bb$abs.res > threshold.res),]$sp.ID) 
        #         idsOUT <- as.numeric (idsOUT) 
        
      } # while
      
      remove <- temp@sp[remove,]
      robs <- temp[remove,drop=F]
      temp <- cv.temp
    }# remove
    
  }
  
  
  
  
  
  nrow=nrow(newdata)
  ncol=ncol(newdata)
 
  
  gg <- newdata
  time <- as.POSIXlt(gg[1,3],"GMT")
    
  df=gg[,4:ncol]
  df <- df[,colSums(is.na(df))<nrow(df)]
  
  gg$tlm= reg.coef[1] + as.matrix(df )  %*%  reg.coef[-1]
  gg$tlm <- as.vector(gg$tlm)
      
  names(gg)[1] = "x"
  names(gg)[2] = "y"
  coordinates(gg) <- ~x+y
    
    
    
  pred.var=1
  if(computeVar==TRUE){ pred.var=2}
      
      
  res= krigeST(as.formula("tres~1"),
               data=temp[,,'tres', drop=F], # [,1:6] # in this case I must to limit for the fist few days
               newdata=STF(as(gg,"SpatialPoints"),
                           time,    # [3]    
                           #srb@data,  # [3,1]
                           time+86400),    # [3]    
               modelList=vgm.model,
               computeVar=computeVar)@data[,pred.var]
      
       
# res= as.vector(res)
  gg$pred<- res + gg$tlm
  gg$pred <- as.vector(gg$pred)

  stfdf <-gg
  
  if (returnList){
    return(list(pred=stfdf,cv =cv.temp, remst =remove, remobs=robs) )
  } else {
    if (computeVar){
      return(data.frame(newdata[,1:3], var=stfdf$pred) )
    } else {
      return(data.frame(newdata[,1:3], pred=stfdf$pred) )
    }
    
  }
}
