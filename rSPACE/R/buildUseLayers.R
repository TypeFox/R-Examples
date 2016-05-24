# Create a use surface from a list of wolverine locations ---------
build.useLayer<-function(map, wolv, Parameters, Example=F){
  nTypes <- length(Parameters$MFratio)
  NotUsed<-1-use_surface(Wolv = wolv[[1]],
              howmuch = Parameters$moveDistQ[1],
              howfar  = Parameters$moveDist[1],
              map     = map,
              trunk   = Parameters$maxDistQ[1])

  if(!Example & nTypes > 1){
    for(ii in 2:nTypes){
      NotUsed<-NotUsed*(1-use_surface(Wolv = wolv[[ii]],
                            howmuch = Parameters$moveDistQ[ii],
                            howfar  = Parameters$moveDist[ii],
                            map     = map,
                            trunk   = Parameters$maxDistQ[ii]))
    }
  }

  if(!is.null(Parameters$repeat.groups))
    if(Parameters$repeat.groups == T)
      NotUsed<-NotUsed^2

  useLayer<-1 - NotUsed
  return(useLayer)
}

# Single group probability of use --------------------------------
use_surface<-function(Wolv, howmuch, howfar, map, trunk = trunk){
  sdXY<-solveSD(howmuch, howfar, map)

  if(trunk > 1 | trunk < 0)
    stop('Invalid value for trunction')
  trunk<-ifelse(trunk==1, 0, qnorm((1 + trunk) / 2))

  USE = .C('use_surface',
              as.double(coordinates(map)[Wolv,1]),   #x_wolv
              as.double(coordinates(map)[Wolv,2]),   #y_wolv
              as.integer(length(Wolv)),              #N_wolv

              as.double(coordinates(map)[,1]),       #x_snow
              as.double(coordinates(map)[,2]),       #y_snow
              use = as.double(getValues(map)),       #snow
              as.integer(length(getValues(map))),    #pixels

              as.double(c(sdXY[,1])),                # double sd_x[]
              as.double(c(sdXY[,2])),                # double sd_y[]
              as.double(c(trunk))                    # double trunc_cutoff[]
              )$use
  return(USE)
  }


# Solve for variance of movement distribution -----------------
solveSD<-function(howmuch, howfar, map){
  isUTM     <- grepl('+proj=utm', proj4string(map))
  isLongLat <- grepl('+proj=longlat', proj4string(map))

  if(isUTM){  # km per map unit conversion
    km_per_x <- 1/1000
    km_per_y <- 1/1000
  } else if(isLongLat){
    mid <- c(mean(c(slot(extent(map), 'xmin'), slot(extent(map),'xmax'))),
             mean(c(slot(extent(map), 'ymin'), slot(extent(map),'ymax'))))
    km_per_x = pointDistance(mid, mid + c(1, 0), lonlat=T) / 1000
    km_per_y = pointDistance(mid, mid + c(0, 1), lonlat=T) / 1000
  }

  cutoff <- qnorm((1+howmuch)/2)     # Cutoff for a standard normal distribution
  howfar_x <- howfar / km_per_x      # Convert km to map units
  howfar_y <- howfar / km_per_y
  sd_x <- howfar_x / cutoff          # Convert cutoff to map scale to get sd
  sd_y <- howfar_y / cutoff

   return(cbind(sd_x,sd_y))
}