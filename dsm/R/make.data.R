#' @importFrom stats aggregate
make.data <- function(response, ddfobject, segdata, obsdata, group,
                      convert.units, availability, strip.width, segment.area){

  # probably want to do something smart here...
  seglength.name<-'Effort'
  segnum.name<-'Sample.Label'
  distance.name<-'distance'
  cluster.name<-'size'

  # Check that observations are between left and right truncation
  # warning only -- observations are excluded below
  # No truncation check for strip transects
  if(!is.null(ddfobject) & !is.null(obsdata[[distance.name]])){
    if(any(obsdata[,distance.name]>ddfobject$meta.data$width)){
      warning("Some observations are outside of detection function truncation!")
    }
  }

  # Estimating group abundance/density
  if(group){
    obsdata[,cluster.name][obsdata[,cluster.name]>0] <- 1
  }

  # if we fitted a detection function
  if (!is.null(ddfobject)){
    # grab the probabilities of detection
    fitted.p <- fitted(ddfobject)

    # remove observations which were not in the detection function
    obsdata <- obsdata[obsdata$object %in% names(fitted.p),]
    # what if there are no matches? Perhaps this is due to the object
    # numbers being wrong? (HINT: yes.)
    if(nrow(obsdata) == 0){
      stop("No observations in detection function matched those in observation table. Check the \"object\" column.")
    }
  }else{
    # strip transects or presence/absence data
    fitted.p <- 1
    if(is.null(strip.width) & (response != "presence")){
      stop("You must specify strip.width for strip transects!")
    }
  }

  # if the ps are all the same (count model) then just grab the 1 unique
  # value
  if(response %in% c("N","abundance","count","n")){
    fitted.p <- unique(fitted.p)
  }

  ## Aggregate response values of the sightings over segments
  if(response %in% c("D","density","Dhat","density.est")){
    responsedata <- aggregate(obsdata[,cluster.name]/(fitted.p*availability),
                                list(obsdata[,segnum.name]), sum)
    off.set <- "none"
  }else if(response %in% c("N","abundance","count","n")){
    responsedata <- aggregate(obsdata[,cluster.name]/availability,
                              list(obsdata[,segnum.name]), sum)
    off.set <- "eff.area"
  }else if(response %in% c("Nhat","abundance.est")){
    responsedata <- aggregate(obsdata[,cluster.name]/(fitted.p*availability),
                              list(obsdata[,segnum.name]), sum)
    off.set<-"area"
  }else if(response == "presence"){
    responsedata <- aggregate(obsdata[,cluster.name],
                                list(obsdata[,segnum.name]), sum)
    responsedata$x[responsedata$x>0] <- 1
    responsedata$x[responsedata$x<1] <- 0
    off.set <- "none"
  }

  ## warn if any observations were not allocated
  responsecheck <- aggregate(obsdata[,cluster.name],
                             list(obsdata[,segnum.name]), sum)
  if(sum(obsdata[,cluster.name]) != sum(responsecheck[,2])){
    message(paste0("Some observations were not allocated to segments!\n",
                   "Check that Sample.Labels match"))
  }

  # name the response data columns
  names(responsedata)<-c(segnum.name,response)

  # Next merge the response variable with the segment records and any
  # response variable that is NA should be assigned 0 because these
  # occur due to 0 sightings
  dat <- merge(segdata,responsedata,by=segnum.name,all.x=T)
  dat[,response][is.na(dat[,response])] <- 0

  if(!is.null(segment.area)){

    # pull the column if segment.area is character
    if(is.character(segment.area)){
      segment.area <- dat[,segment.area]
    }

    dat$off.set <- switch(off.set,
                          eff.area=segment.area*fitted.p,
                          area=segment.area,
                          none=1)

    # if we have density then use that as the response
    if(response %in% c("D","density","Dhat","density.est")){
      dat[,response] <- dat[,response]/(segment.area*convert.units)
    }

    # set the segment area in the data
    dat$segment.area <- segment.area*convert.units

  }else{
    # pull this from the detection function
    if (!is.null(ddfobject)){
      # calculate the "width" of the transect first, make sure we get it right
      # if we are doing left truncation
      width <- ddfobject$meta.data$width
      if(!is.null(ddfobject$meta.data$left)){
        width <- width - ddfobject$meta.data$left
      }
    }else{
    # or use strip.width if we have strip transects
      width <- strip.width
      # note we have to reset off.set
      if(response == "presence"){
        off.set <- "none"
      }else{
        off.set <- "area"
      }
    }

    # check that none of the Effort values are zero
    if(any(dat[,seglength.name]==0)){
      stop(paste0("Effort values for segments: ",
                  paste(which(dat[,seglength.name]==0),collapse=", "),
                  " are 0."))
    }

    # calculate the offset
    #   area we just calculate the area
    #   effective area multiply by p
    #   when density is response, offset should be 1 (and is ignored anyway)

    dat$off.set <- switch(off.set,
                          eff.area=2*dat[,seglength.name]*width*fitted.p,
                          area=2*dat[,seglength.name]*width,
                          none=1)

    # calculate the density (count/area)
    if(response %in% c("D","density","Dhat","density.est")){
      dat[,response] <- dat[,response]/(2*dat[,seglength.name]*
                                        width*convert.units)
    }

    # set the segment area in the data
    dat$segment.area <- 2*dat[,seglength.name]*width*convert.units
  }

  # multiply up by conversion factor
  dat$off.set <- dat$off.set*convert.units

  # Set offset as log of area or effective area
  dat$off.set <- log(dat$off.set)

  return(dat)
}
