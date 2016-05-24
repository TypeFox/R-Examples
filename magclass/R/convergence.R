convergence <- function(origin, aim, start_year=NULL, end_year=NULL,
                        direction=NULL, type="smooth", par=1.5) {

  ### Basic checks ###
  if(!is.magpie(origin)) stop("origin is no magpie object")

  if(is.null(dim(aim))) aim<-as.magpie(array(aim,dim(origin),dimnames(origin)))
  if(!is.magpie(aim)) stop("aim is no magpie object")

  if (all(dimnames(aim)[[1]]!=dimnames(origin)[[1]]))
    stop("regions have to be the same")

  if (ndata(origin)!=1 & !identical(getNames(origin), getNames(aim)))
    stop("If there ist more than one name-column, dimnames have to be the same")

  if(nyears(aim)==1) {
    tmp <- setYears(aim,NULL)
    aim <- origin
    aim[,,] <- tmp
    rm(tmp)
  }
  
  if(nyears(origin)==1) {
    tmp <- setYears(origin,NULL)
    origin <- aim
    origin[,,] <- tmp
    rm(tmp)
  }

  if (any(getYears(aim) != getYears(origin)))
    stop("Objects need the same timesteps, or aim has to have only one timestep")

  if (is.null(start_year)) start_year <- getYears(aim)[1]
  if (is.null(end_year)) end_year <- getYears(aim)[nyears(aim)]

  if(isYear(start_year,with_y=TRUE)) start_year <- substr(start_year,2,5)
  if(isYear(end_year,with_y=TRUE))   end_year   <- substr(end_year,2,5)
  start_year <- as.numeric(start_year)
  end_year   <- as.numeric(end_year)
  if(!isYear(start_year,with_y=FALSE)) stop("wrong year format for convergence aim")
  if(!isYear(end_year,with_y=FALSE)) stop("wrong year format for convergence aim")


  # In the case of direction up or down data should only be manipulated in one
  # direction. Therefor, the aim object is manipulated accordingly
  if(!is.null(direction)) {
    aim <- as.array(aim)
    if(direction=="up")   aim[which(aim<origin)] <- as.array(origin)[which(aim<origin)]
    else if(direction=="down") aim[which(aim>origin)] <- as.array(origin)[which(aim>origin)]
    else stop("Illegal direction setting, only up and down are allowed arguments!")
    aim <- as.magpie(aim)
  }                                                                                            
                        
  years<-new.magpie("GLO",getYears(origin),NULL,getYears(origin,as.integer=TRUE))                                                 
  pos <- (years - start_year)/(end_year - start_year)
  pos[pos<0] <- 0
  pos[pos>1] <- 1
  if (type == "linear") { mix <- pos 
  } else if (type == "s") { mix <- pos^4/(0.07+pos^4)*1.07
  } else if (type == "smooth") {mix <- pos^3/(0.1 + pos^3)
#  } else if (type == "decay") {mix <- pos/(0.5 + pos)*1.5
  } else if (type == "decay") {mix <- pos/(par + pos)*(par+1)
  } else {stop("type does not exist")}
  converged <- aim * mix + origin * (1 - mix)

  return(converged)
}