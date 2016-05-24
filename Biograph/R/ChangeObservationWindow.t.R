ChangeObservationWindow.t <-
function (Bdata,starttime,endtime,covs.dates)
{ #  Check whether Parameters was called
	
   if (missing(covs.dates)) covs.dates <- NULL
  z<- check.par(Bdata)
     namstates <- attr(Bdata,"param")$namstates
  # Check whether starting and ending times are in observation window
  Bdata2 <- Bdata
  locpat <- locpath(Bdata)
  za <- rep(starttime,nrow(Bdata2))
  zb <- rep(endtime,nrow(Bdata2))
  Bdata2 <- subset(Bdata,Bdata$start<zb & Bdata$end>za)
  za <- rep(starttime,nrow(Bdata2))
  zb <- rep(endtime,nrow(Bdata2))
  Bdata2$start <- ifelse (Bdata2$start < za, za,Bdata2$start)
  Bdata2$end <- ifelse (Bdata2$end > zb,zb,Bdata2$end)
  attr(Bdata2,"format.date") <- attr(Bdata,"format.date")
  attr(Bdata2,"format.born") <- attr(Bdata,"format.born")
  attr(Bdata2,"param") <- attr(Bdata,"param")
 # Bdata2$marriage <- ifelse (Bdata2$marriage >= za & Bdata2$marriage <= zb,Bdata2$marriage,0)
  if (is.null(covs.dates)) 
  {for (jj in 1:ncol(Bdata2))
  	{ if (colnames(Bdata2)[jj] %in% covs.dates)
    { # If event after END of observation, date = 0
      Bdata2[,jj]	<- ifelse (Bdata2[,jj] <= zb,Bdata2[,jj],0) }
    } }
   
   for (i in 1:nrow(Bdata2))
  { # state occupied at starttime
  	ns <- nchar(Bdata2$path[i])
    if (ns >1)
    {zx <- c(Bdata2$start[i],Bdata[i,(locpath(Bdata)+1):(locpath(Bdata)+ns-1)],Bdata2$end[i])
     zx <- Bdata[i,(locpath(Bdata)+1):(locpath(Bdata)+ns-1)]
     zx <- unlist (zx)
     zy <- ifelse (zx >= rep(starttime,length(zx)) & zx <= rep(endtime,length(zx)),zx,0)
     # new transition dates
     Bdata2[i,(locpat+1):ncol(Bdata2)] <- NA
     if (sum(zy[zy>0]) > 0) # At least one transition in observation window
      { Bdata2[i,(locpat+1):(locpat+length(zy[zy>0]))] <- zy[zy>0]
      }
 ## Determine state occupied at onset of new observation window
      if (max(zx,na.rm=TRUE) <starttime)
        {ns <- 1
         Bdata2$path[i] <- substr(Bdata$path[i],ns,ns) } else
       { if (length(na.omit(starttime-zx)) > 0)  # vector of NAs
         {ii <- ifelse (starttime > zx[1], min(which(starttime-zx<=0),na.rm=TRUE),NA)} else # starttime after birth
         {ii <- NA}
      if (is.na(ii)) ii <- 1 
      state1 <- substr(Bdata2$path[i],ii,ii)
      pathn <- substr(Bdata2$path[i],ii,ii+length(zy[zy>0]))
     Bdata2$path[i] <- pathn
   }}    
  }
  param <- Parameters(Bdata2)
  attr(Bdata2, "param") <- param
 # attr(Bdata2,"statespace") <- namstates
  attr(Bdata2,"format.date") <- attr(Bdata,"format.date")
  attr(Bdata2,"format.born") <- attr(Bdata,"format.born")
  print("A Biograph object with new observation window is returned.",quote = FALSE)
  return (Bdata =Bdata2)
 }
