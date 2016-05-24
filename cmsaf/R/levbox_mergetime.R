levbox_mergetime <-
function(var,level=1,path,pattern,outfile,lon1=-180,lon2=180,lat1=-90,lat2=90){

  start.time <- Sys.time()

# define standard names of variables and dimensions

   t_name <- "time"
   t_standard_name = "time"
   t_units = "undefined"
   t_calendar = "undefined"

   nb2_units = "1"

   lat_name = "latitude"
   lat_standard_name = "latitude"
   lat_long_name = "latitude"
   lat_units = "degrees_north"
   lat_axis = "Y"

   lon_name = "longitude"
   lon_standard_name = "longitude"
   lon_long_name = "longitude"
   lon_units = "degrees_east"
   lon_axis = "X"

   v_standard_name = "undefined"
   v_long_name = "undefined"
   v_units = "undefined"
   v__FillValue = "undefined"
   v_missing_value = "undefined"

   info = "Created with the CM SAF R toolbox." 
   var_prec="double"

   att_list <- c("standard_name","long_name","units","_FillValue","missing_value","calendar")
   v_att_list <- c("v_standard_name","v_long_name","v_units","v__FillValue","v_missing_value","v_calendar")
  
# get file information

  cat("get file information", "\n")

  filelist <- list.files(path=path, pattern=pattern)
  fdim <- length(filelist)

  file=filelist[1]
  file <- paste(path,"/",file,sep="")
  id <- nc_open(file)

  # get information about dimensions

  dimnames <- names(id$dim)

    # check standard_names of dimensions
      for (i in 1:length(dimnames)){
	sn <- ncatt_get(id,dimnames[i],"standard_name")
	if (length(sn)>0){
	  sn <- sn$value
	  if (sn=="longitude")(lon_name <- dimnames[i])
	  if (sn=="latitude")(lat_name <- dimnames[i])
	  if (sn=="time")(t_name <- dimnames[i])
	}
      }

  for (i in 1:length(dimnames)){
    if (t_name %in% dimnames){
      attnames <- names(id$dim[[i]])
      if ("units" %in% attnames){
	t_units <- ncatt_get(id,t_name,"units")$value}
      if ("calendar" %in% attnames){
	t_calendar <- ncatt_get(id,t_name,"calendar")$value}
    }
  }

  # get information about variables
	
  varnames <- names(id$var)

   if (var %in% varnames){
    for (i in 1:6){
      att_dum <- ncatt_get(id,var,att_list[i])
      if (att_dum$hasatt){
	assign(v_att_list[i],att_dum$value)}
    }

      # get data of first file

	lon <- ncvar_get(id,lon_name)
	lat <- ncvar_get(id,lat_name)
	time1 <- ncvar_get(id,t_name)
	time_len <- length(time1)
	if ("time_bnds" %in% varnames){
	  tbnds1 <- ncvar_get(id,"time_bnds")
	}
	
	lon_limit <- which(lon>=lon1&lon<=lon2)  
	lat_limit <- which(lat>=lat1&lat<=lat2) 

	lon <- lon[lon_limit]
	lat <- lat[lat_limit]
  
	# check for empty lon_limit or lat_limit
	
	if (length(lon_limit)==0|length(lat_limit)==0){
	  nc_close(id)
	  stop("Selected region is outside target area!")
	}

	startx <- min(lon_limit)
	starty <- min(lat_limit)
	countx <- length(lon_limit)
	county <- length(lat_limit)
	countt <- length(time1)

	target <- ncvar_get(id,var,start=c(startx,starty,level,1),count=c(countx,county,1,countt))
   }else{
      nc_close(id)
      stop(cat(paste("Variable ",var," not found! File contains: ",varnames,sep="")),"\n")}

  if (v__FillValue == "undefined"){ 
    v__FillValue = v_missing_value}
  if (v_missing_value == "undefined"){ 
    v_missing_value = v__FillValue}

  nc_close(id)

   if ("time_bnds" %in% varnames){
    time_bnds <- array(NA,dim=c(2,length(time1)))
    time_bnds[,1:length(time1)] <- tbnds1
   }   

  #dum <- max(target,na.rm=T)
  #if (is.integer(dum)){var_prec="short"}

# get time reference

  dt_ref <- get_time(t_units,0)
  unit_ref <- unlist(strsplit(t_units,split=" "))[1]

  # check reference time unit
  if (unit_ref=="minutes"|unit_ref=="Minutes"|unit_ref=="Mins"|unit_ref=="Min"|unit_ref=="min")(unit_ref <- "mins")
  if (unit_ref=="seconds"|unit_ref=="Seconds"|unit_ref=="Secs"|unit_ref=="Sec"|unit_ref=="sec")(unit_ref <- "secs")
  if (unit_ref=="Hours"|unit_ref=="Hour"|unit_ref=="hour")(unit_ref <- "hours")
  if (unit_ref=="Days"|unit_ref=="Day"|unit_ref=="day")(unit_ref <- "days")
  if (unit_ref=="Weeks"|unit_ref=="Week"|unit_ref=="week")(unit_ref <- "weeks")
  if (unit_ref=="Months"|unit_ref=="Month"|unit_ref=="month")(unit_ref <- "months")
  if (unit_ref!="mins"&unit_ref!="secs"&unit_ref!="hours"&unit_ref!="days"&unit_ref!="weeks"&unit_ref!="months")(unit_ref <- "auto")

# create netcdf

  cat("create netcdf", "\n")

    target[is.na(target)] <- v_missing_value
    nb2 <- c(0,1)

    x <- ncdim_def(name="lon",units=lon_units,vals=lon)
    y <- ncdim_def(name="lat",units=lat_units,vals=lat)
    t <- ncdim_def(name="time",units=t_units,vals=time1,unlim=TRUE)
    if ("time_bnds" %in% varnames){
      tb <- ncdim_def(name="nb2",units="1",vals=nb2)
    }

    var1 <- ncvar_def(name=var,units=v_units,dim=list(x,y,t),prec=var_prec)

    if ("time_bnds" %in% varnames){
      var2 <- ncvar_def(name="time_bnds",units="1",dim=list(tb,t),prec=var_prec)
      vars <- list(var1,var2)
      ncnew <- nc_create(outfile,vars)

      ncvar_put(ncnew,var1,target)
      ncvar_put(ncnew,var2,time_bnds)

      ncatt_put(ncnew,var,"standard_name",v_standard_name,prec="text")
      ncatt_put(ncnew,var,"long_name",v_long_name,prec="text")
      ncatt_put(ncnew,var,"_FillValue",v__FillValue,prec=var_prec)
      ncatt_put(ncnew,var,"missing_value",v_missing_value,prec=var_prec)

      ncatt_put(ncnew,"time","standard_name",t_standard_name,prec="text")
      ncatt_put(ncnew,"time","calendar",t_calendar,prec="text")
      ncatt_put(ncnew,"time","bounds","time_bnds",prec="text")

      ncatt_put(ncnew,"lon","standard_name",lon_standard_name,prec="text")
      ncatt_put(ncnew,"lon","long_name",lon_long_name,prec="text")
      ncatt_put(ncnew,"lon","axis",lon_axis,prec="text")

      ncatt_put(ncnew,"lat","standard_name",lat_standard_name,prec="text")
      ncatt_put(ncnew,"lat","long_name",lat_long_name,prec="text")
      ncatt_put(ncnew,"lat","axis",lat_axis,prec="text")

      ncatt_put(ncnew,0,"Info",info,prec="text")

    } else {
      vars <- list(var1)
      ncnew <- nc_create(outfile,vars)

      ncvar_put(ncnew,var1,target)

      ncatt_put(ncnew,var,"standard_name",v_standard_name,prec="text")
      ncatt_put(ncnew,var,"long_name",v_long_name,prec="text")
      ncatt_put(ncnew,var,"_FillValue",v__FillValue,prec=var_prec)
      ncatt_put(ncnew,var,"missing_value",v_missing_value,prec=var_prec)

      ncatt_put(ncnew,"time","standard_name",t_standard_name,prec="text")
      ncatt_put(ncnew,"time","calendar",t_calendar,prec="text")

      ncatt_put(ncnew,"lon","standard_name",lon_standard_name,prec="text")
      ncatt_put(ncnew,"lon","long_name",lon_long_name,prec="text")
      ncatt_put(ncnew,"lon","axis",lon_axis,prec="text")

      ncatt_put(ncnew,"lat","standard_name",lat_standard_name,prec="text")
      ncatt_put(ncnew,"lat","long_name",lat_long_name,prec="text")
      ncatt_put(ncnew,"lat","axis",lat_axis,prec="text")

      ncatt_put(ncnew,0,"Info",info,prec="text")

    }
      
    # get data for specific level and cut desired region

    time_len <- length(time1)

    for (i in 2:fdim){
      cat("\r","loading file ",i," of ",fdim,sep="")
      file=filelist[i]
      file <- paste(path,"/",file,sep="")
      id <- nc_open(file)

      dum_dat <- ncvar_get(id,var,start=c(startx,starty,level,1),count=c(countx,county,1,-1))
      dum_time <- ncvar_get(id,t_name)
      time_len <- time_len+length(dum_time)

      dum_t_units <- ncatt_get(id,t_name,"units")$value
      dt_dum <- get_time(dum_t_units,dum_time)
      if (as.character(dt_ref)=="-4712-01-01 12:00:00"){
        dum_time <- (as.numeric(dt_dum)/86400)+2440587.5
      } else {
        if (unit_ref=="months"){
          dum_time <- round((difftime(dt_dum,dt_ref,units=c("days")))/30.4375)
          dum_time <- as.numeric(dum_time)
        } else {
          dum_time <- difftime(dt_dum,dt_ref,units=c(unit_ref))
        }
      }

      if ("time_bnds" %in% varnames){
	      dum_tb <- ncvar_get(id,"time_bnds",collapse_degen=FALSE)
      }

      nc_close(id)

    dum_dat[is.na(dum_dat)] <- v_missing_value
    startt2 <- time_len-length(dum_time)+1
    countt2 <- length(dum_time)

    if ("time_bnds" %in% varnames){

      ncvar_put(ncnew,var1,dum_dat,start=c(1,1,startt2),count=c(-1,-1,countt2))
      ncvar_put(ncnew,var2,dum_tb,start=c(1,startt2),count=c(-1,countt2))
      ncvar_put(ncnew,t,dum_time, start=startt2, count=countt2 )
      nc_sync(ncnew)

    } else {

      ncvar_put(ncnew,var1,dum_dat,start=c(1,1,startt2),count=c(-1,-1,countt2))
      ncvar_put(ncnew,t,dum_time, start=startt2, count=countt2 )
      nc_sync(ncnew)
    }
  }

    nc_close(ncnew)

end.time <- Sys.time()
cat("\n","processing time: ",round(as.numeric(end.time-start.time,units="secs"),digits=2)," s",sep="", "\n")
}
