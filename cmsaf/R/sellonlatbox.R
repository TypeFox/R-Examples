sellonlatbox <-
function(var,infile,outfile,lon1=-180,lon2=180,lat1=-90,lat2=90){

  start.time <- Sys.time()

# check filename

  filecheck <- checkfile(infile,outfile)

  if (filecheck[[1]]){
    infile <- filecheck[[2]]
    outfile <- filecheck[[3]]  

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

  id <- nc_open(infile)

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

  # get data of first file and cut desired region 

	lon <- ncvar_get(id,lon_name)
	lat <- ncvar_get(id,lat_name)
	time1 <- ncvar_get(id,t_name)
	time_len <- length(time1)
	if ("time_bnds" %in% varnames){
	  tbnds1 <- ncvar_get(id,"time_bnds",collapse_degen=FALSE)
	}
	
	lon_limit <- which(lon>=lon1&lon<=lon2)  
	lat_limit <- which(lat>=lat1&lat<=lat2) 

	lon <- lon[lon_limit]
	lat <- lat[lat_limit]

	startx <- min(lon_limit)
	starty <- min(lat_limit)
	countx <- length(lon_limit)
	county <- length(lat_limit)
	countt <- length(time1)

	data1 <- ncvar_get(id,var,start=c(startx,starty,1),count=c(countx,county,countt))
   }else{
      nc_close(id)
      stop(cat(paste("Variable ",var," not found! File contains: ",varnames,sep="")),"\n")}

  if (v__FillValue == "undefined"){ 
    v__FillValue = v_missing_value}
  if (v_missing_value == "undefined"){ 
    v_missing_value = v__FillValue}

  nc_close(id)

# create netcdf

  cat("create netcdf", "\n")

    if (length(time1)==1){
      dummy <- array(NA,dim=c(dim(data1)[1],dim(data1)[2],1))
      dummy[,,1] <- data1
      data1 <- dummy
    }

    data1[is.na(data1)] <- v_missing_value
    nb2 <- c(0,1)

    x <- ncdim_def(name="lon",units=lon_units,vals=lon)
    y <- ncdim_def(name="lat",units=lat_units,vals=lat)
    t <- ncdim_def(name="time",units=t_units,vals=time1,unlim=TRUE)
    if ("time_bnds" %in% varnames){
      tb <- ncdim_def(name="nb2",units=nb2_units,vals=nb2)
    }

    var1 <- ncvar_def(name=var,units=v_units,dim=list(x,y,t),prec=var_prec)

    if ("time_bnds" %in% varnames){
      var2 <- ncvar_def(name="time_bnds",units="1",dim=list(tb,t),prec="double")
      vars <- list(var1,var2)
      ncnew <- nc_create(outfile,vars)

      ncvar_put(ncnew,var1,data1)
      ncvar_put(ncnew,var2,tbnds1)

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

      ncvar_put(ncnew,var1,data1)

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

    nc_close(ncnew)

end.time <- Sys.time()
cat("processing time: ",round(as.numeric(end.time-start.time,units="secs"),digits=2)," s", sep="","\n")
  } # endif filecheck
}
