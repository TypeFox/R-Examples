selyear <-
function(var,year=c(2000),infile,outfile){

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

      # get details of file

	lon <- ncvar_get(id,lon_name)
	lat <- ncvar_get(id,lat_name)
	time1 <- ncvar_get(id,t_name)
	if ("time_bnds" %in% varnames){
	  tbnds1 <- ncvar_get(id,"time_bnds",collapse_degen=FALSE)
	}
   }else{
      nc_close(id)
      stop(cat(paste("Variable ",var," not found! File contains: ",varnames,sep="")),"\n")}

  if (v__FillValue == "undefined"){ 
    v__FillValue = v_missing_value}
  if (v_missing_value == "undefined"){ 
    v_missing_value = v__FillValue}

  nc_close(id)   

# extract time information

  date.time <- as.Date(get_time(t_units,time1))
  a <- as.character(date.time)
  b <- strsplit(a,"-")
  d <- unlist(b)
  dummy <- length(d)
  dum <- seq(1,dummy,3)
  years <- d[dum]
  years <- as.numeric(years)

 if (length(which(year==years))>=1){

  target <- array(NA,dim=c(length(lon),length(lat),1))
  time_bnds <- array(NA, dim=c(2,1))

 # create netcdf

  cat("create netcdf", "\n")

    target[is.na(target)] <- v_missing_value

    nb2 <- c(0,1)

    x <- ncdim_def(name="lon",units=lon_units,vals=lon)
    y <- ncdim_def(name="lat",units=lat_units,vals=lat)
    t <- ncdim_def(name="time",units=t_units,vals=0,unlim=TRUE)
    if ("time_bnds" %in% varnames){
      tb <- ncdim_def(name="nb2",units="1",vals=nb2)
    }

    var1 <- ncvar_def(name=var,units=v_units,dim=list(x,y,t),prec=var_prec)

    if ("time_bnds" %in% varnames){
      var2 <- ncvar_def(name="time_bnds",units="1",dim=list(tb,t),prec="double")
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

    # extract desired year from infile

      id <- nc_open(infile)
      count <- 1
  
      for (i in 1:length(time1)){
       for (j in 1:length(year)){
	if (years[i]==year[j]){
	  dum_dat <- ncvar_get(id,var,start=c(1,1,i),count=c(-1,-1,1))
	  dum_dat[is.na(dum_dat)] <- v_missing_value
	  ncvar_put(ncnew,var1,dum_dat,start=c(1,1,count),count=c(-1,-1,1))
	  ncvar_put(ncnew,t,time1[i], start=count, count=1)
	  if ("time_bnds" %in% varnames){
	    ncvar_put(ncnew,var2,tbnds1[,i],start=c(1,count),count=c(-1,1))
	  }
	  count <- count+1
	}
       }
      }
      nc_close(id)
      nc_close(ncnew)

 } else { cat("No match. Years are:",years, "\n")}

end.time <- Sys.time()
cat("processing time: ",round(as.numeric(end.time-start.time,units="secs"),digits=2)," s",sep="", "\n")
  } # endif filecheck
}
