year.anomaly <-
function(var,infile,outfile){

  start.time <- Sys.time()

# check filename

  filecheck <- checkfile(infile,outfile)

  if (filecheck[[1]]){
    infile <- filecheck[[2]]
    outfile <- filecheck[[3]]  

# user define section

  limit <- 2601*2601*31	  # This value can be ajusted to avoid RAM overflow  

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
	time_len <- length(time1)
	if ("time_bnds" %in% varnames){
	  tbnds1 <- ncvar_get(id,"time_bnds")
	}
   }else{
      nc_close(id)
      stop(cat(paste("Variable ",var," not found! File contains: ",varnames,sep="")),"\n")}

  if (v__FillValue == "undefined"){ 
    v__FillValue = v_missing_value}
  if (v_missing_value == "undefined"){ 
    v_missing_value = v__FillValue} 

  # check data dimensions and calculate climatology 
  
  if ((length(lon)*length(lat)*time_len)<limit){

    # calculate temporal mean in a sequence depending on limit

    dum_dat <- ncvar_get(id,var,collapse_degen=FALSE)
    cat("get climatology", "\n")
    clim <- rowMeans(dum_dat,dims=2,na.rm=T)
  } else {

   dum1 <- round((limit/length(lon))/length(lat))
   dum2 <- seq(1,time_len,dum1)
   dum3 <- array(dum1,dim=c(length(dum2)))
   cor <- dum1*length(dum2)-time_len
   dum3[length(dum2)] <- dum3[length(dum2)]-cor
  
  sum_data <- array(NA,dim=c(length(lon),length(lat),length(dum2)))
   num <- 0
   cat("get climatology sequential",sep="","\n")
   for (i in 1:length(dum2)){
    dum_dat <- ncvar_get(id,var,start=c(1,1,dum2[i]),count=c(-1,-1,dum3[i]),collapse_degen=FALSE)
    sum_data[,,i] <- rowSums(dum_dat,dims=2,na.rm=T)
    nan <- rowSums(!is.na(dum_dat),dims=2)
    num <- num+nan
   }
    clim <- rowSums(sum_data,dims=2,na.rm=T)/num
  }
   nc_close(id) 

# extract time information

  date.time <- as.Date(get_time(t_units,time1))
  a <- as.character(date.time)
  b <- strsplit(a,"-")
  d <- unlist(b)
  dum <- seq(1,length(d),3)
  year <- as.integer(d[dum])
  yl <- as.integer(levels(factor(year)))
  dummy_vec <- c(1:length(year))

  target <- array(NA,dim=c(length(lon),length(lat),length(yl)))
  time_bnds <- array(NA, dim=c(2,length(yl)))
  count <- 1
  for (i in 1:length(yl)){
    year_dummy <- which(year==yl[i])
    time_bnds[1,count] <- time1[min(year_dummy)]
    time_bnds[2,count] <- time1[max(year_dummy)]
    count <- count+1
   }

# create netcdf

  cat("create netcdf", "\n")

    target[is.na(target)] <- v_missing_value

    nb2 <- c(0,1)
    times <- time_bnds[1,]

    x <- ncdim_def(name="lon",units=lon_units,vals=lon)
    y <- ncdim_def(name="lat",units=lat_units,vals=lat)
    t <- ncdim_def(name="time",units=t_units,vals=times,unlim=TRUE)
    tb <- ncdim_def(name="nb2",units="1",vals=nb2)

    var1 <- ncvar_def(name=var,units=v_units,dim=list(x,y,t),prec=var_prec)
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

  # calculate monthly means

  count <- 1
  for (i in 1:length(yl)){
    year_dummy <- which(year==yl[i])
    startt <- min(dummy_vec[year_dummy])
    countt <- length(year_dummy)
    id <- nc_open(infile)
    dum_dat <- ncvar_get(id,var,start=c(1,1,startt),count=c(-1,-1,countt),collapse_degen=FALSE)
    cat("\r","apply annual mean ",count,sep="")
    mean_data <- rowMeans(dum_dat,dims=2,na.rm=T)
    mean_data <- mean_data-clim	 # get anomaly
    mean_data[is.na(mean_data)] <- v_missing_value
    ncvar_put(ncnew,var1,mean_data,start=c(1,1,count),count=c(-1,-1,1))
    count <- count+1
   }
 nc_close(id)

 nc_close(ncnew)

end.time <- Sys.time()
cat("\n","processing time: ",round(as.numeric(end.time-start.time,units="secs"),digits=2)," s",sep="", "\n")
  } # endif filecheck
}
