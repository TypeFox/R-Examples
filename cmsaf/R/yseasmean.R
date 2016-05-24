yseasmean <-
function(var,infile,outfile){

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
	time_len <- length(time1)
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
  dum <- seq(2,dummy,3)
  mon <- as.integer(d[dum])
  dum <- seq(1,dummy,3)
  year <- as.integer(d[dum])
  yl <- as.integer(levels(factor(year)))
  ml <- as.integer(levels(factor(mon)))
  nmonmeans <- length(yl)*length(ml)
  mul <- year*mon
  dummy_vec <- c(1:length(mon))

  seas <- array(NA,dim=c(4,3,length(yl)))

  for (i in 1:length(yl)){
    win <- which(mon==1&year==yl[i]|mon==2&year==yl[i]|mon==12&year==(yl[i]-1))
    if (length(win)==3){seas[1,,i]<-win}
      spr <- which(mon==3&year==yl[i]|mon==4&year==yl[i]|mon==5&year==yl[i])
    if (length(spr)==3){seas[2,,i]<-spr}
      sum <- which(mon==6&year==yl[i]|mon==7&year==yl[i]|mon==8&year==yl[i])
    if (length(sum)==3){seas[3,,i]<-sum}
      aut <- which(mon==9&year==yl[i]|mon==10&year==yl[i]|mon==11&year==yl[i])
    if (length(aut)==3){seas[4,,i]<-aut}
  }

  target <- array(NA,dim=c(length(lon),length(lat),1))
  tbnds <- array(NA, dim=c(2,1))

# create netcdf

  cat("create netcdf", "\n")

    target[is.na(target)] <- v_missing_value
    nb2 <- c(0,1)

    x <- ncdim_def(name="lon",units=lon_units,vals=lon)
    y <- ncdim_def(name="lat",units=lat_units,vals=lat)
    t <- ncdim_def(name="time",units=t_units,vals=0,unlim=TRUE)
    tb <- ncdim_def(name="nb2",units="1",vals=nb2)

    var1 <- ncvar_def(name=var,units=v_units,dim=list(x,y,t),prec=var_prec)
    var2 <- ncvar_def(name="time_bnds",units="1",dim=list(tb,t),prec=var_prec)
    vars <- list(var1,var2)
    ncnew <- nc_create(outfile,vars)

    ncvar_put(ncnew,var1,target)
    ncvar_put(ncnew,var2,tbnds)

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

  # get data and calculate multi-year seasonal means

  id <- nc_open(infile)
  limit <- (length(yl)*3)*0.25  # 75% of monthly data have to be available
  count <- 1
    for (j in 1:4){
      mon_dummy <- seas[j,,]
	if (sum(is.na(mon_dummy))<=limit){
	  dum_dat <- array(NA,dim=c(length(lon),length(lat),length(mon_dummy)))
	   for (i in 1:length(mon_dummy)){
	    if (!is.na(mon_dummy[i])){
	      dum_dat[,,i] <- ncvar_get(id,var,start=c(1,1,mon_dummy[i]),count=c(-1,-1,1),collapse_degen=FALSE)
	    }
	   }
	  cat("\r","apply multi-year seasonal mean ",count," of 4",sep="")
	  mean_data <- rowMeans(dum_dat,dims=2,na.rm=T)
	  mean_data[is.na(mean_data)] <- v_missing_value
	  tdum <- min(time1[mon_dummy],na.rm=T)
	  tbnds[1,1] <- min(time1[mon_dummy],na.rm=T)
	  tbnds[2,1] <- max(time1[mon_dummy],na.rm=T)
	  ncvar_put(ncnew,var1,mean_data,start=c(1,1,count),count=c(-1,-1,1))
	  ncvar_put(ncnew,t,tdum,start=count,count=1)
	  ncvar_put(ncnew,var2,tbnds,start=c(1,count),count=c(-1,1))
	  count <- count+1
	}
     }

 nc_close(id)

 nc_close(ncnew)

end.time <- Sys.time()
cat("\n","processing time: ",round(as.numeric(end.time-start.time,units="secs"),digits=2)," s",sep="", "\n")
  } # endif filecheck
}
