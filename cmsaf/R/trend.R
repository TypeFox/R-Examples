trend <-
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

   s_name = "sig"
   s_standard_name = "significance"
   s_long_name = "significance based on 95% confidence interval"
   s_units = "1"
   s__FillValue = "undefined"
   s_missing_value = "undefined"
   s_info = "1 = positive significant, 0 = not significant, -1 = negative significant"

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

  target <- array(NA,dim=c(length(lon),length(lat),1))
  target_p <- array(NA,dim=c(length(lon),length(lat),1))
  y <- c(1:time_len)
      cat("fit with linear model...",sep="", "\n")

      for (i in 1:length(lon)){
	prog <- round((100/length(lon))*i)
	cat("\r","progress: ",prog,"%",sep="")
	for (j in 1:length(lat)){
	  dum_dat <- ncvar_get(id,var,start=c(i,j,1),count=c(1,1,-1))
	  if (time_len-(sum(is.na(dum_dat)))>=2){
	    fit <- lm(dum_dat~y,na.action=na.omit)
	    val <- fit$coef[2]*time_len
	    ci <- confint(fit,'y', level=0.95)
	    sig <- 0
	    if (ci[1]*ci[2]<0)(sig <- 0)
	    if (ci[1]<0&ci[2]<0)(sig <- -1)
	    if (ci[1]>0&ci[2]>0)(sig <- 1)
	  } else {
	      val <- NA
	      sig <- NA
	    }
	  target[i,j,1] <- val
	  target_p[i,j,1] <- sig
	}
      }

   nc_close(id)

 # create netcdf

  cat("create netcdf", "\n")

    target[is.na(target)] <- v_missing_value
    target_p[is.na(target)] <- v_missing_value

    time_bnds <- array(NA, dim=c(2,1))
    time_bnds[1,1] <- min(time1)
    time_bnds[2,1] <- max(time1)

    nb2 <- c(0,1)
    times <- time_bnds[1,]

    x <- ncdim_def(name="lon",units=lon_units,vals=lon)
    y <- ncdim_def(name="lat",units=lat_units,vals=lat)
    t <- ncdim_def(name="time",units=t_units,vals=times,unlim=TRUE)
    tb <- ncdim_def(name="nb2",units="1",vals=nb2)

    var1 <- ncvar_def(name=var,units=v_units,dim=list(x,y,t),prec=var_prec)
    var2 <- ncvar_def(name="time_bnds",units="1",dim=list(tb,t),prec="double")
    var3 <- ncvar_def(name=s_name,units=s_units,dim=list(x,y,t),prec="double")

      vars <- list(var1,var2,var3)
      ncnew <- nc_create(outfile,vars)

      ncvar_put(ncnew,var1,target)
      ncvar_put(ncnew,var2,time_bnds)
      ncvar_put(ncnew,var3,target_p)

      ncatt_put(ncnew,var,"standard_name",v_standard_name,prec="text")
      ncatt_put(ncnew,var,"long_name",v_long_name,prec="text")
      ncatt_put(ncnew,var,"_FillValue",v__FillValue,prec=var_prec)
      ncatt_put(ncnew,var,"missing_value",v_missing_value,prec=var_prec)

      ncatt_put(ncnew,s_name,"standard_name",s_standard_name,prec="text")
      ncatt_put(ncnew,s_name,"long_name",s_long_name,prec="text")
      ncatt_put(ncnew,s_name,"description",s_info,prec="text")
      ncatt_put(ncnew,s_name,"_FillValue",s__FillValue,prec="double")
      ncatt_put(ncnew,s_name,"missing_value",s_missing_value,prec="double")

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
     
   nc_close(ncnew)

end.time <- Sys.time()
cat("\n","processing time: ",round(as.numeric(end.time-start.time,units="secs"),digits=2)," s",sep="", "\n")
  } # endif filecheck
}
