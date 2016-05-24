ncinfo <-
function(infile, info="s"){

# define standard names of variables and dimensions

   t_name <- "time"
   t_standard_name = "time"
   t_units = "undefined"

# get file information

  id <- nc_open(infile)

  dimnames <- names(id$dim)
  varnames <- names(id$var)

if (info=="l"){ cat(str(id),"\n")}
if (info=="m"){ print(id)}
if (info=="s"){
     cat("The file:",id$filename,"contains:", "\n")
    if (length(varnames)==1){
      cat("\n","Variable:",sep="", "\n")
      cat(varnames[1], "\n")
    } else {
      cat("\n","Variables:",sep="", "\n")
      for (i in 1:length(varnames)){
	cat(varnames[i], "\n")
      }
    }
    
    # check standard_names of dimensions
    for (i in 1:length(dimnames)){
	    sn <- ncatt_get(id,dimnames[i],"standard_name")
	    if (length(sn)>0){
	      sn <- sn$value
	    if (sn=="time")(t_name <- dimnames[i])
	    }
    }

  cat("\n","With following dimensions:",sep="", "\n")
      for (i in 1:length(dimnames)){
	if (dimnames[i]==t_name){
	    for (j in 1:length(dimnames)){
	      if (t_name %in% dimnames){
	      attnames <- names(id$dim[[i]])
	      if ("units" %in% attnames){
		      t_units <- ncatt_get(id,t_name,"units")$value}
	      }
	    }
	    time1 <- ncvar_get(id,"time")
	    date.time <- as.Date(get_time(t_units,time1))
	    cat("time with length ",length(time1)," (range ",as.character(min(date.time))," to ",as.character(max(date.time)),")",sep="", "\n")
	} else {
	    cat(dimnames[i]," with length ",id$dim[[i]]$len," (range ",min(id$dim[[i]]$vals)," to ",max(id$dim[[i]]$vals),")",sep="", "\n")
	  }
      }
  }

nc_close(id)
}
