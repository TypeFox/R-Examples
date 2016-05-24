date_b <-
function (Bdata,selectday,format.out,covs)
#  e.g format.in = "CMC" and selectday = 1
{  # Determine position of path variable in Biograph object
#   if (is.null(attr(Bdata,"timeunit"))) 
#      {stop ("date.years: attribute timeunit missing from data. Please add attribute.") } else
#      {timeunit <- attr(Bdata,"timeunit")
#       if (timeunit=="year") 
#         {print ("Biograph message: Time unit is already year. Conversion skipped and no new object created.",quote=FALSE); return ()}}  
   if (missing(selectday)) selectday <- 1 
   if (missing(covs)) covs=NULL
    format.in <- attr (Bdata,"format.date")
   format.born <- attr (Bdata,"format.born")
 if (is.null(format.born)) stop ("Error in date_b: add format of birth date as attribute to the data (Parameter file)")
   born <- Bdata$born
   Bdata$start <- date_convert(Bdata$start,format.in,selectday,format.out,born=born,format.born=format.born)
   Bdata$end   <- date_convert(Bdata$end,format.in,selectday,format.out,born=born,format.born=format.born)
 # Determine maximum number of transitions on Biograph object
   locpat <- locpath(Bdata)
   maxtrans <- ncol(Bdata)-locpat
   for (i in (locpat+1):(locpat+maxtrans))
    { Bdata[,i] <- date_convert(Bdata[,i],format.in,selectday,format.out,born=born,format.born=format.born) # MAKE SURE TO SET Bdata = GLHS
    }
    if (!is.null(covs)) 
    { zz <- Bdata[colnames(Bdata) %in% covs]
      for (j in 1:ncol(zz))
      { jj <- as.numeric(which(colnames(Bdata)==colnames(zz)[j]))
      	if (is.factor(Bdata[,jj])) stop ("The date-covariate is a factor variable, not a date.")
      	d <- as.numeric(Bdata[,jj])
      	d[d==0] <- NA
      	Bdata[,jj] <- date_convert (d=Bdata[,jj],format.in,selectday,format.out,born=born,format.born=format.born)
      }} 
   
  # convert birtdate into format.out format (e.g. year)   
  if (format.out!="age")  
    { Bdata$born <- date_convert(Bdata$born,format.born,selectday,format.out=format.out,born=Bdata$born,format.born=format.born)
      attr(Bdata,"format.born") <- format.out
      attr(Bdata,"param")$format.born <- format.out
    }
     attr(Bdata,"format.date") <- format.out
     attr(Bdata,"param")$format.date <- format.out
   return(date=Bdata)
}
