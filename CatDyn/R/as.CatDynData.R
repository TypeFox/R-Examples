as.CatDynData <-
function(x,step,fleet.name,coleff,colcat,colmbw,unitseff,unitscat,unitsmbw,nmult,season.dates)
   {
    if(step != "day" && step != "week" && step != "month")
      {stop("step for the catch dynamics must be 'day', 'week', or 'month' ")}
    if(length(fleet.name) < 1 || length(fleet.name) > 2){stop("The number of fleets must be 1 or 2")}
    if(length(unique(fleet.name)) != length(fleet.name)){stop("Fleets must have unique names")}
    if(any(x[,coleff] < 0) || any(x[,colcat] < 0) || any(x[,colmbw] < 0))
      {stop("Catch, effort, and mean body weight observations must all be non-negative")}
    if(sum(is.na(x[,c(coleff,colcat)])) > 0){stop("NAs are not allowed in the catch, effort series")}
    if(any(coleff > dim(x)[2])){stop("coleff are the columns where the effort data is stored")}
    if(any(colcat > dim(x)[2])){stop("colcat are the columns where the catch data is stored")}
    if(any(colmbw > dim(x)[2])){stop("colmbw are the columns where the mean body weight is stored")}
    if(unitscat != "kg" && unitscat != "ton" && unitscat != "ind") {stop("unitscat must be a character (all fleets same unit) either kilograms ('kg'), tonnes ('ton') or numbers ('ind')")}
    if(unitsmbw != "kg" && unitsmbw != "g" && unitsmbw != "ind") {stop("unitsmbw must be a character (all fleets same unit) either grams ('g'), kilograms ('kg') or numbers ('ind')")}
    if(length(unitscat) > 1 || length(unitsmbw) > 1) {stop("The catches and mean body weight of all fleets must be in the same units, thus this parameter must be of length 1")}
    if(nmult != "bill" && nmult != "mill" && nmult != "thou" && nmult != "hund" && nmult != "ind") {stop("nmult, the catch in numbers multiplier, must be 'bill' (billions), 'mill' (millions), 'thou' (thousands), 'hund' (hundreds), or 'ind' (single individuals)")}
    if(length(nmult) > 1) {stop("Since the catches and mean body weight are in the same unit for both fleets, the nmult parameter must be of length 1")}
    if(paste(substr(season.dates[1],5,5),substr(season.dates[1],8,8),sep="") != "--" || paste(substr(season.dates[2],5,5),substr(season.dates[2],8,8),sep="") != "--")
      {stop("season.dates must be a length-2 character vector (start of season, end of season) with dates in the ISO 8601 international standard, i.e. yyyy-mm-dd")}
    cateff                         <- vector("list",2);
    names(cateff)                  <- c("Properties","Data");
    cateff$Data                    <- vector("list",length(fleet.name));
    names(cateff$Data)             <- fleet.name
    cateff$Properties              <- vector("list",3)
    names(cateff$Properties)       <- c("Units","Fleets","Dates")
    cateff$Properties$Units        <- c(step,unitscat,unitsmbw,nmult)
    names(cateff$Properties$Units) <- c("Time Step","Biomass","Bodymass","NumbersMultiplier")
    cateff$Properties$Fleets       <- data.frame(Fleet=fleet.name,Units=unitseff,stringsAsFactors=FALSE)
    cateff$Properties$Dates        <- season.dates
    names(cateff$Properties$Dates) <- c("StartDate","EndDate")
    if(step == "day")
      {
       if(paste(substr(season.dates[1],5,5),substr(season.dates[1],8,8),sep="") != "--" || paste(substr(season.dates[2],5,5),substr(season.dates[2],8,8),sep="") != "--")
        {stop("dates of weekly catch dynamics must be a length 2 character vector with dates in the ISO 8601 international standard, i.e. yyyy-mm-dd")}
       if(as.numeric(format(as.Date(season.dates[2],"%Y-%m-%d"), "%Y")) - as.numeric(format(as.Date(season.dates[1],"%Y-%m-%d"), "%Y")) > 1)
         {stop("For a daily time step the modeling is seasonal, involving no more than two calendar years")}
       if(as.numeric(format(as.Date(season.dates[1],"%Y-%m-%d"), "%Y")) > as.numeric(format(as.Date(season.dates[2],"%Y-%m-%d"), "%Y")))
         {stop("Dates must be arranged in ascending order")}
       if(as.numeric(format(as.Date(season.dates[1],"%Y-%m-%d"), "%Y")) == as.numeric(format(as.Date(season.dates[2],"%Y-%m-%d"), "%Y")) &&
          as.numeric(format(as.Date(season.dates[1],"%Y-%m-%d"), "%m")) >  as.numeric(format(as.Date(season.dates[2],"%Y-%m-%d"), "%m")))
         {stop("Dates must be arranged in ascending order")}          
       if(as.numeric(format(as.Date(season.dates[1],"%Y-%m-%d"), "%Y")) == as.numeric(format(as.Date(season.dates[2],"%Y-%m-%d"), "%Y")) &&
          as.numeric(format(as.Date(season.dates[1],"%Y-%m-%d"), "%m")) == as.numeric(format(as.Date(season.dates[2],"%Y-%m-%d"), "%m")) &&
          as.numeric(format(as.Date(season.dates[1],"%Y-%m-%d"), "%d")) >  as.numeric(format(as.Date(season.dates[2],"%Y-%m-%d"), "%d")))
         {stop("Dates must be arranged in ascending order")}          
       start <- as.numeric(julian(as.Date(season.dates[1],"%Y-%m-%d"),origin=as.Date(paste(substr(season.dates[1],1,4),"-01-01",sep="")))) + 1;
       for(i in 1:length(fleet.name))
         {
          (if(as.numeric(format(as.Date(season.dates[1],"%Y-%m-%d"), "%Y")) == as.numeric(format(as.Date(season.dates[2],"%Y-%m-%d"), "%Y"))) 
            {end <-  as.numeric(julian(as.Date(season.dates[2],"%Y-%m-%d"),
                                       origin=as.Date(paste(substr(season.dates[2],1,4),"-01-01",sep=""))) + 1)
             cateff$Data[[i]] <- data.frame(time.step=start:end,x[,coleff[i]],x[,colcat[i]],x[,colmbw[i]])}
           else 
            {end <- (as.numeric(julian(as.Date(paste(substr(season.dates[1],1,4),"-12-31",sep="")),
                                       origin=as.Date(paste(substr(season.dates[1],1,4),"-01-01",sep=""))))-start+2) + 
                     as.numeric(julian(as.Date(season.dates[2],"%Y-%m-%d"),origin=as.Date(paste(substr(season.dates[2],1,4),"-01-01",sep="")))) + 1
             cateff$Data[[i]] <- data.frame(time.step=1:end,x[,coleff[i]],x[,colcat[i]],x[,colmbw[i]]);})
#         if(length(start:(start+end)) != dim(x)[1])
#           {stop("Check season.dates, you have a number of days that is not identical to the number of effort/catch/mean body weight observations \n Consider changing the dates")}
#         cateff$Data[[i]]              <- data.frame(time.step=start:end,x[,coleff[i]],x[,colcat[i]],x[,colmbw[i]]);
         names(cateff$Data[[i]])[2:4]  <-c(paste("obseff.",unitseff[i],sep=""),
                                           paste("obscat.",unitscat,sep=""),
                                           paste("obsmbw.",unitsmbw,sep="")); 
         if(nmult == "bill" && unitscat == "ton" && unitsmbw == "g")
          { 
           cateff$Data[[i]]$obscat.bill    <- 1e-9*1e6*cateff$Data[[i]][,3]/cateff$Data[[i]][,4];
          }
         else if(nmult == "bill" && unitscat == "ton" && unitsmbw == "kg")
          {
           cateff$Data[[i]]$obscat.bill    <- 1e-9*1e3*cateff$Data[[i]][,3]/cateff$Data[[i]][,4];
          }
         else if(nmult == "bill" && unitscat == "kg" && unitsmbw == "g")
          {
           cateff$Data[[i]]$obscat.bill    <- 1e-9*1e3*cateff$Data[[i]][,3]/cateff$Data[[i]][,4];
          }
         else if(nmult == "bill" && unitscat == "kg" && unitsmbw == "kg")
          {
           cateff$Data[[i]]$obscat.bill    <- 1e-9*cateff$Data[[i]][,3]/cateff$Data[[i]][,4];
          }
          else if(nmult == "mill" && unitscat == "ton" && unitsmbw == "g")
          { 
           cateff$Data[[i]]$obscat.mill    <- 1e-6*1e6*cateff$Data[[i]][,3]/cateff$Data[[i]][,4];
          }
         else if(nmult == "mill" && unitscat == "ton" && unitsmbw == "kg")
          {
           cateff$Data[[i]]$obscat.mill    <- 1e-6*1e3*cateff$Data[[i]][,3]/cateff$Data[[i]][,4];
          }
         else if(nmult == "mill" && unitscat == "kg" && unitsmbw == "g")
          {
           cateff$Data[[i]]$obscat.mill    <- 1e-6*1e3*cateff$Data[[i]][,3]/cateff$Data[[i]][,4];
          }
         else if(nmult == "mill" && unitscat == "kg" && unitsmbw == "kg")
          {
           cateff$Data[[i]]$obscat.mill    <- 1e-6*cateff$Data[[i]][,3]/cateff$Data[[i]][,4];
          }
         else if(nmult == "thou" && unitscat == "ton" && unitsmbw == "g")
          {
           cateff$Data[[i]]$obscat.thou    <- 1e-3*1e6*cateff$Data[[i]][,3]/cateff$Data[[i]][,4];
          }
         else if(nmult == "thou" && unitscat == "ton" && unitsmbw == "kg")
          {
           cateff$Data[[i]]$obscat.thou    <- 1e-3*1e3*cateff[,3][[i]]/cateff$Data[[i]][,4];
          }
         else if(nmult == "thou" && unitscat == "kg" && unitsmbw == "g")
          {
           cateff$Data[[i]]$obscat.thou    <- 1e-3*1e3*cateff$Data[[i]][,3]/cateff$Data[[i]][,4];
          }
         else if(nmult == "thou" && unitscat == "kg" && unitsmbw == "kg")
          {
           cateff$Data[[i]]$obscat.thou    <- 1e-3*cateff$Data[[i]][,3]/cateff$Data[[i]][,4];
          }
         else if(nmult == "hund" && unitscat == "ton" && unitsmbw == "g")
          {
           cateff$Data[[i]]$obscat.hund    <- 1e-2*1e6*cateff$Data[[i]][,3]/cateff$Data[[i]][,4];
          }
         else if(nmult == "hund" && unitscat == "ton" && unitsmbw == "kg")
          {
           cateff$Data[[i]]$obscat.hund    <- 1e-2*1e3*cateff$Data[[i]][,3]/cateff$Data[[i]][,4];
          }
         else if(nmult == "hund" && unitscat == "kg" && unitsmbw == "g")
          {
           cateff$Data[[i]]$obscat.hund    <- 1e-2*1e3*cateff$Data[[i]][,3]/cateff$Data[[i]][,4];
          }
         else if(nmult == "hund" && unitscat == "kg" && unitsmbw == "kg")
          {
           cateff$Data[[i]]$obscat.hund    <- 1e-2*cateff$Data[[i]][,3]/cateff$Data[[i]][,4];
          }
         else if(nmult == "bill" && unitscat == "ind" && unitsmbw == "ind")
          {
           cateff$Data[[i]]$obscat.bill    <- 1e-9*cateff$Data[[i]][,3]/cateff$Data[[i]][,4];
          }
         else if(nmult == "mill" && unitscat == "ind" && unitsmbw == "ind")
          {
           cateff$Data[[i]]$obscat.mill    <- 1e-6*cateff$Data[[i]][,3]/cateff$Data[[i]][,4];
          }
         else if(nmult == "thou" && unitscat == "ind" && unitsmbw == "ind")
          {
           cateff$Data[[i]]$obscat.thou    <- 1e-3*cateff$Data[[i]][,3]/cateff$Data[[i]][,4];
          }
         else if(nmult == "hund" && unitscat == "ind" && unitsmbw == "ind")
          {
           cateff$Data[[i]]$obscat.hund    <- 1e-2*cateff$Data[[i]][,3]/cateff$Data[[i]][,4];
          }
         else if(nmult == "ind" && unitscat == "ind" && unitsmbw == "ind")
          {
           cateff$Data[[i]]$obscat.ind2    <- cateff$Data[[i]][,3]/cateff$Data[[i]][,4];
          }
         cateff$Data[[i]]$spikecat  <- unlist(10*(cateff$Data[[i]][5]/max(cateff$Data[[i]][5])-cateff$Data[[i]][2]/max(cateff$Data[[i]][2])));
        }
      }
    else if(step == "week")
      {
       if(paste(substr(season.dates[1],5,5),substr(season.dates[1],8,8),sep="") != "--" || paste(substr(season.dates[2],5,5),substr(season.dates[2],8,8),sep="") != "--")
        {stop("dates of weekly catch dynamics must be a length 2 character vector with dates in the ISO 8601 international standard, i.e. yyyy-mm-dd")}
       if(as.numeric(format(as.Date(season.dates[2],"%Y-%m-%d"), "%Y")) - as.numeric(format(as.Date(season.dates[1],"%Y-%m-%d"), "%Y")) > 1)
         {stop("For a weekly time step the modeling is seasonal, involving no more than two calendar years")}
       if(as.numeric(format(as.Date(season.dates[1],"%Y-%m-%d"), "%Y")) > as.numeric(format(as.Date(season.dates[2],"%Y-%m-%d"), "%Y")))
         {stop("Dates must be arranged in ascending order")}
       if(as.numeric(format(as.Date(season.dates[1],"%Y-%m-%d"), "%Y")) == as.numeric(format(as.Date(season.dates[2],"%Y-%m-%d"), "%Y")) &&
          as.numeric(format(as.Date(season.dates[1],"%Y-%m-%d"), "%m")) >  as.numeric(format(as.Date(season.dates[2],"%Y-%m-%d"), "%m")))
         {stop("Dates must be arranged in ascending order")}          
       if(as.numeric(format(as.Date(season.dates[1],"%Y-%m-%d"), "%Y")) == as.numeric(format(as.Date(season.dates[2],"%Y-%m-%d"), "%Y")) &&
          as.numeric(format(as.Date(season.dates[1],"%Y-%m-%d"), "%m")) == as.numeric(format(as.Date(season.dates[2],"%Y-%m-%d"), "%m")) &&
          as.numeric(format(as.Date(season.dates[1],"%Y-%m-%d"), "%d")) >  as.numeric(format(as.Date(season.dates[2],"%Y-%m-%d"), "%d")))
         {stop("Dates must be arranged in ascending order")}          
       start         <- as.numeric(format(as.Date(season.dates[1],"%Y-%m-%d"),"%W"));
       (if(as.numeric(format(as.Date(season.dates[1],"%Y-%m-%d"), "%Y")) == as.numeric(format(as.Date(season.dates[2],"%Y-%m-%d"), "%Y"))) 
          {end <- as.numeric(format(as.Date(season.dates[2],"%Y-%m-%d"),"%W"))}
        else 
          {end <- as.numeric(format(as.Date(paste(substr(season.dates[1],1,4),"-12-31",sep="")),"%W")) + 1 +
                     as.numeric(format(as.Date(season.dates[2],"%Y-%m-%d"),"%W"))})
       for(i in 1:length(fleet.name))
         {
         if(length(start:end) != length(x[,coleff[i]]))
           {stop("Check season.dates, you have a number of weeks that is not identical to the number of effort/catch/mean body weight observations \n Consider changing the dates")}
         cateff$Data[[i]]              <- data.frame(time.step=start:end,x[,coleff[i]],x[,colcat[i]],x[,colmbw[i]]);
         names(cateff$Data[[i]])[2:4]  <-c(paste("obseff.",unitseff[i],sep=""),
                                           paste("obscat.",unitscat,sep=""),
                                           paste("obsmbw.",unitsmbw,sep="")); 
         if(nmult == "bill" && unitscat == "ton" && unitsmbw == "g")
          { 
           cateff$Data[[i]]$obscat.bill    <- 1e-9*1e6*cateff$Data[[i]][,3]/cateff$Data[[i]][,4];
          }
         else if(nmult == "bill" && unitscat == "ton" && unitsmbw == "kg")
          {
           cateff$Data[[i]]$obscat.bill    <- 1e-9*1e3*cateff$Data[[i]][,3]/cateff$Data[[i]][,4];
          }
         else if(nmult == "bill" && unitscat == "kg" && unitsmbw == "g")
          {
           cateff$Data[[i]]$obscat.bill    <- 1e-9*1e3*cateff$Data[[i]][,3]/cateff$Data[[i]][,4];
          }
         else if(nmult == "bill" && unitscat == "kg" && unitsmbw == "kg")
          {
           cateff$Data[[i]]$obscat.bill    <- 1e-9*cateff$Data[[i]][,3]/cateff$Data[[i]][,4];
          }
         else if(nmult == "mill" && unitscat == "ton" && unitsmbw == "g")
          { 
           cateff$Data[[i]]$obscat.mill    <- 1e-6*1e6*cateff$Data[[i]][,3]/cateff$Data[[i]][,4];
          }
         else if(nmult == "mill" && unitscat == "ton" && unitsmbw == "kg")
          {
           cateff$Data[[i]]$obscat.mill    <- 1e-6*1e3*cateff$Data[[i]][,3]/cateff$Data[[i]][,4];
          }
         else if(nmult == "mill" && unitscat == "kg" && unitsmbw == "g")
          {
           cateff$Data[[i]]$obscat.mill    <- 1e-6*1e3*cateff$Data[[i]][,3]/cateff$Data[[i]][,4];
          }
         else if(nmult == "mill" && unitscat == "kg" && unitsmbw == "kg")
          {
           cateff$Data[[i]]$obscat.mill    <- 1e-6*cateff$Data[[i]][,3]/cateff$Data[[i]][,4];
          }
         else if(nmult == "thou" && unitscat == "ton" && unitsmbw == "g")
          {
           cateff$Data[[i]]$obscat.thou    <- 1e-3*1e6*cateff$Data[[i]][,3]/cateff$Data[[i]][,4];
          }
         else if(nmult == "thou" && unitscat == "ton" && unitsmbw == "kg")
          {
           cateff$Data[[i]]$obscat.thou    <- 1e-3*1e3*cateff[,4][[i]]/cateff$Data[[i]][,4];
          }
         else if(nmult == "thou" && unitscat == "kg" && unitsmbw == "g")
          {
           cateff$Data[[i]]$obscat.thou    <- 1e-3*1e3*cateff$Data[[i]][,3]/cateff$Data[[i]][,4];
          }
         else if(nmult == "thou" && unitscat == "kg" && unitsmbw == "kg")
          {
           cateff$Data[[i]]$obscat.thou    <- 1e-3*cateff$Data[[i]][,3]/cateff$Data[[i]][,4];
          }
         else if(nmult == "hund" && unitscat == "ton" && unitsmbw == "g")
          {
           cateff$Data[[i]]$obscat.hund    <- 1e-2*1e6*cateff$Data[[i]][,3]/cateff$Data[[i]][,4];
          }
         else if(nmult == "hund" && unitscat == "ton" && unitsmbw == "kg")
          {
           cateff$Data[[i]]$obscat.hund    <- 1e-2*1e3*cateff$Data[[i]][,3]/cateff$Data[[i]][,4];
          }
         else if(nmult == "hund" && unitscat == "kg" && unitsmbw == "g")
          {
           cateff$Data[[i]]$obscat.hund    <- 1e-2*1e3*cateff$Data[[i]][,3]/cateff$Data[[i]][,4];
          }
         else if(nmult == "hund" && unitscat == "kg" && unitsmbw == "kg")
          {
           cateff$Data[[i]]$obscat.hund    <- 1e-2*cateff$Data[[i]][,3]/cateff$Data[[i]][,4];
          }
         else if(nmult == "bill" && unitscat == "ind" && unitsmbw == "ind")
          {
           cateff$Data[[i]]$obscat.bill    <- 1e-9*cateff$Data[[i]][,3]/cateff$Data[[i]][,4];
          }
         else if(nmult == "mill" && unitscat == "ind" && unitsmbw == "ind")
          {
           cateff$Data[[i]]$obscat.mill    <- 1e-6*cateff$Data[[i]][,3]/cateff$Data[[i]][,4];
          }
         else if(nmult == "thou" && unitscat == "ind" && unitsmbw == "ind")
          {
           cateff$Data[[i]]$obscat.thou    <- 1e-3*cateff$Data[[i]][,3]/cateff$Data[[i]][,4];
          }
         else if(nmult == "hund" && unitscat == "ind" && unitsmbw == "ind")
          {
           cateff$Data[[i]]$obscat.hund    <- 1e-2*cateff$Data[[i]][,3]/cateff$Data[[i]][,4];
          }
         else if(nmult == "ind" && unitscat == "ind" && unitsmbw == "ind")
          {
           cateff$Data[[i]]$obscat.ind2    <- cateff$Data[[i]][,3]/cateff$Data[[i]][,4];
          }
         cateff$Data[[i]]$spikecat  <- unlist(10*(cateff$Data[[i]][5]/max(cateff$Data[[i]][5])-cateff$Data[[i]][2]/max(cateff$Data[[i]][2])));
        }
      }
    else
      {
       if(paste(substr(season.dates[1],5,5),substr(season.dates[1],8,8),sep="") != "--" || paste(substr(season.dates[2],5,5),substr(season.dates[2],8,8),sep="") != "--")
         {stop("dates of monthly catch dynamics must be a length 2 character vector with dates in the ISO 8601 international standard, i.e. yyyy-mm-dd")}
       if(as.numeric(format(as.Date(season.dates[1],"%Y-%m-%d"), "%Y")) > as.numeric(format(as.Date(season.dates[2],"%Y-%m-%d"), "%Y")))
         {stop("Dates must be arranged in ascending order")}
       if(as.numeric(format(as.Date(season.dates[1],"%Y-%m-%d"), "%Y")) == as.numeric(format(as.Date(season.dates[2],"%Y-%m-%d"), "%Y")) &&
          as.numeric(format(as.Date(season.dates[1],"%Y-%m-%d"), "%m")) >  as.numeric(format(as.Date(season.dates[2],"%Y-%m-%d"), "%m")))
         {stop("Dates must be arranged in ascending order")}          
       if(as.numeric(format(as.Date(season.dates[1],"%Y-%m-%d"), "%Y")) == as.numeric(format(as.Date(season.dates[2],"%Y-%m-%d"), "%Y")) &&
          as.numeric(format(as.Date(season.dates[1],"%Y-%m-%d"), "%m")) == as.numeric(format(as.Date(season.dates[2],"%Y-%m-%d"), "%m")) &&
          as.numeric(format(as.Date(season.dates[1],"%Y-%m-%d"), "%d")) >  as.numeric(format(as.Date(season.dates[2],"%Y-%m-%d"), "%d")))
         {stop("Dates must be arranged in ascending order")}          
       start         <- as.numeric(format(as.Date(season.dates[1],"%Y-%m-%d"),"%m"));
       end           <- as.numeric(format(as.Date(season.dates[2],"%Y-%m-%d"),"%m"))*
                       (as.numeric(format(as.Date(season.dates[2],"%Y-%m-%d"),"%Y"))-as.numeric(format(as.Date(season.dates[1],"%Y-%m-%d"),"%Y"))+1);
       if(start >= end) {stop("Dates must be arranged in ascending order")}
       for(i in 1:length(fleet.name))
         {
         cateff$Data[[i]]              <- data.frame(time.step=seq(start,end,1),x[,coleff[i]],x[,colcat[i]],x[,colmbw[i]]);
         names(cateff$Data[[i]])[2:4]  <-c(paste("obseff.",unitseff[i],sep=""),
                                           paste("obscat.",unitscat,sep=""),
                                           paste("obsmbw.",unitsmbw,sep="")); 
         if(nmult == "bill" && unitscat == "ton" && unitsmbw == "g")
          { 
           cateff$Data[[i]]$obscat.bill    <- 1e-9*1e6*cateff$Data[[i]][,3]/cateff$Data[[i]][,4];
          }
         else if(nmult == "bill" && unitscat == "ton" && unitsmbw == "kg")
          {
           cateff$Data[[i]]$obscat.bill    <- 1e-9*1e3*cateff$Data[[i]][,3]/cateff$Data[[i]][,4];
          }
         else if(nmult == "bill" && unitscat == "kg" && unitsmbw == "g")
          {
           cateff$Data[[i]]$obscat.bill    <- 1e-9*1e3*cateff$Data[[i]][,3]/cateff$Data[[i]][,4];
          }
         else if(nmult == "bill" && unitscat == "kg" && unitsmbw == "kg")
          {
           cateff$Data[[i]]$obscat.bill    <- 1e-9*cateff$Data[[i]][,3]/cateff$Data[[i]][,4];
          }
         else if(nmult == "mill" && unitscat == "ton" && unitsmbw == "g")
          { 
           cateff$Data[[i]]$obscat.mill    <- 1e-6*1e6*cateff$Data[[i]][,3]/cateff$Data[[i]][,4];
          }
         else if(nmult == "mill" && unitscat == "ton" && unitsmbw == "kg")
          {
           cateff$Data[[i]]$obscat.mill    <- 1e-6*1e3*cateff$Data[[i]][,3]/cateff$Data[[i]][,4];
          }
         else if(nmult == "mill" && unitscat == "kg" && unitsmbw == "g")
          {
           cateff$Data[[i]]$obscat.mill    <- 1e-6*1e3*cateff$Data[[i]][,3]/cateff$Data[[i]][,4];
          }
         else if(nmult == "mill" && unitscat == "kg" && unitsmbw == "kg")
          {
           cateff$Data[[i]]$obscat.mill    <- 1e-6*cateff$Data[[i]][,3]/cateff$Data[[i]][,4];
          }
         else if(nmult == "thou" && unitscat == "ton" && unitsmbw == "g")
          {
           cateff$Data[[i]]$obscat.thou    <- 1e-3*1e6*cateff$Data[[i]][,3]/cateff$Data[[i]][,4];
          }
         else if(nmult == "thou" && unitscat == "ton" && unitsmbw == "kg")
          {
           cateff$Data[[i]]$obscat.thou    <- 1e-3*1e3*cateff[,3][[i]]/cateff$Data[[i]][,4];
          }
         else if(nmult == "thou" && unitscat == "kg" && unitsmbw == "g")
          {
           cateff$Data[[i]]$obscat.thou    <- 1e-3*1e3*cateff$Data[[i]][,3]/cateff$Data[[i]][,4];
          }
         else if(nmult == "thou" && unitscat == "kg" && unitsmbw == "kg")
          {
           cateff$Data[[i]]$obscat.thou    <- 1e-3*cateff$Data[[i]][,3]/cateff$Data[[i]][,4];
          }
         else if(nmult == "hund" && unitscat == "ton" && unitsmbw == "g")
          {
           cateff$Data[[i]]$obscat.hund    <- 1e-2*1e6*cateff$Data[[i]][,3]/cateff$Data[[i]][,4];
          }
         else if(nmult == "hund" && unitscat == "ton" && unitsmbw == "kg")
          {
           cateff$Data[[i]]$obscat.hund    <- 1e-2*1e3*cateff$Data[[i]][,3]/cateff$Data[[i]][,4];
          }
         else if(nmult == "hund" && unitscat == "kg" && unitsmbw == "g")
          {
           cateff$Data[[i]]$obscat.hund    <- 1e-2*1e3*cateff$Data[[i]][,3]/cateff$Data[[i]][,4];
          }
         else if(nmult == "hund" && unitscat == "kg" && unitsmbw == "kg")
          {
           cateff$Data[[i]]$obscat.hund    <- 1e-2*cateff$Data[[i]][,3]/cateff$Data[[i]][,4];
          }
         else if(nmult == "bill" && unitscat == "ind" && unitsmbw == "ind")
          {
           cateff$Data[[i]]$obscat.bill    <- 1e-9*cateff$Data[[i]][,3]/cateff$Data[[i]][,4];
          }
         else if(nmult == "mill" && unitscat == "ind" && unitsmbw == "ind")
          {
           cateff$Data[[i]]$obscat.mill    <- 1e-6*cateff$Data[[i]][,3]/cateff$Data[[i]][,4];
          }
         else if(nmult == "thou" && unitscat == "ind" && unitsmbw == "ind")
          {
           cateff$Data[[i]]$obscat.thou    <- 1e-3*cateff$Data[[i]][,3]/cateff$Data[[i]][,4];
          }
         else if(nmult == "hund" && unitscat == "ind" && unitsmbw == "ind")
          {
           cateff$Data[[i]]$obscat.hund    <- 1e-2*cateff$Data[[i]][,3]/cateff$Data[[i]][,4];
          }
         else if(nmult == "ind" && unitscat == "ind" && unitsmbw == "ind")
          {
           cateff$Data[[i]]$obscat.ind2    <- cateff$Data[[i]][,3]/cateff$Data[[i]][,4];
          }
         cateff$Data[[i]]$spikecat  <- unlist(10*(cateff$Data[[i]][5]/max(cateff$Data[[i]][5])-cateff$Data[[i]][2]/max(cateff$Data[[i]][2])));
        }
      }
    oldClass(cateff) <- "CatDynData" 
    return(cateff)
   }
