as.GroAgeData <-
function(x, sex, maleskey=NULL, femaleskey=NULL, colsex=NULL, colage, collen, colbw=NULL, colliver=NULL, colgonad=NULL, coldate=NULL, lentype, unitsage, unitslen, unitsbw=NULL, unitsliver=NULL, unitsgonad=NULL, spec)
  {
    sex.set <- c("Females","Males","Pooled","Both","Total")
    units.len.set <- c("mm", "cm", "m")
    units.age.set <- c("Days","Weeks","Months","Years")
    units.wg.set  <- c("g","kg")
    lentype.set   <- c("Total","Standard","Fork")
    if(!sex%in%sex.set) {stop("sex must be either 'F' (females), 'M' (males), 'B' (both), 'P' (pooled), or 'T' (total, included unsexed individuals)")}
    if(sex == "Both" & is.null(maleskey) | sex == "Both" & is.null(colsex)) {stop("When sex='Both' a column number and an integer key has to be provided to identify sexes")}
    if(sex == "Both" & is.null(femaleskey) | sex == "Both" & is.null(colsex)) {stop("When sex='Both' a column number and an integer key has to be provided to identify sexes")}
    if(sex == "Total" & is.null(maleskey) | sex == "Total" & is.null(colsex)) {stop("When sex='Total' (males, females and unsexed) a column number and an integer key has to be provided to identify sexes")}
    if(sex == "Total" & is.null(femaleskey) | sex == "Total" &  is.null(colsex)) {stop("When sex='Total' (males, females and unsexed) a column number and an integer key has to be provided to identify sexes")}
    if(any(unlist(c(x)) < 0, na.rm=TRUE)) {stop("Age and length data must be all non-negative")}
    if(!is.null(colbw) & is.null(unitsbw)) {stop("If the optional body weight data is included, units for it must be provided")}
    if(!is.null(colliver) & is.null(unitsliver)) {stop("If the optional liver weight data is included, units for it must be provided")}
    if(!is.null(colgonad) & is.null(unitsgonad)) {stop("If the optional gonad weight data is included, units for it must be provided")}
    if(is.null(unitsage) | is.null(unitslen) | is.null(spec)) {stop("Units for age and length, as well as species identification, must be provided")}
    if(!unitslen%in%units.len.set) {stop("Units of body length must be from the MKS system of units: 'mm', 'cm', or 'm'")}
    if(!unitsage%in%units.age.set) {stop("Age must be in 'days', 'weeks', 'months', or 'years'")}
    if(!lentype%in%lentype.set) {stop("The type of length measurement must be either 'Total', 'Standard', or 'Fork'")}
    groage            <- vector("list",2);
    names(groage)     <- c("Properties","Data");
    groage$Properties <- data.frame(Species=spec,
                                    Sex=sex,
                                    Units.Age=unitsage,
                                    Units.Length=unitslen,
                                    Length.Type=lentype,
                                    stringsAsFactors=FALSE,
                                    Units.Weight=NA,
                                    Units.Gonad=NA,
                                    Units.Liver=NA,
                                    Start.Month=NA,
                                    End.Month=NA)
    if(!is.null(colbw))
      {
       groage$Properties$Units.Weight <- unitsbw
      }
    if(!is.null(colgonad))
      {
       groage$Properties$Units.Gonad <- unitsgonad
      }
    if(!is.null(colliver))
      {
       groage$Properties$Units.Liver <- unitsliver
      }
    if(!is.null(coldate))
      {
       groage$Properties$Start.Month <- min(as.numeric(format(x[,coldate], "%m")))
       groage$Properties$End.Month   <- max(as.numeric(format(x[,coldate], "%m")))
      }
    if(sex == "Females" | sex == "Males" | sex == "Pooled")
      {
       groage$Data        <- data.frame(1:dim(x)[1],x[,colage],x[,collen])
       names(groage$Data) <-c("Individual","Age","Length");
       if(!is.null(colbw))
         {
          groage$Data$Weight        <- x[,colbw];
          names(groage$Data$Weight) <- paste("Weight",unitsliver,sep=".");
         }
       if(!is.null(colgonad))
         {
          groage$Data$Gonad        <- x[,colgonad];
          names(groage$Data$Gonad) <- paste("Gonad",unitsgonad,sep=".");
         }
       if(!is.null(colliver))
         {
          groage$Data$Liver        <- x[,colliver];
          names(groage$Data$Liver) <- paste("Liver",unitsliver,sep=".");
         }
       if(!is.null(coldate))
         {
          groage$Data$Date    <- x[,coldate];
         }
      }
    else if(sex == "Both")
      {
       groage$Data              <- vector("list",2);
       names(groage$Data)       <- c("Males","Females");
       #
       groage$Data$Males        <- data.frame(1:dim(x[which(x[,colsex]==maleskey),])[1],x[which(x[,colsex]==maleskey),colage],x[which(x[,colsex]==maleskey),collen])
       names(groage$Data$Males) <- c("Individual","Age","Length");
       if(!is.null(colbw))
         {
          groage$Data$Males$Weight        <- x[which(x[,colsex]==maleskey),colbw];
         }
       if(!is.null(colgonad))
         {
          groage$Data$Males$Gonad        <- x[which(x[,colsex]==maleskey),colgonad];
         }
       if(!is.null(colliver))
         {
          groage$Data$Males$Liver        <- x[which(x[,colsex]==maleskey),colliver];
         }
       if(!is.null(coldate))
         {
          groage$Data$Males$Date <- x[which(x[,colsex]==maleskey),coldate];
         }
       groage$Data$Females        <- data.frame(1:dim(x[which(x[,colsex]==femaleskey),])[1],x[which(x[,colsex]==femaleskey),colage],x[which(x[,colsex]==femaleskey),collen])
       names(groage$Data$Females) <- c("Individual","Age","Length");
       if(!is.null(colbw))
         {
          groage$Data$Females$Weight        <- x[which(x[,colsex]==femaleskey),colbw];
         }
       if(!is.null(colgonad))
         {
          groage$Data$Females$Gonad        <- x[which(x[,colsex]==femaleskey),colgonad];
         }
       if(!is.null(colliver))
         {
          groage$Data$Females$Liver        <- x[which(x[,colsex]==femaleskey),colliver];
         }
       if(!is.null(coldate))
         {
          groage$Data$Females$Date <- x[which(x[,colsex]==femaleskey),coldate];
         }
      }
    else
      {
       groage$Data              <- vector("list",3);
       names(groage$Data)       <- c("Males","Females","Unsexed");
       #
       groage$Data$Males        <- data.frame(1:dim(x[which(x[,colsex]==maleskey),])[1],x[which(x[,colsex]==maleskey),colage],x[which(x[,colsex]==maleskey),collen])
       names(groage$Data$Males) <- c("Individual","Age","Length");
       if(!is.null(colbw))
         {
          groage$Data$Males$Weight        <- x[which(x[,colsex]==maleskey),colbw];
         }
       if(!is.null(colgonad))
         {
          groage$Data$Males$Gonad        <- x[which(x[,colsex]==maleskey),colgonad];
         }
       if(!is.null(colliver))
         {
          groage$Data$Males$Liver        <- x[which(x[,colsex]==maleskey),colliver];
         }
       if(!is.null(coldate))
         {
          groage$Data$Males$Date <- x[which(x[,colsex]==maleskey),coldate];
         }
       groage$Data$Females        <- data.frame(1:dim(x[which(x[,colsex]==femaleskey),])[1],x[which(x[,colsex]==femaleskey),colage],x[which(x[,colsex]==femaleskey),collen])
       names(groage$Data$Females) <- c("Individual","Age","Length");
       if(!is.null(colbw))
         {
          groage$Data$Females$Weight        <- x[which(x[,colsex]==femaleskey),colbw];
         }
       if(!is.null(colgonad))
         {
          groage$Data$Females$Gonad        <- x[which(x[,colsex]==femaleskey),colgonad];
         }
       if(!is.null(colliver))
         {
          groage$Data$Females$Liver        <- x[which(x[,colsex]==femaleskey),colliver];
         }
       if(!is.null(coldate))
         {
          groage$Data$Females$Date <- x[which(x[,colsex]==femaleskey),coldate];
         }
       groage$Data$Unsexed      <- data.frame(1:dim(x[-which(x[,colsex] == maleskey | x[,colsex] == femaleskey),])[1],x[-which(x[,colsex] == maleskey | x[,colsex] == femaleskey),colage],x[-which(x[,colsex] == maleskey | x[,colsex] == femaleskey),collen])
       names(groage$Data$Unsexed) <- c("Individual","Age","Length");
       if(!is.null(colbw))
         {
          groage$Data$Unsexed$Weight        <- x[-which(x[,colsex] == maleskey | x[,colsex] == femaleskey),colbw];
         }
       if(!is.null(colgonad))
         {
          groage$Data$Unsexed$Gonad        <- x[-which(x[,colsex] == maleskey | x[,colsex] == femaleskey),colgonad];
         }
       if(!is.null(colliver))
         {
          groage$Data$Unsexed$Liver        <- x[-which(x[,colsex] == maleskey | x[,colsex] == femaleskey),colliver];
         }
       if(!is.null(coldate))
         {
          groage$Data$Unsexed$Date <- x[-which(x[,colsex] == maleskey | x[,colsex] == femaleskey),coldate];
         }
      }
    oldClass(groage) <- "GroAgeData"
    return(groage)
  }
