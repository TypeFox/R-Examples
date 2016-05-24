addDataValues <- function(DataZoo=NULL, Date=NULL, Value=NULL, ValueAccuracy=rep(NA, NCOL(DataZoo)), Site, Variable, Offset=rep(NA, NCOL(DataZoo)), OffsetType=rep('No', NCOL(DataZoo)), CensorCode = rep("nc", NCOL(DataZoo)), Qualifier=rep("No", NCOL(DataZoo)), Method=rep('No', NCOL(DataZoo)), Source, Sample=rep("No", NCOL(DataZoo)), DerivedFrom=NULL, QualityControlLevel, tolerance=0){
	#mandatory: DataValue,  Site, Variable, CensorCode, Method, Source, QualityControlLevel
	todo("UnitConversion")
	if(is.null(DataZoo) & (is.null(Date) | is.null(Value))) stop("You must provide information either as DataZoo or Date and Value")
	if(!is.null(Value)){
		#stopifnot(is.numeric(Value))
	}

	if(is.null(Date)) Date <- index(DataZoo)
	stopifnot("POSIXt" %in%  class(Date))
	if(NROW(Date)==1){Date <- rep(Date, NROW(Value))}
	Date.by.col <- NROW(Date)==NCOL(Value) 

	if(is.null(Value)) Value <- coredata(DataZoo)
	if(is.null(dim(Value))) dim(Value)  <-  c(NROW(Value),NCOL(Value) )

	#Replace empty strings by NA and issue warning
	emptyValue <- Value == ""
	emptyValue[is.na(emptyValue)]  <- FALSE
	if(any(emptyValue)){
		Value[emptyValue] <- NA
		warning("Some data values included empty strings. These are treated as NA.")
	}
	
	#check variables with foreign keys
	SiteID <- expandVar(Site, nrow=NROW(Value), ncol=NCOL(Value), checkID=TRUE, table="Site")
	VariableID <- expandVar(Variable, nrow=NROW(Value), ncol=NCOL(Value), checkID=TRUE, table="Variable")
	OffsetTypeID <- expandVar(OffsetType, nrow=NROW(Value), ncol=NCOL(Value), checkID=TRUE, table="OffsetType")
	QualifierID <- expandVar(Qualifier, nrow=NROW(Value), ncol=NCOL(Value), checkID=TRUE, table="Qualifier")
	MethodID <- expandVar(Method, nrow=NROW(Value), ncol=NCOL(Value), checkID=TRUE, table="Method")
	CensorCodeNum <- expandVar(CensorCode, nrow=NROW(Value), ncol=NCOL(Value), checkID=TRUE, table="CensorCode")
	SourceID <- expandVar(Source, nrow=NROW(Value), ncol=NCOL(Value), checkID=TRUE, table="Source")
	SampleID <- expandVar(Sample, nrow=NROW(Value), ncol=NCOL(Value), checkID=TRUE, table="Sample")
	QualityControlLevelID <- expandVar(QualityControlLevel, nrow=NROW(Value), ncol=NCOL(Value), checkID=TRUE, table="QualityControlLevel")

	#Adjust Dimension for Values without foreign key
	Offset <- expandVar(Offset, nrow=NROW(Value), ncol=NCOL(Value))
	# expandVar does not work for POSIX
	#Date <- expandVar(Date, nrow=NROW(Value), ncol=NCOL(Value))
	ValueAccuracy <- expandVar(ValueAccuracy, nrow=NROW(Value), ncol=NCOL(Value))


	
	for(column in 1:NCOL(Value)){
		cat("Importing column ", column, "out of", NCOL(Value), "\n")
		#Determine all cases of unique IDs
		metadata <- data.frame( Site=SiteID[,column], Variable=VariableID[,column], 
				Offset=Offset[,column], OffsetType=OffsetTypeID[,column], 
				CensorCode=CensorCodeNum[,column], Qualifier=QualifierID[,column],
			       	Method=MethodID[,column], Source=SourceID[,column], 
				Sample=SampleID[,column],QualityControlLevel=QualityControlLevelID[,column])
		c.m <- unique(metadata) #cases.metadata
		metadata$rownr <- 1:NROW(metadata)
		#import unique cases of meta data separately (e.g. Worldbank)
		for(ca in 1:NROW(c.m)){ #each case

			#select corresponding data from DataZoo
			row.sel <- sort(merge(c.m[ca,], metadata, all.y=FALSE)$rownr)
			if(Date.by.col){
			    import.dates.check <- Date[column]
			    if(length(row.sel)>1){
				cat("Conflicting data in import data set: multiple values with same metadata and time stamp\n")
				cat("Dates affected: ", strftime(import.dates.check), "\n")
				cat("ID of metadata:\n")
				print(c.m[ca,])
				stop()

			    }
			} else {
			    import.dates.check <- Date[row.sel]
			}
			idc.duplicates <- duplicated(import.dates.check)
			if(any(idc.duplicates)){
				cat("Conflicting data in import data set: multiple values with same metadata and time stamp\n")
				cat("Dates affected: ", strftime(unique(import.dates.check[idc.duplicates])), "\n")
				cat("ID of metadata:\n")
				print(c.m[ca,])
				stop()
			}

			#check for existing entries
			database.entries  <- NULL
			if(Date.by.col){
				time.range <- range(Date[column])
				order.date <- Date[column]

			} else {
				time.range <- range(Date)
				order.date <- Date[row.sel]
			}
      
			database.entries <- getDataValues(from=time.range[1], to=time.range[2], tz="global", Site=c.m$Site[ca], Variable=c.m$Variable[ca], Offset=c.m$Offset[ca], OffsetType=c.m$OffsetType[ca], CensorCode=c.m$CensorCodeNum[ca], Qualifier=c.m$Qualifier[ca], Method=c.m$Method[ca], Source=c.m$Source[ca], Sample=c.m$Sample[ca],QualityControlLevel=c.m$QualityControlLevel[ca], show.deleted=TRUE, all.ID=TRUE)
		  if(NROW(database.entries)>0){
			   
    entriesXtsValues = xts(database.entries@data, order.by=index(database.entries@time))
    to.test <- merge(entriesXtsValues, xts(Value[row.sel,column], order.by=order.date), join="right")
        
				stopifnot(NROW(to.test) == NROW(row.sel))
				if(NCOL(to.test)!=2){
					cat("Error while checking for duplicates. Contact the maintainer\n")

					cat("database.entries:\n")
					cat("NROW(database.entries):", NROW(database.entries), "\n")

					print(database.entries)
					str(database.entries)
					cat("\n\nValue[row.sel,column]:\n")
					print(Value[row.sel,column])
					stop("Error while checking for duplicates. Contact the maintainer\n")
				}
				names(to.test) <- c("inDatabase", "toImport")
				if(any(different <- (abs(to.test$inDatabase - to.test$toImport) > tolerance), na.rm=TRUE)){
					if(interactive()){
						plot.zoo(to.test$inDatabase, col="green", main=paste("Import of column", column), ylim=range(coredata(to.test), na.rm=TRUE), type="b")
						points(to.test$toImport[different], col="red")
						legend("top", c("database records", "differing (to import)"), col=c("green", "red"), lty=c(1,NA), pch=c(NA,1))
						cat("Data to import is not matching data in database for", sum(different, na.rm=TRUE), "values (See plot)\nWhat shall I do?\n  1) Discard data to import and import remaining, missing data\n  2) Overwrite values in database with new values\n	0) Stop and let you modify the data to import before another attempt\nEnter a number (0-2) for your choise.\n")
						choice <- "impossible"
						attempts <- 0
						if(!interactive()) {
							choice=1
							warning("non-interactive session. Not replacing data in database")
						}
						while(choice == "impossible"){
							next.step <- readline("What is your choice? ")
							choice <- switch(next.step, "1"=1,"2"=2,"0"=0, "impossible")
							attempts <- attempts + 1
							if(attempts == 10) choice <- 0
						} 
					} else {
						 choice <- 1 #Don't change too much
					}

					if(choice==0){
						cat("returning a spacetime object with the data in the database and the data to be imported\n")
						return(to.test)
					} else if (choice==2){	
						differentUpdate <- different
						differentUpdate[is.na(differentUpdate)] <- FALSE
						database.entries@data[differentUpdate,] <- as.numeric(to.test$toImport[differentUpdate])
						updateDataValues(database.entries, paste("Replacement upon import on", date()))
					}
					
					#nothing to do for choice 1 because is.na(different) is false for differing values, so data will not be imported
				}
    
				the.missing <- is.na(different)
				do.import <- !is.na(to.test$toImport) & the.missing
				to.import <- to.test$toImport[do.import]
			} else {
				to.import <- Value[row.sel,column]
				do.import <- rep(FALSE, NROW(Value))
				do.import[row.sel] <- TRUE
			}
			if(any(do.import)){
				if(Date.by.col){
					theDate <- rep(Date[column], sum(do.import))
				} else {
					theDate <- Date[do.import]
				}
				# strftime(theDate, "%z") funktioniert. As.numeric produced a warning: NAs introduced by coercion
				thetz <- strftime(theDate, "%z")
				todo("Correct functioning of tz")

				IaddDataValues(getOption("odm.handler"),localDateTime=theDate, values=to.import, TZ=thetz, SiteID=SiteID[do.import,column], VariableID=VariableID[do.import,column], Offset=Offset[do.import,column], OffsetTypeID=OffsetTypeID[do.import,column], CensorCode=CensorCodeNum[do.import], QualifierID=QualifierID[do.import,column], MethodID=MethodID[do.import,column], SourceID=SourceID[do.import,column], SampleID=SampleID[do.import,column],QualityControlLevelID=QualityControlLevelID[do.import,column], valueAccuracy=ValueAccuracy[do.import,column])
			}
					
			#import data
		}
	}



}
