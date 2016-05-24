drFit <-
function(gcFitData, control=grofit.control())
{

# /// check input values and initialize some variables
if (is(control)!="grofit.control") stop("control must be of class grofit.control!")

EC50.table  <- NULL
all.EC50    <- NA
table.tests <- table((gcFitData[,1])[which((gcFitData[,4]==TRUE)&(is.na(gcFitData[,control$parameter])==FALSE))])
distinct    <- names(table.tests)
EC50        <- list()
EC50.boot   <- list()

validdata   <- cbind(as.character(distinct), table.tests)
colnames(validdata) <- c("TestID", "Number")
rownames(validdata) <- rep("     ", times = dim(validdata)[1])

# /// print informations about submitted data to screen
if (control$suppress.messages==FALSE){
	cat("\n")
	cat("=== EC 50 Estimation ==============================\n")
	cat("---------------------------------------------------\n")
	cat("--> Checking data ...\n")
	cat(paste("--> Number of distinct tests found:", as.character(length(distinct))),"\n")
	cat("--> Valid datasets per test: \n")
	print(validdata, quote = FALSE)
}

if (TRUE%in%(table.tests<control$have.atleast)){
	cat(paste("Warning: following tests have not enough ( <", as.character(control$have.atleast-1),") datasets:\n"))
	cat(distinct[(table.tests<control$have.atleast)])
	cat("These tests will not be regarded\n")
	distinct <- distinct[table.tests>=control$have.atleast]
}

if ( (length(distinct)) == 0 ){
    	cat(paste("There are no tests having enough ( >", as.character(control$have.atleast-1), ") datasets!\n"))
}
else{
	# /// Loop over all tests to obtain DR curve
	for (i in 1:length(distinct)){
		# /// choose current data
		conc <- (gcFitData[,3])[which(gcFitData[,1]==distinct[i])]
		test <- (gcFitData[,control$parameter])[which(gcFitData[,1]==distinct[i])]
		drID <- distinct[i]

		# /// fit spline
		EC50[[i]] <- drFitSpline(conc, test, drID, control)
		
		# /// start bootstrapping
		if (control$nboot.dr>0){
		    EC50.boot[[i]] <- drBootSpline(conc, test, drID, control)
		}
		else{
		# /// create empty object
		    EC50.boot[[i]]        <- list(raw.x=conc, raw.y=test, drID=drID, boot.x= NA, boot.y=NA, boot.drSpline = NA, ec50.boot=NA, bootFlag=FALSE, control=control)
		    class(EC50.boot[[i]]) <- "drBootSpline"
		}
		description <- data.frame("Test"=distinct[i], "log.x"=control$log.x.dr, "log.y"=control$log.y.dr, "Samples"=control$nboot.dr)
		out.row     <- cbind(description, summary(EC50[[i]]), summary(EC50.boot[[i]]))
		EC50.table  <- rbind(EC50.table, out.row)
	} # end of for(i in...)
} # end of else of if ( (length(distinct))...

drFit        <- list(raw.data=gcFitData, drTable=EC50.table,  drBootSplines = EC50.boot, drFittedSplines=EC50, control=control)
class(drFit) <- "drFit"
drFit

}

