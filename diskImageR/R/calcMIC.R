#' Convert RAD to MIC based on built-in or provided parameters or datasets

#' @description Used to convert RAD into MIC. In all cases the linear relationship is log2(MIC) regressed onto RAD^2. This conversion can be based on a) existing built-in data from a number of species/drug combinations, b) a user-supplied slope and intercept of log2(MIC) regressed on RAD^2 for the species/drug combination of interest, c) a user supplied file containing MIC information from lines previously analyzed by diskImageR for RAD, or d) a user supplied file containing both RAD and MIC information. Note that user-supplied data should be MIC and RAD, not log2(MIC) and RAD^2 (the function will do this authomatically).

#' @param projectName the short name you have been using for the project.
#' @inheritParams twoParamPlot
#' @inheritParams oneParamPlot
#' @inheritParams plotRaw
#' @param RAD a numeric value the the critical level of the radius of inhibition (i.e., resistance) parameter to use for MIC. Currently only \code{RAD} = "80" (80\% reduction in growth), \code{RAD} = "50" (50\% reduction in growth), and \code{RAD} = "20" (20\% reduction in growth) are supported [Default = "20"]
#' @param addBreakpoints Indicates whether to add breakpoint lines to the standard curve plot (if the user has supplied data to generate a standard curve)

#' @return In all cases the function will return an updated .csv file that contains the MIC values that correspond to calculated RAD values in the directory "parameter_files" in the main project directory. If the user has supplied their own MIC data the function will also save the calculated model parameters into a separate file and will plot the linear relationship and line of best fit.

#' @export

calcMIC <- function(projectName, type="df", RAD="20", height = 4, width = 6, addBreakpoints = TRUE, savePDF = TRUE, popUp = TRUE){
	ZOIvalue <- RAD
	if(.Platform$OS.type=="windows"){
	    drugFile <- file.path(.libPaths(), "diskImageR", "knownMIC-RAD.csv")[1]
		drugFile <- gsub("Program Files", "progra~1", drugFile)
		knownSppDrug <- read.csv(drugFile)
	}
	else{
		knownSppDrug <- read.csv(file.path(.libPaths(), "diskImageR", "knownMIC-RAD.csv"))
	}
	# knownSppDrug <- read.csv("knownMIC-RAD.csv")
	head(knownSppDrug[,1:2])	
	#Check whether file exists in environment or prompt user to load it
	if(type == "df"){
		if(paste(projectName, ".df", sep="") %in% ls(globalenv())) dataframe <- eval(parse(text=paste(projectName, ".df", sep="")))
		else stop(paste(projectName, ".df not found in working environment. Please load with function 'readExistingDF'", sep=""))
		}
	if(type == "ag"){
		if(paste(projectName, ".ag", sep="") %in% ls(globalenv())) dataframe <- eval(parse(text=paste(projectName, ".ag", sep="")))
		else stop(paste(projectName, ".ag not found in working environment. Please load with function 'readExistingAG'", sep=""))	
		}	
	useBuiltIn <- readline("Do you want to use built-in data for existing species/drug combinations? [y/n] ")
	if(useBuiltIn=="y"){
		print(knownSppDrug[,1:2])
		sppDrug <- readline("Please choose the number that corresponds to the species/drug combination you wish to use: ")
		curvePars <- c(knownSppDrug[sppDrug, 3], knownSppDrug[sppDrug, 4])
		cat("_____________________________________________________________________")
		cat(paste("\nintercept: ", curvePars[1], "\tslope: ", curvePars[2], sep=""))		
		cat(paste("\nReference: ", knownSppDrug[sppDrug,8], "\nThis was based on ", knownSppDrug[sppDrug, 7], " isolates. The log2(MIC), RAD relationship was ", knownSppDrug[sppDrug, 5], " with an R^2 value of ", knownSppDrug[sppDrug, 6], "\n", sep=""))		
		relat <- knownSppDrug[sppDrug, 5]
	}
	else{
		haveSI<-readline("Do you know the slope and intercept of log2(MIC) regressed on RAD? [y/n]  ")
		if(haveSI == "y"){
	   		intercept <- readline("Enter the intercept: ")
	   		intercept <- as.numeric(intercept)
	   		slope <- readline("Enter the slope: ")
	   		slope <- as.numeric(slope)
	   		curvePars <- c(intercept, slope)
	   		relation <- readline("Is this from a linear (L: RAD) or quadratic (Q: RAD^2) relationship? [L/Q] ")
	   		relat <- ifelse(relation == "L", "linear", "quadratic")
			cat(paste("\nintercept: ", curvePars[1], "\tslope: ", curvePars[2], "\n", sep=""))		
	   	}
		if(haveSI =="n"){
	  		exFile <- readline("Do you have a standard curve file that provides MIC and RAD data?  [y/n] ")   
	  		if(exFile == "n") stop("You must use data from built-in combinations, provide the required parameters or a standard curve file to proceed.")
		  	else{
 		 		sameData <- readline("Does your MIC data correspond to the same strains as your current project? [y/n] ")
				if(sameData == "y"){
					MICFile <- tcltk::tk_choose.files(caption = "Select the MIC standard curve file (containing MIC and RAD data in a text file, comma delimited)") 
					MICdata<- read.csv(MICFile, header=TRUE,sep=",") 
					if (ncol(MICdata)<2) stop("Wrong data format: the file must contain two columns, the first containing the line name, and the second wtih corresponding MIC values \n")
					else{
						MIC_names <- MICdata[,1]
						MIC<-0
						N1 <- dataframe[,1]
						for (i in 1:length(N1)){  
							ind <- grep(N1[i],MIC_names);
							if (length(ind)>0) MIC[i]<-mean(MICdata[ind,2])
							else MIC[i]<-NA
						}
						RAD <- subset(dataframe, names %in% MIC_names)[paste("RAD", RAD, sep="")]
						if(length(RAD) != length(MIC) & type=="df") stop(paste("The length of the MIC file does not match the length of ", projectName, ".df. If you have replicates for RAD please run 'aggregateData' and the rerun 'calcMIC' with 'type=\"ag\"'. If you have replicates for MIC please average in the standard curve file before proceeding", sep=""))
						if(length(RAD) != length(MIC) & type=="df") stop(paste("The length of the MIC file does not match the length of ", projectName, ".ag.", sep=""))						
					}
					fitL<-lm(log2(MIC)~RAD, na.action=na.exclude)
					fitQ<-lm(log2(MIC)~I(RAD^2), na.action=na.exclude)	
					if(summary(fitL)$adj.r.squared > summary(fitQ)$adj.r.squared){
						relat <- "linear"
						A <- summary(fitL)$coefficients[1]
						B <- summary(fitL)$coefficients[2]
						R2 <- summary(fitL)$adj.r.squared
						}
					else{
						relat <- "quadratic"
						A <- summary(fitQ)$coefficients[1]
						B <- summary(fitQ)$coefficients[2]
						R2 <- summary(fitQ)$adj.r.squared
						}
					curvePars <-c(A, B)
					cat(paste("\nThe best fit relationship is ", relat, " (R^2 = ", round(R2, 2), ") \nwith intercept: ", round(curvePars[1],2), "\tand slope: ", round(curvePars[2],2), "\n", sep=""))
				}
				if(sameData == "n"){
					MICFile <- tcltk::tk_choose.files(caption = "Select the MIC standard curve file (containing MIC and RAD data in a text file, comma delimited)") 		
					MICdata <- read.csv(MICFile, header=TRUE)
					if (ncol(MICdata)<2) stop("Wrong data format: the file must contain at least two columns, one containing RAD values ('RAD'), and one with corresponding MIC values ('MIC') \n")
					else{
						RAD <- MICdata$RAD
						MIC <- MICdata$MIC
						ddnv <- split(MICdata, MIC)
						meanRAD <- unlist(lapply(ddnv, function(x) mean(x$RAD)))
						meanMIC <- as.numeric(names(ddnv))
						fitL<-lm(log2(MIC)~RAD, na.action=na.exclude)
						fitQ<-lm(log2(MIC)~I(RAD^2), na.action=na.exclude)					
						if(summary(fitL)$adj.r.squared > summary(fitQ)$adj.r.squared){
							relat <- "linear"
							A <- summary(fitL)$coefficients[1]
							B <- summary(fitL)$coefficients[2]
							R2 <- summary(fitL)$adj.r.squared
							}
						else{
							relat <- "quadratic"
							A <- summary(fitQ)$coefficients[1]
							B <- summary(fitQ)$coefficients[2]
							R2 <- summary(fitQ)$adj.r.squared
							}
						curvePars <-c(A, B)
					cat(paste("\nThe best fit relationship is ", relat, " (R^2 = ", round(R2, 2), ") \nwith intercept: ", round(curvePars[1],2), "\tand slope: ", round(curvePars[2],2), "\n", sep=""))
						}
					}
				}			
		if(savePDF){
			t <- file.path(getwd(), "figures", projectName,  "RAD-MIC_standardCurve.pdf")					
			pdf(t, width=width, height=height)
			par(oma=c(1, 4, 1, 1))
			if(relat == "quadratic"){
				plot(RAD^2, log2(MIC), xlim=c(0, max(RAD^2)+20), yaxt="n", xaxt="n", xlab="", ylab="")
				abline(fitQ, col="red")
				mtext(expression(paste(RAD^2, " (mm)", sep="")), side=1, outer=FALSE, line=3) 		
				}
			else{
				plot(RAD, log2(MIC), xlim=c(0, max(RAD)+5), yaxt="n", xaxt="n", xlab="", ylab="")
				abline(fitL, col="red")
				mtext(expression(paste(RAD, " (mm)", sep="")), side=1, outer=FALSE, line=3) 		
				}
			axis(1, cex.axis=0.8)
			yVal <- seq(min(log2(c(unique(MIC)))), max(log2(c(unique(MIC)))), by=1)			
			axis(2, at=yVal, las=2, labels=2^yVal, cex.axis=0.8)
			if(addBreakpoints){
				abline(v=(19-6)^2/2, col="grey")
				abline(v=(14.5-6)^2/2, col="grey")
				abline(h=log2((32+64)/2), col="grey")
				abline(h=log2((8+16)/2), col="grey")
				}			
			title <- as.list(expression(paste(log[2], "(MIC)", sep=""), paste("actual values indicated, ", mu, "g/mL)")))
			mtext(do.call(expression, title), side=2, cex=0.8, line = c(3.5, 2.5))		
			dev.off()
			cat(paste("\nFigure saved: ", t, sep=""))
			if(popUp){
				tt <- paste("open ",t)
				system(tt)
			}
	}	
			paramName <- file.path(getwd(), "parameter_files", projectName, paste("RAD-MIC_parameters.csv", sep=""))
			params <- data.frame(intercept = curvePars[1], slope = curvePars[2])
			write.csv(params, file=paramName, row.names=FALSE)	
			cat(paste("\n\nThe calculated parameters have been saved here: \n", paramName, "\t and can be used in the future for the same species/drug combination\n", sep=""))
			}
		}

	if(relat == "quadratic")	MIC <- c(round(2^curvePars[1]*2^(curvePars[2]*dataframe[paste("RAD", ZOIvalue, sep="")]^2), 2))
	else MIC <- c(round(2^curvePars[1]*2^(curvePars[2]*dataframe[paste("RAD", ZOIvalue, sep="")]), 2))
	# The three lines below would round the MIC to the typical values that are tested:
	# typicalMIC <- c(0.03, 0.0625, 0.12, 0.25, 0.5, 1, 2, 4, 8, 16, 32, 64, 128)	
	# place <- apply(t(MIC), 1, function(x) which.min(typicalMIC - x))
	# MIC <- typicalMIC[place]
	# upDataframe <- data.frame(dataframe, MIC = roundMIC)
	upDataframe <- data.frame(dataframe, MIC = MIC)	
	dfName <- paste(projectName, ".df", sep="")	
	filename <- file.path(getwd(), "parameter_files", projectName, paste(projectName, "_df.csv", sep=""))
	write.csv(upDataframe, file=filename, row.names=FALSE)	
	cat(paste("\n", dfName, " has been updated and written to the global environment", sep=""))
	cat(paste("\n\nSaving file: ", filename, sep=""))
	# assign(dfName, upDataframe, envir=globalenv())
	assign(dfName, upDataframe, inherits=TRUE)
	}