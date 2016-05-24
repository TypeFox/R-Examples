#' Dataframe creation

#' @description \code{createDataframe} saves the calculated resistance, perseverence and sensitivity estimates

#' @inheritParams maxLik
#' @param nameVector either a logial value or a character vector. Supported values are \code{nameVector} = "TRUE" to assign the photograph name to the 'name' column, \code{nameVector} = "FALSE" to assign th photograph number to the 'name' column, or \code{nameVector} = a vector the same length as the number of photographs indicating the desired names.
#' @param typeVector a logical value. \code{typeVector} = "TRUE" will add a 'type' vector to the dataframe using values found in the \code{typePlace} position of the photograph names (see \code{\link{IJMacro}} for more details) while \code{typeVector} = "FALSE" will not add a type column.
#' @param typePlace a number that indicates the position of the photograph name to be stored as the 'type' vector'. Defaults to 2. For more details see \code{\link{IJMacro}}
#' @param typeName a character string that indicates what to name the typeVector. Defaults to "type".
#' @param removeClear a logical value that indicates whether to remove the clear halo picture from the dataset (i.e., is this picture an experimental picture, or one solely included to use as a clear halo). Defaults to FALSE.

#' @details A dataframe with 11 columns:
#' \itemize{
#' 		\item\bold{name:} determined by \code{nameVector}, either photograph names, photograph numbers, or a user-supplied list of names
#'	 	\item\bold{line:} the first components of the \code{namesVector}; everything that comes before the first "_" in the photograph name
#' 		\item\bold{type:} the location within the \code{name} of the phograph type is supplied by \code{typePlace}. Use \code{\link{addType}} if more than one type column are desired.
#' 		\item\bold{RAD80, RAD50, RAD20:} resistance parameters, coresponding to the distance in mm of 80\%, 50\% and 20\% reduction in growth
#' 		\item\bold{FoG80, FoG50, FoG20:} perseverence parameters, coresponding to the fraction of growth achieved above the 80\%, 50\% and 20\% reduction in growth points
#' 		\item\bold{slope:} sensitivity parameter, indicating the slope at the midpoint, i.e., how rapidly the population changes from low growth to full growth
#'	}
	
#' @return A dataframe "projectName.df" is saved to the global environment and a .csv file "projectName_df.csv" is exported to the "parameter_files" directory. 

#' @examples 
#' \dontrun{
#' createDataframe("myProject", clearHalo=1)
#' createDataframe("myProject", clearHalo=1, removeClear = TRUE, typeName = "drugAmt")
#' }

#' @export


createDataframe <- function(projectName, clearHalo, diskDiam = 6, maxDist = 30, standardLoc = 2.5, removeClear = FALSE, nameVector=TRUE, typeVector=TRUE, typePlace=2, typeName = "type"){
	if(!(hasArg(clearHalo))){
		cont <- readline(paste("Please specify photograph number with a clear halo ", sep=""))
		clearHalo <- as.numeric(cont)
	}
	data <- eval(parse(text=projectName))
	df <- data.frame()
	dotedge <- diskDiam/2 + 0.4
	newdir <- file.path(getwd(), "parameter_files")
	newdir2 <- file.path(getwd(), "parameter_files", projectName)
	newdir3 <- file.path(getwd(), "figures", projectName)
	
	filename <- file.path(getwd(), "parameter_files", projectName, paste(projectName, "_df.csv", sep=""))

	if (!file.exists(newdir)){		
		dir.create(newdir, showWarnings = FALSE)
		cat(paste("\n\tCreating new directory: ", newdir), sep="")
		}
	if (!file.exists(newdir2)){		
		dir.create(newdir2, showWarnings = FALSE)
		cat(paste("\nCreating new directory: ", newdir2), sep="")
		}
	if (!file.exists(newdir3)){		
		dir.create(newdir3, showWarnings = FALSE)
		cat(paste("\nCreating new directory: ", newdir3), sep="")
		}
	df <- data.frame(row.names = seq(1, length(data)))

	ML <- paste(projectName, ".ML", sep="")	
	ML2 <- paste(projectName, ".ML2", sep="")
	ML <- eval(parse(text=ML))
	ML2 <- eval(parse(text=ML2))	

	dotMax <- max(sapply(data, function(x) {x[which(x[,1] > standardLoc)[1], 2]})) 		
	stand <-c( sapply(data, function(x) {dotMax-x[which(x[,1] > standardLoc)[1], 2]}))
	clearHaloData <- data[[clearHalo]]
	startX <- which(clearHaloData[,1] > dotedge+0.5)[1]
	stopX <- which(clearHaloData[,1] > maxDist - 0.5)[1]
	clearHaloData <- clearHaloData[startX:stopX, 1:2]
	clearHaloData$x <- clearHaloData$x + stand[clearHalo] 
	clearHaloData$distance <- clearHaloData$distance - (dotedge+0.5)
	clearHaloStand <- clearHaloData[1,2]

	slope <- sapply(c(1:length(data)), .findSlope, data=data, ML=ML, stand = stand, dotedge = dotedge, maxDist = maxDist, clearHaloStand = clearHaloStand)

	FoG.df <-  sapply(c(1:length(data)), .findFoG, data=data, ML=ML, ML2 = ML2, stand = stand, dotedge = dotedge,  maxDist = maxDist, clearHaloStand = clearHaloStand, standardLoc = standardLoc)

	x80 <- unlist(FoG.df[1,])	
	x50 <- unlist(FoG.df[2,])	
	x20 <- unlist(FoG.df[3,])
	FoG80 <- unlist(FoG.df[4,])	
	FoG50 <- unlist(FoG.df[5,])	
	FoG20 <- unlist(FoG.df[6,])	
	maxFoG <- unlist(FoG.df[7,])		
	maxFoG80 <- unlist(FoG.df[8,])		
	maxFoG50 <- unlist(FoG.df[9,])		
	maxFoG20 <- unlist(FoG.df[10,])				
	
	FoG80[slope < 5] <- NA
	FoG50[slope < 5] <- NA
	FoG20[slope < 5] <- NA
	x80[slope < 5] <- 1
	x50[slope < 5] <- 1
	x20[slope < 5] <- 1

	aveFoG80 <- FoG80/x80
	aveFoG50 <- FoG50/x50
	aveFoG20 <- FoG20/x20	

	param <- data.frame(RAD80 =round(x80, digits=0), RAD50 = round(x50, digits=0), RAD20 = round(x20, digits=0), FoG80 = round(FoG80/maxFoG80, digits=2), FoG50 = round(FoG50/maxFoG50, digits=2), FoG20 = round(FoG20/maxFoG20, digits=2), slope=round(slope, digits=1))
	
	if (is.logical(nameVector)){
		if (nameVector){
			line <- unlist(lapply(names(data), function(x) strsplit(x, "_")[[1]][1]))
			df <- data.frame(name = names(data), line)
			}
			
		if (!nameVector){
			line <- seq(1, length(data))
			df <- data.frame(name = names(data), line, df)	
		}
	}
	if (!is.logical(nameVector)){
		line <- nameVector
		names <- unlist(lapply(names(data), function(x) strsplit(x, "_")[[1]][1]))
		df <- data.frame(names=names, line=line, df)	
		}

	if (typeVector){	
			type <- unlist(lapply(names(data), function(x) strsplit(x, "_")[[1]][typePlace]))
			df <- data.frame(df, type, param)
		}
	else {
			df$type <- 1
			df <- data.frame(df, param)
		}
		
	names(df)[3] <- typeName

	df <- df[order(df$line),] 
	df$FoG80[df$FoG80 >1] <- 1
	df$FoG50[df$FoG50 >1] <- 1
	df$FoG20[df$FoG20 >1] <- 1	
	df$FoG80[df$RAD80 == 0] <- NA
	df$FoG50[df$RAD50 == 0] <- NA
	df$FoG20[df$RAD20 == 0] <- NA
	df$FoG80[df$RAD80 == 1] <- NA
	df$FoG50[df$RAD50 == 1] <- NA
	df$FoG20[df$RAD20 == 1] <- NA

	if (removeClear)	df <- df[-clearHalo,]
	 	
	write.csv(df, file=filename, row.names=FALSE)	
	
	dfName <- paste(projectName, ".df", sep="")
	cat(paste("\n", dfName, " has been written to the global environment", sep=""))
	cat(paste("\nSaving file: ", filename,  sep=""))
	cat(paste("\n", projectName, "_df.csv can be opened in MS Excel.",  sep=""))
	# assign(dfName, df, envir=globalenv())
	assign(dfName, df, inherits=TRUE)
	}

#Determine the slope

.findSlope <- function(data, ML, i, stand, clearHaloStand, dotedge = 3.4,  maxDist = 35){
	startX <- which(data[[i]][,1] > dotedge+0.5)[1]
	stopX <- which(data[[i]][,1] > maxDist - 0.5)[1]
	data[[i]] <- data[[i]][startX:stopX, 1:2]
	data[[i]]$x <- data[[i]]$x + stand[i] - clearHaloStand 
	data[[i]]$distance <- data[[i]]$distance - (dotedge+0.5)
	xx <- seq(log(data[[i]]$distance[1]), log(max(data[[i]][,1])), length=200)
	yy<- .curve(ML[[i]]['par'][1]$par[1], ML[[i]]['par'][1]$par[2], ML[[i]]['par'][1]$par[3],xx)
	yycor <- (yy+min(data[[i]]$x))
	xcross <- exp(ML[[i]]['par'][1]$par[2])
	xxmid <- which.max(exp(xx) > xcross)
	if ((xxmid-10) > 1){
		xxSlope <- xx[(xxmid-10):(xxmid+10)]
		yySlope <- yy[(xxmid-10):(xxmid+10)]
		}
	else {
		xxSlope <- xx[1:(xxmid+10)]
		yySlope <- yy[1:(xxmid+10)]
	}
	slope <- lm(yySlope ~ xxSlope)$coefficients[2]
	return(slope)
}

.findFoG <- function(data, ML, ML2, stand, clearHaloStand, dotedge = 3.4, maxDist = 35, standardLoc = 2.5, i){	
	startX <- which(data[[i]][,1] > dotedge+0.5)[1]
	stopX <- which(data[[i]][,1] > maxDist - 0.5)[1]
	data[[i]] <- data[[i]][startX:stopX, 1:2]
	data[[i]]$x <- data[[i]]$x + stand[i] - clearHaloStand 
	data[[i]]$distance <- data[[i]]$distance - (dotedge+0.5)
	xx <- seq(log(data[[i]]$distance[1]), log(max(data[[i]][,1])), length=200) 
	yy<- .curve2(ML2[[i]]$par[1], ML2[[i]]$par[2], ML2[[i]]$par[3], ML2[[i]]$par[5], ML2[[i]]$par[6], ML2[[i]]$par[7], xx) 
	ploty <- data[[i]]$x
	ploty[ploty < 0] <-0
	asym <- (ML[[i]]$par[1]+min(data[[i]]$x))

	xx <- seq(log(data[[i]]$distance[1]), log(max(data[[i]][,1])), length=200) 				
	yy <- (yy+min(data[[i]]$x))
	yy[yy < 0] <- 0		
	x80 <- xx[which.max(yy> asym * 0.2)]
	x50 <- xx[which.max(yy> asym * 0.5)]
	x20 <- xx[which.max(yy> asym * 0.8)]
	if (x20 < x50) x20 <- xx[which.max(yy> yy[length(yy)] * 0.8)]

	if(exp(x80)>1) xx80 <- seq(log(data[[i]]$distance[1]), log(round(exp(x80))), length=200)
	else xx80 <- seq(log(data[[i]]$distance[1]), log(data[[i]]$distance[2]), length=200)

	if(exp(x50)>1) xx50 <- seq(log(data[[i]]$distance[1]), log(round(exp(x50))), length=200)	
	else xx50 <- seq(log(data[[i]]$distance[1]), log(data[[i]]$distance[2]), length=200)

	if(exp(x20)>1) xx20 <- seq(log(data[[i]]$distance[1]), log(round(exp(x20))), length=200)
	else xx20 <- seq(log(data[[i]]$distance[1]), log(data[[i]]$distance[2]), length=200)

	yy <- .curve2(ML2[[i]]$par[1], ML2[[i]]$par[2], ML2[[i]]$par[3], ML2[[i]]$par[5], ML2[[i]]$par[6], ML2[[i]]$par[7], xx)		
	yy80 <- .curve2(ML2[[i]]$par[1], ML2[[i]]$par[2], ML2[[i]]$par[3], ML2[[i]]$par[5], ML2[[i]]$par[6], ML2[[i]]$par[7], xx80)	
	yy50<- .curve2(ML2[[i]]$par[1], ML2[[i]]$par[2], ML2[[i]]$par[3], ML2[[i]]$par[5], ML2[[i]]$par[6], ML2[[i]]$par[7], xx50)
	yy20<- .curve2(ML2[[i]]$par[1], ML2[[i]]$par[2], ML2[[i]]$par[3], ML2[[i]]$par[5], ML2[[i]]$par[6], ML2[[i]]$par[7], xx20) 
		
	yy <- (yy+min(data[[i]]$x))	
	yy[yy < 0] <- 0.1
	yy80 <- (yy80+min(data[[i]]$x))
	yy80[yy80 < 0] <- 0.1
	yy50 <- (yy50+min(data[[i]]$x))
	yy50[yy50 < 0] <- 0.1
	yy20 <- (yy20+min(data[[i]]$x))
	yy20[yy20 < 0] <- 0.1
			
	id <- order(xx)
	id80 <- order(xx80)
	id50 <- order(xx50)
	id20 <- order(xx20)		

	maxFoG <- sum(diff(xx[id])*zoo::rollmean(yy[id], 2))	
	maxFoG80 <- exp(x80)*max(yy80)
	maxFoG50 <- exp(x50)*max(yy50)	
	maxFoG20 <- exp(x20)*max(yy20)	
		
	FoG80 <- sum(diff(exp(xx80[id80]))*zoo::rollmean(yy80[id80], 2))		
	FoG50 <- sum(diff(exp(xx50[id50]))*zoo::rollmean(yy50[id50], 2))		
	FoG20 <- sum(diff(exp(xx20[id20]))*zoo::rollmean(yy20[id20], 2))		
	
	 param <- data.frame(x80 = round(exp(x80), digits=0), x50 = round(exp(x50), digits=2), x20 = round(exp(x20), digits=0) , FoG80 = round(FoG80, digits=0), FoG50= round(FoG50, digits=0), FoG20= round(FoG20, digits=0), maxFoG = round(maxFoG, digits=0), maxFoG80 = round(maxFoG80, digits=0), maxFoG50 = round(maxFoG50, digits=0), maxFoG20 = round(maxFoG20, digits=0))	
	 
	 if (exp(param$x80)<1) 	param$x80 <- 1
	 if (exp(param$x50)<1)	param$x50 <- 1	 
	 if (exp(param$x20)<1)	param$x20 <- 1	  
	 return(param)	
	}
