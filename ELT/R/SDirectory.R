

## ------------------------------------------------------------------------ ##
##  Script R for  "Constructing Entity Specific Prospective Mortality Table ##
##                 Adjustment to a reference"                               ##
## ------------------------------------------------------------------------ ##
##  Script        SDirectory.R                                              ##
## ------------------------------------------------------------------------ ##
##  Description   Functions to manage the directories                       ##
## ------------------------------------------------------------------------ ##
##  Authors       Tomas Julien, Frederic Planchet and Wassim Youssef        ##
##                julien.tomas@univ-lyon1.fr                                ##
##                frederic.planchet@univ-lyon1.fr                           ##
##                wassim.g.youssef@gmail.com                                ##
## ------------------------------------------------------------------------ ##
##  Version       01 - 2013/11/06                                           ##
## ------------------------------------------------------------------------ ##

## ------------------------------------------------------------------------ ##
##  Definition of the functions                                             ##
## ------------------------------------------------------------------------ ##


.WarningInvalidAge = function(d, l, x1, x2, t1){
		if(min(x1)<min(x2)){
				stop("The minimum of the age range selected is lower than the age range of the reference. Please, select a valid age range.", call. = F)
			}
		if(max(x1)>max(x2)){
				stop("The maximum of the age range selected is greater than the age range of the reference. Please, select a valid age range.", call. = F)
			}
		xx <- as.numeric(rownames(d))
		if (length(d[d[x1 - min(xx) + 1, ] > l[x1 - min(xx) + 1, ]]) > 0){
				SelMat <- which(d[x1 - min(xx) + 1, ] > l[x1 - min(xx) + 1, ], arr.ind = T)
				StopMat <- matrix(, nrow(SelMat), ncol(SelMat))
				StopMat[, 1] <- SelMat[, 1] + min(x1) - 1
				StopMat[, 2] <- t1[SelMat[, 2]]
				colnames(StopMat) <- c("Age", "Year")
				print("The observed deaths are larger than the observed number of individuals for the following year(s) and age(s):")
				print(StopMat)
				stop("Please, select a valid age range.", call. = F)
			}
	}

## ------------------------------------------------------------------------ ##
##  .CreateDirectory function                                                ##
## ------------------------------------------------------------------------ ##

.CreateDirectory = function(Path){
	CurrentWorkplace <- getwd()
	MyPath <- paste(CurrentWorkplace, Path , sep="")
	if(!file.exists(paste(".",Path,sep=""))){
		print(paste("Create the directory", MyPath))
		dir.create(MyPath, recursive = T)
	}
}
