

## ------------------------------------------------------------------------ ##
##  Script R for  "Constructing Entity Specific Prospective Mortality Table ##
##                 Adjustment to a reference"                               ##
## ------------------------------------------------------------------------ ##
##  Script        S12.R                                                     ##
## ------------------------------------------------------------------------ ##
##  Description   Comparison between the methods                            ##
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

## ------------------------------------------------------------------------ ##
##  Comparison Level 1                                                      ##
## ------------------------------------------------------------------------ ##

.CompLevel1 = function(ListValidationLevel1){
	print("Comparison 1st Level of validation")
	ColNames <- vector(,length(ListValidationLevel1))
	for(i in 1:length(ListValidationLevel1)){
		ColNames[i] <- ListValidationLevel1[[i]][[1]][[1]]$NameMethod
		}
	ValLevel1 <- vector("list", length(ListValidationLevel1[[1]][[1]]))
	names(ValLevel1) <- names(ListValidationLevel1[[1]][[1]])
	for(i in 1:length(ListValidationLevel1[[1]][[1]])){
		ValLevel1[[i]] <- vector("list",5)
		names(ValLevel1[[i]]) <- names(ListValidationLevel1[[1]][[1]][[1]])[1:length(ValLevel1[[i]])]
		for (j in 1:length(ValLevel1[[i]])){
			RowNames <- rownames(ListValidationLevel1[[1]][[1]][[i]][[j]])
			ValLevel1[[i]][[j]] <- matrix(,length(RowNames), length(ColNames))
			rownames(ValLevel1[[i]][[j]]) <- RowNames
			colnames(ValLevel1[[i]][[j]]) <- ColNames
			for(k in 1:length(ListValidationLevel1)){
				ValLevel1[[i]][[j]][,k] <- ListValidationLevel1[[k]][[1]][[i]][[j]]
				}
			}
		}
	print(ValLevel1)
	return(ValLevel1)
	}
	
## ------------------------------------------------------------------------ ##
##  Comparison Level 2                                                      ##
## ------------------------------------------------------------------------ ##

.CompLevel2 = function(ListValidationLevel2){
	print("Comparison 2nd Level of validation")
	ColNames <- vector(,length(ListValidationLevel2))
	for(i in 1:length(ListValidationLevel2)){
		ColNames[i] <- ListValidationLevel2[[i]][[1]]$NameMethod
		}
	ValLevel2 <- vector("list", length(ListValidationLevel2[[1]]))
	names(ValLevel2) <- names(ListValidationLevel2[[1]])
	for(i in 1:length(ListValidationLevel2[[1]])){
		ValLevel2[[i]] <- vector("list",2)
		names(ValLevel2[[i]]) <- names(ListValidationLevel2[[1]][[1]])[1:length(ValLevel2[[i]])]
		for (j in 1:length(ValLevel2[[i]])){
			RowNames <- rownames(ListValidationLevel2[[1]][[i]][[j]])
			ValLevel2[[i]][[j]] <- matrix(,length(RowNames), length(ColNames))
			rownames(ValLevel2[[i]][[j]]) <- RowNames
			colnames(ValLevel2[[i]][[j]]) <- ColNames
			for(k in 1:length(ListValidationLevel2)){
				ValLevel2[[i]][[j]][,k] <- ListValidationLevel2[[k]][[i]][[j]]
				}
			}
		}
	print(ValLevel2)
	return(ValLevel2)
	}

## ------------------------------------------------------------------------ ##
##  Comparison Level 3                                                      ##
## ------------------------------------------------------------------------ ##

.CompLevel3 = function(ListValidationLevel3){
	print("Comparison 3rd Level of validation")
	ColNames <- vector(,length(ListValidationLevel3))
	for(i in 1:length(ListValidationLevel3)){
		ColNames[i] <- ListValidationLevel3[[i]][[1]][[1]]$NameMethod
		}
	ValLevel3 <- vector("list", length(ListValidationLevel3[[1]][[1]]))
	names(ValLevel3) <- names(ListValidationLevel3[[1]][[1]])
	for(i in 1:length(ListValidationLevel3[[1]][[1]])){
		ValLevel3[[i]] <- vector("list",3)
		names(ValLevel3[[i]]) <- names(ListValidationLevel3[[1]][[1]][[1]])[1:length(ValLevel3[[i]])]
		for (j in 1:length(ValLevel3[[i]])){
			RowNames <- rownames(ListValidationLevel3[[1]][[1]][[i]][[j]])
			ValLevel3[[i]][[j]] <- matrix(,length(RowNames), length(ColNames))
			rownames(ValLevel3[[i]][[j]]) <- RowNames
			colnames(ValLevel3[[i]][[j]]) <- ColNames
			for(k in 1:length(ListValidationLevel3)){
				ValLevel3[[i]][[j]][,k] <- ListValidationLevel3[[k]][[1]][[i]][[j]]
				}
			}
		}
	print(ValLevel3)
	return(ValLevel3)
	}

## ------------------------------------------------------------------------ ##
##  Comparison of the methods                                               ##
## ------------------------------------------------------------------------ ##



ComparisonMethods = function(ListOutputs, ListValidationLevel1, ListValidationLevel2, ListValidationLevel3, MyData = MyData, Plot = F, ColorComp = c("#FF6590","#309BFF","#AD79FC","#3CAB5F"), LtyComp=rep(1,4), AgeCrit){
	print("Comparisons of the output produced by the 3 levels of validation")
	CompLevel <- vector("list", 3)
	names(CompLevel) <- c("Level 1", "Level 2", "Level 3")
	CompLevel[[1]] <- .CompLevel1(ListValidationLevel1)
	CompLevel[[2]] <- .CompLevel2(ListValidationLevel2)
	CompLevel[[3]] <- .CompLevel3(ListValidationLevel3)
	if(Plot == T){
		Path <- "Results/Graphics/Comparison"
		.CreateDirectory(paste("/",Path,sep=""))
		print(paste("Create the plot comparing the fits in .../",Path," ...", sep=""))
		for (i in 1:(length(MyData)-1)){
			for(j in MyData[[i]]$YearCom){
				png(filename=paste(Path,"/FitsLog-",j,"-",names(MyData)[i],".png",sep=""), width  = 2100, height = 2100, res=300, pointsize= 12)
				.ComparisonFitsMethodsLog(ListOutputs, j, MyData[[i]], AgeCrit, names(MyData)[i], ColorComp, LtyComp)
				dev.off()
				png(filename=paste(Path,"/Fits-",j,"-",names(MyData)[i],".png",sep=""), width  = 2100, height = 2100, res=300, pointsize= 12)
				.ComparisonFitsMethods(ListOutputs, j, MyData[[i]], AgeCrit, names(MyData)[i], ColorComp, LtyComp)
				dev.off()
				}	
			}
		print(paste("Create the plot comparing the residuals in .../", Path," ...", sep=""))
		for (i in 1:(length(MyData)-1)){
			for(j in MyData[[i]]$YearCom){
				png(filename=paste(Path,"/Residuals-",j,"-",names(MyData)[i],".png",sep=""), width  = 3400, height = 1000 * length(ListValidationLevel1), res=300, pointsize= 12)
				.ComparisonResidualsMethods(ListValidationLevel1, j, AgeCrit, names(MyData)[i], ColorComp)
				dev.off()
				}
			}			
		print(paste("Create the plot comparing the fitted deaths in .../", Path," ...", sep=""))
		for (i in 1:(length(MyData)-1)){
			for(j in MyData[[i]]$YearCom){	
				png(filename=paste(Path,"/FittedDeaths-",j,"-",names(MyData)[i],".png",sep=""), width  = 2100, height = 2100, res=300, pointsize= 12)
				.ComparisonFittedDeathsMethods(ListValidationLevel1, j, MyData[[i]], AgeCrit, names(MyData)[i], ColorComp, LtyComp)
				dev.off()
				}
			}
		AgeCoh <- c(55, 70, 85)
		AgeComp <- as.numeric(rownames(ListValidationLevel3[[1]][["PerLifeExpFitted"]][[i]]$PerLifeExp))
		print(paste("Create the plot comparing the coherence of the fitted periodic life expectancies (until ", max(AgeComp) ," years) for age ", AgeCoh[AgeCoh <= max(AgeComp) & AgeCoh >= min(AgeComp)]," in .../", Path," ...", sep=""))
		for(k in 1:length(AgeCoh)){
			if(AgeCoh[k] <= max(AgeComp) & AgeCoh[k] >= min(AgeComp)){
				for (i in 1 : (length(MyData)-1)){
					png(filename=paste(Path,"/CoherencePerLifeExp-", AgeCoh[k], "-", names(MyData)[i],".png", sep=""), width  = 2100, height = 2100, res=300, pointsize= 12)		
					.ComparisonTrendsMethods(ListValidationLevel3, AgeCoh[k], names(MyData)[i], ColorComp, LtyComp)
					dev.off()
					}
				}
			}
		}
	return(CompLevel)
	}
