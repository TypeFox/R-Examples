

## ------------------------------------------------------------------------ ##
##  Script R for  "Constructing Entity Specific Prospective Mortality Table ##
##                 Adjustment to a reference"                               ##
## ------------------------------------------------------------------------ ##
##  Script         SExcel.R                                                 ##
## ------------------------------------------------------------------------ ##
##  Description   Functions to store results in Excel files                 ##
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
##  ExportObsStatInExcel function                                           ##
## ------------------------------------------------------------------------ ##

.ExportHistoryInExcel = function(MyData){
	Path <- "Results/Excel"
	.CreateDirectory(paste("/",Path,sep=""))
	ExcelPath <- paste(Path, "/History.xlsx", sep = "")
	print(paste("Export the history in .../", ExcelPath, " ...", sep=""))
	if(!file.exists(paste("./",ExcelPath,sep=""))){
		WorkBook <- createWorkbook()
		SheetExpo <- createSheet(WorkBook, sheetName=paste("Exposition", names(MyData)[1]))
		saveWorkbook(WorkBook, file = ExcelPath)
		}
	if(file.exists(paste("./",ExcelPath,sep=""))){
		WorkBook <- loadWorkbook(file = ExcelPath)
		}
	for (i in 1 : (length(MyData)-1)){	
		Deaths <- cbind(0:120, MyData[[i]]$Deaths)
		Expo <- cbind(0:120, MyData[[i]]$Expo)
		Indi <- cbind(0:120, MyData[[i]]$Indi)
		colnames(Deaths) <- colnames(Expo) <- colnames(Indi) <- c("Age", colnames(MyData[[i]]$Deaths))
		NamesSheets <- names(getSheets(WorkBook))
		if(length(NamesSheets[NamesSheets==paste("Exposition", names(MyData)[i])])==1){
			SheetExpo <- getSheets(WorkBook)[[paste("Exposition", names(MyData)[i])]]
			}
		if(length(NamesSheets[NamesSheets==paste("Exposition", names(MyData)[i])])==0){
			SheetExpo <- createSheet(WorkBook, sheetName=paste("Exposition", names(MyData)[i]))
			}
		addDataFrame(Expo, SheetExpo, row.names=F, col.names=T, startRow=1)
		if(length(NamesSheets[NamesSheets==paste("Individuals", names(MyData)[i])])==1){
			SheetIndi <- getSheets(WorkBook)[[paste("Individuals", names(MyData)[i])]]
			}
		if(length(NamesSheets[NamesSheets==paste("Individuals", names(MyData)[i])])==0){
			SheetIndi <- createSheet(WorkBook, sheetName=paste("Individuals", names(MyData)[i]))
			}
		addDataFrame(Indi, SheetIndi, row.names=F, col.names=T, startRow=1)
		if(length(NamesSheets[NamesSheets==paste("Deaths", names(MyData)[i])])==1){
			SheetDeaths <- getSheets(WorkBook)[[paste("Deaths", names(MyData)[i])]]
			}
		if(length(NamesSheets[NamesSheets==paste("Deaths", names(MyData)[i])])==0){
			SheetDeaths <- createSheet(WorkBook, sheetName=paste("Deaths", names(MyData)[i]))
			}
		addDataFrame(Deaths, SheetDeaths, row.names=F, col.names=T, startRow=1)	
		}
	saveWorkbook(WorkBook, file = ExcelPath)	
	}

## ------------------------------------------------------------------------ ##
##  .ExportFinalTablesInExcel function                                       ##
## ------------------------------------------------------------------------ ##

.ExportFinalTablesInExcel = function(ModCompletion, MyData, AgeMethod){
	Path <- "Results/Excel"
	.CreateDirectory(paste("/",Path,sep=""))
	ExcelPath <- paste(Path, "/FinalTables.xlsx", sep = "")
	print(paste("Export the completed tables in .../", ExcelPath, " ...", sep=""))
		if(!file.exists(paste("./",ExcelPath,sep=""))){
		WorkBook <- createWorkbook()
		SheetExpo <- createSheet(WorkBook, sheetName=paste(ModCompletion[[1]]$NameMethod, names(MyData)[1]))
		saveWorkbook(WorkBook, file = ExcelPath)
		}
	if(file.exists(paste("./",ExcelPath,sep=""))){
		WorkBook <- loadWorkbook(file = ExcelPath)
		}
	for (i in 1 : (length(MyData)-1)){	
		Final <- cbind(min(AgeMethod):130, ModCompletion[[i]]$QxtFinal)
#		Final <- as.data.frame(Final, col.names = c("Age", MyData[[i]]$YearRef))
		Final <- as.data.frame(Final)
		colnames(Final) <- c("Age", colnames(ModCompletion[[i]]$QxtFinal))
		NamesSheets <- names(getSheets(WorkBook))
		if(length(NamesSheets[NamesSheets==paste(ModCompletion[[i]]$NameMethod, names(MyData)[i])])==1){
			SheetFinal <- getSheets(WorkBook)[[paste(ModCompletion[[i]]$NameMethod, names(MyData)[i])]]
			}
		if(length(NamesSheets[NamesSheets==paste(ModCompletion[[i]]$NameMethod, names(MyData)[i])])==0){
			SheetFinal <- createSheet(WorkBook, sheetName=paste(ModCompletion[[i]]$NameMethod, names(MyData)[i]))
			}
		addDataFrame(Final, SheetFinal, row.names=F, col.names=T, startRow=1)
		}
	saveWorkbook(WorkBook, file = ExcelPath)	
	}	

## ------------------------------------------------------------------------ ##
##  .ExportValidationL1inExcel function                                      ##
## ------------------------------------------------------------------------ ##

.ExportValidationL1inExcel = function(CritLevel1){
	if(CritLevel1[[1]]$NameMethod == "Method1"){NberCol <- 0}
	if(CritLevel1[[1]]$NameMethod == "Method2"){NberCol <- 1}
	if(CritLevel1[[1]]$NameMethod == "Method3"){NberCol <- 2}
	if(CritLevel1[[1]]$NameMethod == "Method4"){NberCol <- 3}
	Path <- "Results/Excel"
	.CreateDirectory(paste("/",Path,sep=""))
	ExcelPath <- paste(Path, "/Validation.xlsx", sep = "")
	print(paste("Export of the results of the validation first level in .../", ExcelPath, " ...", sep=""))
	if(!file.exists(paste("./",ExcelPath,sep=""))){
		WorkBook <- createWorkbook()
		ValLevel1 <- createSheet(WorkBook, sheetName=paste("Validation Level 1", names(CritLevel1)[1]))
		saveWorkbook(WorkBook, file = ExcelPath)
		}
	if(file.exists(paste("./",ExcelPath,sep=""))){
		WorkBook <- loadWorkbook(file = ExcelPath)
		}
	for (i in 1 : length(CritLevel1)){	
		NamesSheets <- names(getSheets(WorkBook))
		if(length(NamesSheets[NamesSheets==paste("Validation Level 1", names(CritLevel1)[i])])==1){
			ValLevel1 <- getSheets(WorkBook)[[paste("Validation Level 1", names(CritLevel1)[i])]]
			}
		if(length(NamesSheets[NamesSheets==paste("Validation Level 1", names(CritLevel1)[i])])==0){
			ValLevel1 <- createSheet(WorkBook, sheetName=paste("Validation Level 1", names(CritLevel1)[i]))
			}
		addDataFrame(names(CritLevel1[[i]])[1], ValLevel1, row.names=F, col.names=F)
		addDataFrame(rownames(CritLevel1[[i]][1][[1]]), ValLevel1, row.names=F, col.names=F, startRow=4)
		addDataFrame(CritLevel1[[i]][1], ValLevel1, row.names=F, col.names=T, startRow=3, startColumn = 2 + NberCol)
		addDataFrame(names(CritLevel1[[i]])[2], ValLevel1, row.names=F, col.names=F, startRow=9)
		addDataFrame(rownames(CritLevel1[[i]][2][[1]]), ValLevel1, row.names=F, col.names=F,startRow=12)
		addDataFrame(CritLevel1[[i]][2], ValLevel1, row.names=F, col.names=T,startRow=11, startColumn = 2 + NberCol)
		addDataFrame(names(CritLevel1[[i]])[3], ValLevel1, row.names=F, col.names=F, startRow=18)
		addDataFrame(rownames(CritLevel1[[i]][3][[1]]), ValLevel1, row.names=F, col.names=F,startRow=21)
		addDataFrame(CritLevel1[[i]][3], ValLevel1, row.names=F, col.names=T,startRow=20, startColumn = 2 + NberCol)
		addDataFrame(names(CritLevel1[[i]])[4], ValLevel1, row.names=F, col.names=F, startRow=27)
		addDataFrame(rownames(CritLevel1[[i]][4][[1]]), ValLevel1, row.names=F, col.names=F,startRow=30)
		addDataFrame(CritLevel1[[i]][4], ValLevel1, row.names=F, col.names=T,startRow=29, startColumn = 2 + NberCol)
		addDataFrame(names(CritLevel1[[i]])[5], ValLevel1, row.names=F, col.names=F, startRow=33)
		addDataFrame(rownames(CritLevel1[[i]][5][[1]]), ValLevel1, row.names=F, col.names=F, startRow=36)
		addDataFrame(CritLevel1[[i]][5], ValLevel1, row.names=F, col.names=T, startRow=35, startColumn = 2 + NberCol)
		}
	saveWorkbook(WorkBook, file = ExcelPath)	
	}

## ------------------------------------------------------------------------ ##
##  .ExportValidationL2inExcel function                                      ##
## ------------------------------------------------------------------------ ##

.ExportValidationL2inExcel = function(CritLevel2){
	if(CritLevel2[[1]]$NameMethod == "Method1"){NberCol <- 0}
	if(CritLevel2[[1]]$NameMethod == "Method2"){NberCol <- 1}
	if(CritLevel2[[1]]$NameMethod == "Method3"){NberCol <- 2}
	if(CritLevel2[[1]]$NameMethod == "Method4"){NberCol <- 3}
	Path <- "Results/Excel"
	.CreateDirectory(paste("/",Path,sep=""))
	ExcelPath <- paste(Path, "/Validation.xlsx", sep = "")
	print(paste("Export of the results of the validation first level in .../", ExcelPath, " ...", sep=""))
	if(!file.exists(paste("./",ExcelPath,sep=""))){
		WorkBook <- createWorkbook()
		ValLevel2 <- createSheet(WorkBook, sheetName=paste("Validation Level 2", names(CritLevel2)[1]))
		saveWorkbook(WorkBook, file = ExcelPath)
		}
	if(file.exists(paste("./",ExcelPath,sep=""))){
		WorkBook <- loadWorkbook(file = ExcelPath)
		}
	for (i in 1 : length(CritLevel2)){	
		NamesSheets <- names(getSheets(WorkBook))
		if(length(NamesSheets[NamesSheets==paste("Validation Level 2", names(CritLevel2)[i])])==1){
			ValLevel2 <- getSheets(WorkBook)[[paste("Validation Level 2", names(CritLevel2)[i])]]
			}
		if(length(NamesSheets[NamesSheets==paste("Validation Level 2", names(CritLevel2)[i])])==0){
			ValLevel2 <- createSheet(WorkBook, sheetName=paste("Validation Level 2", names(CritLevel2)[i]))
			}
		addDataFrame(names(CritLevel2[[i]])[1], ValLevel2, row.names=F, col.names=F)
		addDataFrame(rownames(CritLevel2[[i]][1][[1]]), ValLevel2, row.names=F, col.names=F, startRow=4)
		addDataFrame(CritLevel2[[i]][1], ValLevel2, row.names=F, col.names=T, startRow=3, startColumn = 2 + NberCol)
		addDataFrame(names(CritLevel2[[i]])[2], ValLevel2, row.names=F, col.names=F, startRow=13)
		addDataFrame(rownames(CritLevel2[[i]][2][[1]]), ValLevel2, row.names=F, col.names=F,startRow=16)
		addDataFrame(CritLevel2[[i]][2], ValLevel2, row.names=F, col.names=T,startRow=15, startColumn = 2 + NberCol)
		}
	saveWorkbook(WorkBook, file = ExcelPath)	
	}

## ------------------------------------------------------------------------ ##
##  .ExportSingleIndiciesinExcel function                                    ##
## ------------------------------------------------------------------------ ##

.ExportSingleIndiciesinExcel = function(SI, MyData){
	if(SI[[1]]$NameMethod == "Method1"){NberCol <- 0}
	if(SI[[1]]$NameMethod == "Method2"){NberCol <- 1}
	if(SI[[1]]$NameMethod == "Method3"){NberCol <- 2}
	if(SI[[1]]$NameMethod == "Method4"){NberCol <- 3}
	Path <- "Results/Excel"
	.CreateDirectory(paste("/",Path,sep=""))
	ExcelPath <- paste(Path, "/Validation.xlsx", sep = "")
	print(paste("Export of the singles indices summarizing the lifetime probability distribution in .../", ExcelPath, " ...", sep=""))
	if(!file.exists(paste("./",ExcelPath,sep=""))){
		WorkBook <- createWorkbook()
		ValLevel3 <- createSheet(WorkBook, sheetName=paste("Validation Level 3", names(MyData)[1]))
		saveWorkbook(WorkBook, file = ExcelPath)
		}
	if(file.exists(paste("./",ExcelPath,sep=""))){
		WorkBook <- loadWorkbook(file = ExcelPath)
		}
	for (i in 1 : (length(MyData)-1)){	
		NamesSheets <- names(getSheets(WorkBook))
		if(length(NamesSheets[NamesSheets==paste("Validation Level 3", names(MyData)[i])])==1){
			ValLevel3 <- getSheets(WorkBook)[[paste("Validation Level 3", names(MyData)[i])]]
			}
		if(length(NamesSheets[NamesSheets==paste("Validation Level 3", names(MyData)[i])])==0){
			ValLevel3 <- createSheet(WorkBook, sheetName=paste("Validation Level 3", names(MyData)[i]))
			}
		addDataFrame("Singles indices summarizing the life time probability distribution", ValLevel3, row.names=F, col.names=F)
		addDataFrame("Median age at death", ValLevel3, row.names=F, col.names=F,startRow=3)
		addDataFrame(as.matrix(rownames(SI[[i]]$MedianAge)), ValLevel3, row.names=F, col.names=F, startRow=6)
		addDataFrame(SI[[i]]$MedianAge, ValLevel3, row.names=F, col.names=T, startRow=5,startColumn=2 + NberCol)
		addDataFrame("Entropy", ValLevel3, row.names=F, col.names=F,startRow=13)
		addDataFrame(as.matrix(rownames(SI[[i]]$Entropy)), ValLevel3, row.names=F, col.names=F, startRow=16)
		addDataFrame(SI[[i]]$Entropy, ValLevel3, row.names=F, col.names=T, startRow=15,startColumn=2 + NberCol)
		addDataFrame(paste("Cohort life expectancy for cohort in",min(MyData[[i]]$YearCom)), ValLevel3, row.names=F, col.names=F,startRow=23)
		addDataFrame(as.matrix(rownames(SI[[i]]$CohortLifeExp)), ValLevel3, row.names=F, col.names=F, startRow=26)
		addDataFrame(SI[[i]]$CohortLifeExp, ValLevel3, row.names=F, col.names=T, startRow=25,startColumn=2 + NberCol)
		}
	saveWorkbook(WorkBook, file = ExcelPath)	
	}

## ------------------------------------------------------------------------ ##
##  .ExportCohortLifeExp5inExcel function                                    ##
## ------------------------------------------------------------------------ ##

.ExportCohortLifeExp5inExcel = function(CohortLifeExp5){
	Path <- "Results/Excel"
	.CreateDirectory(paste("/",Path,sep=""))
	ExcelPath <- paste(Path, "/Validation.xlsx", sep = "")
	print(paste("Export the cohort life expectancies over 5 years in .../", ExcelPath, " ...", sep=""))
	if(!file.exists(paste("./",ExcelPath,sep=""))){
		WorkBook <- createWorkbook()
		EspGenPart5 <- createSheet(WorkBook, sheetName=paste(CohortLifeExp5[[1]]$NameMethod,"- CLE 5 -", names(CohortLifeExp5)[1]))
		saveWorkbook(WorkBook, file = ExcelPath)
		}
	if(file.exists(paste("./",ExcelPath,sep=""))){
		WorkBook <- loadWorkbook(file = ExcelPath)
		}
	for (i in 1 : length(CohortLifeExp5)){	
		NamesSheets <- names(getSheets(WorkBook))
		if(length(NamesSheets[NamesSheets==paste(CohortLifeExp5[[i]]$NameMethod,"- CLE 5 -", names(CohortLifeExp5)[i])])==1){
			EspGenPart5 <- getSheets(WorkBook)[[paste(CohortLifeExp5[[i]]$NameMethod,"- CLE 5 -", names(CohortLifeExp5)[i])]]
			}
		if(length(NamesSheets[NamesSheets==paste(CohortLifeExp5[[i]]$NameMethod,"- CLE 5 -", names(CohortLifeExp5)[i])])==0){
			EspGenPart5 <- createSheet(WorkBook, sheetName=paste(CohortLifeExp5[[i]]$NameMethod,"- CLE 5 -", names(CohortLifeExp5)[i]))
			}
		addDataFrame(paste("Cohort life expectancy over 5 years,",names(CohortLifeExp5)[i],"population,",CohortLifeExp5[[i]]$NameMethod), EspGenPart5, row.names=F, col.names=F)
		addDataFrame(as.data.frame(CohortLifeExp5[[i]]$CohortLifeExp5), EspGenPart5, row.names=T, col.names=T, startRow=3)
		}
	saveWorkbook(WorkBook, file = ExcelPath)	
	}

## ------------------------------------------------------------------------ ##
##  .ExportPeriodicLifeExpinExcel function                                   ##
## ------------------------------------------------------------------------ ##

.ExportPeriodicLifeExpinExcel = function(PerLifeExpFitted, AgeComp){
	Path <- "Results/Excel"
	.CreateDirectory(paste("/",Path,sep=""))
	ExcelPath <- paste(Path, "/Validation.xlsx", sep = "")
	print(paste("Export the periodic life expectancies in .../", ExcelPath, " ...", sep=""))
	if(!file.exists(paste("./",ExcelPath,sep=""))){
		WorkBook <- createWorkbook()
		PerLifeExpRes <- createSheet(WorkBook, sheetName=paste(PerLifeExpFitted[[1]]$NameMethod,"- PLE -", names(PerLifeExpFitted)[1]))
		saveWorkbook(WorkBook, file = ExcelPath)
		}
	if(file.exists(paste("./",ExcelPath,sep=""))){
		WorkBook <- loadWorkbook(file = ExcelPath)
		}
	for (i in 1 : length(PerLifeExpFitted)){	
		NamesSheets <- names(getSheets(WorkBook))
		if(length(NamesSheets[NamesSheets==paste(PerLifeExpFitted[[i]]$NameMethod,"- PLE -", names(PerLifeExpFitted)[i])])==1){
			PerLifeExpRes <- getSheets(WorkBook)[[paste(PerLifeExpFitted[[i]]$NameMethod,"- PLE -", names(PerLifeExpFitted)[i])]]
			}
		if(length(NamesSheets[NamesSheets==paste(PerLifeExpFitted[[i]]$NameMethod,"- PLE -", names(PerLifeExpFitted)[i])])==0){
			PerLifeExpRes <- createSheet(WorkBook, sheetName=paste(PerLifeExpFitted[[i]]$NameMethod,"- PLE -", names(PerLifeExpFitted)[i]))
			}
		addDataFrame(paste("Periodic life expectancies (until",max(AgeComp),"years),",names(PerLifeExpFitted)[i],"population,",PerLifeExpFitted[[i]]$NameMethod), PerLifeExpRes, row.names=F, col.names=F)
		addDataFrame(as.data.frame(PerLifeExpFitted[[i]]$PerLifeExp), PerLifeExpRes, row.names=T, col.names=T, startRow=3)
		}
	saveWorkbook(WorkBook, file = ExcelPath)	
	}
