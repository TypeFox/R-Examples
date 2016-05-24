### R code from vignette source 'fscaret.Rnw'

###################################################
### code chunk number 1: install (eval = FALSE)
###################################################
## install.packages("fscaret", dependencies = c("Depends", "Suggests"))


###################################################
### code chunk number 2: install (eval = FALSE)
###################################################
## install.packages("fscaret", dependencies = c("Depends"))


###################################################
### code chunk number 3: load_data (eval = FALSE)
###################################################
## basename_file <- "My_database"
## file_name <- paste(basename_file,".csv",sep="")


###################################################
### code chunk number 4: load_data (eval = FALSE)
###################################################
## matrixTrain <- read.csv(file_name,header=TRUE,sep="\t",
## 		strip.white = TRUE, na.strings = c("NA",""))


###################################################
### code chunk number 5: load_data (eval = FALSE)
###################################################
## matrixTrain <- as.data.frame(matrixTrain)


###################################################
### code chunk number 6: funcRegPred_all
###################################################
library(fscaret)
data(funcRegPred)
funcRegPred


###################################################
### code chunk number 7: funcClassPred_all
###################################################
library(fscaret)
data(funcClassPred)
funcClassPred


###################################################
### code chunk number 8: fscaret_example (eval = FALSE)
###################################################
## my_res_foba <- myFS$VarImp$model$foba
## my_res_foba <- structure(my_res_foba,class="train")


###################################################
### code chunk number 9: fscaret_example (eval = FALSE)
###################################################
## 
## library(fscaret)
## data(dataset.train)
## data(dataset.test)
## 
## trainDF <- dataset.train
## testDF <- dataset.test
## 
## myFS<-fscaret(trainDF, testDF, myTimeLimit = 5, preprocessData=TRUE,
## 	      Used.funcRegPred=c("pcr","pls"), with.labels=TRUE,
## 	      supress.output=TRUE, no.cores=1)
## myRES_tab <- myFS$VarImp$matrixVarImp.MSE[1:10,]
## myRES_tab <- subset(myRES_tab, select=c("pcr","pls","SUM%","ImpGrad","Input_no"))
## myRES_rawMSE <- myFS$VarImp$rawMSE
## myRES_PPlabels <- myFS$PPlabels


###################################################
### code chunk number 10: fscaret_example
###################################################
library(fscaret)

# if((Sys.info()['sysname'])=="SunOS"){
myRES_tab <- data.frame(pcr = c(5.862841e+01, 1.567799e+01, 1.916511e+01, 2.519981e-01, 1.872058e-02, 
				2.880832e-04, 5.880416e-04, 7.190168e-05, 1.570926e-06, 1.081909e-06),
                        pls = c(5.227714e+01, 2.741963e+01, 1.995465e+01, 3.161112e-01, 3.079973e-02, 
				1.324904e-03, 2.781880e-04, 6.894892e-05, 2.697715e-06, 1.078743e-06),
                        "SUM" = c(1.000000e+02, 3.885975e+01, 3.527303e+01, 5.122461e-01, 4.465089e-02, 
				1.454379e-03, 7.810516e-04, 1.270005e-04, 3.848898e-06, 1.948191e-06),
			ImpGrad=c(0.000000, 61.140253, 9.229896, 98.547769, 91.283313, 96.742777, 46.296556, 
				  83.739807, 96.969384, 49.383143),
			Input_no=c("4","5","22","23","2","13","9","1","17","21"))
names(myRES_tab)[length(myRES_tab)-2]<-"SUM%"			
			
myRES_rawMSE <- data.frame(pcr = c(716.6597),
			  pls = c(671.8195))
			  
myRES_PPlabels <- data.frame("Orig Input No"=c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 
						14, 15, 16, 17, 18, 19, 20, 22, 23, 24, 25, 26, 
						27, 29),
			      Labels = c("Balaban.index", "Dreiding.energy", 
					  "Fused.aromatic.ring.count", "Hyper.wiener.index", 
					  "Szeged.index", "Ring.count.of.atom", 
					  "pI", "Quaternary_structure", 
					  "PLGA_Mw", "La_to_Gly", 
					  "PVA_conc_inner_phase", "PVA_conc_outer_phase", 
					  "PVA_Mw", "Inner_phase_volume", 
					  "Encaps_rate", "PLGA_conc", 
					  "PLGA_to_Placticizer", "diss_pH", 
					  "diss_add", "Prod_method", 
					  "Asymmetric.atom.count.1", "Hyper.wiener.index.1", 
					  "Szeged.index.1", "count", 
					  "pH_14_logd", "bpKa2", 
					  "Cyclomatic.number.2"))
			            
# } else {

# data(dataset.train)
# data(dataset.test)

# trainDF <- dataset.train
# testDF <- dataset.test

# myFS<-fscaret(trainDF, testDF, myTimeLimit = 5, preprocessData=TRUE,regPred=TRUE,
#	      Used.funcRegPred=c("pcr","pls"), with.labels=TRUE,
#	      supress.output=TRUE, no.cores=1, saveModel=FALSE)

# myRES_tab <- myFS$VarImp$matrixVarImp.MSE[1:10,]

# myRES_rawMSE <- myFS$VarImp$rawMSE
# myRES_PPlabels <- myFS$PPlabels

# }

myRES_tab <- subset(myRES_tab, select=c("pcr","pls","SUM%","ImpGrad","Input_no"))



###################################################
### code chunk number 11: fscaret_example_class
###################################################

# library(MASS)
# 
# # make testing set
# data(Pima.te)
# 
# Pima.te[,8] <- as.numeric(Pima.te[,8])-1
# 
# myDF <- Pima.te
# 
# myFS.class<-fscaret(myDF, myDF, myTimeLimit = 20, preprocessData=FALSE, with.labels=TRUE, classPred=TRUE, regPred=FALSE, Used.funcClassPred=c("knn","rpart"),
# 	      supress.output=FALSE, no.cores=1)
# 
# print(myFS.class)
# myRES.class_tab <- myFS.class$VarImp$matrixVarImp.MeasureError[,]
# myRES.class_tab <- subset(myRES.class_tab, select=c("knn","rpart","SUM%","ImpGrad","Input_no"))
# myRES.class_rawError <- myFS.class$VarImp$rawMeasureError


###################################################
### code chunk number 12: fscaret_example_class (eval = FALSE)
###################################################
## library(MASS)
## 
## # make testing set
## data(Pima.te)
## 
## Pima.te[,8] <- as.numeric(Pima.te[,8])-1
## 
## myDF <- Pima.te
## 
## myFS.class<-fscaret(myDF, myDF, myTimeLimit = 5, preprocessData=FALSE,
## 		    with.labels=TRUE, classPred=TRUE,regPred=FALSE, 
## 		    Used.funcClassPred=c("knn","rpart"), supress.output=TRUE, no.cores=1)
## myRES.class_tab <- myFS.class$VarImp$matrixVarImp.MeasureError
## myRES.class_tab <- subset(myRES.class_tab, select=c("knn","rpart","SUM%","ImpGrad","Input_no"))
## myRES.class_rawError <- myFS.class$VarImp$rawMeasureError


###################################################
### code chunk number 13: fscaret_example
###################################################
# Print out the Variable importance results for MSE scaling
print(myRES_tab)


###################################################
### code chunk number 14: fscaret_example
###################################################
# Print out the generalization error for models
print(myRES_rawMSE)


###################################################
### code chunk number 15: fscaret_example
###################################################
# Print out the reduced number of inputs after preprocessing
print(myRES_PPlabels)


###################################################
### code chunk number 16: barPlot
###################################################

# Present variable importance on barplot
a=0.9
b=0.7
c=2

# if((Sys.info()['sysname'])=="SunOS"){
myFS <- NULL
myFS$VarImp$matrixVarImp.MSE <- myRES_tab

# }


lk_row.mse=nrow(myFS$VarImp$matrixVarImp.MSE)

setEPS()

barplot1 <- barplot(myFS$VarImp$matrixVarImp.MSE$"SUM%"[1:(a*lk_row.mse)],
	    cex.names=b, las = c, xlab="Variables", ylab="Importance Sum%",
	    names.arg=c(myFS$VarImp$matrixVarImp.MSE$Input_no[1:(a*lk_row.mse)]))
	    
lines(x = barplot1, y = myFS$VarImp$matrixVarImp.MSE$"SUM%"[1:(a*lk_row.mse)])
points(x = barplot1, y = myFS$VarImp$matrixVarImp.MSE$"SUM%"[1:(a*lk_row.mse)])



###################################################
### code chunk number 17: fscaret_example_class (eval = FALSE)
###################################################
## # Print out the Variable importance results for F-measure scaling
## print(myRES.class_tab)


###################################################
### code chunk number 18: fscaret_example_class (eval = FALSE)
###################################################
## # Print out the generalization error for models
## print(myRES.class_rawError)


###################################################
### code chunk number 19: fscaret_issue (eval = FALSE)
###################################################
## library(fscaret)
## myFuncRegPred <- funcRegPred[which(funcRegPred!="partDSA")]
## 
## print(funcRegPred)
## 
## myFS<-fscaret(trainDF, testDF, myTimeLimit = 12*60*60, preprocessData=TRUE,regPred=TRUE,
##         Used.funcRegPred=myFuncRegPred, with.labels=TRUE,
## 	      supress.output=TRUE, no.cores=NULL, saveModel=FALSE)


