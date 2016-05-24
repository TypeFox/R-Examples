# file NominalLogisticBiplot/R/auxLibrary.R
# copyright (C) 2012-2013 J.C. Hernandez and J.L. Vicente-Villardon
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 or 3 of the License
#  (at your option).
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/
#


#This function plot the category points for the variable we have choosed as the first parameter
#, and the corresponding tesselation
#----------------------Parameters
  #voronoiVariable: object with the information needed about the variable we are treating.
  #LabelVar:boolean variable to choose if we want to see the labels of the variables
  #LabelInd:boolean variable to choose if we want to see the labels of the individuals
  #CexInd: parameter to specify the size of the individual points
  #CexVar: parameter to specify the size of the variable centers
  #ColorInd: parameter to specify the color of the individual points
  #ColorVar: parameter to specify the color of the variable centers
  #PchInd: parameter to specify the symbol of the individual points
  #PchVar: parameter to specify the symbol of the variable centers
  #AtLeastR2 <- 0.01   #It establishes the cut value to plot a variable attending to its R2 Nagelkerke value.
  #lines: boolean variable that specify if lines of the tesselation will be plotted.
  #LabValVar: this paremeter is an array with the values of the labels for the variable which will be plotted 
plot.voronoiprob <- function(voronoiVariable,LabelVar,CexVar,ColorVar,PchVar,AtLeastR2=0.01,lines = TRUE,LabValVar=NULL) {
    vorprob = voronoiVariable$vorprob
    regVisible = voronoiVariable$regVisible
    nameVar = voronoiVariable$nameVar
    R2 = voronoiVariable$model$Nagelkerke

   if(R2 >= AtLeastR2){
    	CoorReal = cbind(vorprob$x, vorprob$y)
    	CoorDummy = cbind(vorprob$dummy.x, vorprob$dummy.y)

   points(vorprob$Centers[, 1], vorprob$Centers[, 2], col = ColorVar, cex = CexVar,pch = PchVar)
   if(LabelVar){
      if(is.null(LabValVar)){
         text(vorprob$Centers[, 1], vorprob$Centers[, 2], paste(regVisible,"_",nameVar ,sep="") , col = ColorVar, cex = CexVar,pos=1,offset=0.2)
      }else{
         text(vorprob$Centers[, 1], vorprob$Centers[, 2], LabValVar[regVisible] , col = ColorVar, cex = CexVar,pos=1,offset=0.2)      
      }
   } 

   if(lines){
    	for (i in 1:dim(CoorReal)[1]) {
    		if (vorprob$n1[i] > 0)
    			segments(CoorReal[i, 1], CoorReal[i, 2], CoorReal[vorprob$n1[i], 1], CoorReal[vorprob$n1[i], 2], col = ColorVar)
    		else segments(CoorReal[i, 1], CoorReal[i, 2], CoorDummy[(-1 * vorprob$n1[i]), 1], CoorDummy[(-1 * vorprob$n1[i]), 2], col = ColorVar)
    		if (vorprob$n2[i] > 0)
    			segments(CoorReal[i, 1], CoorReal[i, 2], CoorReal[vorprob$n2[i], 1], CoorReal[vorprob$n2[i], 2], col = ColorVar)
    		else segments(CoorReal[i, 1], CoorReal[i, 2], CoorDummy[(-1 * vorprob$n2[i]), 1], CoorDummy[(-1 * vorprob$n2[i]), 2], col = ColorVar)
    		if (vorprob$n3[i] > 0)
    			segments(CoorReal[i, 1], CoorReal[i, 2], CoorReal[vorprob$n3[i], 1], CoorReal[vorprob$n3[i], 2], col = ColorVar)
    		else segments(CoorReal[i, 1], CoorReal[i, 2], CoorDummy[(-1 * vorprob$n3[i]), 1], CoorDummy[(-1 * vorprob$n3[i]), 2], col = ColorVar)
    	}
	 }
  }
}

#This function plot the category points for the variable we have choosed as the first parameter
#, and the corresponding tesselation
#----------------------Parameters
  #voronoiVariable: object with the information needed about the variable we are treating.
  #AtLeastR2 <- 0.01   #It establishes the cut value to plot a variable attending to its R2 Nagelkerke value.
  #x: row coordinates for the plane we have choosed for the study
  #line: boolean variable that specify if lines of the tesselation will be plotted.
  #LabelVar:boolean variable to choose if we want to see the labels of the variables
  #LabelInd:boolean variable to choose if we want to see the labels of the individuals
  #CexInd: parameter to specify the size of the individual points
  #CexVar: parameter to specify the size of the variable centers
  #ColorInd: parameter to specify the color of the individual points
  #ColorVar: parameter to specify the color of the variable centers
  #PchInd: parameter to specify the symbol of the individual points
  #PchVar: parameter to specify the symbol of the variable centers
  #LabValVar: this paremeter is an array with the values of the labels for the variable which will be plotted 
plot2CategLine <- function(voronoiVariable,x,AtLeastR2,line,LabelVar,CexVar,ColorVar,PchVar,LabValVar=NULL){
      beta = voronoiVariable$model$beta
      nameVar = voronoiVariable$nameVar
      R2 = voronoiVariable$model$Nagelkerke
      equivFit = voronoiVariable$equivFit

     if(R2 >= AtLeastR2){
        if(is.null(nrow(equivFit))){
           numValDif = 2
         }else{
           numValDif = max(equivFit[1,])
         }
         nRowsdata = nrow(x)
         centerg = matrix(0,1,2)
         numP = 0
         for (i in 1:nRowsdata){
            if(beta[1,1] + beta[1,2]*x[i,1] + beta[1,3]*x[i,2] > 0){
                numP = numP + 1
                centerg[1,1] = centerg[1,1] + x[i,1]
                centerg[1,2] = centerg[1,2] + x[i,2]
            }
         }
         if(numP == 0){  
              centerg[1,1] = sum(x[,1])
              centerg[1,2] = sum(x[,2])
              numP = nRowsdata
         }
         centerg[1,1] = centerg[1,1]/numP
         centerg[1,2] = centerg[1,2]/numP

         PosCatFit = matrix(0,1,2)
         if(is.null(nrow(equivFit))){
           PosCatFit[1,1]=1
           PosCatFit[1,2]=2
         }else{
           n=1
           for(k in 1:numValDif){
              if(equivFit[2,k] > 0){
                PosCatFit[1,n] = k
                n = n + 1
              }
           }
         }

         if(line){
            abline(-beta[1,1]/beta[1,3],-beta[1,2]/beta[1,3])
         }
         points(centerg[1,1], centerg[1,2],pch=PchVar,cex=CexVar,col=ColorVar)
         categF1 = 1/(1 + exp(beta[1,1] + beta[1,2]*centerg[1,1] + beta[1,3]*centerg[1,2]))
         categF2 = exp(beta[1,1] + beta[1,2]*centerg[1,1] + beta[1,3]*centerg[1,2])/(1 + exp(beta[1,1] + beta[1,2]*centerg[1,1] + beta[1,3]*centerg[1,2]))
         firstCateg = FALSE
         if(categF1 > categF2){
            labelF = PosCatFit[1,1]
            firstCateg = TRUE
         }else{
            labelF = PosCatFit[1,2]
         }
         if(LabelVar){
            if(is.null(LabValVar)){
              text(centerg[1,1], centerg[1,2], paste(labelF,"_",nameVar,sep="") ,col=ColorVar, cex = CexVar,pos=1,offset=0.2)
            }else{
              text(centerg[1,1], centerg[1,2], LabValVar[labelF] ,col=ColorVar, cex = CexVar,pos=1,offset=0.2)
            }         
         }
         
         xPointMed = (-beta[1,1]/beta[1,3] + centerg[1,1]*(beta[1,3]/beta[1,2]) - centerg[1,2])/(beta[1,3]/beta[1,2] + beta[1,2]/beta[1,3])
         yPointMed = (-beta[1,2]/beta[1,3])* xPointMed -(beta[1,1]/beta[1,3])

         xSimet = 2*xPointMed - centerg[1,1]
         ySimet = 2*yPointMed - centerg[1,2]

         points(xSimet, ySimet,pch=PchVar,cex=CexVar,col=ColorVar)

         if(firstCateg == TRUE){
            labelT = PosCatFit[1,2]
         }else{
            labelT = PosCatFit[1,1]
         }
         
         if(LabelVar){
            if(is.null(LabValVar)){
               text(xSimet, ySimet, paste(labelT,"_",nameVar,sep="") ,col=ColorVar, cex = CexVar,pos=1,offset=0.2)
            }else{
               text(xSimet, ySimet, LabValVar[labelT] ,col=ColorVar, cex = CexVar,pos=1,offset=0.2)            
            }         
         
         }

     }
} 

mvvSingleVariable <- function(nameVar,numcateg,beta,varstudyC,rowCoords,planex,planey,
                              numFactors,QuitNotPredicted,penalization = 0.2, cte = TRUE,
                              tol = 1e-04, maxiter = 200,ShowResults=TRUE){
    nRowsdata <- nrow(varstudyC)
    numValDif = numcateg

    varNomBiplot <- list()
    varNomBiplot$model <- list()
    varNomBiplot$model$beta = beta
    varNomBiplot$model$Nagelkerke = 0.02        
    varNomBiplot$nameVar = nameVar
    varNomBiplot$vorprob = 0   
    varNomBiplot$regVisible = 0
    varNomBiplot$numFit = 0
    varNomBiplot$equivFit = 0

    coords = cbind(rowCoords[,planex],rowCoords[,planey])            
    
    varNomBiplot = mvvObjectFill(varNomBiplot,numValDif,beta,coords,varstudyC,QuitNotPredicted,penalization=penalization,cte=cte,tol=tol,maxiter=maxiter,ShowResults=ShowResults)

    class(varNomBiplot)='variableNominalBiplot'
    return(varNomBiplot)

}

#This function print on screen principal characteristics of the nominal biplot object calculated for the data
#----------------------Parameters
  #BLM: object of class type group.variables.tesselations
WriteMultinomialLogisticBiplot <- function(BLM){
      print(paste("Coordinates of the rows for the plane:",BLM$planex,"-",BLM$planey,"\n",sep=""))
      print(BLM$rowCoords)
      for(i in 1:BLM$numVar){
        print(paste("Variable:",BLM$biplotNom[,i]$nameVar,sep=""))
        print(paste("PCC:",BLM$biplotNom[,i]$PCC,sep=""))
        print("Coefficients from adjusted model")
        print(BLM$biplotNom[,i]$model$beta)
        print(paste("Pseudo R2 Nagelkerke:",BLM$biplotNom[,i]$model$Nagelkerke,sep=""))
        print("Tesselation information:")
        print(BLM$biplotNom[,i]$vorprob)        
        print("-----------------")
      }
}


#This function check the data set to study keeping its information in a class
#----------------------Parameters
  #datanom: it could be a data.frame or a matrix with the nominal data
#Para que funcione bien las categorias no deben estar codificadas con ceros. Solo se corrige en el 
#caso de dos categorias, pero no en el resto.  
CheckDataSet <- function(datanom){
    
   typeDataFrame = FALSE     
   nRowInit = nrow(datanom)
   datanom <- na.omit(datanom)
    
    nRow = nrow(datanom)
    if(nRow < nRowInit){
      print("Be careful. Some rows have been deleted because they present NA values. Check your data.")
    }
    
    RowNames = NULL
    ColNames = NULL
    LevelNames = list()
    contLevel = 1      
        
    notNumeric = NULL
    if(is.data.frame(datanom)){
      typeDataFrame = TRUE
      RowNames = rownames(datanom)
      ColNames = colnames(datanom)
      for(i in 1:dim(datanom)[2]){        
          if(is.factor(datanom[,i])){
              LevelNames[[contLevel]] = levels(datanom[,i])
              contLevel = contLevel + 1
              datanom[,i] = as.numeric(datanom[,i])
          }else if(is.numeric(datanom[,i])){
                    if(!is.integer(datanom[,i])){
                       #print(paste("Variable ",i," will be transformed to integer.",sep=""))
                       datanom[,i] = as.integer(datanom[,i])
                    }
                    LevelNames[[contLevel]] = c(1:max(datanom[,i]))
                    contLevel = contLevel + 1 
                }else{
                  print(paste("Variable ",i," will be omitted because it is not nominal",sep=""))
                  notNumeric = cbind(notNumeric,i)
                }
      }
      notNumeric = as.vector(notNumeric)
      dataSet = NULL
      for(j in 1:dim(datanom)[2]){
        if(length(which(notNumeric==j))==0){
          dataSet = cbind(dataSet,datanom[,j])
        }
      }
    }else if(is.matrix(datanom)){
              RowNames = dimnames(datanom)[[1]][1:nrow(datanom)]
              ColNames = dimnames(datanom)[[2]][1:ncol(datanom)]
              datanom <- as.data.frame(datanom)
              for(i in 1:dim(datanom)[2]){
                 if(!is.null(levels(datanom[,i]))){
                   LevelNames[[contLevel]] = levels(datanom[,i])
                 }else{
                   LevelNames[[contLevel]] = c(1:max(datanom[,i])) 
                 }
                 contLevel = contLevel + 1
                 datanom[,i] = as.numeric(datanom[,i])                 
                 if(!is.integer(datanom[,i])){
                   #print(paste("Variable ",i," will be transformed to integer.",sep=""))
                   datanom[,i] = as.integer(datanom[,i])
                 }        
              }
              dataSet = datanom
          }else{
            stop("Data set should be a data frame or a matrix")
          }
    #In case that RowNames or ColNames are NULL we fix the name of the variables and rows
   	if (is.null(RowNames)){
  		RowNames <- rownames(datanom, do.NULL = FALSE, prefix = "I")
  		dimnames(datanom)[[1]] = RowNames
    }

   	if (is.null(ColNames)){
  		ColNames <- colnames(datanom, do.NULL = FALSE, prefix = "V")
  		dimnames(datanom)[[2]] = ColNames
    }

      numVarDef <- ncol(dataSet)    
      datanomcats = apply(dataSet[,1:numVarDef], 2, function(x) nlevels(as.factor(x)))
      for(i in 1:numVarDef){        
        if(max(dataSet[,i]) > datanomcats[i]){
          print(paste("Please, it would be desirable that you recode variable ", i ,
                " from the data set because its maximum value exceeds the distinct values it presents,
                so there are some categories that are not present",sep=""))
          columValuesOrd = sort(unique(dataSet[,i]))
          print(paste("Variable ",dimnames(datanom)[[1]][i]," only take the values:",sep=""))
          print(columValuesOrd)
          ActLevelNames = NULL
          newColdataSet = dataSet[,i]      
          for(j in 1:datanomcats[i]){
              newColdataSet = replace(newColdataSet,newColdataSet==columValuesOrd[j],j)
              if(typeDataFrame){
                if(length(LevelNames[[i]]) > 0){
                    ActLevelNames <- c(ActLevelNames,LevelNames[[i]][columValuesOrd[j]])
                }
              }else{
                  ActLevelNames = columValuesOrd
              }
          }
          dataSet[,i] = newColdataSet
          LevelNames[[i]] = ActLevelNames
        }else{
          if((datanomcats[i] == 2)&(max(dataSet[,i]) < datanomcats[i])){
            dataSet[,i] = dataSet[,i] + 1
            ActLevelNames = NULL
            columValuesOrd = sort(unique(dataSet[,i]))
            for(j in 1:datanomcats[i]){
              if(typeDataFrame){
                if(length(LevelNames[[i]]) > 0){
                    ActLevelNames <- c(ActLevelNames,LevelNames[[i]][columValuesOrd[j]])
                }
              }else{
                  ActLevelNames = columValuesOrd
              }
            }
            LevelNames[[i]] = ActLevelNames
          }
        }
      }#end for
    
    dimnames(dataSet)[[1]] = RowNames
    dimnames(dataSet)[[2]] = ColNames
     
    data.nominal<-list()
    data.nominal$datanom = dataSet
    data.nominal$RowNames = RowNames
    data.nominal$ColumNames = ColNames
    data.nominal$LevelNames = LevelNames       
    
    class(data.nominal)='data.nominal'
    return(data.nominal)

}



#This function plot the text beside the individual points at a concrete positicion.
#, depending from the positive or negative coordinate in the first axis.
#----------------------Parameters
  #A: matrix with the row coordinates.
  #CexPoints:parameter to specify the size of the text for each individual points
  #ColorPoints: parameter to specify the color of the text for each individual points
textsmart <- function(A, CexPoints, ColorPoints) {
	n = dim(A)[1]

	if (length(CexPoints == 1))
		CexPoints = rep(CexPoints, n)

	if (length(ColorPoints == 1))
		ColorPoints = rep(ColorPoints, n)
	for (i in 1:n) {
		if (A[i, 1] > 0)
			markerpos = 4
		else markerpos = 2
		text(A[i, 1], A[i, 2], rownames(A)[i], col = ColorPoints[i], pos = markerpos, offset = 0.2)
	}
}


#Function that calculates the variable models using the estimation of individuals by different methods.
  #----------------------Parameters--------------
  #datanom: matrix with the data to do the analysis
  #indRedCoords: estimation of the coordinates for the individuals in complete reduced dimension
  #tol = 1e-04, maxiter = 100, penalization = 0.1,: items used in RidgeMultinomialRegression procedure to calculate the parameters.
  #showIter = boolean parameter to decide if we want to show iteration information in RidgeMultinomialRegression
#This function returns an object with so many columns as variables and it contains for each of them the 
#ridge regression model with all the information.
CalculateVariableModels <- function(datanom,indRedCoords, penalization, tol, maxiter,showIter){
  nColsdata <- dim(datanom)[2]
  VariableModels = 0
	for (j in 1:nColsdata){
	  #print(paste("CalculateVariableModels: column ",j,sep="")) 
  	model = RidgeMultinomialRegression(datanom[, j], indRedCoords, penalization = penalization, tol = tol, maxiter = maxiter,showIter = showIter)  	
  	varEstimation = model
  	varEstimation$name = dimnames(datanom)[[2]][j]
  	VariableModels=cbind(VariableModels,varEstimation)
  }                                                              
  VariableModels=VariableModels[,2:(nColsdata+1)]
  return(VariableModels)
}



#Function that builds the complete tesselations object with the information of all the variables.
  #----------------------Parameters--------------
  #nlbo: object with the estimation of the rows and the parameters used in this estimation
  #planex,planey: these two parameters keep the plane we choose for the representation
  #QuitNotPredicted: boolean parameter that indicates if we will be represent on the graph the categories not predicted
  #ReestimateInFocusPlane: boolean parameter to choose if we will reestimate the models in the concrete plane or we will choose the columms of our plane in beta matrix
  #ShowResults: boolean parameter to show results from the process of estimation
#This function returns an object with the number of variables, the coordinates on the plane selected and 
#a list with objects that contain the information about the tesselation and fitting for each variable.
BuildTesselationsObject <- function(nlbo,planex=1,planey=2,QuitNotPredicted,ReestimateInFocusPlane,ShowResults){

  nRowsdata <- dim(nlbo$dataSet$datanom)[1]
  nColsdata <- dim(nlbo$dataSet$datanom)[2]
  numVar <- ncol(nlbo$dataSet$datanom)    
  datanomcats = apply(nlbo$dataSet$datanom[,1:numVar], 2, function(x) nlevels(as.factor(x)))
  numFactors <- nlbo$numFactors
  fittingVars <- matrix(0,4,numVar)  
  dimnames(fittingVars)[[2]]= dimnames(nlbo$dataSet$datanom)[[2]][1:numVar]
  dimnames(fittingVars)[[1]]= c("PCC","CoxSnell","Macfaden","Nagelkerke")
  rowscoor = nlbo$RowsCoords
  
  x <- cbind(rowscoor[,planex],rowscoor[,planey])

  variableTess = 0

  for(l in 1:numVar){
    if(ShowResults) print(paste("Iteracion.Variable l=",l,sep=""))
    nameVar = nlbo$dataSet$ColumNames[l]
    numValDif = datanomcats[[l]]
    y<-matrix(0,nRowsdata,1)
    for(i in 1:nRowsdata){
      y[i,1] = nlbo$dataSet$datanom[i,l]
    }

    betaRidge = nlbo$VariableModels[,l]
    if(ReestimateInFocusPlane == TRUE){
      betaRidge = RidgeMultinomialRegression(y,x,nlbo$penalization,nlbo$cte,nlbo$tol,nlbo$maxiter,ShowResults)
    }else{
      betaRidge$beta = cbind(betaRidge$beta[,1],betaRidge$beta[,planex + 1],betaRidge$beta[,planey + 1])      
    }
    mvv = MakeVoronoiVariable(l,nlbo,planex,planey,betaRidge,QuitNotPredicted,ShowResults)

    fittingVars[1,l]= round(mvv$PCC,2)    
    fittingVars[2,l]= round(mvv$model$CoxSnell,2)
    fittingVars[3,l]= round(mvv$model$MacFaden,2)                                   
    fittingVars[4,l]= round(mvv$model$Nagelkerke,2)
       
    variableTess=cbind(variableTess,mvv)
  }
  variableTess=variableTess[,2:(nColsdata+1)]

  if(ShowResults) print("Fitting adjustment for the variables")
  if(ShowResults) print(fittingVars)

  modelgvt<-list()
  modelgvt$rowCoords = x
  modelgvt$planex = planex
  modelgvt$planey = planey
  modelgvt$biplotNom = variableTess
  modelgvt$numVar = numVar

  class(modelgvt)='group.variables.tesselations'
  return(modelgvt)

}


#Function that creates for a concrete variable an object with all the information about the tesselation and 
#fitting for it.
  #------------------Parameters-----------------
  #ivar: number of the column in the dataset for the variable choosed for the analysis
  #nlbo: object with the estimation of the rows and the parameters used in this estimation
  #planex,planey: plane selected for the analysis of the variable
  #betaRidge: parameters estimated for this variable in a first step of the procedure
  #ShowResults: boolean parameter to show results from the process of estimation
#This function returns an object with the information of the model fitted for this variable, with the tesselation
#calculated for it,with the percentage of correct classifications, with the number of visible categories
#and the equivalence with the original values of the variable.
MakeVoronoiVariable <- function(ivar,nlbo,planex,planey,betaRidge,QuitNotPredicted,ShowResults){
    nRowsdata <- dim(nlbo$dataSet$datanom)[1]
    nColsdata <- dim(nlbo$dataSet$datanom)[2]
    numVar <- ncol(nlbo$dataSet$datanom)    
    datanomcats = apply(nlbo$dataSet$datanom[,1:numVar], 2, function(x) nlevels(as.factor(x)))
    numFactors = nlbo$numFactors
    l=ivar
    nameVar = dimnames(nlbo$dataSet$datanom)[[2]][l]
    numValDif = datanomcats[[l]]

    varNomBiplot <- list()
    varNomBiplot$model = betaRidge
    varNomBiplot$nameVar = nameVar
    varNomBiplot$vorprob = 0   
    varNomBiplot$regVisible = 0
    varNomBiplot$numFit = 0
    varNomBiplot$equivFit = 0

    beta = betaRidge$beta
    varstudyC = nlbo$dataSet$datanom[1:nRowsdata,l]
    coords = cbind(nlbo$RowsCoords[,planex],nlbo$RowsCoords[,planey])           
    
    varNomBiplot = mvvObjectFill(varNomBiplot,numValDif,beta,coords,varstudyC,QuitNotPredicted,nlbo$penalization,nlbo$cte,nlbo$tol,nlbo$maxiter,nlbo$show)

    class(varNomBiplot)='variableNominalBiplot'
    return(varNomBiplot)

}

mvvObjectFill <- function(varNomBiplot,numValDif,beta,coords,varstudyC,QuitNotPredicted=TRUE,penalization = 0.2, cte = TRUE, tol = 1e-04, maxiter = 200,ShowResults=TRUE){
                          
    nRowsdata = nrow(coords)
    varFitPCC = AdjustFitting(coords,varstudyC,beta,showTable=ShowResults)
    varNomBiplot$PCC = varFitPCC$PCC

    if(nrow(beta) == 1){ 
        varNomBiplot$vorprob = 0
        varNomBiplot$regVisible = 0
        varNomBiplot$numFit = 2
        varNomBiplot$equivFit = 0
        numValFitNow = apply(varFitPCC$varfit, 2, function(x) nlevels(as.factor(x)))
         if(numValFitNow == 1){
          varNomBiplot$numFit = 1
          varNomBiplot$equivFit = varFitPCC$varfit[1,1]
         }
    }else{ 
        numValFitNow = apply(varFitPCC$varfit, 2, function(x) nlevels(as.factor(x)))
        if(numValFitNow == 1){
          varNomBiplot$numFit = 1
         varNomBiplot$equivFit = varFitPCC$varfit[1,1]
        }else{  
            vorprob=Generators(beta)
         	  ngrup = nrow(beta) + 1
            regVisible = matrix(0,1,ngrup - sum(vorprob$hideCat))
            contregVisible = 1
            for(i in 1:ngrup){
              if(vorprob$equivRegiones[2,i] > 0){
                regVisible[1,contregVisible]= vorprob$equivRegiones[1,i]
                contregVisible = contregVisible + 1
              }
            }
            numValFitAnt = numValDif
            if(numValFitNow < numValFitAnt){
                equivFit<-matrix(0,3,numValDif)

                for(j in 1:numValDif){
                  equivFit[1,j] = j        
                }
                for(i in 1:nRowsdata){
                  if(sum(equivFit[2,]) < numValDif){
                    equivFit[2,varFitPCC$varfit[i,1]] = 1 
                  }
                }
                numFit = 0
                for(j in 1:numValDif){
                  if(equivFit[2,j] > 0){
                    numFit = numFit + 1
                    equivFit[3,j] = numFit  
                  }
                }
                yp<-matrix(0,nRowsdata,1)
                for(i in 1:nRowsdata){
                  yp[i,1] = equivFit[3,varFitPCC$varfit[i,1]]
                }
                x <- coords

                if(QuitNotPredicted == TRUE){
                   while(numValFitNow < numValFitAnt){
                      betaFRidge = RidgeMultinomialRegression(yp, x, penalization, cte, tol, maxiter)
                      betaF = betaFRidge$beta
                      if(nrow(betaF) == 1){
                        varNomBiplot$model$beta = betaF  
                        varNomBiplot$vorprob = 0
                        varNomBiplot$regVisible = 0
                        varNomBiplot$numFit = numValFitNow
                        varNomBiplot$equivFit = equivFit
                        break
                      }else{
                          vorprobF = Generators(betaF)
    
                        	ngrupF = nrow(betaF) + 1
                          regVisibleF = matrix(0,1,ngrupF - sum(vorprobF$hideCat))
                          contregVisibleF = 1
                          for(i in 1:ngrupF){
                            indexL = which(equivFit[3,] == i)
                            if(vorprobF$hideCat[i] == 0){
                              regVisibleF[1,contregVisibleF]= equivFit[1,indexL]
                              contregVisibleF = contregVisibleF + 1
                            }
                          }
                          varfit<-matrix(0,nRowsdata,1)
                          varstudyCF = yp
                          varFitPCCF = AdjustFitting(coords,varstudyCF,betaF,showTable=ShowResults) 
                     }
                      yp<-matrix(0,nRowsdata,1)
                      maxVarFit = max(varFitPCCF$varfit)
                      varfitValues=matrix(0,1,maxVarFit)
                      for(i in 1:nRowsdata){
                        varfitValues[1,varFitPCCF$varfit[i,1]]=1
                      }
                      for(i in 1:nRowsdata){
                        yp[i,1] = sum(varfitValues[,1:varFitPCCF$varfit[i,1]])
                      }
    
                      equivFit[2,]=0
                      for(i in 1:nRowsdata){
                          equivFit[2,which(equivFit[3,]==varFitPCCF$varfit[i,1])] = 1 
                      }
                      numFit = 0
                      equivFit[3,]=0
                      for(j in 1:numValDif){
                        if(equivFit[2,j] > 0){
                          numFit = numFit + 1
                          equivFit[3,j] = numFit  
                        }
                      }
    
                      numValFitAnt = numValFitNow
                      numValFitNow = apply(varFitPCCF$varfit, 2, function(x) nlevels(as.factor(x)))
                   }#Fin del while
    
                   if(nrow(betaF) > 1){
                      varNomBiplot$vorprob = vorprobF
                      varNomBiplot$regVisible = regVisibleF
                      varNomBiplot$numFit = numValFitNow
                      varNomBiplot$equivFit = equivFit
                   }
                }else{
                      varNomBiplot$vorprob = vorprob
                      varNomBiplot$regVisible = regVisible
                      varNomBiplot$numFit = numValFitNow
                      varNomBiplot$equivFit = equivFit
                 }
            }else{
                  varNomBiplot$vorprob = vorprob
                  varNomBiplot$regVisible = regVisible
                  varNomBiplot$numFit = numValFitNow
                  varNomBiplot$equivFit = ngrup
            }
        }
    }#fin del if-else
    
    return(varNomBiplot)
}


#Function that calculates category points result from invert the tesselation given by a variable
#----------------------Parameters--------------
  #beta: parameters of the multinomial logistic regression for the variable chosen
  #Borders: matrix(nborders x 2) with the original borders of the tesselation
  #BordersWH: matrix(nborders x 2) with the borders of the tesselation considering hidden categories, so that
            #Borders has been renumbered quitting hidden categories
  #nborders: number of borders of the tesselation
  #ngrupvisible: number of visible categories for the variable  
InvertTesselation <- function(beta,Borders,BordersWH,nborders,ngrupvisible) {

	ngrup = nrow(beta) + 1
	a = matrix(0, ngrup, ngrup)
	b = matrix(0, ngrup, ngrup)

	# Calculate the straight lines separating each pair of categories

	for (i in 1:(ngrup - 1)) {
		a[1, (i + 1)] = -1 * beta[i, 1]/beta[i, 3]
		b[1, (i + 1)] = -1 * beta[i, 2]/beta[i, 3]
	}

	if (ngrup > 2) {
		for (i in 1:(ngrup - 2)) for (j in (i + 1):(ngrup - 1)) {
			a[(i + 1), (j + 1)] = (beta[j, 1] - beta[i, 1])/(beta[i, 3] - beta[j, 3])
			b[(i + 1), (j + 1)] = (beta[j, 2] - beta[i, 2])/(beta[i, 3] - beta[j, 3])
		}
	}

  A = matrix(0, nborders, 2 * ngrupvisible)
	B = matrix(0, nborders, 2 * ngrupvisible)
	d = matrix(0, nborders, 1)

	for (i in 1:nborders) {
	  l = Borders[i, 1]
		m = Borders[i, 2]

		lwh = BordersWH[i, 1]
		mwh = BordersWH[i, 2]

		A[i, (2 * lwh - 1)] = b[l, m]
		A[i, (2 * lwh)] = -1
		A[i, (2 * mwh - 1)] = b[l, m]
		A[i, (2 * mwh)] = -1

		B[i, (2 * lwh - 1)] = -1/b[l, m]
		B[i, (2 * lwh)] = -1
		B[i, (2 * mwh - 1)] = 1/b[l, m]
		B[i, (2 * mwh)] = 1
		d[i, 1] = -2 * a[l, m]
	}

	G = rbind(A, B)
	e = rbind(d, matrix(0, dim(d)))
	Coord = ginv(t(G) %*% G) %*% t(G) %*% e
	Coord
}

#Function that calculates a matrix with the hidden categories for a set of real points.
#----------------------Parameters--------------
  #IndReal: matrix with the indices for each of the real points from a tesselation
  #ngrup: number of categories of the variable
  #nreal: number of real points in the tesselation
#The final matrix has ones in positions corresponding to hiden categories
HideCategories <- function(IndReal,ngrup,nreal){
  hc = matrix(1,ngrup,1)
  for(i in 1:nreal){
    hc[IndReal[i,1],1] = 0
    hc[IndReal[i,2],1] = 0
    hc[IndReal[i,3],1] = 0
  }
  hc
}

#Function that calculates the solutions of a second grade equation
#----------------------Parameters--------------
  #a,b,c: coefficients of x^2, x and independent term.
Eq2gSolve <- function(a,b,c){
    solns<-matrix(rep(NA,length(a)))
    solns<-cbind(solns,solns)
    colnames(solns)<-c("soln 1","soln 2")
    if(a==0 && b!=0){
      solns[1,1]= (-1)*c/b
      solns[1,2]= (-1)*c/b;
    }else if(a==0 && b==0){
      print("Coefficients a and b of the polynomial are zero")
    }else if((b^2 -4*a*c) < 0){
      print("The polynomial has complex roots")
    }else{
      solns[1,1]<-((-1)*b + sqrt(b^2 - 4*a*c))/(2*a)
      solns[1,2]<-((-1)*b - sqrt(b^2 - 4*a*c))/(2*a)
    }
    solns
}

#Function that calculates the percentage of correct classifications in a variable. It makes a 
#contingency table to analyse that variable.
#----------------------Parameters--------------
  #coords: coordinates x and y for the variable on the plane we have choosed.
  #varstudyC: values of the variable for all the individuals
  #beta: parameters estimated by multinomial logistic regression.
  #showTable: boolean parameter to indicate if we want table results will be showed
AdjustFitting <- function(coords,varstudyC,beta,showTable) {
	ngrup = nrow(beta) + 1
  xp = coords[,1]
	yp = coords[,2]

  coorrows = nrow(coords)
  coorcols = ncol(coords)


  xn = cbind(matrix(1, coorrows, 1), xp, yp)
	probab = matrix(0, coorrows, ngrup)

	for (i in 1:coorrows) {
		suma = 1
		for (j in 1:(ngrup - 1)) suma = suma + exp(sum(beta[j, ] * xn[i, ]))
		for (j in 1:(ngrup - 1)) probab[i, (j + 1)] = exp(sum(beta[j, ] * xn[i, ]))/suma
		probab[i, 1] = 1/suma
	}
  varfitted <- matrix(0,coorrows,1)
	for (i in 1:coorrows) {
	   varfitted[i] = which.max(probab[i,])
	}

  if(showTable){
      cTable <- CrossTable(varfitted,varstudyC,expected=TRUE,prop.r=TRUE,prop.c=TRUE,
                prop.t=TRUE,prop.chisq=FALSE,chisq=FALSE,fisher=FALSE,
                resid=TRUE,sresid=TRUE)
   }

  goodClasif = 0
 	for (i in 1:coorrows) {
 	  if(varfitted[i]==as.integer(varstudyC[i])){
      goodClasif = goodClasif + 1
    }
  }
  PCC <- round(goodClasif*100/coorrows,2)
  if(showTable) print(paste("Number of observations:",coorrows,",Good Clasified:",goodClasif,sep=""))
  if(showTable) print(paste("Percentage of good Clasifications:",round(goodClasif*100/coorrows,2),"%",sep=""))

  adjustFit = list()

	adjustFit$varfit = varfitted
	adjustFit$PCC = PCC
	class(adjustFit) <- "fitting"
	return(adjustFit)
 }

#Auxiliar function used in IntersectSegments function. It sees what is the relative position
#of a point in relation with a segment.
#----------------------Parameters--------------
  #segP1x,segP1y,segP2x,segP2y: coordinates x and y for the extremes of the segment.
  #px,py: x and y coordinates of a point
#It calculates the vectorial product to see the z- coordinate
# > 0 to the right, < 0 left, = 0 over the segment
Right_left <- function(segP1x,segP1y,segP2x,segP2y,px,py){
  v1x = px - segP1x
  v1y = py - segP1y
  v2x = segP2x - segP1x
  v2y = segP2y - segP1y
  vectorialProd = ((v1x*v2y)-(v2x*v1y))
  return(vectorialProd)
}

#Function that verifies if two segments are cut.
#When I ask ((border_p1*border_p2) <= 0) I verify if both are at the same side, and this is because
#both are positive or negative at the same time. 
#----------------------Parameters--------------
  #8 coordinates of the points that forms the two segments (2 points of each segment)
IntersectSegments <- function(seg1P1x,seg1P1y,seg1P2x,seg1P2y,seg2P1x,seg2P1y,seg2P2x,seg2P2y){
  border_p1 = Right_left(seg1P1x,seg1P1y,seg1P2x,seg1P2y,seg2P1x,seg2P1y)
  border_p2 = Right_left(seg1P1x,seg1P1y,seg1P2x,seg1P2y,seg2P2x,seg2P2y)
  Intersect = FALSE
  if((border_p1*border_p2) <= 0){
    border_p1 = Right_left(seg2P1x,seg2P1y,seg2P2x,seg2P2y,seg1P1x,seg1P1y)
    border_p2 = Right_left(seg2P1x,seg2P1y,seg2P2x,seg2P2y,seg1P2x,seg1P2y)
    if((border_p1*border_p2) <= 0){
      Intersect = TRUE
    }
  }
  return(Intersect)
}

#Function that builds the diagonal matrix with the values passed as parameter.
#----------------------Parameters--------------
  #d: value or vector that will be presented in the diagonal matrix
diagonal <- function(d){
    n=length(d)
    D=diag(1,n,n)
    diag(D)<-d
    D
}


patterns_eq <- function(nnodos, dims) {
	I = matrix(1:nnodos)
	for (i in 2:dims) {
		nf = dim(I)[1]
		nc = dim(I)[2]
		I2 = cbind(kronecker(matrix(I[1, ], 1, nc), matrix(1, nnodos, 1)), (1:nnodos))
		for (j in 2:nf)
     I2 = rbind(I2, cbind(kronecker(matrix(I[j, ], 1, nc), matrix(1, nnodos, 1)), (1:nnodos)))
		I = I2
	}
	return(I)
}

logit <- function(p) {
	logit = log(p/(1 - p))
	return(logit)
}

ColMax <- function (X)
{
  dimens = dim(X)
  n = dimens[1]
  p = dimens[2]
  Maxs = matrix(0, p, 1)
  for (j in (1:p))
    Maxs[j] = max(X[,j])

  return(Maxs)
}

EvalPolylogist <- function(X, par, Ncats){
	MaxCat=dim(par)[1]
	dims=dim(par)[2]-1
	nitems=dim(par)[3]
	nnodos=dim(X)[1]
	Numcats=sum(Ncats)
	CumCats=cumsum(Ncats)

	P=matrix(0,nnodos,Ncats[1])
  	if(is.vector(par[,,1])){
  	    z = exp(cbind(matrix(1, nnodos, 1), X) %*% as.matrix(par[,,1])[,1:(Ncats[1]-1)])
  	}else{
  	    z = exp(cbind(matrix(1, nnodos, 1), X) %*% t(par[,,1])[,1:(Ncats[1]-1)])
  	}
 	tot= 1+rowSums(z)
	P[,1]=1/tot
	for (k in 1:(Ncats[1]-1))
  	P[,(k+1)]=z[,k]/tot

	PT=P

	for (j in 2:nitems){
		P=matrix(0,nnodos,Ncats[j])
    	if(is.vector(par[,,1])){
     	    z = exp(cbind(matrix(1, nnodos, 1), X) %*% as.matrix(par[,,j])[,1:(Ncats[j]-1)])
    	}else{
    	    z = exp(cbind(matrix(1, nnodos, 1), X) %*% t(par[,,j])[,1:(Ncats[j]-1)])
    	}
  	tot= 1+rowSums(z)
  	P[,1]=1/tot
  	for (k in 1:(Ncats[j]-1))
     	P[,(k+1)]=z[,k]/tot

  	PT=cbind(PT,P)
	}

	return(PT)
}



