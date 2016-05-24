# file OrdinalLogisticBiplot/R/OrdinalLogisticBiplot.R
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

OrdinalLogisticBiplot <- function(datanom,sFormula=NULL,numFactors=2,method="EM",rotation="varimax",
                          metfsco="EAP",nnodos = 10, tol = 1e-04, maxiter = 100, penalization = 0.1,cte=TRUE,
                          show=FALSE,ItemCurves = FALSE,initial=1,alfa=1){

  #We have to check if datanom is a data frame or a matrix to tackle with sFormula parameter
  if(!(is.null(sFormula))){
    if(is.data.frame(datanom)){
      datanom = model.frame(formula=sFormula,data=datanom)
    }else if(is.matrix(datanom)){
              datanom = model.frame(formula=sFormula,data=as.data.frame(datanom))
              datanom = as.matrix(datanom)
          }else{
            print("It is not posible to use the formula passed as parameter. Data are not a data frame nor matrix")
          }
  }
  if(ncol(datanom) <= numFactors){
    stop("It is not posible to reduce dimension because number of Factors is lower than variables in our data set")
  }
  dataSet = CheckDataSet(datanom)
  #We have now in dataSet structure all our data, the names for the row and colums and de levels if exists
  datanom = dataSet$datanom

  nRowsdata <- dim(datanom)[1]
  nColsdata <- dim(datanom)[2]
  numVar <- ncol(datanom)    
  datanomcats = apply(datanom[,1:numVar], 2, function(x) nlevels(as.factor(x)))

  x <- matrix(0,nRowsdata,numFactors)
  
  if (method == "EM"){
    olb = OrdinalLogBiplotEM(datanom,dim = numFactors, nnodos = nnodos, tol = tol, maxiter = maxiter,
                               penalization = penalization,initial=initial,show=show,alfa=alfa)
    print(olb$ColumnParameters$fit)
 
    coefs = olb$ColumnParameters$coefficients
    thresholds = olb$ColumnParameters$thresholds
    x = as.matrix(olb$RowCoordinates)
    Fitting = olb$ColumnParameters$fit
    LogLik = olb$LogLikelihood        
    FactorLoadings = olb$loadings
    h2 = olb$r2
    rotation=" "
    metfsco=" "

  }else if (method == "MIRT"){ 
      rows = EstimationRowsMIRT(datanom,numFactors = numFactors,metfsco=metfsco,rotation=rotation,maxiter=maxiter)  #Se estima una vez y luego se elige el plano que queremos estudiar en GetOrdinalBiplotObjectMIRT
      fittedMirtModels = CalculateFittingIndicatorsMirt(rows)
      x=as.matrix(rows$estimRows)                   
      coefs = as.matrix(rows$sepCoefMirt$xCoeffic)
      thresholds = rows$sepCoefMirt$indCoeffic
      fittedInd = matrix(0,nColsdata,8)                 
      dimnames(fittedInd)[[1]]= dimnames(datanom)[[2]][1:dim(datanom)[2]]        
      dimnames(fittedInd)[[2]]= c("logLik","Deviance","df","p-value","PCC","CoxSnell","Macfaden","Nagelkerke")  
      sumLogLik = 0
      for(i in 1:nColsdata){
         fittedInd[i,]=c(fittedMirtModels[,i]$logLik,fittedMirtModels[,i]$Deviance,
                            fittedMirtModels[,i]$df,fittedMirtModels[,i]$pval,
                            fittedMirtModels[,i]$PercentClasif,fittedMirtModels[,i]$CoxSnell,
                            fittedMirtModels[,i]$MacFaden,fittedMirtModels[,i]$Nagelkerke) 
         sumLogLik = sumLogLik + fittedMirtModels[,i]$logLik
      }
      Fitting = fittedInd                                             
      LogLik = sumLogLik
      if(numFactors > 1){
         FactorLoadings = rows$summ$rotF
         h2 = rows$summ$h2
      }else{     
         FactorLoadings = rows$summ$F
         h2 = ""   
      }

      rotation=rotation
      metfsco=metfsco
      nnodos=" "
      olb = rows
  }                                                                                                               
  
   if((numFactors > 1)|(method == "EM")){
     FactorLoadingsComm = cbind(FactorLoadings,h2)
     dimnames(FactorLoadingsComm)[[2]] = c(paste(c(rep("F_",numFactors)),c(1:numFactors),sep=""),"Communalities")   
   }else{
        FactorLoadingsComm = FactorLoadings
        dimnames(FactorLoadingsComm)[[2]] = "F_1"
   }
   dimnames(FactorLoadingsComm)[[1]] = dataSet$ColumNames


   dimnames(x)[[1]]=dimnames(datanom)[[1]]    
   dimnames(x)[[2]]=c(1:numFactors)
   dimnames(coefs)[[1]] = dataSet$ColumNames
   dimnames(coefs)[[2]] = paste(c(rep("Dim",numFactors)),c(1:numFactors))
   dimnames(thresholds)[[1]] = dataSet$ColumNames
   dimnames(thresholds)[[2]] = c(1:ncol(thresholds))

    ordinal.logistic.biplot<-list()
    ordinal.logistic.biplot$dataSet = dataSet
    ordinal.logistic.biplot$RowCoords = x
    ordinal.logistic.biplot$Ncats = datanomcats
    ordinal.logistic.biplot$estimObject = olb
    ordinal.logistic.biplot$Fitting = Fitting 
    ordinal.logistic.biplot$coefs = coefs  
    ordinal.logistic.biplot$thresholds = thresholds   
    ordinal.logistic.biplot$NumFactors = numFactors
    ordinal.logistic.biplot$Coordinates = method
    ordinal.logistic.biplot$Rotation = rotation
    ordinal.logistic.biplot$Methodfscores = metfsco
    ordinal.logistic.biplot$NumNodos = nnodos
    ordinal.logistic.biplot$tol = tol
    ordinal.logistic.biplot$maxiter = maxiter
    ordinal.logistic.biplot$penalization = penalization
    ordinal.logistic.biplot$cte = cte
    ordinal.logistic.biplot$show = show
    ordinal.logistic.biplot$ItemCurves =  ItemCurves
    ordinal.logistic.biplot$LogLik =  LogLik
    ordinal.logistic.biplot$FactorLoadingsComm = FactorLoadingsComm    

    class(ordinal.logistic.biplot)='ordinal.logistic.biplot'
    return(ordinal.logistic.biplot)
}

#This function shows a summary of the principal results from the estimation for rows and variables.
#----------------------Parameters
  #x: object with the information needed about all the variables and individuals
summary.ordinal.logistic.biplot <- function(object,data = FALSE,rowCoords = FALSE,coefs = FALSE,loadCommun = FALSE,...) {
      x=object
      if(x$Coordinates == "EM"){
         	cat(paste(" Ordinal Logistic Biplot Estimation ", "with Ridge Penalization :", x$penalization, ", EM algorithm and logit link"), "\n")        
         	cat("\n Percentage of correct classifications,Pseudo R-squared measures and other indicators: \n")
         	print(x$Fitting)
      }else{
         	cat(paste(" Ordinal Logistic Biplot Estimation ", "using MIRT method"), "\n")
          cat(paste("Rotation", x$Rotation), "\n")
          cat(paste("Method of fscores calculation:", x$Methodfscores), "\n")
         	cat("\n Percentage of correct classifications,Pseudo R-squared measures and other indicators: \n")          
         	print(x$Fitting)
      }
      cat(paste("\n Number of factors for the reduced solution:",x$NumFactors,"\n",sep=""))
      cat("\n Number of categories of the variables:\n")
     	print(x$Ncats)

      if(data){                                                     
    	    cat("\n n: ", nrow(x$dataSet$datanom), "\n")
    	}
    	if(rowCoords){
    	   cat("\nCoordinates for the individuals: ","\n")    	
         print(x$RowCoords)
      }
      if(coefs){                                                                    
       	cat("\n Coefficients:\n")
       	print(x$coefs)
       	cat("\n Thresholds:\n")
       	print(x$thresholds)
     	}
     	if(loadCommun){
       	cat("\n Factor Loadings and Communalities:\n")
       	print(x$FactorLoadingsComm)
    	}
}

plot.ordinal.logistic.biplot <- function(x,planex=1,planey=2,AtLeastR2 = 0.01,
        xlimi=-1.5,xlimu=1.5,ylimi=-1.5,ylimu=1.5,margin = 0,
        ShowAxis = TRUE, PlotVars = TRUE, PlotInd = TRUE, LabelVar = TRUE,
        LabelInd = TRUE,CexInd = NULL, CexVar = NULL, ColorInd = NULL, ColorVar = NULL,
        PchInd = NULL, PchVar = NULL,showIIC=FALSE,iicxi=-1.5,iicxu=1.5,
        legendPlot = FALSE,PlotClus = FALSE,Clusters=NULL,chulls = TRUE,centers = TRUE,
        colorCluster = NULL,ConfidentLevel=NULL,addToExistingPlot=FALSE,...) {
  olbo = x
  if(olbo$NumFactors == 1){
    stop("There is only one factor and it is not posible to represent the biplot on two dimensions")
  }
  
  if(planex == planey){
    stop("The plane of the biplot is not correct. It should be selected different factors.")
  }
  
  if((planex > olbo$NumFactors) | (planey > olbo$NumFactors)){
    stop(paste("With ",olbo$NumFactors, " factors it is not posible to analyze the plane ",planex,"-",planey,". Please,
                 select another plane.",sep=""))
  }

  n = nrow(olbo$dataSet$datanom)
	p = ncol(olbo$dataSet$datanom)
	RowNames = olbo$dataSet$RowNames
	VarNames = olbo$dataSet$ColumNames

	DimNames = "Dim 1"
	for (i in 2:olbo$NumFactors)
     DimNames = c(DimNames, paste("Dim", i))

  # Determining sizes and colors of the points
	if (is.null(CexInd)){
		CexInd = rep(0.5, n)
	}else{
     if (length(CexInd) == 1){
       CexInd = rep(CexInd, n)
     }else if(length(CexInd) < n){
             CexInd = rep(0.5, n)
           }else{
             CexInd = CexInd[1:n]
           }
	}

  if (is.null(ColorInd)){
  		ColorInd = rep("black", n)
	}else{
     if (length(ColorInd) == 1){
       ColorInd = rep(ColorInd, n)
     }else if(length(ColorInd) < n){
             ColorInd = rep("black", n)
           }else{
             ColorInd = ColorInd[1:n]
           }
	}

 	if (is.null(PchInd)){
		PchInd = rep(1, n)
	}else{
     if (length(PchInd) == 1){
       PchInd = rep(PchInd, n)
     }else if(length(PchInd) < n){
             PchInd = rep(1, n)
           }else{
             PchInd = PchInd[1:n]
           }
	}

	if (is.null(CexVar)){
		CexVar = rep(0.8, p)
	}else{
     if (length(CexVar) == 1){
       CexVar = rep(CexVar, p)
     }else if(length(CexVar) < p){
             print("It has been specified lower cex values for the variables than variables")
             CexVar = rep(0.8, p)
           }else{
             CexVar = CexVar[1:p]
           }
  }

	if (is.null(PchVar)){
    PchVar = c(0:(p-1))
	}else{
    if (length(PchVar) == 1){
       PchVar = c(0:(p-1))
     }else if(length(PchVar) < p){
             print("It has been specified lower pch values for the variables than variables")
             PchVar = c(0:(p-1))
           }else{
             PchVar = PchVar[1:p]
           }
  }
  if (is.null(ColorVar)){
		ColorVar = colors()[20 + 2*c(1:p)]
	}else{
     if (length(ColorVar) == 1){
       ColorVar = colors()[20 + 2*c(1:p)]
     }else if(length(ColorVar) < p){
             print("It has been specified lower color values for the variables than variables")
             ColorVar = colors()[20 + 2*c(1:p)]
           }else{
             ColorVar = ColorVar[1:p]
           }
  }

	if (ShowAxis) {
		xaxt = "s"
		yaxt = "s"
	} else {
		xaxt = "n"
		yaxt = "n"
	}

  if ((margin < 0) | (margin > 0.3))
		margin = 0

  xmin= xlimi - (xlimu - xlimi) * margin
  xmax= xlimu + (xlimu - xlimi) * margin
  ymin= ylimi - (ylimu - ylimi) * margin
  ymax= ylimu + (ylimu - ylimi) * margin

  if(PlotInd == TRUE){
    if(addToExistingPlot == TRUE){
      points(olbo$RowCoords[,planex], olbo$RowCoords[,planey], pch = PchInd, col=ColorInd, cex = CexInd)
    }else{
      dev.new()
      plot(olbo$RowCoords[,planex], olbo$RowCoords[,planey], cex = CexInd, col=ColorInd, pch = PchInd, asp=1, xaxt = xaxt, yaxt = yaxt ,xlim=c(xmin,xmax),ylim=c(ymin,ymax),
        main="Ordinal Logistic Biplot", xlab=paste("Axis ",planex,sep=""), ylab=paste("Axis ",planey,sep=""))        
    }
    
    if(LabelInd == TRUE){
         text(olbo$RowCoords[,planex], olbo$RowCoords[,planey],row.names(RowNames), cex = CexInd,col=ColorInd,pos=1,offset=0.1)
    }
  }else{
    if(!addToExistingPlot){
      plot(olbo$RowCoords[,planex], olbo$RowCoords[,planey], cex = 0,asp=1, xaxt = xaxt, yaxt = yaxt ,xlim=c(xmin,xmax),ylim=c(ymin,ymax),
        main="Ordinal Logistic Biplot", xlab=paste("Axis ",planex,sep=""), ylab=paste("Axis ",planey,sep=""))
    }
  }

  if(olbo$Coordinates == "EM"){
      olb = olbo$estimObject
      catOrdBiplotPenal = GetOrdinalBiplotObjectPenal(olbo$dataSet$ColumNames,olb,planex,planey)
  }else{
      if(olbo$Coordinates == "MIRT"){
           rows = olbo$estimObject
           catOrdBiplot = GetOrdinalBiplotObjectMIRT(rows,planex,planey)
      }else{
          stop("Coordinates for the items has not been specified.")
      }
  }

  if(PlotVars){
      if(olbo$Coordinates == "EM"){
        D = 1
        plot.ordinalBiplotPenal(catOrdBiplotPenal,olbo$NumFactors,D,planex,planey,xi=xlimi,xu=xlimu,
              yi=ylimi,yu=ylimu,margin = margin, CexVar = CexVar,ColorVar = ColorVar,
              PchVar = PchVar,levelsVar = olbo$dataSet$LevelNames)
      }else{
        D = 1.702
        plot.ordinalBiplotPenal(catOrdBiplot,olbo$NumFactors,D,planex,planey,xi=xlimi,xu=xlimu,
              yi=ylimi,yu=ylimu,margin = margin, CexVar = CexVar,ColorVar = ColorVar,
              PchVar = PchVar,levelsVar = olbo$dataSet$LevelNames)
      }
      if(legendPlot){
        legend("bottomright", legend=VarNames, col= ColorVar,pch=PchVar,cex=0.7)
      }
  }
  
  if (PlotClus) {
    if (is.null(Clusters)){
      Clusters=as.factor(ones(c(n,1)))
    }else{
      Clusters=Clusters
    }
    A = cbind(olbo$RowCoords[,planex],olbo$RowCoords[,planey])
    if(!is.null(colorCluster)){
        col = colorCluster
    }else{
        col = colors()[20 + 2*c(1:length(levels(Clusters)))]
    }
    if(!is.null(ConfidentLevel)){
      if((ConfidentLevel < 0) | (ConfidentLevel > 1)){
         stop("The ConfidentLevel parameter has not a valid value. It should be between 0 and 1") 
      }
    }
    PlotClusters(A, Clusters, colors = col, chulls = chulls,centers = centers,ConfidentLevel=ConfidentLevel)
  }
  
  if(showIIC){
    for(v in 1:ncol(olbo$dataSet$datanom)){
      nameVariable = VarNames[v]
      if(olbo$Coordinates == "EM"){
          coeffic = catOrdBiplotPenal$matBiplot[,v]$coef
          slopeort = catOrdBiplotPenal$matBiplot[,v]$slope
          D = 1
      }else{
          if(olbo$Coordinates == "MIRT"){
              coeffic = catOrdBiplot$matBiplot[,v]$coef
              slopeort = catOrdBiplot$matBiplot[,v]$slope
              D = 1.702
          }else{
              stop("Coordinates for the items has not been specified.")
          }        
      }
      numcat = max(olbo$dataSet$datanom[,v])
      xi = iicxi
      xu = iicxu
      plotCurvesCategoriesVariable(coeffic,slopeort,D,numcat,nameVariable,xi,xu,planex,planey)
    }
  }

}
