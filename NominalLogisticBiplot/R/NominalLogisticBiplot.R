# file NominalLogisticBiplot/R/NominalLogisticBiplot.R
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

#Function that calculate the object for the Nominal Logistic Biplot with the 
#estimation for the rows and columns depending on the parameters selected.
  #------------------Parameters-----------------
  #datanom: matrix with the data to do the analysis
  #sFormula: S Formula for selecting variables from a data frame o a matrix
  #numFactors: number of dimensions retained
  #method: Valid values: "EM", "ACM", "MIRT" o "PCOA". This is the method to calculate the row coordinates
  #rotation: "varimax". Type of rotation if the coordinates election is MIRT
  #metfsco : it can be: "EAP", "ML","MAP", o "WLE". Type of method to calculate the MIRT coordinates for the rows.
  #nnodos = 10: number of nodes in the quadrature procedure with EM coordinates method.
  #tol = 1e-04, maxiter = 100, penalization = 0.1, cte, show: items used in the alternated method used in EM procedure to calculate the parameters estimates.
#This function returns a list with the original data matrix and the coordinates for the individuals
#with so many colums like factors.
NominalLogisticBiplot <- function(datanom,sFormula=NULL,numFactors=2,method="EM",
          rotation="varimax",metfsco="EAP",nnodos = 10, tol = 1e-04, maxiter = 100,
          penalization = 0.1,cte=TRUE, initial=1,alfa=1,show=FALSE){       

  #We have to check if datanom is a data frame or a matrix to tackle with sFormula parameter
  if(!(is.null(sFormula))){
    if(is.data.frame(datanom)){
      datanom = model.frame(formula=sFormula,data=datanom)
    }else if(is.matrix(datanom)){
              datanom = model.frame(formula=sFormula,data=as.data.frame(datanom))
              datanom = as.matrix(datanom)
          }else{
            print("It is not posible to use the formula passed as parameter. Data should be a data frame or a matrix")
          }
  }
  if(ncol(datanom) <= numFactors){
    stop("It is not posible to reduce dimension because number of factors is lower than variables in our data set")
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
    xEM = NominalLogBiplotEM(datanom, dim = numFactors, nnodos = nnodos, tol = tol,
                   maxiter = maxiter, penalization = penalization,initial=initial,showResults=show)
    x = xEM$RowCoordinates
    VariableModels = xEM$ColumnModels
    rotation=" "
    metfsco=" "

  }else if (method == "MIRT"){
      varStudy = matrix(0,nRowsdata,numVar)
      for(i in 1:numVar){
        varStudy[,i]= as.numeric(datanom[,i])
      }
      dimnames(varStudy)[[2]]=dimnames(datanom)[[2]][1:numVar]     
      technical = list()
      technical$NCYCLES = maxiter
      nominalHighLow = matrix(0,2,numVar)
      for(i in 1:numVar){
      nominalHighLow[1,i] = as.numeric(max(as.matrix(datanom)[,i]))
      nominalHighLow[2,i] = as.numeric(min(as.matrix(datanom)[,i]))
      }                                                                           
      mod = mirt(varStudy,numFactors,"nominal",rotate = rotation,nominal.highlow=nominalHighLow,technical=technical)           
      coef(mod)
      tabscores<-fscores(mod,full.scores=TRUE,method=metfsco)
      x=as.matrix(tabscores[,1:numFactors])
      VariableModels = CalculateVariableModels(datanom,x, penalization, tol, maxiter,show)
      nnodos=" "
      
  }else if (method == "ACM"){ 
      G = NominalMatrix2Binary(datanom)
    	corr=afc(G,dim=numFactors)
    	x=as.matrix(corr$RowCoordinates[,1:numFactors]*10)       
      VariableModels = CalculateVariableModels(datanom,x, penalization, tol, maxiter,show)
      rotation=" "
      metfsco=" "
      nnodos=" "
   }else{
      datanom.hamm <- NominalDistances(datanom[1:nRowsdata,1:numVar])
     datanom.pco <- PCoA(datanom.hamm, r = numFactors)
      x = as.matrix(datanom.pco$RowCoordinates)
      VariableModels = CalculateVariableModels(datanom,x, penalization, tol, maxiter,show)
      rotation=" "
      metfsco=" "
      nnodos=" "

   }
   dimnames(x)[[1]]=dimnames(datanom)[[1]]    
   dimnames(x)[[2]]=c(1:numFactors)
   
    nominal.logistic.biplot<-list()
    nominal.logistic.biplot$dataSet = dataSet
    nominal.logistic.biplot$RowsCoords = x
    nominal.logistic.biplot$VariableModels = VariableModels
    nominal.logistic.biplot$NumFactors = numFactors
    nominal.logistic.biplot$Method = method
    nominal.logistic.biplot$Rotation = rotation
    nominal.logistic.biplot$Methodfscores = metfsco
    nominal.logistic.biplot$NumNodos = nnodos
    nominal.logistic.biplot$tol = tol
    nominal.logistic.biplot$maxiter = maxiter
    nominal.logistic.biplot$penalization = penalization
    nominal.logistic.biplot$cte = cte
    nominal.logistic.biplot$show = show
    
    class(nominal.logistic.biplot)='nominal.logistic.biplot'
    return(nominal.logistic.biplot)
}

#This function shows a summary of the principal results from the estimation for rows and variables.
#----------------------Parameters
  #object: object with the information needed about all the variables and individuals
summary.nominal.logistic.biplot <- function(object,completeEstim = FALSE,coorInd = FALSE, ...) {
      x = object
    	cat(paste("Nominal Logistic Biplot Estimation", "with Ridge Penalizations", x$penalization, " and logit link"), "\n")
    	cat("n: ", nrow(x$dataSet$datanom), "\n")
    	aic=0
    	bic=0
    	for(i in 1:ncol(x$dataSet$datanom)){
    	     aic = aic + x$VariableModels[,i]$AIC
    	     bic = bic + x$VariableModels[,i]$BIC    	     
    	}
    	cat("AIC: ", aic, "\n")
    	cat("BIC: ", bic, "\n")
    	
	    smmr = matrix(0,8,ncol(x$dataSet$datanom))
      dimnames(smmr)[[1]] = c("loglikelihood","Deviance","AIC","BIC","CoxSnell","Nagelkerke","McFaden","PCC") 
      dimnames(smmr)[[2]] = x$dataSet$ColumNames
    	cat("\nGoodness-of-Fit Statistics for the variables:\n")
 	    for(j in 1:ncol(x$dataSet$datanom)){
    	    smmr[1,j] = x$VariableModels[,j]$logLik
    	    smmr[2,j] = x$VariableModels[,j]$Deviance
    	    smmr[3,j] = x$VariableModels[,j]$AIC
    	    smmr[4,j] = x$VariableModels[,j]$BIC
    	    smmr[5,j] = x$VariableModels[,j]$CoxSnell
    	    smmr[6,j] = x$VariableModels[,j]$Nagelkerke
    	    smmr[7,j] = x$VariableModels[,j]$MacFaden
    	    smmr[8,j] = x$VariableModels[,j]$PercentCorrect
	    }
	    print(smmr)

    	if(coorInd){
    	   cat("\nCoordinates for the individuals: \n\n")
       	 print(x$RowsCoords)
    	}

    	if(completeEstim){
    	    cat("\nCoefficients of the estimated variables: \n")
    	    dimNamesCols = "ind"
    	    dimNamesCols = c(dimNamesCols, "x1")
          if(x$NumFactors > 1){    	    
        	   for (i in 2:x$NumFactors)
                dimNamesCols = c(dimNamesCols, paste("x", i,sep=""))
          }
          dimNamesRows = " "
        	for(i in 1:ncol(x$dataSet$datanom)){
              cat("Variable: ",x$dataSet$ColumNames[i],"\n")
            	coefs <- array(0, c(nrow(x$VariableModels[,i]$beta), x$NumFactors + 1, 4))
            	dimnames(coefs) = list(dimNamesRows,dimNamesCols, c("Estimate coefficients", "Std. Error", "z value", "Pr(>|z|)"))
            	coefs[,,1] <- x$VariableModels[,i]$beta
            	coefs[,,2] <- x$VariableModels[,i]$stderr
            	coefs[,,3] <- x$VariableModels[,i]$beta/x$VariableModels[,i]$stderr
            	coefs[,,4] <- 2 * (1 - pnorm(abs(coefs[,,3])))
            	print(coefs)      	
            	cat("\n\n ----------------------------------------------- \n") 
        	}
    	}
}

#This function plot the row points with the category points according to the options for all the variables.
#----------------------Parameters
  #nlbo: object with the information needed about all the variables and individuals
  #planex,planey: these two parameters keep the plane we choose for the representation  
  #QuitNotPredicted: boolean parameter that indicates if we will be represent on the graph the categories not predicted
  #ReestimateInFocusPlane: boolean parameter to choose if we will reestimate the models in the concrete plane or we will choose the columms of our plane in beta matrix
  #proofMode: boolean variable to plot each variable with its tesselation separate from the others.
              #Otherwise(false value) it shows all the variables in a single window with the legend for identifies each variable.
  #AtLeastR2Nagel <- 0.01   #It establishes the cut value to plot a variable attending to its R2 Nagelkerke value.
  #xlimi,xlimu,ylimi,ylimu: Limits to plot the graph
  #linesVoronoi: boolean variable that specify if lines of each tesselation will be plotted.
  #ShowAxis:boolean variable to choose if we want to see the axis in the plot
  #PlotVars:boolean variable to choose if we want to see the position of the variables
  #PlotInd:boolean variable to choose if we want to see the position of the individuals
  #LabelVar:boolean variable to choose if we want to see the labels of the variables
  #LabelInd:boolean variable to choose if we want to see the labels of the individuals
  #CexInd: parameter to specify the size of the individual points. It could be an array with the cex information for each row
  #CexVar: parameter to specify the size of the variable centers.It could be an array with the cex information for each variable
  #ColorInd: parameter to specify the color of the individual points.It could be an array with the color information for each row
  #ColorVar: parameter to specify the color of the variable centers.It could be an array with the color information for each variable
  #SmartLabels: boolean variable to plot the label text for the individuals to the right or left of the points            
  #PchInd: parameter to specify the symbol of the individual points.It could be an array with the pch information for each row
  #PchVar: parameter to specify the symbol of the variable centers.It could be an array with the pch information for each variable
  #LabelValuesVar: list with the labels text that will be plotted for all the variables
  #ShowResults: boolean parameter to show results from the process of estimation  
plot.nominal.logistic.biplot <- function(x,planex=1,planey=2,QuitNotPredicted=TRUE,
        ReestimateInFocusPlane=TRUE,proofMode=FALSE,AtLeastR2 = 0.01,
        xlimi=-1.5,xlimu=1.5,ylimi=-1.5,ylimu=1.5,linesVoronoi = FALSE,
        ShowAxis = TRUE, PlotVars = TRUE, PlotInd = TRUE, LabelVar = TRUE,
        LabelInd = TRUE,CexInd = NULL, CexVar = NULL, ColorInd = NULL, ColorVar = NULL,
        SmartLabels = FALSE, PchInd = NULL, PchVar = NULL,LabelValuesVar=NULL,ShowResults=FALSE,...) {
  nlbo = x
  
  if(nlbo$NumFactors == 1){
    stop("There is only one factor and it is not posible to represent the biplot on two dimensions.")
  }
  
  if(planex == planey){
    stop("The plane of the biplot is not correct. It should be selected different factors.")
  }
  
  if((planex > nlbo$NumFactors) | (planey > nlbo$NumFactors)){
    stop(paste("With ",nlbo$NumFactors, " factors it is not posible to analyze the plane ",planex,"-",planey,". Please,
                 select another plane.",sep=""))
  }
  
  BLM = BuildTesselationsObject(nlbo,planex,planey,QuitNotPredicted,ReestimateInFocusPlane,ShowResults)
  if(ShowResults){
    WriteMultinomialLogisticBiplot(BLM)
  }
  
  n = nrow(nlbo$dataSet$datanom)
	p = ncol(nlbo$dataSet$datanom)
	RowNames = nlbo$dataSet$RowNames
	VarNames = nlbo$dataSet$ColumNames
	
	DimNames = "Dim 1"
	for (i in 2:nlbo$NumFactors)
     DimNames = c(DimNames, paste("Dim", i))

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
  
	if (!is.null(LabelValuesVar)){
     if (length(LabelValuesVar) > p){
       LabelValuesVar = LabelValuesVar[1:p]    	   
   	   for(i in 1:p){
   	      if(!is.null(nlbo$dataSet$LevelNames[[i]])){
     	      if(length(LabelValuesVar[[i]]) < length(nlbo$dataSet$LevelNames[[i]])){
     	          LabelValuesVar[[i]] = c(1:max(nlbo$dataSet$datanom[,i]))
     	      }else{
          	    LabelValuesVar[[i]] = LabelValuesVar[[i]][1:max(nlbo$dataSet$datanom[,i])]    	      
     	      }
   	      }else{
   	           LabelValuesVar[[i]] = c(1:max(nlbo$dataSet$datanom[,i]))
          }
       }
     }else{
       print("Labels for the values of the variables probably will not be showed as desired because it exists
             variables not specified in the parameter") 
       for(i in 1:length(LabelValuesVar)){
          if(!is.null(nlbo$dataSet$LevelNames[[i]])){
               if(length(LabelValuesVar[[i]]) < length(nlbo$dataSet$LevelNames[[i]])){
         	          LabelValuesVar[[i]] = c(1:max(nlbo$dataSet$datanom[,i]))
         	     }else{
              	    LabelValuesVar[[i]] = LabelValuesVar[[i]][1:max(nlbo$dataSet$datanom[,i])]    	      
         	     }
    	     }else{
    	        LabelValuesVar[[i]] = c(1:max(nlbo$dataSet$datanom[,i]))
    	     }
       }
       if ((p - length(LabelValuesVar)) > 0){
           LabelValuesVarAdd = list()
           for(j in 1:(p - length(LabelValuesVar))){
              if(is.null(nlbo$dataSet$LevelNames[[length(LabelValuesVar) + j]])){
          	    LabelValuesVarAdd[[j]] = c(1:max(nlbo$dataSet$datanom[,length(LabelValuesVar) + j])) 
         	    }else{
         	      LabelValuesVarAdd[[j]] = nlbo$dataSet$LevelNames[[length(LabelValuesVar) + j]]
         	    }       
           }
           LabelValuesVar = c(LabelValuesVar,LabelValuesVarAdd) 
       }
     }
  }else{
  	 LabelValuesVar = list()  
     for(i in 1:p){
        if(is.null(nlbo$dataSet$LevelNames[[i]])){
    	    LabelValuesVar[[i]] = c(1:max(nlbo$dataSet$datanom[,i])) 
   	    }else{
   	      LabelValuesVar[[i]] = nlbo$dataSet$LevelNames[[i]]
   	    }
     }
  }
   
	if (ShowAxis) {
		xaxt = "s"
		yaxt = "s"
	} else {
		xaxt = "n"
		yaxt = "n"
	}
 
  if(proofMode == FALSE){
    dev.new()
    if(PlotInd == TRUE){
      plot(BLM$rowCoords[,1], BLM$rowCoords[,2], cex = CexInd,, col=ColorInd, pch = PchInd, asp=1, xaxt = xaxt, yaxt = yaxt ,xlim=c(xlimi,xlimu),ylim=c(ylimi,ylimu),
          main="Nominal Logistic Biplot", xlab=paste("Axis ",BLM$planex,sep=""), ylab=paste("Axis ",BLM$planey,sep=""))
      if(LabelInd == TRUE){
        if(SmartLabels){
  			   textsmart(cbind(BLM$rowCoords[,1], BLM$rowCoords[,2]), CexPoints = CexInd, ColorPoints = ColorInd)
			  }else{
           text(BLM$rowCoords[,1], BLM$rowCoords[,2],row.names(BLM$rowCoords), cex = CexInd,col=ColorInd,pos=1,offset=0.1)
        }
      }
    }else{
       plot(BLM$rowCoords[,1], BLM$rowCoords[,2], cex = 0,asp=1, xaxt = xaxt, yaxt = yaxt ,xlim=c(xlimi,xlimu),ylim=c(ylimi,ylimu),
          main="Nominal Logistic Biplot", xlab=paste("Axis ",BLM$planex,sep=""), ylab=paste("Axis ",BLM$planey,sep=""))
    }
    if(PlotVars){
      legend("bottomright", legend=VarNames, col= ColorVar,pch=PchVar,cex=0.5)
    }
  }

  if(PlotVars){
    for(l in 1:(BLM$numVar)){
      if(proofMode == TRUE){
        dev.new()
        if(PlotInd == TRUE){
          plot(BLM$rowCoords[,1], BLM$rowCoords[,2], cex = CexInd, col=ColorInd, pch = PchInd,asp=1, xaxt = xaxt, yaxt = yaxt ,xlim=c(xlimi,xlimu),ylim=c(ylimi,ylimu),
              main=paste("Nominal Logistic Biplot for ",BLM$biplotNom[,l]$nameVar,sep=""), xlab=paste("Axis ",BLM$planex,sep=""), ylab=paste("Axis ",BLM$planey,sep=""))
          if((LabelInd) & (SmartLabels)){
      			 textsmart(cbind(BLM$rowCoords[,1], BLM$rowCoords[,2]), CexPoints = CexInd, ColorPoints = ColorInd)
          }else if((LabelInd) & (SmartLabels==FALSE)){
             text(BLM$rowCoords[,1], BLM$rowCoords[,2],row.names(BLM$rowCoords), cex = CexInd,col=ColorInd,pos=1,offset=0.1)             
          } 
        }else{
           plot(BLM$rowCoords[,1], BLM$rowCoords[,2], cex = 0,asp=1, xaxt = xaxt, yaxt = yaxt ,xlim=c(xlimi,xlimu),ylim=c(ylimi,ylimu),
              main="Nominal Logistic Biplot", xlab=paste("Axis ",BLM$planex,sep=""), ylab=paste("Axis ",BLM$planey,sep=""))
        }
      }
      
      if(BLM$biplotNom[,l]$numFit==1){
         Barx=sum(BLM$rowCoords[,1])/nrow(BLM$rowCoords)
          Bary=sum(BLM$rowCoords[,2])/nrow(BLM$rowCoords)
          points(Barx,Bary,pch=PchVar,cex=CexVar,col=ColorVar)
          if(LabelVar){
              if(is.null(LabelValuesVar)){
                	text(Barx,Bary, paste(BLM$biplotNom[,l]$equivFit,"_",BLM$biplotNom[,l]$nameVar ,sep="") , col = ColorVar, cex = CexVar,pos=1,offset=0.1)
              }else{
               	  text(Barx,Bary, LabelValuesVar[[l]][BLM$biplotNom[,l]$equivFit] , col = ColorVar, cex = CexVar,pos=1,offset=0.1)              
              }
           } 
      }else if(BLM$biplotNom[,l]$numFit==2){
              if(max(nlbo$dataSet$datanom[,l]) == 2){
                  x = cbind(BLM$rowCoords[,1],BLM$rowCoords[,2])
                  plot2CategLine(BLM$biplotNom[,l],x,AtLeastR2,line=linesVoronoi,LabelVar=LabelVar,CexVar=CexVar,ColorVar=ColorVar[l],PchVar=PchVar[l],LabValVar=LabelValuesVar[[l]])
               }else if(QuitNotPredicted == TRUE){
                         x = cbind(BLM$rowCoords[,1],BLM$rowCoords[,2])
                         plot2CategLine(BLM$biplotNom[,l],x,AtLeastR2,line=linesVoronoi,LabelVar=LabelVar,CexVar=CexVar,ColorVar=ColorVar[l],PchVar=PchVar[l],LabValVar=LabelValuesVar[[l]])                 
                     }else{
                         plot.voronoiprob(BLM$biplotNom[,l],LabelVar=LabelVar,CexVar=CexVar,ColorVar=ColorVar[l],PchVar=PchVar[l],AtLeastR2=AtLeastR2,lines=linesVoronoi,LabValVar=LabelValuesVar[[l]])
                     }
             }else{
                plot.voronoiprob(BLM$biplotNom[,l],LabelVar=LabelVar,CexVar=CexVar,ColorVar=ColorVar[l],PchVar=PchVar[l],AtLeastR2=AtLeastR2,lines=linesVoronoi,LabValVar=LabelValuesVar[[l]])
             }
    }#end for
  }#end if(PlotVars)
}

