# file OrdinalLogisticBiplot/R/auxLibrary_OLB.R
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

GetOrdinalBiplotObjectMIRT <- function(modelMirt,planex=1,planey=2){

    nColsdataF <- ncol(modelMirt$dataFactor)
    varStudy = modelMirt$dataFactor

    matBiplot=0
    for(nvarColum in 1:nColsdataF){
      numcat = max(varStudy[,nvarColum])
      coeffic=modelMirt$coefMirt[[nvarColum]]
      numFactors = length(coeffic)- (numcat - 1)
      coeffic = c(modelMirt$coefMirt[[nvarColum]][1:numFactors],(-1)*modelMirt$coefMirt[[nvarColum]][(numFactors+1):length(coeffic)])
      D=1.702
      nameVariable = dimnames(varStudy)[[2]][nvarColum]
      ordBip = CalculateOrdinalBiplotGeneral(nameVariable,numcat,coeffic,planex,planey,modelMirt$numFactors,D)
      matBiplot=cbind(matBiplot,ordBip)
    }
    matBiplot=matBiplot[,2:(nColsdataF+1)]

    result = list()
    result$coefMirt = modelMirt$coefMirt
    result$sepCoefMirt = modelMirt$sepCoefMirt
    result$numFactors = modelMirt$numFactors
    result$rotation = modelMirt$rotation
    result$metfsco = modelMirt$metfsco
    result$planex = planex
    result$planey = planey
  	result$scores = modelMirt$estimRows
    result$matBiplot = matBiplot
    result$summaryMirt = modelMirt$summ
  	class(result) <- "CategOrdBiplot"
  	return(result)
}

GetOrdinalBiplotObjectPenal <- function(ColumNames,olb,planex=1,planey=2){

    nColsdataF <- nrow(olb$Ncats)
    numFactors = ncol(olb$RowCoordinates)

    matBiplot=0
    for(nvarColum in 1:nColsdataF){
      numcat = olb$Ncats[nvarColum,1]
      coeffic = ExtractCoefficientsPenal(olb,nvarColum,numFactors,numcat)
      D = 1
      nameVariable = ColumNames[nvarColum]
      ordBipVar = CalculateOrdinalBiplotGeneral(nameVariable,numcat,coeffic,planex,planey,numFactors,D)

      matBiplot=cbind(matBiplot,ordBipVar)
    }
    matBiplot=matBiplot[,2:(nColsdataF+1)]

    result = list()
    result$models = olb
  	result$matBiplot = matBiplot
  	class(result) <- "CategOrdBiplotPen"
  	return(result)
}

EstimationRowsMIRT <- function(dataFile,numFactors=2,metfsco="EAP",rotation="varimax",maxiter=100)
{
  nRowsdataF <- dim(dataFile)[1]
  nColsdataF <- dim(dataFile)[2]
  varStudy = matrix(0,nRowsdataF,nColsdataF)
  for(i in 1:nColsdataF){
    varStudy[,i]= factor(dataFile[,i])
  }
  dimnames(varStudy)[[2]]=dimnames(dataFile)[[2]]
  datanomcats = apply(varStudy, 2, function(x) nlevels(as.factor(x)))
  
  technical = list()
  technical$NCYCLES = maxiter
  
  (modF<-mirt(varStudy,numFactors,rotate=rotation,,technical=technical))
  summ = summary(modF)
  
  tabscores<-fscores(modF,full.scores=TRUE,method=metfsco)
  xSubi = as.matrix(tabscores[,1:numFactors])
  catsVarMax = max(datanomcats)
  sepCoefMirt = ExtractMirtCoefficients(datanomcats,cmodF = coef(modF),numFactors = numFactors)
  
  result = list()
  result$estimRows = xSubi
  result$dataFactor = varStudy
  result$rotation = rotation
  result$metfsco = metfsco
  result$numFactors = numFactors
  result$coefMirt = coef(modF)
  result$sepCoefMirt = sepCoefMirt
  result$summ = summ
  class(result) <- "EstimationRowsMIRT"
  return(result)
  
}



zeros <- function(dims){
   if(is.vector(dims)){ 
     if(length(dims) < 2){
        stop("Zeros function parameter with only one dimension.")
     }else{
        if(length(dims) > 2){
            print("Too much dimensions passed to zeros function. It will choose the first two")
            dims = dims[1:2]
            ret = matrix(0,dims[1],dims[2])
        }else{ 
            ret = matrix(0,dims[1],dims[2])
        }
     }
   }else{
      stop("Zeros function must receive a vector parameter")
   }
   return (ret)
}

ones <- function(dims){
   if(is.vector(dims)){ 
     if(length(dims) < 2){
        stop("Ones function parameter with only one dimension.")
     }else{
        if(length(dims) > 2){
            print("Too much dimensions passed to Ones function. It will choose the first two")
            dims = dims[1:2]
            ret = matrix(1,dims[1],dims[2])
        }else{ 
            ret = matrix(1,dims[1],dims[2])
        }
     }
   }else{
      stop("Ones function must receive a vector parameter")
   }
   return (ret)
}


plot.ordinalBiplotPenal <- function(catOrdBiplotPenal,numFactors,D=1,planex=1,planey=2,xi=-3.5,xu=3.5,yi=-3.5,yu=3.5,
       margin = 0, CexVar,ColorVar, PchVar,levelsVar) {

    nColsdataF <- ncol(catOrdBiplotPenal$matBiplot)

    for(nvarColum in 1:nColsdataF){
      ordBipVar = catOrdBiplotPenal$matBiplot[,nvarColum]
      plot.ordBipVariable(ordBipVar,D,planex,planey,xi,xu,yi,yu,margin,numFactors,CexVar[nvarColum],
                ColorVar[nvarColum],PchVar[nvarColum],levelsVar = levelsVar[[nvarColum]])        
    }
}

plot.ordBipVariable <- function(ordBipVar,D,planex,planey,xi,xu,yi,yu,margin,
                                numFactors,CexVar,ColorVar,PchVar,levelsVar = NULL){

   orderUp = TRUE
   numcat = ordBipVar$numcat
   abline(0,ordBipVar$slope,lty=1,col=ColorVar)

   pointsBoxVarName = pointsIntersectLineWindow(ordBipVar$slope,xi,xu,yi,yu)

   if(ordBipVar$order == TRUE){
      xi = min(ordBipVar$pointsc[,1]-1)
      xu = max(ordBipVar$pointsc[,1]+1)
      yi = min(ordBipVar$pointsc[,2]-1)
      yu = max(ordBipVar$pointsc[,2]+1)
   }else{
      xi = min(ordBipVar$pointprob[,2]-1)
      xu = max(ordBipVar$pointprob[,2]+1)
      yi = min(ordBipVar$pointprob[,3]-1)
      yu = max(ordBipVar$pointprob[,3]+1)
   }
   pointsBox = pointsIntersectLineWindow(ordBipVar$slope,xi,xu,yi,yu)
   ang = atan(ordBipVar$slope) * 180/pi

   if(ordBipVar$order==TRUE){
     orderUp = isOrderCurvesAscent(ordBipVar,D,planex,planey)
     for(k in 1:(numcat -1)){
        points(ordBipVar$pointsc[k,1],ordBipVar$pointsc[k,2],pch=PchVar,col=ColorVar,cex=CexVar)
     }

     pointsJoin = rbind(pointsBox,ordBipVar$pointsc)
     orderPointsJoin = pointsJoin[order(pointsJoin[,1]),]
     MedPoints = CalculateMediaBetPoints(orderPointsJoin)

     for(k in 1:numcat){
        if(!is.null(levelsVar)){
          if(orderUp == TRUE){
            labelText = levelsVar[k]
          }else{
            labelText = levelsVar[numcat-k+1]
          }
        }else{
          if(orderUp == TRUE){
            labelText = k
          }else{
            labelText = numcat-k+1
          }
        }
        if((ordBipVar$slope < 1) &
                 (ordBipVar$slope > (-1))){
          text(MedPoints[k,1],MedPoints[k,2] + 0.02,labelText,pos=3,offset=0.1, srt = ang,cex=CexVar,col=ColorVar)
          if(k==numcat){
            text(pointsBoxVarName[1,1],pointsBoxVarName[1,2],ordBipVar$var,pos=1,offset=0.1, srt = ang,cex=CexVar,col=ColorVar)
          }
        }else{
          if(ordBipVar$slope > 1){
            markerpos = 2
          }else{
            markerpos = 4
          }
          text(MedPoints[k,1],MedPoints[k,2],pos=markerpos,offset=0.4, srt = ang,labelText,cex=CexVar,col=ColorVar)
          if(k==numcat){
            text(pointsBoxVarName[1,1],pointsBoxVarName[1,2],ordBipVar$var,pos=1,offset=0.1, srt = ang,cex=CexVar,col=ColorVar)
          }
        }
     }#end for

   }else{
       for(m in 1:(numcat-1)){
         if(ordBipVar$pointprob[m,1] > 0){
           numpoints = m
         }
       }
       pointprobJoin = rbind(pointsBox,ordBipVar$pointprob[1:numpoints,2:3])
       orderPointProbJoin = pointprobJoin[order(pointprobJoin[,1]),]
       orderPointProbJoinCateg = cbind(orderPointProbJoin,0,0)
       for(s in 2:(nrow(orderPointProbJoin)-1)){
          orderPointProbJoinCateg[s,3] = ordBipVar$pointprob[which(ordBipVar$pointprob[,2] == orderPointProbJoin[s,1]),4]
          orderPointProbJoinCateg[s,4] = ordBipVar$pointprob[which(ordBipVar$pointprob[,2] == orderPointProbJoin[s,1]),5]
       }

       Categories = matrix(0,nrow(orderPointProbJoin)-1,1)
       if(orderPointProbJoinCateg[2,3] == 1){
          Categories[1,1] = 1
          for(n in 2:(nrow(orderPointProbJoin)-1)){
            Categories[n,1] = orderPointProbJoinCateg[n,4]
          }
       }else{
          Categories[1,1] = orderPointProbJoinCateg[2,4]
          for(n in 2:(nrow(orderPointProbJoin)-1)){
            Categories[n,1] = orderPointProbJoinCateg[n,3]
          }
       }

       MedPoints = CalculateMediaBetPoints(orderPointProbJoin)
       MedPointsCateg = cbind(MedPoints,Categories)
       for(n in 1:numpoints){
         points(ordBipVar$pointprob[n,2],ordBipVar$pointprob[n,3],pch=PchVar,col=ColorVar,cex=CexVar)
       }
       for(k in 1:(numpoints+1)){
          if((ordBipVar$slope < 1) &
                 (ordBipVar$slope > (-1))){
            if(!is.null(levelsVar)){
              text(MedPointsCateg[k,1],MedPointsCateg[k,2]+0.02,levelsVar[MedPointsCateg[k,3]],pos=3,offset=0.1, srt = ang,cex=CexVar,col=ColorVar)
            }else{
              text(MedPointsCateg[k,1],MedPointsCateg[k,2]+0.02,MedPointsCateg[k,3],pos=3,offset=0.1, srt = ang,cex=CexVar,col=ColorVar)
            }
            if(k==(numpoints+1)){
              text(pointsBoxVarName[1,1],pointsBoxVarName[1,2],ordBipVar$var,pos=1,offset=0.1, srt = ang,cex=CexVar,col=ColorVar)
            }
          }else{
            if(ordBipVar$slope > 1){
              markerpos = 2
            }else{
              markerpos = 4
            }
            if(!is.null(levelsVar)){
              text(MedPointsCateg[k,1],MedPointsCateg[k,2],levelsVar[MedPointsCateg[k,3]],pos=markerpos,offset=0.4, srt = ang,cex=CexVar,col=ColorVar)
            }else{
              text(MedPointsCateg[k,1],MedPointsCateg[k,2],MedPointsCateg[k,3],pos=markerpos,offset=0.4, srt = ang,cex=CexVar,col=ColorVar)
            }
            if(k==(numpoints+1)){
              text(pointsBoxVarName[1,1],pointsBoxVarName[1,2],ordBipVar$var,pos=1,offset=0.1, srt = ang,cex=CexVar,col=ColorVar)
            }
          }
       }
   }
}

isOrderCurvesAscent <- function(ordBipVar,D=1,planex=1,planey=2){
    numFactors = length(ordBipVar$coef) - ordBipVar$numcat + 1
    numcat = ordBipVar$numcat
    coeffic = ordBipVar$coef
    if(ordBipVar$order == TRUE){
      x = max(ordBipVar$pointsc[,1]) + 1
    }else{
      x = max(ordBipVar$pointprob[,2]) + 1
    }
    yFirstCurve = 1/(1 + exp(D*(coeffic[planex]*x + coeffic[planey]*(ordBipVar$slope*x) - coeffic[numFactors + 1])))
    yLastCurve = 1/(1 + exp(-D*(coeffic[planex]*x + coeffic[planey]*(ordBipVar$slope*x) - coeffic[numFactors + numcat - 1])))
    if(yFirstCurve > yLastCurve){
        result = FALSE
     }else{
        result = TRUE
    }
    return (result)
}

plotCurvesCategoriesVariable <- function(coeffic,slopeort,D,numcat,nameVariable,xi,xu,planex,planey){
    dev.new()
    x<-seq(xi,xu,length=1000)
    y = slopeort * x
    z = sqrt(x^2 + y^2)
    escProd = x*coeffic[planex] + y*coeffic[planey]
    signPos = as.numeric(escProd >= 0)
    signNeg = as.numeric(escProd < 0)
    signZ = signPos - signNeg
    numFactors = length(coeffic) - numcat + 1
    abcissa = z*signZ
    fx = 1/(1 + exp(D*(coeffic[planex]*x + coeffic[planey]*(slopeort*x) - coeffic[numFactors + 1])))
    plot(abcissa,fx,type="l",cex=0.1,xlim=c(xi,xu),ylim=c(0,1),main=paste("Curves of the variable:", nameVariable,sep=""))

    if(numcat > 2){
      for(s in 2:(numcat - 1)){
          y = (1/(1 + exp(D*(coeffic[planex]*x + coeffic[planey]*(slopeort*x) - coeffic[numFactors + s])))) -
              (1/(1 + exp(D*(coeffic[planex]*x + coeffic[planey]*(slopeort*x) - coeffic[numFactors + s - 1]))))
          points(z*signZ,y,type="l")
      }
    }

    ylast = 1- 1/(1 + exp(D*(coeffic[planex]*x + coeffic[planey]*(slopeort*x) - coeffic[numFactors + numcat - 1])))
    points(z*signZ,ylast,type="l")

}

testit <- function(x){
  p1 <- proc.time()
  Sys.sleep(x)
  proc.time() - p1
}

CalculateMediaBetPoints <- function(orderPointsJoin){
    rows = nrow(orderPointsJoin)
    mediaP = matrix(0,rows-1,2)
    for(i in 1:(rows -1)){
      mediaP[i,1] = (orderPointsJoin[i,1] + orderPointsJoin[i+1,1])/2
      mediaP[i,2] = (orderPointsJoin[i,2] + orderPointsJoin[i+1,2])/2
    }

    return(mediaP)
}

pointsIntersectLineWindow <- function(slope,xi,xu,yi,yu){
    pointsBox = matrix(0,2,2)

    if(slope == 0){
      pointsBox[1,1] = xu
      pointsBox[1,2] = 0
      pointsBox[2,1] = xi
      pointsBox[2,2] = 0
    }else{
      cutMatrix = matrix(0,4,2)
      cutMatrix[1,1] = xi
      cutMatrix[1,2] = xi * slope
      cutMatrix[2,1] = xu
      cutMatrix[2,2] = xu * slope
      cutMatrix[3,1] = yi/slope
      cutMatrix[3,2] = yi
      cutMatrix[4,1] = yu/slope
      cutMatrix[4,2] = yu
      oCutMatrix = cutMatrix[order(cutMatrix[,1]),]
      pointsBox = oCutMatrix[2:3,]
    }
    return (pointsBox)
}

Eq2gSolve <- function(a,b,c){
    solns<-matrix(rep(NA,length(a)))
    solns<-cbind(solns,solns)
    colnames(solns)<-c("soln 1","soln 2")
    if(a==0 && b!=0){
      solns[1,1]= (-1)*c/b
      solns[1,2]= (-1)*c/b;
    }else if(a==0 && b==0){
      print("Coefficients a and b of the polinomial are zero")
    }else if((b^2 -4*a*c) < 0){
      print("Polinomial has complex roots")
    }else{
      solns[1,1]<-((-1)*b + sqrt(b^2 - 4*a*c))/(2*a)
      solns[1,2]<-((-1)*b - sqrt(b^2 - 4*a*c))/(2*a)
    }
    solns
}


ExtractCoefficientsPenal <- function(olb,nvarColum,numFactors,numcat)
{
   coeffic = matrix(0,1,numFactors + numcat - 1)
   for(j in 1:numFactors){
      coeffic[j] = olb$ColumnParameters$coefficients[nvarColum,j]
   }
   for(k in 1:(numcat - 1)){
      coeffic[numFactors + k] = (1) * olb$ColumnParameters$thresholds[nvarColum,k]
   }

   return(coeffic)
}

ExtractMirtCoefficients <- function(datanomcats,cmodF,numFactors)
{
   numVar = length(cmodF) - 1
   catsVarMax = max(datanomcats)
   
   matrixMirt = matrix(0,numVar,numFactors + catsVarMax - 1)
   for(i in 1:(length(cmodF) - 1)){
      if(datanomcats[i]==2){                 
          matrixMirt[i,] = cmodF[[i]][1:(length(cmodF[[i]])-2)]
      }else{
          if(length(cmodF[[i]]) < length(matrixMirt[i,])){
             matrixMirt[i,1:length(cmodF[[i]])] = cmodF[[i]]
          }else{
             matrixMirt[i,] = cmodF[[i]]
          }
      }
   }
   xCoeffic = matrixMirt[,1:numFactors]
   indCoeffic = matrixMirt[,(numFactors + 1):ncol(matrixMirt)]

   result = list()
   result$xCoeffic = xCoeffic
   result$indCoeffic = as.matrix(indCoeffic)
   return(result)

}

CalculateOrdinalBiplotGeneral <- function(nameVariable,numcat,coeffic,planex,planey,numFactors,D){
   orderPoints=TRUE     
   if((coeffic[planex] == 0) & (coeffic[planey] == 0)){
      stop(paste("Coefficients estimated for the variable ",nameVariable," on the plane ",planex,"-",planey," are zero and it is not posible to calculate the biplot",sep=""))
   }

   if(coeffic[planey] == 0){
       slope = NA 
   }else{
       slope = (-1)*coeffic[planex]/coeffic[planey]
   }
   ordorig=matrix(0,1,numcat-1)

   if(numcat > 2){
     if((exp((-1)*D*coeffic[numFactors+1])- 2*exp((-1)*D*coeffic[numFactors+2])) <= 0){
        ordorig[1]= NA
     }else{
        ordorig[1] = ((-1)/D)*log(exp((-1)*D*coeffic[numFactors+1])- 2*exp((-1)*D*coeffic[numFactors+2]))     
     }
     if((numcat-2) > 1){
       for(j in 2:(numcat-2)){
            num= exp((-1)*D*coeffic[j-1+numFactors])-2*exp((-1)*D*coeffic[j+numFactors])+
                    exp((-1)*D*coeffic[j+1+numFactors])
           denom = exp((-1)*D*(coeffic[j-1+numFactors]+coeffic[j+numFactors])) +
                        (-2)*exp((-1)*D*(coeffic[j-1+numFactors]+coeffic[j+1+numFactors]))+
                          exp((-1)*D*(coeffic[j+1+numFactors]+coeffic[j+numFactors]))
           if((num/denom) < 0){
              ordorig[j]= NA  
           }else{
              ordorig[j]= (1/D)*log(num/denom)
           }                          
       }
     }
     num= exp((-1)*D*coeffic[numcat -2 +numFactors])-2*exp((-1)*D*coeffic[numcat -1 +numFactors]);
     denom = exp((-1)*D*(coeffic[numcat -2+numFactors]+coeffic[numcat -1+numFactors]));
     if((num/denom) < 0){
        ordorig[numcat -1]= NA   
     }else{
        ordorig[numcat -1]= (1/D)*log(num/denom)     
     }                               
   }else if(numcat == 2){
            ordorig[1] = (coeffic[1+numFactors]+coeffic[numcat-1+numFactors])/2
         }else{
            stop("There is a variable with the same value for all the items. Revise the data set.")
         }
   if(!is.na(slope)){
      if(slope == 0){
          slopeort = NA
      }else{
          slopeort = -1/slope
      }      
   }else{
      slopeort = 0
   }

   pointprob=matrix(0,numcat-1,5)

   pointsc=matrix(0,numcat -1,2)
   for(k in 1:(numcat -1)){
      if(!is.na(ordorig[k])){
        if(coeffic[planey] == 0){
            pointsc[k,1] = ordorig[k]/coeffic[planex]
            pointsc[k,2] = 0            
        }else{
          if(coeffic[planex] == 0){
              pointsc[k,1] = 0
              pointsc[k,2] = ordorig[k]/coeffic[planey]            
          }else{
              pointsc[k,1]= ordorig[k]/((slopeort-slope)*(coeffic[planey]))
              pointsc[k,2]= slopeort*pointsc[k,1]     
          }
        }
      }else{
        pointsc[k,1] = NA
        pointsc[k,2] = NA
      }
   }

   if(length(sort(ordorig)) != length(ordorig)){
       pointprobFull = PointsCurveIntersect(coeffic,numcat,planex,planey,D)
       pointprob = PointsBiplotCurves(pointprobFull,numcat)
       orderPoints=FALSE
   }else{
       orddescendent = TRUE  
       if(coeffic[planex] == 0){
          compareColumn = 2
       }else{
          compareColumn = 1
       }
       ordPointsc = sort(pointsc[,compareColumn])
       for(p in 1:(numcat -1)){
          if(pointsc[p,compareColumn] !=  ordPointsc[numcat-p]){
              orddescendent = FALSE
          }
       }
       if((all(ordPointsc == t(pointsc[,compareColumn]))==FALSE)
                                        &&(orddescendent == FALSE)){
           orderPoints=FALSE
           pointprobFull = PointsCurveIntersect(coeffic,numcat,planex,planey,D)
           pointprob = PointsBiplotCurves(pointprobFull,numcat)
       }
   }
   
   vectcos = matrix(0,1,numFactors)
   denomCosines = 0
   for(n in 1:numFactors){
      denomCosines = denomCosines + coeffic[n]^2
   }
   for(n in 1:numFactors){
      vectcos[n]=abs(coeffic[n])/sqrt(denomCosines)
   }

	result = list()
	result$var = nameVariable
	result$cosines = vectcos
	result$numcat = numcat
  result$coef = coeffic
	result$slope = slopeort
	result$order = orderPoints
	result$pointsc = pointsc
	result$pointprob = pointprob
	class(result) <- "ordinalBiplotPOLR"
	return(result)
}

PointsCurveIntersect <- function(coeffic,numcat,planex,planey,D){
  numFactors = length(coeffic) - (numcat - 1)
  pointprobFull = matrix(0,numcat*(numcat-1)/2,5)
  posInsert = 0
  
   for(i in 1:(numcat -1)){
    for(j in (i + 1):numcat){
         if(i == 1){
            if(j == numcat){
                #case 1-n
                 posInsert = posInsert + 1
                 pointprobFull[posInsert,1]= 1/(1+exp(D*((coeffic[numcat-1+numFactors] - coeffic[1+numFactors])/2)))
                 pointprobFull[posInsert,2]= (coeffic[planex]*(coeffic[numcat-1+numFactors]+coeffic[1+numFactors]))/(2*(coeffic[planex]^2+coeffic[planey]^2))
                 pointprobFull[posInsert,3]= (coeffic[planey]*(coeffic[numcat-1+numFactors]+coeffic[1+numFactors]))/(2*(coeffic[planex]^2+coeffic[planey]^2))
                 pointprobFull[posInsert,4] = 1
                 pointprobFull[posInsert,5] = numcat
            }else{
                 if(j==2){    
                      #case 1-2
                      if((exp((-1)*D*coeffic[numFactors+1])
                                    - 2*exp((-1)*D*coeffic[numFactors+2]))>0){
                        posInsert = posInsert + 1 
                        oo12 = ((-1)/D)*log(exp((-1)*D*coeffic[numFactors+1])
                                    - 2*exp((-1)*D*coeffic[numFactors+2]))
                        corteProb12 = 1/(1+exp(D*(oo12 - coeffic[1+numFactors])))
                        pointprobFull[posInsert,1]= corteProb12
                        pointprobFull[posInsert,2]= (coeffic[planex]*oo12)/(coeffic[planex]^2+coeffic[planey]^2)
                        pointprobFull[posInsert,3]= (coeffic[planey]*oo12)/(coeffic[planex]^2+coeffic[planey]^2)
                        pointprobFull[posInsert,4]= 1
                        pointprobFull[posInsert,5]= 2                                             
                      }else{
                        print(paste("Curves 1-",j," do not intersect at any point.",sep=""))                                      
                      }
                 }else{
                    #case 1-j, with j#numcat, and j > 1 
                     categ = j                     
                      #(e(-D(d1+di-1))-e(-D(d1+di))-e(-D(di-1+di)))x^2-2e(-D*di)x -1=0
                      a= exp((-1)*D*(coeffic[1+numFactors]+coeffic[categ-1+numFactors]))-exp((-1)*D*(coeffic[1+numFactors]+coeffic[categ+numFactors]))-exp((-1)*D*(coeffic[categ-1+numFactors]+coeffic[categ+numFactors]))
                      b= (-2)*exp((-1)*D*coeffic[categ+numFactors])
                      c= -1
                      xsolvect1i=matrix(-99,1,2)
                      xsolvect1i= Eq2gSolve(a,b,c)
                      if(xsolvect1i[1,1]>0){
                        xsolve=xsolvect1i[1,1]
                      }else if(xsolvect1i[1,2]>0){
                        xsolve=xsolvect1i[1,2]
                      }else{
                        print(paste("Curves 1-",categ," do not intersect at any point",sep=""))              
                        xsolve = -1
                      }
                       if(xsolve > 0){
                          posInsert = posInsert + 1
                          oo1i=(1/D)*log(xsolve)
                          corteProb = 1/(1+(exp(D*(oo1i-coeffic[1+numFactors])))) 
                          pointprobFull[posInsert,1]= corteProb
                          pointprobFull[posInsert,2]= (coeffic[planex]*oo1i)/(coeffic[planex]^2+coeffic[planey]^2)
                          pointprobFull[posInsert,3]= (coeffic[planey]*oo1i)/(coeffic[planex]^2+coeffic[planey]^2)
                          pointprobFull[posInsert,4]= 1
                          pointprobFull[posInsert,5]= categ
                      }
                 }
            }
         }else{
            if(j == numcat){
              if(i == (numcat -1)){
                #case (n-1)-n                
                oonm1n = (exp((-1)*D*coeffic[numcat - 2 + numFactors])- 2*exp((-1)*D*coeffic[numcat - 1 +numFactors]))/
                              (exp((-1)*D*(coeffic[numcat - 1 + numFactors] + coeffic[numcat - 2 + numFactors])))
                if(oonm1n > 0){
                    posInsert = posInsert + 1
                    corteProbnm1n = 1- 1/(1+exp((-1)*D*coeffic[numcat - 1+numFactors])*oonm1n)                
                    pointprobFull[posInsert,1]= corteProbnm1n
                    pointprobFull[posInsert,2]= (coeffic[planex]*(1/D)*log(oonm1n))/(coeffic[planex]^2+coeffic[planey]^2)
                    pointprobFull[posInsert,3]= (coeffic[planey]*(1/D)*log(oonm1n))/(coeffic[planex]^2+coeffic[planey]^2)
                    pointprobFull[posInsert,4]= numcat - 1
                    pointprobFull[posInsert,5]= numcat
                }else{
                    print(paste("Curves ",numcat-1,"-",numcat," do not intersect at any point",sep=""))
                }
              }else{
                #case i-n
                categ  = i                
                #(e(-D(di-1+di+dn-1))x^2+2e(-D*(di+dn-1))x+(e(-Ddn-1)+e(-Ddi)-e(-Ddi-1)))=0
                a= exp((-1)*D*(coeffic[categ-1+numFactors]+coeffic[categ+numFactors]+
                              coeffic[numcat-1+numFactors]))
                b= 2*exp((-1)*D*(coeffic[categ+numFactors]+coeffic[numcat-1+numFactors]))
                c= exp((-1)*D*(coeffic[categ+numFactors])) +
                      exp((-1)*D*(coeffic[numcat-1+numFactors])) -
                        exp((-1)*D*(coeffic[categ-1+numFactors]))
                xsolvectin=matrix(-99,1,2)
                xsolvectin= Eq2gSolve(a,b,c)
                if(xsolvectin[1,1]>0){
                  xsolvein=xsolvectin[1,1]
                }else if(xsolvectin[1,2]>0){
                  xsolvein=xsolvectin[1,2]
                }else{
                  print(paste("Curves ",i,"-",numcat," do not intersect at any point",sep=""))
                  xsolvein = -1
                }
                if(xsolvein > 0){
                    posInsert = posInsert + 1
                    ooin=(1/D)*log(xsolvein)
                    corteProb= 1- 1/(1+exp((-1)*D*coeffic[numcat - 1+numFactors])*xsolvein)                                
                    pointprobFull[posInsert,1]= corteProb
                    pointprobFull[posInsert,2]= (coeffic[planex]*ooin)/(coeffic[planex]^2+coeffic[planey]^2)
                    pointprobFull[posInsert,3]= (coeffic[planey]*ooin)/(coeffic[planex]^2+coeffic[planey]^2)
                    pointprobFull[posInsert,4]= categ
                    pointprobFull[posInsert,5]= numcat
                }
              }
            }else{
                if(j==(i+1)){
                  #Case i-(i+1)
                   categ = i + 1                   
                   num= -2*exp((-1)*D*coeffic[categ -1 +numFactors])+
                      exp((-1)*D*coeffic[categ-2+numFactors])+
                          exp((-1)*D*coeffic[categ+numFactors])
                   denom = -2*exp((-1)*D*(coeffic[categ-2+numFactors]+coeffic[categ+numFactors]))+
                          exp((-1)*D*(coeffic[categ-2+numFactors]+coeffic[categ-1+numFactors]))+
                          exp((-1)*D*(coeffic[categ+numFactors]+coeffic[categ-1+numFactors]))
                   expDA=(num)/(denom)
                   if(expDA > 0){
                       posInsert = posInsert + 1
                       corteProb=(expDA*(exp((-1)*D*coeffic[categ-1+numFactors])-exp((-1)*D*coeffic[categ+numFactors])))/
                        ((1+expDA*exp((-1)*D*coeffic[categ-1+numFactors]))*(1+expDA*exp((-1)*D*coeffic[categ+numFactors])))                                       
                       pointprobFull[posInsert,1]= corteProb
                       pointprobFull[posInsert,2]= (coeffic[planex]*((1/D)*log(expDA)))/(coeffic[planex]^2+coeffic[planey]^2)
                       pointprobFull[posInsert,3]= (coeffic[planey]*((1/D)*log(expDA)))/(coeffic[planex]^2+coeffic[planey]^2)
                       pointprobFull[posInsert,4]= categ -1
                       pointprobFull[posInsert,5]= categ
                   }else{
                      print(paste("Curves ",i,"-",i+1," do not intersect at any point",sep=""))
                   }
                }else{
                    #case i-j, with i#j and j>(i+1)
                    categ = j
                    jcomp = i
                    a= exp((-1)*D*(coeffic[categ-1+numFactors]+coeffic[jcomp-1+numFactors]+coeffic[jcomp+numFactors]))-
                          exp((-1)*D*(coeffic[categ+numFactors]+coeffic[jcomp-1+numFactors]+coeffic[jcomp+numFactors]))-
                             exp((-1)*D*(coeffic[categ-1+numFactors]+coeffic[jcomp-1+numFactors]+coeffic[categ+numFactors]))+
                                 exp((-1)*D*(coeffic[categ-1+numFactors]+coeffic[jcomp+numFactors]+coeffic[categ+numFactors]))
                    b= 2*(exp((-1)*D*(coeffic[categ-1+numFactors]+coeffic[jcomp+numFactors]))-
                              exp((-1)*D*(coeffic[categ+numFactors]+coeffic[jcomp-1+numFactors])))
                    c= exp((-1)*D*(coeffic[categ-1+numFactors]))-exp((-1)*D*(coeffic[categ+numFactors]))-
                          exp((-1)*D*(coeffic[jcomp-1+numFactors]))+exp((-1)*D*(coeffic[jcomp+numFactors]))
                    xsolvectjim=matrix(-99,1,2)
                    xsolvectjim= Eq2gSolve(a,b,c)
                    if(xsolvectjim[1,1]>0){
                      xsolvejim=xsolvectjim[1,1]
                    }else if(xsolvectjim[1,2]>0){
                      xsolvejim=xsolvectjim[1,2]
                    }else{
                      print(paste("Curves ",i,"-",j," do not intersect at any point",sep=""))
                      xsolvejim = -1
                    }
                    if(xsolvejim > 0){
                        posInsert = posInsert + 1
                        oojim1=(1/D)*log(xsolvejim)
                        corteProb= (1/(1+(exp(-D*coeffic[categ+numFactors])*xsolvejim)))-(1/(1+(exp(-D*coeffic[categ-1+numFactors])*xsolvejim)))                                          
                        pointprobFull[posInsert,1]= corteProb
                        pointprobFull[posInsert,2]= (coeffic[planex]*oojim1)/(coeffic[planex]^2+coeffic[planey]^2)
                        pointprobFull[posInsert,3]= (coeffic[planey]*oojim1)/(coeffic[planex]^2+coeffic[planey]^2)
                        pointprobFull[posInsert,4]= jcomp
                        pointprobFull[posInsert,5]= categ
                    }
                }
            }
         }
    }
  }
  return(pointprobFull)
}


SeekRowCompleteProb <- function(pointprobFull,coef1,coef2){
    row = 0
    for(r in 1:nrow(pointprobFull)){
        if((pointprobFull[r,4] == coef1) && (pointprobFull[r,5] == coef2)){
            row = r
            break
        }
    }
    if(row == 0){
        print(paste("Row with categories:",coef1,"-",coef2," doesn't exist.",sep=""))
    }
    return(row)
}

PointsBiplotCurves <- function(pointprobFull,numcat){

      axisOrthogonal = FALSE
      if((all(as.integer(pointprobFull[,2]))==0)==TRUE){
          axisOrthogonal = TRUE
      }
  
      pointprob = matrix(0,nrow(pointprobFull),5)
      rowpointprob = 1  
      catref = 1
      ordyprobref = 0
      for(s in 1:(numcat*(numcat -1)/2)){
         if((pointprobFull[s,4] == 1) &&
            (pointprobFull[s,1] > ordyprobref)){
            ordyprobref = pointprobFull[s,1]
            fila = s
         }
      }
     pointprob[1,1] = pointprobFull[fila,1]
     pointprob[1,2] = pointprobFull[fila,2]
     pointprob[1,3] = pointprobFull[fila,3]
     pointprob[1,4] = pointprobFull[fila,4]
     pointprob[1,5] = pointprobFull[fila,5]

     catref = pointprob[1,5]
     if(axisOrthogonal){
        xcatref = pointprob[1,3]  
     }else{
        xcatref = pointprob[1,2]     
     }
     
     if(catref == (numcat - 1)){
        srow = SeekRowCompleteProb(pointprobFull,numcat - 1,numcat)
        if(srow > 0){
          rowpointprob = rowpointprob + 1
          pointprob[rowpointprob,1] = pointprobFull[srow,1]
          pointprob[rowpointprob,2] = pointprobFull[srow,2]
          pointprob[rowpointprob,3] = pointprobFull[srow,3]
          pointprob[rowpointprob,4] = pointprobFull[srow,4]
          pointprob[rowpointprob,5] = pointprobFull[srow,5]
          catref = numcat
        }        
     }else{
        while(catref < numcat){
            if(catref == (numcat - 1)){
                srow = SeekRowCompleteProb(pointprobFull,numcat - 1,numcat)
                if(srow > 0){
                    rowpointprob = rowpointprob + 1
                    pointprob[rowpointprob,1] = pointprobFull[srow,1]
                    pointprob[rowpointprob,2] = pointprobFull[srow,2]
                    pointprob[rowpointprob,3] = pointprobFull[srow,3]
                    pointprob[rowpointprob,4] = pointprobFull[srow,4]
                    pointprob[rowpointprob,5] = pointprobFull[srow,5]
                    catref = numcat
                }
            }else{
                xValues = matrix(0,2,1) 
                for(j in (catref + 1):numcat){
                   findRow = SeekRowCompleteProb(pointprobFull,catref,j)
                   if(findRow > 0){
                      if(axisOrthogonal){
                        columnAdded = c(pointprobFull[findRow,3],j) 
                      }else{
                        columnAdded = c(pointprobFull[findRow,2],j) 
                      }
                      xValues = cbind(xValues,columnAdded)                      
                   }                   
                }
                if(ncol(xValues) == 1){
                    stop(paste("Reference category ",catref," does not intersect with any other upper category",sep=""))
                }else{
                    xValues = as.matrix(xValues[,2:ncol(xValues)])                    
                    if(xcatref > xValues[1,1]){
                        col = which.max(xValues[1,])
                        cunion = xValues[2,col]
                    }else{
                        col = which.min(xValues[1,])
                        cunion = xValues[2,col]
                    }
                        
                     srow = SeekRowCompleteProb(pointprobFull,catref,cunion)
                    if(srow > 0){
                        rowpointprob = rowpointprob + 1
                        pointprob[rowpointprob,1] = pointprobFull[srow,1]
                        pointprob[rowpointprob,2] = pointprobFull[srow,2]
                        pointprob[rowpointprob,3] = pointprobFull[srow,3]
                        pointprob[rowpointprob,4] = catref
                        pointprob[rowpointprob,5] = cunion
                         catref = cunion
                        if(axisOrthogonal){
                            xcatref = pointprob[rowpointprob,3]
                        }else{
                            xcatref = pointprob[rowpointprob,2]
                        }
                     }else{
                        stop(paste("An error ocurred in PointsBiplotCurves function.Categories:",catref,"-",union,sep=""))                    
                    }                     
                }                             
            }
        }
     }

    return(pointprob)
}                                


CalculateFittingIndicatorsMirt <- function(x) {
  
  n <- nrow(x$estimRows)
  p <- ncol(x$estimRows)
  
  fitModelsMirt=0
  for(r in 1:ncol(x$dataFactor)){                 
    J = max(x$dataFactor[,r])  
    y = x$dataFactor[,r]
    Y = matrix(0, n, J)
    for (i in 1:n){
      if (y[i] > 0){
        Y[i, y[i]] = 1
      }
    }
    R = matrix(0, n, J)
    for (i in 1:n){                                                                                           
      R[i, ] = cumsum(Y[i, ])
    }
    A = 1.702*(-1)*x$sepCoefMirt$indCoeffic[r,1:(J-1)]                                      	    
    eta = x$estimRows %*% (1.702*(as.matrix(x$sepCoefMirt$xCoeffic)[r,]))    
    ETA = matrix(1, n, 1) %*% A - eta %*% matrix(1, 1, (J - 1))
    PIA = exp(ETA)/(1 + exp(ETA))
    PIA = cbind(PIA, matrix(1, n, 1))
    PI = matrix(0, n, J)
    PI[, 1] = PIA[, 1]
    PI[, 2:J] = PIA[, 2:J] - PIA[, 1:(J - 1)]
    Rho = log(PIA[, 1:(J - 1)]/(PIA[, 2:J] - PIA[, 1:(J - 1)]))
    gRho = log(1 + exp(Rho))
    L = sum(R[, 1:(J - 1)] * Rho - R[, 2:J] * gRho)
    
    Deviance = -2 * sum(L)
    
    model <- list()
    model$nameVar = dimnames(x$dataFactor)[[2]][r]
    model$nobs=n
    model$J=J
    model$nvar=p
    model$fitted.values = PI
    model$pred = matrix(max.col(PI), n, 1)
    model$clasif = table(y, model$pred)
    model$PercentClasif = sum(y == model$pred)/n
    model$coefficients = as.matrix(x$sepCoefMirt$xCoeffic)[r,]
    model$thresholds = x$sepCoefMirt$indCoeffic[r,]
    model$logLik = L
    model$Deviance = Deviance
    
    # Null Model ---------------------------------
    Beta = matrix(0, p, 1)
    eta = x$estimRows %*% Beta    
    ETA = matrix(1, n, 1) %*% A - eta %*% matrix(1, 1, (J - 1))
    PIA = exp(ETA)/(1 + exp(ETA))
    PIA = cbind(PIA, matrix(1, n, 1))
    PI = matrix(0, n, J)
    PI[, 1] = PIA[, 1]
    PI[, 2:J] = PIA[, 2:J] - PIA[, 1:(J - 1)]
    Rho = log(PIA[, 1:(J - 1)]/(PIA[, 2:J] - PIA[, 1:(J - 1)]))
    gRho = log(1 + exp(Rho))
    
    model$DevianceNull = -2 * sum(R[, 1:(J - 1)] * Rho - R[, 2:J] * gRho)
    model$Dif=(model$DevianceNull - Deviance)
    model$df=p
    model$pval=1-pchisq(model$Dif, df =  model$df)
    model$CoxSnell=1-exp(-1*model$Dif/n)
    model$Nagelkerke= model$CoxSnell/(1-exp((model$DevianceNull/(-2)))^(2/n))
    model$MacFaden=1-(model$Deviance/model$DevianceNull)
    #model$iter=iter
    
    fitModelsMirt=cbind(fitModelsMirt,model)
  }
  fitModelsMirt=fitModelsMirt[,2:(ncol(x$dataFactor)+1)]
  
  return(fitModelsMirt)
}


logit <- function(p) {
	logit = log(p/(1 - p))
	return(logit)
}

ColMax <- function(X) {
	dimens = dim(X)
	n = dimens[1]
	p = dimens[2]
	Maxs = matrix(0, p, 1)
	for (j in (1:p)) Maxs[j] = max(X[, j])
	return(Maxs)
}


EvalOrdlogist <- function(X, par, Ncats) {
	MaxCat = max(Ncats)
	dims = dim(par$coefficients)[2]
	nitems = dim(par$coefficients)[1]
	nnodos = dim(X)[1]
	Numcats = sum(Ncats)
	CumCats = cumsum(Ncats)
	eta = X %*% par$coefficients[1, ]
  ETA = matrix(1, nnodos, 1) %*% t(par$thresholds[1,1:(Ncats[1]-1)]) - eta %*% matrix(1, 1, (Ncats[1] - 1))	
	PIA = exp(ETA)/(1 + exp(ETA))
	PIA = cbind(PIA, matrix(1, nnodos, 1))
	PI = matrix(0, nnodos, Ncats[1])
	PI[, 1] = PIA[, 1]
	PI[, 2:Ncats[1]] = PIA[, 2:Ncats[1]] - PIA[, 1:(Ncats[1] - 1)]
	PT = PI
	for (j in 2:nitems) {
		eta = X %*% par$coefficients[j, ]
    ETA = matrix(1, nnodos, 1) %*% t(par$thresholds[j,1:(Ncats[j]-1)]) - eta %*% matrix(1, 1, (Ncats[j] - 1))			
		PIA = exp(ETA)/(1 + exp(ETA))
		PIA = cbind(PIA, matrix(1, nnodos, 1))
		PI = matrix(0, nnodos, Ncats[j])
		PI[, 1] = PIA[, 1]
		PI[, 2:Ncats[j]] = PIA[, 2:Ncats[j]] - PIA[, 1:(Ncats[j] - 1)]
		PT = cbind(PT, PI)
	}
	return(PT)
}

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
