# PMML (Predictive Model Markup Language) Transformations 
#
# Copyright (c) 2015 Zementis, Inc.
#
# This file is part of the pmmlTransformations package 
#
# The pmmlTransformations package is free: you can redistribute it and/or 
# modify it under the terms of the GNU General Public License as published 
# by the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# The pmmlTransformations package is distributed in the hope that it will 
# be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. Please see the
# GNU General Public License for details (http://www.gnu.org/licenses/).
############################################################################
#
# Author: Tridivesh Jena
#
#---------------------------------------------------------------------------

MapXform <- function(boxdata,xformInfo,table=NA,defaultValue=NA,mapMissingTo=NA,...)
{
	newrow <- NULL
	colnamesGiven <- FALSE
        j <- 0
	origFieldName <- NULL
	derivedFieldName <- NA
        sampleMin <- NA
        sampleMax<- NA
        xformedMin <- NA
        xformedMax <- NA
        centers <- NA
        scales <- NA
	missingValue <- NA
  
  functionXform <- NA
  
	dataMatrix <- NULL
	default <- NA

	newBoxData <- Initialize(boxdata)
  
	if(!is.list(xformInfo))
 	{
# EXTRACT DATA FROM xformInfo 
#[a,b->c][d,d->s]
        minf <- as.character(xformInfo)

        #separate variable names and data types 
        split0 <- strsplit(minf,"]\\[")[[1]]
        split0[1] <- gsub("\\[","",split0[1])

# is dataTypes given?
        given <- length(split0)

# mapVal: a,b->c
# datType: d,d->s
        if(given == 2)
        {
         mapVal <- split0[1]
         split0[2] <- gsub("]","",split0[2])
         datType <- split0[2]
        } else
        {
         split0[1] <- gsub("]","",split0[1]) 
         mapVal <- split0[1]
         datType <- NA
        }

# mapVal: a,b-> c
        mapVal <- gsub("^[ ]*","",mapVal)
        mapVal <- gsub("[ $]*","",mapVal)

# split the variable and dataType strings
        if(grepl("[^-]->",mapVal))
        {
                st <- "->"
        } else
        {
                st <- "-->"
        }

        val <- strsplit(mapVal,st)[[1]]

# inval: a,b
        inVal <- val[1]
        valsplit <- strsplit(inVal,",")[[1]]

# inVals: "a" "b"
        inVals <- valsplit
        for(i in 1:length(inVals))
        {
         inVals[i] <- gsub("^[ ]*","",inVals[i])
         inVals[i] <- gsub("[ $]*","",inVals[i])
        }

# outVal: "c"
        outVal <- val[2]
        outVal <- gsub("^[ ]*","",outVal)
        outVal <- gsub("[ $]*","",outVal)

# if data types provided
        inDats <- NULL
        outDat <- NULL
        if(!is.na(datType))
        {
         if(grepl("[^-]->",datType))
         {
                 st <- "->"
         } else
         {
                 st <- "-->"
         }
         datsplit <- strsplit(datType,st)[[1]]
         inDat <- datsplit[1]

# inDats: "d" "d"
         inDats <- strsplit(inDat,",")[[1]]
         for(i in 1:length(inDats))
         {
          inDats[i] <- gsub("^[ ]*","",inDats[i])
          inDats[i] <- gsub("[ $]*","",inDats[i])
         }

# outDat: "s"
         outDat <- datsplit[2]
         outDat <- gsub("^[ ]*","",outDat)
         outDat <- gsub("[ $]*","",outDat)
        }

# convert double, integer to numeric for fieldData
        if(!is.null(outDat))
        {
         if((outDat == "double") || (outDat == "integer"))
         {
          outDat <- "numeric"
         }
        } else
        {
          outDat <- "string"
        }

# make variable name and data type rows to bind on data map later
        vnames <- NULL
        vdtype <- NULL
        for(i in 1:length(inVals))
        {
         vnames <- c(vnames,inVals[i])
         if(is.null(inDats[i]))
         {
          vdtype <- c(vdtype,"string")
         } else
         {
          vdtype <- c(vdtype,inDats[i])
         }
        }
        vnames <- c(vnames,outVal)

        if(is.na(outDat))
        {
         vdtype <- c(vdtype,"string")
        } else
        {
          vdtype <- c(vdtype,outDat)
        }

#EXTRACT DATA FROM table
# read data from csv file
	tabl <- as.character(table) 
        file <- scan(tabl,what=character(0),sep=",",quiet=T)
        ndat <- length(file)
        nrows <- length(scan(tabl,what=character(0),sep="\n",quiet=T))
        numcols <- ndat/nrows
        dataMatrix <- matrix(file,nrow=nrows,byrow=TRUE)

#defaultValue=f,mapMissingTo=g
	if(!is.na(defaultValue))
	{
	  default <- as.character(defaultValue)
	}
	if(!is.na(mapMissingTo))
	{
	  missingValue <- as.character(mapMissingTo)
	}

# add variable info to the data matrix
	top <- rbind(vnames,vdtype)
	dataMatrix <- rbind(top,dataMatrix)
	rownames(dataMatrix) <- NULL

	type <- "derived"
##the following origFieldName assignment is replaced by one line below ( origFieldName <- paste(inVals,collapse=",") )
# 	oFN <- NULL
# 	for(i in 1:length(inVals))
# 	{
# 	oFN <- c(oFN,inVals[i])
# 	}
# 	origFieldName <- list(oFN)


#         origfieldName <- paste(inVals,collapse=",") #variable name is misspelled?
  origFieldName <- paste(inVals,collapse=",") #variable name corrected; works for first documentation example but not for second

  if(is.null(outDat))
	{
	 dataType <- "string"
	} else
	{
	 dataType <- outDat
	}
	fieldsMap <- list(dataMatrix)
	transform <- "MapValues"
	derivedFieldName <- outVal


# 	suppressWarnings(newrow <- data.frame(type,dataType,I(origFieldName),sampleMin,sampleMax,xformedMin,xformedMax,centers,scales,I(fieldsMap),transform,default,missingValue,row.names=derivedFieldName,check.names=FALSE))
  suppressWarnings(newrow <- data.frame(type,dataType,I(origFieldName),sampleMin,sampleMax,xformedMin,xformedMax,centers,scales,I(fieldsMap),transform,default,missingValue,functionXform,row.names=derivedFieldName,check.names=FALSE))
  suppressWarnings(newBoxData$fieldData <- rbind(newBoxData$fieldData,newrow))

	newcol <- NULL
	newmatrixcol <- NULL

         if(!is.na(default))
         {
           if(outDat == "numeric")
           {
            newcol <- rep(as.numeric(default),nrow(newBoxData$data))
	    newmatrixcol <- rep(as.numeric(default),nrow(newBoxData$data))
           } else if(outDat == "boolean"){
	    newcol <- rep(as.logical(default),nrow(newBoxData$data))
	   } 
           else
           {
            newcol <- rep(default,nrow(newBoxData$data))
	    newmatrixcol <- rep(default,nrow(newBoxData$data))
           }
         } else
         {
          newcol <- rep(NA,nrow(newBoxData$data))
	  newmatrixcol <- rep(NA,nrow(newBoxData$data))
         }

        if(!is.na(missingValue))
        {
          for(i in 1:length(inVals))
          {
            na <- which(is.na(newBoxData$data[,inVals[i]])) 
            for(j in na)
            {
              if(outDat == "numeric")
              {
                newcol[j] <- rep(as.numeric(missingValue),nrow(newBoxData$data))
                newmatrixcol[j] <- rep(as.numeric(missingValue),nrow(newBoxData$data))
              } else if(outDat == "boolean")
	      {
                newcol[j] <- rep(as.logical(missingValue),nrow(newBoxData$data))
                newmatrixcol[j] <- rep(as.logical(missingValue),nrow(newBoxData$data))
              } else
              {
                newcol[j] <- rep(missingValue,nrow(newBoxData$data))
                newmatrixcol[j] <- rep(missingValue,nrow(newBoxData$data))
              }
            }
          }
        }

	 # for each mapvalue row given except the top 2 (var name and dataType)
         for( j in 3:nrow(dataMatrix))
         {
	  if(outDat == "numeric")
	  {
# for each column of dataMatrix (ie, each input map value): find if input data has that value,
# (do this by creating a matrix with the input values from dataMatrix repeated as many times as the number of 
# data input variables; so each input value has the same matrix row to compare with)
# result is a set of True and False. Apply function 'all' by row to see if all values in a row are true
# resulting rows are the rows which match all the input map values  
newcol[apply(as.matrix(newBoxData$data[,dataMatrix[1,1:(ncol(dataMatrix)-1)]]==matrix(rep(dataMatrix[j,1:(ncol(dataMatrix)-1)],nrow(newBoxData$data)),nrow=nrow(newBoxData$data),byrow=T)),1,all)]<-as.numeric(dataMatrix[j,ncol(dataMatrix)])
newmatrixcol[apply(as.matrix(newBoxData$data[,dataMatrix[1,1:(ncol(dataMatrix)-1)]]==matrix(rep(dataMatrix[j,1:(ncol(dataMatrix)-1)],nrow(newBoxData$data)),nrow=nrow(newBoxData$data),byrow=T)),1,all)]<-as.numeric(dataMatrix[j,ncol(dataMatrix)])
	  } else if(outDat == "boolean"){
newcol[apply(as.matrix(newBoxData$data[,dataMatrix[1,1:(ncol(dataMatrix)-1)]]==matrix(rep(dataMatrix[j,1:(ncol(dataMatrix)-1)],nrow(newBoxData$data)),nrow=nrow(newBoxData$data),byrow=T)),1,all)]<-as.logical(dataMatrix[j,ncol(dataMatrix)])
newmatrixcol[apply(as.matrix(newBoxData$data[,dataMatrix[1,1:(ncol(dataMatrix)-1)]]==matrix(rep(dataMatrix[j,1:(ncol(dataMatrix)-1)],nrow(newBoxData$data)),nrow=nrow(newBoxData$data),byrow=T)),1,all)]<-as.logical(dataMatrix[j,ncol(dataMatrix)])
	  } else
	  {
newcol[apply(as.matrix(newBoxData$data[,dataMatrix[1,1:(ncol(dataMatrix)-1)]]==matrix(rep(dataMatrix[j,1:(ncol(dataMatrix)-1)],nrow(newBoxData$data)),nrow=nrow(newBoxData$data),byrow=T)),1,all)]<-dataMatrix[j,ncol(dataMatrix)]
newmatrixcol[apply(as.matrix(newBoxData$data[,dataMatrix[1,1:(ncol(dataMatrix)-1)]]==matrix(rep(dataMatrix[j,1:(ncol(dataMatrix)-1)],nrow(newBoxData$data)),nrow=nrow(newBoxData$data),byrow=T)),1,all)]<-dataMatrix[j,ncol(dataMatrix)]
	  }
         }

     col <- as.matrix(newcol)
     newBoxData$data <- data.frame(newBoxData$data,col,check.names=FALSE)
     if(outDat=="string")
         newBoxData$data[,ncol(newBoxData$data)] <- as.factor(newBoxData$data[,ncol(newBoxData$data)])
     colnames(newBoxData$data)[ncol(newBoxData$data)] <- dataMatrix[1,ncol(dataMatrix)]
     rownames(newBoxData$data) <- NULL

     if(!is.null(newBoxData$matrixData))
     {
      matrixcol <- as.matrix(newmatrixcol) 
      newBoxData$matrixData <- cbind(newBoxData$matrixData,matrixcol)
      colnames(newBoxData$matrixData) <- colnames(newBoxData$data)
      rownames(newBoxData$matrixData) <- NULL
     }
     } else 
     {
      for(k in 1:length(xformInfo))
      {
	ifelse(is.list(xformInfo),xform<-xformInfo[[k]],xform<-table[[k]])
	datatypes<-c("string","String","double","Double","boolean","Boolean","integer","Integer","float","Float")
#	if(!(TRUE %in% (xform %in% datatypes)))
#	if(!("mapMissingTo" %in% rownames(xform)))
# check if data types given in 2nd row
# Warning: if the map itself involves terms like "double"; have to make sure 2nd row is datatype; 
# dont depend on default
#        if(!((xform[2,1] %in% datatypes) && (xform[2,2] %in% datatypes)))
        if(!(xform[2,ncol(xform)] %in% datatypes))
	{
	  datype<-rep("string",ncol(xform))
          for(l in 1:ncol(xform)){xform[[l]]<-as.character(xform[[l]])}
	  xform <- rbind(xform[1,],datype,xform[-1,])
	}
	orig <- as.matrix(xform[1,-ncol(xform)])
	colnames(orig) <- NULL
	rownames(orig) <- NULL
	origFieldName <- list(orig) #this line is present in the current CRAN package
#         origFieldName <- paste(xform,collapse=",") #pmml() function will break if using this line
 
	derivedFieldName <- as.character(xform[1,ncol(xform)])
	default <- defaultValue[k]
	missingValue <- mapMissingTo[k]
        if((as.character(xform[2,ncol(xform)]) == "double") || (as.character(xform[2,ncol(xform)]) == "integer"))
        {
         outDat <- "numeric"
        } else if(as.character(xform[2,ncol(xform)]) == "double")
	{
	  outDat <- "boolean"
	} else
	{
	  outDat <- "string"
	} 

	if(derivedFieldName %in% rownames(newBoxData$fieldData)[newBoxData$fieldData[,"type"]=="derived"])
	{
	  newBoxData$fieldData <- newBoxData$fieldData[-which(rownames(newBoxData$fieldData)==derivedFieldName),]
	}

# 	suppressWarnings(newrow <- data.frame("derived",outDat,I(origFieldName),sampleMin,sampleMax,xformedMin,xformedMax,centers,scales,I(list(as.matrix(xform))),"MapValues",defaultValue[k],mapMissingTo[k],row.names=derivedFieldName,check.names=FALSE))
#   colnames(newrow)<-c("type","dataType","origFieldName","sampleMin","sampleMax","xformedMin","xformedMax","centers","scales","fieldsMap","transform","default","missingValue")

  suppressWarnings(newrow <- data.frame("derived",outDat,I(origFieldName),sampleMin,sampleMax,xformedMin,xformedMax,centers,scales,I(list(as.matrix(xform))),"MapValues",defaultValue[k],mapMissingTo[k],functionXform,row.names=derivedFieldName,check.names=FALSE))
  colnames(newrow)<-c("type","dataType","origFieldName","sampleMin","sampleMax","xformedMin","xformedMax","centers","scales","fieldsMap","transform","default","missingValue","functionXform")



  suppressWarnings(newBoxData$fieldData <- rbind(newBoxData$fieldData,newrow))
#?need to convert xform to a matrix? 
#?need to remove rownames?
	dataMatrix <- as.matrix(xform)

        newcol <- NULL
        newmatrixcol <- NULL
        if(!is.na(default))
        {
          if(outDat == "numeric")
          {
           newcol <- rep(as.numeric(default),nrow(newBoxData$data))
           newmatrixcol <- rep(as.numeric(default),nrow(newBoxData$data))
          } else if(outDat == "boolean"){
           newcol <- rep(as.logical(default),nrow(newBoxData$data))
          }
          else
          {
           newcol <- rep(default,nrow(newBoxData$data))
           newmatrixcol <- rep(default,nrow(newBoxData$data))
          }
        } 

        if(!is.na(missingValue))
        {
          for(i in 1:length(orig[1,]))
          {
            na <- which(is.na(newBoxData$data[,orig[1,i]])) 
            for(j in na)
            {
              if(outDat == "numeric")
              {
                newcol[j] <- rep(as.numeric(missingValue),nrow(newBoxData$data))
                newmatrixcol[j] <- rep(as.numeric(missingValue),nrow(newBoxData$data))
              } else if(outDat == "boolean"){
                newcol[j] <- rep(as.logical(missingValue),nrow(newBoxData$data))
                newmatrixcol[j] <- rep(as.logical(missingValue),nrow(newBoxData$data))
              } else
              {
                newcol[j] <- rep(missingValue,nrow(newBoxData$data))
                newmatrixcol[j] <- rep(missingValue,nrow(newBoxData$data))
              }
            }
          }
        }

         # for each mapvalue row given except the top 2 (var name and dataType)
         for( j in 3:nrow(dataMatrix))
         {
          if(outDat == "numeric")
          {
# for each column of dataMatrix (ie, each input map value): find if input data has that value,
# (do this by creating a matrix with the input values from dataMatrix repeated as many times as the number of
# data input variables; so each input value has the same matrix row to compare with)
# result is a set of True and False. Apply function 'all' by row to see if all values in a row are true
# resulting rows are the rows which match all the input map values
newcol[apply(as.matrix(newBoxData$data[,dataMatrix[1,1:(ncol(dataMatrix)-1)]]==matrix(rep(dataMatrix[j,1:(ncol(dataMatrix)-1)],nrow(newBoxData$data)),nrow=nrow(newBoxData$data),byrow=T)),1,all)]<-as.numeric(dataMatrix[j,ncol(dataMatrix)])
newmatrixcol[apply(as.matrix(newBoxData$data[,dataMatrix[1,1:(ncol(dataMatrix)-1)]]==matrix(rep(dataMatrix[j,1:(ncol(dataMatrix)-1)],nrow(newBoxData$data)),nrow=nrow(newBoxData$data),byrow=T)),1,all)]<-as.numeric(dataMatrix[j,ncol(dataMatrix)])
          } else if(outDat == "boolean"){
newcol[apply(as.matrix(newBoxData$data[,dataMatrix[1,1:(ncol(dataMatrix)-1)]]==matrix(rep(dataMatrix[j,1:(ncol(dataMatrix)-1)],nrow(newBoxData$data)),nrow=nrow(newBoxData$data),byrow=T)),1,all)]<-as.logical(dataMatrix[j,ncol(dataMatrix)])
newmatrixcol[apply(as.matrix(newBoxData$data[,dataMatrix[1,1:(ncol(dataMatrix)-1)]]==matrix(rep(dataMatrix[j,1:(ncol(dataMatrix)-1)],nrow(newBoxData$data)),nrow=nrow(newBoxData$data),byrow=T)),1,all)]<-as.logical(dataMatrix[j,ncol(dataMatrix)])
          } else
          {
newcol[apply(as.matrix(newBoxData$data[,dataMatrix[1,1:(ncol(dataMatrix)-1)]]==matrix(rep(dataMatrix[j,1:(ncol(dataMatrix)-1)],nrow(newBoxData$data)),nrow=nrow(newBoxData$data),byrow=T)),1,all)]<-dataMatrix[j,ncol(dataMatrix)]
newmatrixcol[apply(as.matrix(newBoxData$data[,dataMatrix[1,1:(ncol(dataMatrix)-1)]]==matrix(rep(dataMatrix[j,1:(ncol(dataMatrix)-1)],nrow(newBoxData$data)),nrow=nrow(newBoxData$data),byrow=T)),1,all)]<-dataMatrix[j,ncol(dataMatrix)]
          }
         }

#         if(!is.na(missingValue))
#         {
#          mis <- which(is.na(newcol))
#          newcol[mis] <- missingValue
#         }

      	col <- as.matrix(newcol)
      	matrixcol <- as.matrix(newmatrixcol)

     	newBoxData$data <- data.frame(newBoxData$data,col,check.names=FALSE)
        colnames(newBoxData$data)[ncol(newBoxData$data)] <- dataMatrix[1,ncol(dataMatrix)]
        rownames(newBoxData$data) <- NULL

	if(!is.null(newBoxData$matrixData))
	{
     	 newBoxData$matrixData <- cbind(newBoxData$matrixData,matrixcol)
         colnames(newBoxData$matrixData) <- colnames(newBoxData$data)
         rownames(newBoxData$matrixData) <- NULL
        }
     }
    }
     return(newBoxData)
}
