# PMML (Predictive Model Markup Language) Transformations 
#
# Copyright (c) 2015 Zementis, Inc.
#
# This file is part of the pmmlTransformations package. 
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

DiscretizeXform <-
function(boxdata,xformInfo,table,defaultValue=NA,mapMissingTo=NA,...)
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
        fieldsMap <- NA
	dataMatrix <- NULL
	default <- NA
	missingValue <- NA
        sep <- ":"
	functionXform <- NA
  
	newBoxData <- Initialize(boxdata)
#	if(!is.numeric(newBoxData$data))
#	  stop("Non-numeric matrices not yet supported for transformations")

#if command given in older format 
 	if(!is.list(xformInfo))
        {
# EXTRACT DATA FROM xformInfo 
#[a->c][d->s]
        minf <- as.character(xformInfo)

        #separate variable names and data types 
        split0 <- strsplit(minf,"]\\[")[[1]]
        split0[1] <- gsub("\\[","",split0[1])

# is dataTypes given?
        given <- length(split0)

# discretize: a->c
# dataTypes: d->s
        if(given == 2)
        {
         discretize <- split0[1]
         split0[2] <- gsub("]","",split0[2])
         datType <- split0[2]
        } else
        {
         split0[1] <- gsub("]","",split0[1]) 
         discretize <- split0[1]
         datType <- NA
        }

# discretize: a->c
# split the variable and dataType strings
        if(grepl("[^-]->",discretize))
        {
                st <- "->"
        } else 
        {
                st <- "-->"
        }

        val <- strsplit(discretize,st)[[1]]

# inVal: a
        inVal <- gsub(" ","",val[1])

# outVal: "c"
        outVal <- gsub(" ","",val[2])

# if data types provided
        inDat <- NULL
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
         inDat <- gsub(" ","",datsplit[1])

# outDat: "s"
         outDat <- gsub(" ","",datsplit[2])
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
	vname <- NULL	
        vdtype <- NULL
        for(i in 1:length(inVal))
        {
         vname <- c(vname,inVal)
         if(is.null(inDat))
         {
          vdtype <- c(vdtype,"numeric")
         } else
         {
          vdtype <- c(vdtype,inDat)
         }
        }
        vname <- c(vname,outVal)

        if(is.na(outDat))
        {
         vdtype <- c(vdtype,"string")
        } else
        {
          vdtype <- c(vdtype,outDat)
        }

# placeholder for interval and values
	vname <- c(vname,"leftInterval","leftValue","rightInterval","rightValue")
        vdtype <- c(vdtype,"string","double","string","double")

#EXTRACT DATA FROM table
# read interval data from csv file
	tabl <- as.character(table) 
        file <- scan(tabl,what=character(0),sep=",",quiet=T)
        ndat <- length(file)
        nrows <- length(scan(tabl,what=character(0),sep="\n",quiet=T))
        numcols <- ndat/nrows
        dataMatrix <- matrix(file,nrow=nrows,byrow=TRUE)

# add columns for left/right intervalType/value
        dataMatrix <- cbind(dataMatrix,NA)
        dataMatrix <- cbind(dataMatrix,NA)
        dataMatrix <- cbind(dataMatrix,NA)
        dataMatrix <- cbind(dataMatrix,NA)

#defaultValue=f,mapMissingTo=g
	if(!is.na(defaultValue))
	{
	  default <- as.character(defaultValue)
#          if(!((defaultValue==1) || (defaultValue==0) || toupper(defaultValue)==TRUE || toupper(defaultValue)==FALSE))
#          {
#            stop("defaultValue must be a proper boolean value")
#          }
	}

	if(!is.na(mapMissingTo))
	{
	  missingValue <- as.character(mapMissingTo)
#          if(!((missingValue==1) || (missingValue==0) || toupper(missingValue)==TRUE || toupper(missingValue)==FALSE))
#          {
#            stop("missingValue must be a proper boolean value")
#          }
	}


# add variable info to the data matrix
	top <- rbind(vname,vdtype)
	dataMatrix <- rbind(top,dataMatrix)
	rownames(dataMatrix) <- NULL

	type <- "derived"
	origFieldName <- inVal

	if(is.null(outDat))
	{
	 dataType <- "string"
	} else
	{
	 dataType <- outDat
	}

	transform <- "discretize"
	derivedFieldName <- outVal

         # for each field map row given except the top 2 (name and type)
         for( j in 3:nrow(dataMatrix))
         {
           leftValue <- NA
           rightValue <- NA
           leftInterval <- NA
           rightInterval <- NA

           if(grepl(sep,dataMatrix[j,1]))
           {
            range <- strsplit(dataMatrix[j,1],sep)[[1]]

            if(grepl("^\\[",range[1]))
            {
             leftValue <- gsub("\\[","",range[1])
             leftInterval <- "closed"
            } else if(grepl("^\\(",range[1]))
            {
             leftValue <- gsub("\\(","",range[1])
             leftInterval <- "open"
            }

            if(grepl("\\]$",range[2]))
            {
             rightValue <- gsub("\\]","",range[2])
             rightInterval <- "Closed"
            } else if(grepl("\\)$",range[2]))
            {
             rightValue <- gsub("\\)","",range[2])
             rightInterval <- "Open"
            }
# end if both left and right limits given
           } else
           {
            range <- dataMatrix[j,1]
            if(grepl("^\\[",range[1]))
            {
             leftValue <- gsub("\\[","",range[1])
             leftInterval <- "closed"
             rightInterval <- "Open"
            } else if(grepl("^\\(",range[1]))
            {
             leftValue <- gsub("\\(","",range[1])
             leftInterval <- "open"
             rightInterval <- "Open"
            } else if(grepl("\\]$",range[1]))
            {
             rightValue <- gsub("\\]","",range[1])
             leftInterval <- "open"
             rightInterval <- "Closed"
            } else if(grepl("\\)$",range[1]))
            {
             rightValue <- gsub("\\)","",range[1])
             leftInterval <- "open"
             rightInterval <- "Open"
            }
           }

	   dataMatrix[j,3] <- leftInterval
	   dataMatrix[j,4] <- leftValue
           dataMatrix[j,5] <- rightInterval
           dataMatrix[j,6] <- rightValue
          }

	colnames(dataMatrix) <- dataMatrix[1,]
        fieldsMap <- list(dataMatrix)
#         suppressWarnings(newrow <- data.frame(type,dataType,I(origFieldName),sampleMin,sampleMax,xformedMin,xformedMax,centers,scales,I(fieldsMap),transform,default,missingValue,row.names=derivedFieldName,check.names=FALSE))
        suppressWarnings(newrow <- data.frame(type,dataType,I(origFieldName),sampleMin,sampleMax,xformedMin,xformedMax,centers,scales,I(fieldsMap),transform,default,missingValue,functionXform,row.names=derivedFieldName,check.names=FALSE))

        suppressWarnings(newBoxData$fieldData <- rbind(newBoxData$fieldData,newrow))

	newBoxData <- .xformData(dataMatrix,default,missingValue,dataType,newBoxData)
     } else
     {
#for each xform indicated
      for(k in 1:length(xformInfo))
      {
	dataMatrix<-as.matrix(xformInfo[[k]])
        colnames(dataMatrix) <- NULL
	datatypes<-c("string","String","double","Double","boolean","Boolean","integer","Integer","float","Float")

	for(l in 1:ncol(dataMatrix)){dataMatrix[[l]]<-as.character(dataMatrix[[l]])}

	if(!(dataMatrix[2,1] %in% datatypes))
	  dataMatrix<-rbind(c("string","numeric","string","double","string","double"),dataMatrix)

        fieldsMap <- list(dataMatrix)
        dataType <- dataMatrix[2,2]
	colnames(dataMatrix) <- dataMatrix[1,]

	origFieldName <- dataMatrix[1,1]
	derivedFieldName <-dataMatrix[1,2]
	default <- defaultValue[k]
	missingValue <- mapMissingTo[k]

# 	suppressWarnings(newrow <- data.frame(type="derived",dataType=dataType,I(origFieldName),sampleMin,sampleMax,xformedMin,xformedMax,centers,scales,I(fieldsMap),transform="discretize",default,missingValue,row.names=derivedFieldName,check.names=FALSE))
	suppressWarnings(newrow <- data.frame(type="derived",dataType=dataType,I(origFieldName),sampleMin,sampleMax,xformedMin,xformedMax,centers,scales,I(fieldsMap),transform="discretize",default,missingValue,functionXform,row.names=derivedFieldName,check.names=FALSE))
	
  suppressWarnings(newBoxData$fieldData <- rbind(newBoxData$fieldData,newrow))

	colnames(dataMatrix) <- dataMatrix[1,]

	newBoxData <- .xformData(dataMatrix,default,missingValue,dataType,newBoxData)
       }
     }
     return(newBoxData)
     }
      .xformData <- function(dataMatrix,default,missingValue,dataType,newBoxData)
      {
        type <- NULL

        dataMatrix <- tolower(dataMatrix)

        for(j in 3:nrow(dataMatrix))
	{
	  if(is.na(dataMatrix[j,"leftValue"]) && (dataMatrix[j,"rightInterval"] == "closed"))
	  {
	    type[j] <- 1
	    next
	  }
          if(is.na(dataMatrix[j,"leftValue"]) && (dataMatrix[j,"rightInterval"] == "open"))
          {
            type[j] <- 2
	    next
          }
          if((dataMatrix[j,"leftInterval"] == "closed") && is.na(dataMatrix[j,"rightValue"]))
          {
            type[j] <- 7
	    next
          }
          if((dataMatrix[j,"leftInterval"] == "open") && is.na(dataMatrix[j,"rightValue"]))
          {
            type[j] <- 8
	    next
          }
	  if((dataMatrix[j,"leftInterval"] == "closed") && (dataMatrix[j,"rightInterval"] == "closed"))
          {
            type[j] <- 3
	    next
          }
          if((dataMatrix[j,"leftInterval"] == "closed") && (dataMatrix[j,"rightInterval"] == "open"))
          {
            type[j] <- 4
	    next
          }
          if((dataMatrix[j,"leftInterval"] == "open") && (dataMatrix[j,"rightInterval"] == "closed"))
          {
            type[j] <- 5
	    next
          }
          if((dataMatrix[j,"leftInterval"] == "open") && (dataMatrix[j,"rightInterval"] == "open"))
          {
            type[j] <- 6
	    next
          }
	}

	newcol <- NULL
#	origName <- as.character(dataMatrix[1,1])
#	derivedName <- as.character(dataMatrix[1,2])
        origName <- colnames(dataMatrix)[1]
        derivedName <- colnames(dataMatrix)[2]

         if(!is.numeric(newBoxData$data[1,1]))
           stop("Non-numeric matrices not yet supported for transformations")

# Tridi: 9/20/13: Decided that values should be initialized to defaultt rather than missing

        if(!is.na(default))
        {
          if(dataType == "numeric")
          {
           newcol <- rep(as.numeric(default),nrow(newBoxData$data))
           newmatrixcol <- rep(as.numeric(default),nrow(newBoxData$data))
          } else if(dataType == "boolean"){
           newcol <- rep(as.logical(default),nrow(newBoxData$data))
          }
          else
          {
           newcol <- rep(default,nrow(newBoxData$data))
           newmatrixcol <- rep(default,nrow(newBoxData$data))
          }
        } else if(!is.na(missingValue))
        {
          if(dataType == "numeric")
          {
           newcol <- rep(as.numeric(missingValue),nrow(newBoxData$data))
           newmatrixcol <- rep(as.numeric(missingValue),nrow(newBoxData$data))
          } else if(missingValue == "boolean"){
           newcol <- rep(as.logical(missingValue),nrow(newBoxData$data))
          }
          else
          {
           newcol <- rep(missingValue,nrow(newBoxData$data))
           newmatrixcol <- rep(missingValue,nrow(newBoxData$data))
          }
        } else
        {
         newcol <- rep(NA,nrow(newBoxData$data))
         newmatrixcol <- rep(NA,nrow(newBoxData$data))
        }


	origFieldName <- colnames(dataMatrix)[1]
        derivedName <- colnames(dataMatrix)[2]

	 # for each field map row given except the top 2 (name and type)
         for( j in 3:nrow(dataMatrix))
         {
#print("ROW BEGIN")
#print(proc.time())
          if(type[j] == 1)
          {
           newcol[newBoxData$data[,origName]<=as.numeric(dataMatrix[j,"rightValue"])] <- dataMatrix[j,derivedName] 
          }
          if(type[j] == 2)
          {
           newcol[newBoxData$data[,origName]<as.numeric(dataMatrix[j,"rightValue"])] <- dataMatrix[j,derivedName]                                                                 
          }
	  if(type[j] == 3)
	  {
	   newcol[(as.numeric(dataMatrix[j,"leftValue"])<=newBoxData$data[,origName]) & (newBoxData$data[,origFieldName]<=as.numeric(dataMatrix[j,"rightValue"]))] <- dataMatrix[j,derivedName]
	  }
          if(type[j] == 4)
          {
           newcol[(as.numeric(dataMatrix[j,"leftValue"])<=newBoxData$data[,origName]) & (newBoxData$data[,origFieldName]<as.numeric(dataMatrix[j,"rightValue"]))] <- dataMatrix[j,derivedName]
          }
          if(type[j] == 5)
          {
           newcol[(as.numeric(dataMatrix[j,"leftValue"])<newBoxData$data[,origName]) & (newBoxData$data[,origFieldName]<=as.numeric(dataMatrix[j,"rightValue"]))] <- dataMatrix[j,derivedName]
          }
          if(type[j] == 6)
          {
           newcol[(as.numeric(dataMatrix[j,"leftValue"])<newBoxData$data[,origName]) & (newBoxData$data[,origFieldName]<as.numeric(dataMatrix[j,"rightValue"]))] <- dataMatrix[j,derivedName]
          }
          if(type[j] == 7)
          {
           newcol[as.numeric(dataMatrix[j,"leftValue"])<=newBoxData$data[,origName]] <- dataMatrix[j,derivedName]
          }
          if(type[j] == 8)
          {
           newcol[as.numeric(dataMatrix[j,"leftValue"])<newBoxData$data[,origName]] <- dataMatrix[j,derivedName]
          }
#print("ROW END")
#print(proc.time())
	}

      col <- as.matrix(newcol)
      colnames(col) <- colnames(dataMatrix)[2]
      rownames(col) <- NULL

     if(dataType == "numeric")
     {
       newBoxData$data <- data.frame(newBoxData$data,as.numeric(col),check.names=FALSE)
       colnames(newBoxData$data)[ncol(newBoxData$data)] <- derivedName 
     } else
     {
       newBoxData$data <- data.frame(newBoxData$data,col,check.names=FALSE)
       newBoxData$data[,dim(newBoxData$data)[2]] <- as.factor(newBoxData$data[,dim(newBoxData$data)[2]]) 
     }
#new
     if(!is.null(newBoxData$matrixData))
     {
      if(dataType == "numeric")
      {
       newBoxData$matrixData <- cbind(newBoxData$matrixData,as.numeric(col))
       colnames(newBoxData$matrixData) <- colnames(newBoxData$data)
      } else
      {
       newBoxData$matrixData <- cbind(newBoxData$matrixData,col)
      }
     }

     return(newBoxData)
    }
