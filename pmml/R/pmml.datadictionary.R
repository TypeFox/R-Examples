# PMML: Predictive Model Markup Language
#
# Copyright (c) 2009-2015, some parts by Togaware Pty Ltd and other by Zementis, Inc. 
#
# This file is part of the PMML package for R.
#
# The PMML package is free software: you can redistribute it and/or 
# modify it under the terms of the GNU General Public License as
# published by the Free Software Foundation, either version 2 of 
# the License, or (at your option) any later version.
#
# The PMML package is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. Please see the
# GNU General Public License for details (http://www.gnu.org/licenses/).
######################################################################################


.pmmlDataDictionary <- function(field, dataset=NULL, weights=NULL, transformed=NULL)
{
  # 090806 Generate and return a DataDictionary element that includes
  # each supplied field.
  #
  # field$name is a vector of strings, and includes target
  # field$class is indexed by fields$name
  # field$levels is indexed by fields$name
  #
  # 091003 If the dataset is supplied then also include an Interval
  # element within the DataField for each numeric variable.

  number.of.fields <- length(field$name)

  if(field$name[1] == "ZementisClusterIDPlaceHolder" || field$name[1] == "ZementisHiddenTargetField")
  {
   begin<-2
  } else 
  {
   begin <- 1
  }

  namelist <- list()
  dnamelist <- list() 
  optypelist <- list()
  datypelist <- NULL
  fname <- NULL
  data.fields <- list()
  # discrete place holder variable name as set in pmml.naiveBayes.R
  DPL1 <- "DiscretePlaceHolder"
  DPL2 <- "Temp"
  DPL3 <- "predictedScore"

  if(!is.null(transformed))
  {
   for(i in 1:nrow(transformed$fieldData))
   {
     # Determine the operation type

      type <- as.character(transformed$fieldData[i,"dataType"])
      if(type == "numeric")
      {
        datypelist[[row.names(transformed$fieldData)[i]]] <- "double"
      } else if(type == "logical")
      {
        datypelist[[row.names(transformed$fieldData)[i]]] <- "boolean"
      } else
      {
        datypelist[[row.names(transformed$fieldData)[i]]] <- "string"
      }

      if(type == "numeric")
      {
        optypelist[[row.names(transformed$fieldData)[i]]] <- "continuous"
      } else
      {
	optypelist[[row.names(transformed$fieldData)[i]]] <- "categorical"
      }

   }
   if(field$name[1] == "survival")
   {
    datypelist[[field$name[1]]] <- "double"
    optypelist[[field$name[1]]] <- "continuous"
   }
   if(DPL1 %in% field$name)
   {
     datypelist[[DPL1]] <- "string"
     optypelist[[DPL1]] <- "categorical"
   }
   if(DPL2 %in% field$name)
   {
     datypelist[[DPL2]] <- "string"
     optypelist[[DPL2]] <- "categorical"
   }
   if(DPL3 %in% field$name)
   {
     datypelist[[DPL3]] <- "double"
     optypelist[[DPL3]] <- "continuous"
   } 
  } else
  {
   for(i in begin:number.of.fields)
   {
   fname <- field$name[i]
   if(length(grep("as\\.factor\\(",field$name[i])) == 1)
   {
        fname <- gsub("as.factor\\((\\w*)\\)","\\1", field$name[i], perl=TRUE)
   }
     # Determine the operation type

    optype <- "UNKNOWN"
    datype <- "UNKNOWN"
    values <- NULL

    if (field$class[[field$name[i]]] == "numeric")
    {
      optypelist[[fname]] <- "continuous"
      datypelist[[fname]] <- "double"
    }
    else if (field$class[[field$name[i]]] == "logical")
    {
      optypelist[[fname]] <- "categorical"
      datypelist[[fname]] <- "boolean"
    }
    else if (field$class[[field$name[i]]] == "factor")
    {
      optypelist[[fname]] <- "categorical"
      datypelist[[fname]] <- "string"
    }
   }
  }

  for (i in begin:number.of.fields)
  {
    # DataDictionary -> DataField
    if(!is.null(transformed) && i!=1)
    {
       if(transformed$fieldData[field$name[i],"type"] == "original")
       {
         if(!(.removeAsFactor(field$name[i]) %in% namelist))
         {
          namelist <- c(namelist,.removeAsFactor(field$name[i]))
         }
       
       } else
       {
         ofnames <- strsplit(transformed$fieldData[field$name[i],"origFieldName"][[1]],",")[[1]]
         for(j in 1:length(ofnames))
         {
          ofname <- gsub("^\\s+|\\s+$","",ofnames[j])
          hname <- transformed$fieldData[ofname,"origFieldName"]
          ancestorField <- ofname
          while(!is.na(hname))
          {
           ancestorField <- hname
           hname <- transformed$fieldData[hname,"origFieldName"]
          }
          fname <- .removeAsFactor(ancestorField)
          if((!(fname %in% namelist)) && (!(fname %in% dnamelist)))
          {
           namelist <- c(namelist,fname)
           if(!(.removeAsFactor(field$name[i]) %in% dnamelist)) 
             dnamelist <- c(dnamelist, .removeAsFactor(field$name[i]))
          }
        }
       } 
    } else
    {
      fName <- field$name[i]
      if(!is.na(field$class[fName]) && field$class[fName] == "factor")
        optypelist[[fName]] <- "categorical"
      
      if(length(grep("as\\.factor\\(",field$name[i])) == 1)
        fName <- gsub("as.factor\\((\\w*)\\)","\\1", field$name[i], perl=TRUE)

      if(!is.na(field$class[fName]) && field$class[fName] == "factor")
        optypelist[[fName]] <- "categorical"

      if(!(fName %in% namelist) && fName != "ZementisClusterIDPlaceHolder")
        namelist <- c(namelist,fName)
    }
  }

  # DataDictionary
  data.dictionary <- xmlNode("DataDictionary",
                             attrs=c(numberOfFields=length(namelist)))

  if (! is.null(weights) && length(weights))
    data.dictionary <-append.XMLNode(data.dictionary, xmlNode("Extension",
                                                              attrs=c(name="Weights",
                                                              value=weights, extender="Rattle")))

  nmbr <- 1
  for(ndf2 in 1:length(namelist))
  {
    optype <- optypelist[[namelist[ndf2][[1]]]]
    datype <- datypelist[[namelist[ndf2][[1]]]]
    data.fields[[nmbr]] <- xmlNode("DataField", attrs=c(name=namelist[ndf2],
                                             optype=optype, dataType=datype))

   # DataDictionary -> DataField -> Interval
   fname <- namelist[ndf2][[1]]
   if (optypelist[[fname]] == "continuous" && !is.null(dataset) && fname != "survival")
   {
    dataval <- NULL
    for(j in 1:length(dataset[[fname]]))
    {
     dataval<-c(dataval,as.numeric(dataset[[fname]][j]))
    }

    interval <-  xmlNode("Interval",
                           attrs=c(closure="closedClosed", 
			     leftMargin=min(dataval, na.rm=TRUE), # 091025 Handle missing values
                             rightMargin=max(dataval, na.rm=TRUE))) # 091025 Handle missing values
    data.fields[[nmbr]] <- append.XMLNode(data.fields[[nmbr]], interval)
   }

   # DataDictionary -> DataField -> Value
   name <- namelist[nmbr][[1]]
   if (optypelist[[name]] == "categorical")
   {
     if(is.null(field$levels[[name]]) && !is.null(transformed))
     {
       lev <- levels(as.list(unique(transformed$data[name]))[[1]])
       for (j in seq_along(lev))
       {
         data.fields[[nmbr]][[j]] <- xmlNode("Value",
                          attrs=c(value=.markupSpecials(lev[j])))
       }
     } else
     {
       for (j in seq_along(field$levels[[namelist[nmbr][[1]]]]))
       {
         data.fields[[nmbr]][[j]] <- xmlNode("Value",
                          attrs=c(value=field$levels[[namelist[nmbr][[1]]]][j])) 
                          # attrs=c(value=.markupSpecials(field$levels[[namelist[nmbr][[1]]]][j])))
       }
     }
   }

   data.dictionary <- append.XMLNode(data.dictionary, data.fields[[nmbr]])
   nmbr <- nmbr + 1
  }

  return(data.dictionary)
}

.pmmlDataDictionarySurv <- function(field, timeName, statusName, dataset=NULL, weights=NULL, transformed=NULL)
{
  # Tridi 012712
  # modify for a survival model. Survival forests do not typically have
  # a predicted field. Add a generic predicted field.If a predicted
  # field is included, this field will just be ignored by the model. 
  # 090806 Generate and return a DataDictionary element that incldues
  # each supplied field.
  #
  # field$name is a vector of strings, and includes target
  # field$class is indexed by fields$name
  # field$levels is indexed by fields$name
  #
  # 091003 If the dataset is supplied then also include an Interval
  # element within the DataField for each numeric variable.

  number.of.fields <- length(field$name)
  ii<-0

  optypelist <- list()
  namelist <- list()
  datypelist <- list()
  data.fields <- list()

  if(field$name[1] == "ZementisClusterIDPlaceHolder")
  {
   begin <- 2
  } else 
  {
   begin <- 1
  }

  if(!is.null(transformed))
  {
   for(i in 1:nrow(transformed$fieldData))
   {
     # Determine the operation type

      type <- as.character(transformed$fieldData[i,"dataType"])
      if(type == "numeric")
      {
        datypelist[[row.names(transformed$fieldData)[i]]] <- "double"
      } else if(type == "logical")
      {
        datypelist[[row.names(transformed$fieldData)[i]]] <- "boolean"
      } else 
      {
        datypelist[[row.names(transformed$fieldData)[i]]] <- "categorical"
      }

      if(type == "numeric")
      {
        optypelist[[row.names(transformed$fieldData)[i]]] <- "continuous"
      } else
      {
        optypelist[[row.names(transformed$fieldData)[i]]] <- "categorical"
      }
   }
   if(field$name[1] == "survival")
   {
    datypelist[[field$name[1]]] <- "double"
    optypelist[[field$name[1]]] <- "continuous"
   }
  } else
  {
   for(i in begin:number.of.fields)
   {
     # Determine the operation type

    optype <- "UNKNOWN"
    datype <- "UNKNOWN"
    values <- NULL

    if (field$class[[field$name[i]]] == "numeric")
    {
      optypelist[[field$name[i]]] <- "continuous"
      datypelist[[field$name[i]]] <- "double"
    } else if (field$class[[field$name[i]]] == "logical")
    {
      optypelist[[field$name[i]]] <- "categorical"
      datypelist[[field$name[i]]] <- "boolean"
    }
    else if (field$class[[field$name[i]]] == "factor")
    {
      optypelist[[field$name[i]]] <- "categorical"
      datypelist[[field$name[i]]] <- "string"
    }
   }
  }

  for (i in 1:number.of.fields)
  {
    if(length(grep(":",field$name[i])) == 1){
    } else 
    {
    ii<-ii+1 

    # DataDictionary -> DataField
     if(!is.null(transformed))
     {
       if(is.na(transformed$fieldData[field$name[i],"origFieldName"]))
       {
        if(is.na(transformed$fieldData[field$name[i],"transform"]))
        {
         if(!(field$name[i] %in% namelist))
         {
          namelist <- c(namelist,field$name[i])
         }
        }
       } else
       {

           ofname <- transformed$fieldData[field$name[i],"origFieldName"][[1]]
# following for loop not needed as multiple parents via mapvalues dealt with above
           for(j in 1:length(ofname))
           {
            fname <- ofname[j]
            while(!is.na(ofname[j]))
            {
             fname <- ofname[j]
             xvalue <- transformed$fieldData[fname,"transform"]
             if(!is.na(xvalue) && xvalue == "MapValues")
             {
              parents <- transformed$fieldData[fname,"origFieldName"][[1]]
              for(j in 1:length(parents))
              {
               if(!(parents[j] %in% namelist))
               {
                namelist <- c(namelist,parents[j])
               }
              }
              fname <- NA
              break
             }
             ofname[j] <- transformed$fieldData[ofname[j],"origFieldName"][[1]]
            }
            if(!(fname %in% namelist))
            {
             namelist <- c(namelist,fname)
            }
           }
       }
     } else
# if no transforms
     {
      if(length(grep("as\\.factor\\(",field$name[ii])) == 1)
        fName <- gsub("as.factor\\((\\w*)\\)","\\1", field$name[ii], perl=TRUE)
      else
        fName <- field$name[ii]

      data.fields[[ii]] <- xmlNode("DataField", attrs=c(name=fName,
                                                optype=optypelist[[fName]],
                                                dataType=datypelist[[fName]]))

      if(!(fName %in% namelist))
      {
       namelist <- c(namelist,fName)
      }
    
     }
    }
   }

    # DataDictionary -> DataField -> Interval
    nmbr <- 1
    for(ndf2 in 1:length(namelist))
    {
     fname <- namelist[[ndf2]]
     if (optypelist[[fname]] == "continuous" && ! is.null(dataset))
     {
      interval <-  xmlNode("Interval",
                           attrs=c(closure="closedClosed",
                             leftMargin=min(dataset[[namelist[[ndf2]]]],
                               na.rm=TRUE), # 091025 Handle missing values
                             rightMargin=max(dataset[[namelist[[nmbr]]]],
                               na.rm=TRUE))) # 091025 Handle missing values
      data.fields[[ii]] <- append.XMLNode(data.fields[[ii]], interval)
     }
    
    # DataDictionary -> DataField -> Value

    if (optypelist[[fname]] == "categorical")
    {
      if(is.null(field$levels[[fname]]) && !is.null(transformed))
      {
        lev <- levels(as.list(unique(transformed$data[fname]))[[1]])
        for (j in seq_along(lev))
        {
          data.fields[[ii]][[j]] <- xmlNode("Value",
                          attrs=c(value=.markupSpecials(lev[j])))
        }
      } else
      {
        for (j in seq_along(field$levels[[fname]]))
        {
          data.fields[[ii]][[j]] <- xmlNode("Value",
                          attrs=c(value=.markupSpecials(field$levels[[fname]][j])))
        }
      }
    }

#    if (optypelist[[fname]] == "categorical")
#      for (j in seq_along(field$levels[[fname]]))
#        data.fields[[ii]][j] <- xmlNode("Value",
#                                         attrs=c(value=
#                                           .markupSpecials(field$levels[[fname]][j])))
    }

  if (! is.null(weights) && length(weights))
    data.dictionary <-append.XMLNode(data.dictionary, xmlNode("Extension",
                                                              attrs=c(name="Weights",
                                                                value=weights,
                                                                extender="Rattle")))

  data.fields[[ii+1]] <- xmlNode("DataField", attrs=c(name=statusName,
                                                optype="continuous",dataType="double"))
  data.fields[[ii+2]] <- xmlNode("DataField", attrs=c(name=timeName,
                                                optype="continuous",dataType="double"))
  data.fields[[ii+3]] <- xmlNode("DataField", attrs=c(name="cumulativeHazard",
                                                optype="continuous",dataType="double"))

  data.dictionary <- xmlNode("DataDictionary",
                             attrs=c(numberOfFields=length(namelist)+3))
  data.dictionary <- append.XMLNode(data.dictionary, data.fields)


  return(data.dictionary)
}

