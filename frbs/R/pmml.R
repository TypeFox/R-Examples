# PMML: Predictive Model Markup Language
#
# Copyright (c) 2009-2013, some parts by Togaware Pty Ltd and other by Zementis, Inc. 
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

#######################################################################
# MAIN PMML FUNCTION
#
# modified 081513 to add implementation of transformations generator
# and transformations element addition;  tridivesh.jena@zementis.com
# 
pmml <- function(model=NULL,
                 model.name="frbs_Model",
                 app.name="frbsPMML",
                 description=NULL,
                 copyright=NULL,
                 transforms=NULL,
                 ...)
{
  if(is.null(model) && !is.null(transforms))
  {
    field <- NULL
    field$name <- names(transforms$fieldData)
    field$class <- transforms$fieldData[,"dataType"]
    names(field$class) <- row.names(transforms$fieldData)

    return(.pmmlLocalTransformations(field, transforms, NULL))
  }
  else
  {
    UseMethod("pmml")
  }
}

########################################################################
# UTILITY FUNCTIONS

.markupSpecials <- function(x)
  gsub("<", "&lt;", gsub(">", "&gt;", gsub("&", "&amp;", x)))

.generateCopyright <- function()
{
  return(paste("Copyright (c)", format(Sys.time(), "%Y"), Sys.info()["user"]))
}

.sdecimal2binary <- function(x)
{
  return(rev(.sdecimal2binary.smallEndian(x)))
}

.sdecimal2binary.smallEndian <- function(x)
{
  if (x==0) return(0)
  if (x<0) stop("Sorry, the input must be positive")
  dec <- x

  n <- floor(log(x)/log(2))
  bin <- c(1)
  dec <- dec - 2 ^ n

  while(n > 0)
  {
    if (dec >= 2 ^ (n-1)) {bin <- c(bin,1); dec <- dec - 2 ^ (n-1)}
    else bin <- c(bin,0)
    n <- n - 1
  }
  return(bin)
}

# Function .pmmlRootNode

.pmmlRootNode <- function(version)
{
  requireNamespace("XML", quietly = TRUE)
 
 if (version == "1.0")
    node <- XML::xmlNode("frbsPMML",
                    attrs=c(version="1.0",
                      xmlns="http://sci2s.ugr.es/dicits/",
                      "xmlns:xsi"="http://www.w3.org/2001/XMLSchema-instance", 
                      "xsi:schemaLocation"=paste("http://sci2s.ugr.es/dicits/",
                        "http://sci2s.ugr.es/dicits/")))
  else
    node <- XML::xmlNode("frbsPMML",
                    attrs=c(version="4.0",
                      xmlns="http://www.dmg.org/PMML-4_0",
                      "xmlns:xsi"="http://www.w3.org/2001/XMLSchema-instance",
                      "xsi:schemaLocation"=paste("http://www.dmg.org/PMML-4_0",
                        "http://dicits.ugr.es/software/frbsJpmml/")))

  return(node)
}

.pmmlHeader <- function(description, copyright, app.name)
{
	requireNamespace("XML", quietly = TRUE)
 
	if (is.null(copyright)) copyright <- .generateCopyright()
    
    # Header Node
    header <- XML::xmlNode("Header", attrs=c(copyright=copyright, description=description))
    
    # Header -> Extension for user info
    header <- XML::append.XMLNode(header, XML::xmlNode("Extension", attrs=c(name="user",value=sprintf("%s", Sys.info()["user"]), extender=app.name)))
    
    # Header -> Application
    VERSION <- "1.4"
    header <- XML::append.XMLNode(header, XML::xmlNode("Application", attrs=c(name=app.name, version=VERSION)))
    
    # Header -> Timestamp
    header <- XML::append.XMLNode(header, XML::xmlNode("Timestamp", sprintf("%s", Sys.time())))
    
    return(header)
}


.pmmlLocalTransformations <- function(field, transforms=NULL, LTelement=NULL)
{

  # 090806 Generate and return a LocalTransformations element that incldues
  # each supplied field.
  #
  # field$name is a vector of strings, and includes target
  # field$class is indexed by fields$name

  # LocalTransformations
  requireNamespace("XML", quietly = TRUE)
  
  if(is.null(LTelement))
  {
   local.transformations <- XML::xmlNode("LocalTransformations")
  }
  target <- field$name[1]

  if(!is.null(transforms))
  {
    inputs <- transforms$fieldData

# list of all fields derived from the target field
   targetDL <- NULL
   targetDL <- c(targetDL,target)

# not used code to make sure list of unique elements 
#       flist <- flist[!duplicate(flist)]

# code to output all fields, possibly to allow user to output any derived fields via OutputField element
   for(i in 1:nrow(inputs))
   {

    if(inputs[i,"origFieldName"] %in% targetDL)
    {
     targetDL <- c(targetDL,rownames(inputs)[i])
     if(rownames(inputs)[i] %in% field$name[-1])
     {
       stop("Target variable and derivations are not allowed to be used as input variables.")
     }
    } else 
    {
     fname <- rownames(inputs)[i]

     if(inputs[fname,"type"] == "derived" && fname != target)
     {
      if(inputs[fname,"transform"] == "zxform")
      {
       origName <- inputs[fname,"origFieldName"]
       missing <- inputs[fname,"missingValue"]
       dfNode <- XML::xmlNode("DerivedField",attrs=c(name=fname,dataType="double",optype="continuous"))
       if(!is.na(missing))
       {
         ncNode <- XML::xmlNode("NormContinuous",attrs=c(mapMissingTo=missing,field=origName))
       } else
       {
         ncNode <- XML::xmlNode("NormContinuous",attrs=c(field=origName))     
       }

       o1 <- as.numeric(inputs[fname,"centers"])
       o2 <- as.numeric(inputs[fname,"centers"]) + as.numeric(inputs[fname,"scales"])
       lnNode1 <- XML::xmlNode("LinearNorm",attrs=c(orig=o1,norm="0"))
       lnNode2 <- XML::xmlNode("LinearNorm",attrs=c(orig=o2,norm="1"))

       ncNode <- XML::append.XMLNode(ncNode, lnNode1)
       ncNode <- XML::append.XMLNode(ncNode, lnNode2)
       dfNode <- XML::append.XMLNode(dfNode, ncNode)
#       local.transformations <- append.XMLNode(local.transformations, dfNode)
      } else if(inputs[fname,"transform"] == "minmax")
      {
       origName <- inputs[fname,"origFieldName"]
       missing <- inputs[fname,"missingValue"]
       dfNode <- XML::xmlNode("DerivedField",attrs=c(name=fname,dataType="double",optype="continuous"))
       if(!is.na(missing))
       {
         ncNode <- XML::xmlNode("NormContinuous",attrs=c(mapMissingTo=missing,field=origName)) 
       } else
       {
         ncNode <- XML::xmlNode("NormContinuous",attrs=c(field=origName)) 
       }

       o1 <- inputs[fname,"sampleMin"]
       n1 <- inputs[fname,"xformedMin"]
       o2 <- inputs[fname,"sampleMax"]
       n2 <- inputs[fname,"xformedMax"]
       lnNode1 <- XML::xmlNode("LinearNorm",attrs=c(orig=o1,norm=n1))
       lnNode2 <- XML::xmlNode("LinearNorm",attrs=c(orig=o2,norm=n2))
       ncNode <- XML::append.XMLNode(ncNode, lnNode1)
       ncNode <- XML::append.XMLNode(ncNode, lnNode2)
       dfNode <- XML::append.XMLNode(dfNode, ncNode)
#       local.transformations <- append.XMLNode(local.transformations, dfNode)
      } else if(inputs[fname,"transform"] == "MapValues")
      {
       map <- inputs[fname,"fieldsMap"][[1]]

       dtype <- map[2,ncol(map)]
       if(dtype == "numeric")
       {
	dtype <- "double"
	otype <- "continuous"
       } else if(dtype == "boolean")
       {
        dtype <- "boolean"
        otype <- "categorical" 
       } else 
       {
        dtype <- "string"
	otype <- "categorical"
       }

       dfNode <- XML::xmlNode("DerivedField",attrs=c(name=fname,dataType=as.character(dtype),optype=otype))
       default <- inputs[fname,"default"]
       missing <- inputs[fname,"missingValue"]
       if(dtype == "boolean")
       {
	if((default==1) || (toupper(default)==TRUE))
	{
	 default="true"
	} else
	{
	 default="false"
	}
        if((missing==1) || (toupper(missing)==TRUE))
        {
         missing="true"
        } else
        {
         missing="false"
        }
       }

       if(!is.na(default) && !is.na(missing))
       {
	mapvNode <- XML::xmlNode("MapValues",attrs=c(mapMissingTo=missing,defaultValue=default,outputColumn="output"))
       } else if(!is.na(default) && is.na(missing))
       {
	mapvNode <- XML::xmlNode("MapValues",attrs=c(defaultValue=default,outputColumn="output"))
       } else if(is.na(default) && !is.na(missing))
       {
	mapvNode <- XML::xmlNode("MapValues",attrs=c(mapMissingTo=missing,outputColumn="output"))
       } else
       {
	mapvNode <- XML::xmlNode("MapValues",attrs=c(outputColumn="out"))
       } 

       for(j in 1:(ncol(map)  - 1))
       {
	colname <- paste("input",j,sep="")
	val <- as.character(map[1,j])
	fcpNode <- XML::xmlNode("FieldColumnPair",attrs=c(field=val,column=colname))
	mapvNode <- XML::append.XMLNode(mapvNode,fcpNode)
       }

       inline <- XML::xmlNode("InlineTable")
       for(j in 3:nrow(map))
       {
	row <- XML::xmlNode("row")
	for(k in 1:(ncol(map) - 1))
	{
	 initNode <- XML::xmlNode(paste("input",k,sep=""),value=as.character(map[j,k]))
         row <- XML::append.XMLNode(row, initNode)
	}
	out <- XML::xmlNode("output",value=as.character(map[j,ncol(map)]))
        row <- XML::append.XMLNode(row, out)
	inline <- XML::append.XMLNode(inline,row)
       }

       mapvNode <- XML::append.XMLNode(mapvNode,inline)
       dfNode <- XML::append.XMLNode(dfNode,mapvNode)
      } else if(inputs[fname,"transform"] == "NormDiscrete")
      {
       map <- inputs[fname,"fieldsMap"][[1]]
       dfName <- row.names(inputs)[i]
       missing <- inputs[fname,"missingValue"]
 
       dfNode <- XML::xmlNode("DerivedField",attrs=c(name=dfName,dataType="double",optype="continuous"))
       if(!is.na(missing))
       {
         normNode <- XML::xmlNode("NormDiscrete",attrs=c(field=as.character(inputs[fname,"origFieldName"]),value=as.character(map[1]),mapMissingTo=missing))
       } else
       {
         normNode <- XML::xmlNode("NormDiscrete",attrs=c(field=as.character(inputs[fname,"origFieldName"]),value=as.character(map[1])))
       }
       dfNode <- XML::append.XMLNode(dfNode,normNode)
      } else if(inputs[fname,"transform"] == "discretize")
      {
	maps <- inputs[fname,"fieldsMap"][[1]]
	missingVal <- inputs[fname,"missingValue"]
	defVal <- inputs[fname,"default"]

	origName <- inputs[fname,"origFieldName"] 
	map <- maps[c(-1,-2),]
	dtype <- as.character(inputs[fname,"dataType"])
        if(dtype == "numeric")
        {
         dtype <- "double"
         otype <- "continuous"
        }

	# The following doesnt work as there seems to be no way in PMML to have predicates
	# which indicate if a boolean variable is true or false; and this issue comes up
	# when derived fields of type boolean are used in Tree models
	# We cannot use the operator "isIn" as  in 
	#  <... booleanOperator="isIn> <Array type="string>"TRUE"</Attay
	# as then we need an array of type boolean
	# and that is not allowed; and ADAPA complains if boolean vaiable is tested as being 
	# contained in an Arraya of type string  
	# else if(dtype == "boolean")

        # {
        # dtype <- "boolean"
        # otype <- "categorical"
        # }

	else
        {
         dtype <- "string"
         otype <- "categorical"
        }

        if(dtype == "boolean")
        {
         if((default==1) || (toupper(default)==TRUE))
         {
          default="true"
         } else
         {
          default="false"
         }
         if((missing==1) || (toupper(missing)==TRUE))
         {
          missing="true"
         } else
         {
          missing="false"
         }
        }

	dfNode <- XML::xmlNode("DerivedField",attrs=c(name=fname,dataType=dtype,optype=otype))
        if(!is.na(defVal) && !is.na(missingVal))
        {
         discNode <- XML::xmlNode("Discretize",attrs=c(field=origName,mapMissingTo=missingVal,defaultValue=defVal)) 
        } else if(!is.na(defVal) && is.na(missingVal))
        {
         discNode <- XML::xmlNode("Discretize",attrs=c(field=origName,defaultValue=defVal))
        } else if(is.na(defVal) && !is.na(missingVal))
        {
         discNode <- XML::xmlNode("Discretize",attrs=c(field=origName,mapMissingTo=missingVal))
        } else
        {
         discNode <- XML::xmlNode("Discretize",attrs=c(field=origName))
        }

	for(i in 1:nrow(map))
	{ 
	 dbinNode <- XML::xmlNode("DiscretizeBin",attrs=c(binValue=map[i,2]))
	 clsr <- paste(map[1,3],map[i,5],sep="")
	 if(!is.na(map[i,4]))
  	 {
	  if(!is.na(map[i,6]))
	  {
	   intrNode <- XML::xmlNode("Interval",attrs=c(closure=clsr,leftMargin=map[i,4],rightMargin=map[i,6]))
	  } else
	  {
	   intrNode <- XML::xmlNode("Interval",attrs=c(closure=clsr,leftMargin=map[i,4]))
	  }
	 } else
	 {
	  intrNode <- XML::xmlNode("Interval",attrs=c(closure=clsr,rightMargin=map[i,6]))
	 }
	 dbinNode <- XML::append.XMLNode(dbinNode,intrNode)
	 discNode <- XML::append.XMLNode(discNode,dbinNode)
	}
	dfNode <- XML::append.XMLNode(dfNode,discNode) 

      }

      if(is.null(LTelement))
      {
       local.transformations <- XML::append.XMLNode(local.transformations, dfNode)
      } else
      {
       LTelement <- XML::append.XMLNode(LTelement, dfNode)
      }
     }
     }
    }
  }

  if(is.null(LTelement))
  {
   return(local.transformations)
  } else
  {
   return(LTelement)
  }
}

#####################################################################
# PMML Output element

.pmmlOutput <- function(field, target=NULL, optype=NULL)
{
  requireNamespace("XML", quietly = TRUE)
  
  number.of.fields <- length(field$name)

  output <- XML::xmlNode("Output")
  output.fields <- list()

  for (i in 1:number.of.fields)
  {
    if (field$name[i]==target)
    {
      targetout = target
      if(length(grep("as\\.factor\\(",targetout)) == 1)
      {
        targetout <- gsub("as.factor\\((\\w*)\\)","\\1", targetout, perl=TRUE)
      }

      if (is.null(optype))
        output.fields[[1]] <- XML::xmlNode("OutputField",
                                      attrs=c(name=gsub(" ","",paste("Predicted_",targetout)),
                                        feature="predictedValue"))
      else
        output.fields[[1]] <- XML::xmlNode("OutputField",
                                      attrs=c(name=gsub(" ","",paste("Predicted_",targetout)),
                                        optype=optype,
                                        dataType=ifelse(optype=="continuous",
                                          "double", "string"),
                                        feature="predictedValue"))

     {
      for (j in seq_along(field$levels[[field$name[i]]]))
        output.fields[[j+1]] <- XML::xmlNode("OutputField",
                                        attrs=c(name=paste("Probability_",
                                                  field$levels[[field$name[i]]][j],
                                                  sep=""),
                                          optype="continuous",
                                          dataType = "double",
                                          feature="probability",
                                          value= field$levels[[field$name[i]]][j]))
     }
    }
  }
  
  output$children <- output.fields
  return(output)
}

#######################################################################
# MAIN PMML FUNCTION
#
# modified 081513 to add implementation of transformations generator
# and transformations element addition;  tridivesh.jena@zementis.com
# 

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
  requireNamespace("XML", quietly = TRUE)
  
  number.of.fields <- length(field$name)

  if(field$name[1] == "ZementisClusterIDPlaceHolder" || field$name[1] == "ZementisHiddenTargetField")
  {
   begin<-2
  } else 
  {
   begin <- 1
  }

  namelist <- list()
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
     {

      fName <- field$name[i]
      if(!is.na(field$class[fName]) && field$class[fName] == "factor")
      {
       optypelist[[fName]] <- "categorical"
      }

      if(length(grep("as\\.factor\\(",field$name[i])) == 1)
        fName <- gsub("as.factor\\((\\w*)\\)","\\1", field$name[i], perl=TRUE)

      if(!is.na(field$class[fName]) && field$class[fName] == "factor")
      {
       optypelist[[fName]] <- "categorical"
      }


      if(!(fName %in% namelist) && fName != "ZementisClusterIDPlaceHolder")
      {
       namelist <- c(namelist,fName)
      }

     }
  }

  # DataDictionary

  data.dictionary <- XML::xmlNode("DataDictionary",
                             attrs=c(numberOfFields=length(namelist)))

  if (! is.null(weights) && length(weights))
    data.dictionary <- XML::append.XMLNode(data.dictionary, XML::xmlNode("Extension",
                                                              attrs=c(name="Weights",
                                                              value=weights, extender="Rattle")))

  nmbr <- 1

  for(ndf2 in 1:length(namelist))
  {
    optype <- optypelist[[namelist[ndf2][[1]]]]

    datype <- datypelist[[namelist[ndf2][[1]]]]
    data.fields[[nmbr]] <- XML::xmlNode("DataField", attrs=c(name=namelist[ndf2],
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

    interval <-  XML::xmlNode("Interval",
                           attrs=c(closure="closedClosed", 
			     leftMargin=min(dataval, na.rm=TRUE), # 091025 Handle missing values
                             rightMargin=max(dataval, na.rm=TRUE))) # 091025 Handle missing values
    data.fields[[nmbr]] <- XML::append.XMLNode(data.fields[[nmbr]], interval)
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
         data.fields[[nmbr]][[j]] <- XML::xmlNode("Value",
                          attrs=c(value=.markupSpecials(lev[j])))
       }
     } else
     {
       for (j in seq_along(field$levels[[namelist[nmbr][[1]]]]))
       {
         data.fields[[nmbr]][[j]] <- XML::xmlNode("Value",
                          attrs=c(value=.markupSpecials(field$levels[[namelist[nmbr][[1]]]][j])))
       }
     }
   }

   data.dictionary <- XML::append.XMLNode(data.dictionary, data.fields[[nmbr]])
   nmbr <- nmbr + 1
  }

  return(data.dictionary)
}

# MAIN PMML FUNCTION
#
# modified 081513 to add implementation of transformations generator
# and transformations element addition;  tridivesh.jena@zementis.com
# 
.pmmlMiningSchema <- function(field, target=NULL, transformed=NULL)
{
	requireNamespace("XML", quietly = TRUE)
	
    namelist <- .origFieldList(field, transformed)

    mining.schema <- XML::xmlNode("MiningSchema")
    target <- .removeAsFactor(target)  

    for(j in 1:length(namelist)) {
  
        if(!is.na(namelist[j])) {        
            usage <- ifelse(namelist[j] == target, "predicted", "active")
     
            if(namelist[j]=="Temp" || namelist[j]=="DiscretePlaceHolder") {
            # If field name is the naive bayes categorical field place holder, add missingValueReplacement
                if(length(field$levels[[namelist[j]]])==1) {
                     mf <- XML::xmlNode("MiningField", attrs=c(name=namelist[j],
                            usageType=usage,missingValueReplacement=field$levels[[namelist[j]]]))
                }
            } else {
                mf <- XML::xmlNode("MiningField", attrs=c(name=namelist[j], usageType=usage))
            }
            
            mining.schema <- XML::append.XMLNode(mining.schema, mf)
         }
    }

    return(mining.schema)
}

.origFieldList <- function(field, transformed=NULL)
{
  requireNamespace("XML", quietly = TRUE)
  
  # Create a list of original field names from which any input fields may be derived from 
  number.of.fields <- length(field$name)
  mining.fields <- list()

  if(field$name[1] == "ZementisClusterIDPlaceHolder" || field$name[1] == "ZementisHiddenTargetField")
  {
   begin <- 2
  } else
  {
   begin <- 1
  }

  mining.schema <- XML::xmlNode("MiningSchema")
  namelist <- NULL
  for (i in begin:number.of.fields)
  {
    if(!is.null(transformed) && i!=1)
    {
       if(is.na(transformed$fieldData[field$name[i],"origFieldName"]))
       {
         if(!(field$name[i] %in% namelist))
         {
          namelist <- c(namelist,field$name[i])
         }
       
       } else
       {
           ofname <- transformed$fieldData[field$name[i],"origFieldName"][[1]]
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

              for(k in 1:length(parents))
              {
               if(!(parents[k] %in% namelist))
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
    } else {

        fName <- .removeAsFactor(field$name[i])  
     
        if(!(fName %in% namelist) && fName != "ZementisClusterIDPlaceHolder") {
            namelist <- c(namelist,fName)
        }
    }
    
  }

  return(namelist)
}


.removeAsFactor <- function(fieldName)
{
    if(length(grep("as\\.factor\\(",fieldName)) == 1)
    {
        fieldName <- gsub("as.factor\\((\\w*)\\)","\\1", fieldName, perl=TRUE)
    } 
    return (fieldName)
}

.getNamespace <- function(x)
{
  return(paste("http://sci2s.ugr.es/dicits/-",x,sep=""))
}


