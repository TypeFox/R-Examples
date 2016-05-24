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

########################################################################
# Linear Model PMML exporter
#
# Implemented: 091015 Graham Williams based on pmml.lm

pmml.coxph <- function(model,
                       model.name="CoxPH_Survival_Regression_Model",
                       app.name="Rattle/PMML",
                       description="CoxPH Survival Regression Model",
                       copyright=NULL,
                       transforms=NULL,
		       unknownValue=NULL,
                       ...)
{
  if (! inherits(model, "coxph")) stop("Not a legitimate coxph object")

  # Collect the required information.

  # Tridi Zementis: detect if special terms exist which are not supported 
  vars <- names(attributes(model$terms)$dataClasses)
  coefs <- names(coefficients(model))
  for(i in 1:length(vars))
  {
    if(grepl("cluster\\(",vars[i]) || grepl("tt\\(",vars[i]))
      stop("Special model equation terms 'cluster' and 'tt' not yet supported in PMML")
  }
  for(i in 1:length(coefs))
  {
    if(grepl(":strata\\(",coefs[i]))
      stop("Multiplicative strata variables not yet supported in PMML")
  }
 
  # For a regression, all variables will have been used except those
  # with a NA coefficient indicating singularities. We mark
  # singularities as inactive shortly.

  terms <- attributes(model$terms)
  
  field <- NULL
  field2<-NULL
  numFields <- length(terms$dataClasses)
  numFields2 <- 0
  for(i in 1:numFields)
  {
    if(!grepl("strata",names(terms$dataClasses)[i]) && !grepl("Surv",names(terms$dataClasses)[i]) )
    {
      field$name[i] <- names(terms$dataClasses)[i]
      field$class[i] <- terms$dataClasses[i]
      names(field$class)[i] <- field$name[i]

      field2$name[i] <- names(terms$dataClasses)[i]
      field2$class[i] <- terms$dataClasses[i]
      names(field2$class)[i] <- field2$name[i]
      numFields2 <- numFields2 + 1
    }
  } 
  # 091020 The target field is actually "risk" of the event, and it is
  # not an actual supplied variable. Notice that in pmml.lm we get the
  # target from field$name[1], which in our case here is "Surv(time,
  # status)". We could get the risk score as status - since risk is
  # the probability of the event occuring in comparison to the
  # population. But let's introduce a new variable, "risk" as numeric,
  # as the predicted variable.
  # 11/13/12 Since the output calculated is survival, call the predicted field survival 
  # field$name[1] <- sub(')', '', sub('Surv\\([^,]*, *', '', field$name[1]))
  field$name[1] <- 'survival'
  field$class[1] <- "numeric"
  names(field$class)[1] <- field$name[1]
  field2$name[1] <- 'survival'
  field2$class[1] <- "numeric"
  names(field2$class)[1] <- field2$name[1]
  numFields2 <- numFields2 + 1

  # Tridi Zementis: Include startTime, endTime, status and strata variable if present
  isStrata <- FALSE
  starttimeVar <- FALSE
  endtimeVar <- "" 
  statusVar <- "" 
  strataVar <- "" 
  model.type <- "coxph"

  survObject <- names(terms$dataClasses)[1]
  survObject0 <- gsub("Surv\\(","",survObject)
  survObject1 <- gsub("\\)","",survObject0) 
  survList <- strsplit(survObject1,",")
  if(length(survList[[1]]) == 2)
  {
    endtimeVar <- gsub(" ","",survList[[1]][1])
    statusVar <- gsub(" ","",survList[[1]][2])
  } else 
  {
    starttimeVar <- gsub(" ","",survList[[1]][1])
    endtimeVar <- gsub(" ","",survList[[1]][2])
    statusVar <- gsub(" ","",survList[[1]][3])
  }

  for(i in 1:length(attributes(model$terms)$term.labels))
  {
    termVar <- attributes(model$terms)$term.labels[i]
    if(grepl("strata",termVar) != 0)
    {
      isStrata <- TRUE
      strataVar0 <- gsub("strata\\(","",termVar)
      strataVar <- gsub("\\)","",strataVar0)
  #if someone used as.factor to emphasize the strata is a factor, remove that
      if (length(grep( "^as.factor\\(", strataVar) ))
      {
        strataVar <- sub("^as.factor\\((.*)\\)", "\\1", strataVar)
      }

    }
  }

  if((endtimeVar %in% field$name) || (starttimeVar %in% field$name) || 
	(statusVar %in% field$name) || (strataVar %in% field$name))
  {
    stop("Predictor fields cannot be status, time or strata variables") 
  }

  # 090103 Support transforms if available.
  
  orig.names <- field$name
  orig.class <- field$class

  number.of.fields <- length(field$name)

  target <- field$name[1]

  # 091020 TODO Do we want to do this for coxph? 090501 Identify those
  # who are singularities. For numerics, this is easy since the names
  # are just the variable names. For categorics this gets tricky
  # because the names include the levels. So we need to keep in
  # inactive the actual variable name, if all coefficients for that
  # variable are NAs.
  
  inactive <- names(which(is.na(coef(model))))
  active <- names(which(!is.na(coef(model))))

  # These are the actual variable names.
  
  tmp <- sapply(sapply(field$name, grep, inactive), length)
  inactive.vars <- names(tmp[tmp>0])
  tmp <- sapply(sapply(field$name, grep, active), length)
  active.vars <- names(tmp[tmp>0])
  
  # Now remove any which have any non-NA levels. This final list is
  # passed on as the definitive list of nonactive variables

  inactive <- setdiff(inactive.vars, active.vars)

  # 091020 Do we need to modify this for coxph?
  numfac <- 0 
  for (i in 1:number.of.fields)
  {
    if (field$class[[field$name[i]]] == "factor")
    {
      numfac <- numfac + 1
      # 081004 gjw Test if the data is available in the model, as it
      # would be for a glm (but not an lm), since if the target
      # variable is categoric then the levels are not recorded in
      # xlevels for the target variable, so we will need to get the
      # levels from the data itself.
      if (is.null(model$data))
        field$levels[[field$name[i]]] <- model$xlevels[[field$name[i]]]
      else
        field$levels[[field$name[i]]] <- levels(model$data[[field$name[i]]])

      #Tridi 11/1/12: remove any 'as.factor' from field names
      if (length(grep("^as.factor\\(", field$name[i])))
      {
        field$name[i] <- sub("^as.factor\\((.*)\\)", "\\1", field$name[i])
        names(field$class)[i] <- sub("^as.factor\\((.*)\\)", "\\1", names(field$class)[i])
        names(field$levels)[numfac] <- sub("^as.factor\\((.*)\\)", "\\1", names(field$levels)[numfac])
      }
    }
  }

  # PMML

  pmml <- .pmmlRootNode("4.2")

  # PMML -> Header

  pmml <- append.XMLNode(pmml, .pmmlHeader(description, copyright, app.name))

  # PMML -> RegressionModel

  # Zementis: Different start node depending on existence of startTimeVariable 
  # or baselineStrataVariable attributes 
  if(isStrata)
  {
    if(starttimeVar==FALSE)
    { 
      the.model <- xmlNode("GeneralRegressionModel",
                       attrs=c(modelType="CoxRegression",
                         modelName=model.name,
                         functionName="regression",
                         algorithmName="coxph",
                         endTimeVariable=endtimeVar,
                         statusVariable=statusVar,
			 baselineStrataVariable=strataVar))
        field2$name[numFields2+1] <- endtimeVar
        field2$class[numFields2+1] <- "numeric"
        names(field2$class)[numFields2+1] <- field2$name[numFields2+1]
        field2$name[numFields2+2] <- statusVar
        field2$class[numFields2+2] <- "numeric"
        names(field2$class)[numFields2+2] <- field2$name[numFields2+2]
        field2$name[numFields2+3] <- strataVar
        field2$class[numFields2+3] <- "factor"
        names(field2$class)[numFields2+3] <- field2$name[numFields2+3]
    } else
    {
      the.model <- xmlNode("GeneralRegressionModel",
                       attrs=c(modelType="CoxRegression",
                         modelName=model.name,
                         functionName="regression",
                         algorithmName="coxph",
                         endTimeVariable=endtimeVar,
			 startTimeVariable=starttimeVar,
                         statusVariable=statusVar,
                         baselineStrataVariable=strataVar))
        field2$name[numFields2+1] <- endtimeVar
        field2$class[numFields2+1] <- "numeric"
        names(field2$class)[numFields2+1] <- field2$name[numFields2+1]
        field2$name[numFields2+2] <- starttimeVar
        field2$class[numFields2+2] <- "numeric"
        names(field2$class)[numFields2+2] <- field2$name[numFields2+2]
        field2$name[numFields2+3] <- statusVar
        field2$class[numFields2+3] <- "numeric"
        names(field2$class)[numFields2+3] <- field2$name[numFields2+3]
        field2$name[numFields2+4] <- strataVar
        field2$class[numFields2+4] <- "factor"
        names(field2$class)[numFields2+4] <- field2$name[numFields2+4]

    }
  } else
  {
    if(starttimeVar==FALSE)
    {
      the.model <- xmlNode("GeneralRegressionModel",
                       attrs=c(modelType="CoxRegression",
                         modelName=model.name,
                         functionName="regression",
                         algorithmName="coxph",
                         endTimeVariable=endtimeVar,
                         statusVariable=statusVar))
        field2$name[numFields2+1] <- endtimeVar
        field2$class[numFields2+1] <- "numeric"
        names(field2$class)[numFields2+1] <- field2$name[numFields2+1]
        field2$name[numFields2+2] <- statusVar
        field2$class[numFields2+2] <- "numeric"
        names(field2$class)[numFields2+2] <- field2$name[numFields2+2]

    } else
    {
      the.model <- xmlNode("GeneralRegressionModel",
                       attrs=c(modelType="CoxRegression",
                         modelName=model.name,
                         functionName="regression",
                         algorithmName="coxph",
			 startTimeVariable=starttimeVar,
                         endTimeVariable=endtimeVar,
                         statusVariable=statusVar))
        field2$name[numFields2+1] <- starttimeVar
        field2$class[numFields2+1] <- "numeric"
        names(field2$class)[numFields2+1] <- field2$name[numFields2+1]
        field2$name[numFields2+2] <- endtimeVar
        field2$class[numFields2+2] <- "numeric"
        names(field2$class)[numFields2+2] <- field2$name[numFields2+2]
        field2$name[numFields2+3] <- statusVar
        field2$class[numFields2+3] <- "numeric"
        names(field2$class)[numFields2+3] <- field2$name[numFields2+3]

    }
  }

  # PMML -> DataDictionary

  pmml <- append.XMLNode(pmml, .pmmlDataDictionary(field2,transformed=transforms))

  # PMML -> RegressionModel -> MiningSchema

  the.model <- append.XMLNode(the.model, .pmmlMiningSchema(field2, target, transforms, unknownValue=unknownValue))

  # Tridi Zementis: Add output fields to output both hazard and survival

  output <- xmlNode("Output")

  out1 <- xmlNode("OutputField",attrs=c(name="Predicted_hazard", feature="predictedValue"))
  output <- append.XMLNode(output, out1)

  out2 <- xmlNode("OutputField",attrs=c(name="SurvivalProbability",feature="transformedValue"))
  applyExp <- xmlNode("Apply",attrs=c("function"="exp"))
  const <- xmlNode("Constant","-1.0")
  applyMult <- xmlNode("Apply",attrs=c("function"="*"))
  fieldref <- xmlNode("FieldRef",attrs=c(field="Predicted_hazard"))

  applyMult <- append.XMLNode(applyMult, const)
  applyMult <- append.XMLNode(applyMult, fieldref)
  applyExp <- append.XMLNode(applyExp, applyMult)
  out2 <- append.xmlNode(out2, applyExp)

  output <- append.XMLNode(output, out2)
  the.model <- append.XMLNode(the.model, output)
 
  # PMML -> TreeModel -> LocalTransforms

# Wen ... commented to see if it break anything ... august 2013
#  if (.supportTransformExport(transforms))
#    the.model <- append.XMLNode(the.model, .gen.transforms(transforms))

  # test of Zementis xform functions
  if(!is.null(transforms))
  {
    the.model <- append.XMLNode(the.model, .pmmlLocalTransformations(field, transforms))
  }

 plNode <- xmlNode("ParameterList")
 num <- 0
 for(i in 1:length(names(coefficients(model)))){
   pname <- paste("p",num)
   pname <- gsub(" ","",pname)
   num <- num + 1
   pnode <- xmlNode("Parameter",attrs=c(name=pname,label=names(coefficients(model))[i],
					referencePoint=model$means[[i]]))
   plNode <- append.XMLNode(plNode,pnode)
 }

 the.model <- append.XMLNode(the.model,plNode)

  flNode <- xmlNode("FactorList")
  for(i in 2:number.of.fields)
  {
    if(field$class[i] == "factor")
    {
      pdNode <- xmlNode("Predictor",attrs=c(name=field$name[i]))
      flNode <- append.XMLNode(flNode,pdNode)
    } 
  } 

  the.model <- append.XMLNode(the.model,flNode)

  cvNode <- xmlNode("CovariateList")
  for(i in 2:number.of.fields)
  {
    if(field$class[i] == "numeric")
    {
      pdNode <- xmlNode("Predictor",attrs=c(name=field$name[i]))
      cvNode <- append.XMLNode(cvNode,pdNode)
    }
  }

  the.model <- append.XMLNode(the.model,cvNode)

  ppm <- xmlNode("PPMatrix")
  # for each variable name 
  for(j in 1:length(model$coefficients)){
    # interactive terms 
    if(length(grep(".+:.+",names(coefficients(model))[j])) == 1){
       # for each multiplicative variable name 
       for(k in 1:length(strsplit(names(coefficients(model))[j],":")[[1]])){
         # go through all the fields and find the one matching the field name
         for(f in 2:number.of.fields){
           if(field$class[[field$name[f]]] == "factor"){
             if(length(grep(field$name[f],strsplit(names(coefficients(model))[j],":")[[1]][k])) == 1){
               modfield <- gsub(field$name[f],"",strsplit(names(coefficients(model))[j],":")[[1]][k])
               ppcell <- xmlNode("PPCell",attrs=c(value=modfield,predictorName=field$name[f],
                                                             parameterName=gsub(" ","",paste("p",j-1))))
               ppm <- append.XMLNode(ppm,ppcell) 
             } 
           } else{
              if(length(grep(field$name[f],strsplit(names(coefficients(model))[j],":")[[1]][k])) == 1){
              ppcell <- xmlNode("PPCell",attrs=c(value="1",predictorName=field$name[f],
                                                             parameterName=gsub(" ","",paste("p",j-1))))
              ppm <- append.XMLNode(ppm,ppcell)
           }
         }
       }
     }
   } else {
# categorical terms
       for(f in 2:number.of.fields){

         if(field$class[[field$name[f]]] == "factor"){
           if(length(grep(field$name[f],names(coefficients(model))[j])) == 1){
            if (length(grep("^as.factor\\(", names(coefficients(model))[j]))==1)
            {
             modfield <- sub("^as.factor\\((.*)\\)", "\\1", names(coefficients(model))[j])
             modfield <- gsub(field$name[f],"",modfield)
            } else {
              modfield <- gsub(field$name[f],"",names(coefficients(model))[j])
            }
            ppcell <- xmlNode("PPCell",attrs=c(value=modfield,predictorName=field$name[f],
                                                          parameterName=gsub(" ","",paste("p",j-1))))
            ppm <- append.XMLNode(ppm,ppcell)
           }
         } else{
# numerical terms

             if(length(grep(field$name[f],names(coefficients(model))[j])) == 1){
             ppcell <- xmlNode("PPCell",attrs=c(value="1",predictorName=field$name[f],
                                                           parameterName=gsub(" ","",paste("p",j-1))))
             ppm <- append.XMLNode(ppm,ppcell)
         }
       }
     }
   }
  }
  the.model <- append.XMLNode(the.model,ppm)

  pmNode <- xmlNode("ParamMatrix")
  for(i in 1:length(model$coefficients))
  {
   if( !is.na(coefficients(model)[i]) )
   {
     pcNode <- xmlNode("PCell",attrs=c(parameterName=gsub(" ","",paste("p",i-1)), df="1",
                                        beta=as.numeric(coefficients(model)[i])))
     pmNode <- append.XMLNode(pmNode,pcNode)
   } else
   {
     stop("Model coefficients did not converge and resulted in singular values")
   }
  }
  the.model <- append.XMLNode(the.model,pmNode)

  # Zementis: value which represents an event taking place is not given in the R object, so cannot 
  # make the EventValue node 
  #eventNode<-xmlCommentNode("EventValues expected by R are either [0,1] or [1,2] or [TRUE,FALSE]")
  #the.model <- append.XMLNode(the.model,eventNode)

  # Zementis: call the R basehaz function to get the baselineHazard values
#  CumHazard <- survival:::basehaz(model)
  sfit<-survival::survfit(model)
  H<- -log(sfit$surv)

  strata<-sfit$strata
  if (!is.null(strata))
      strata<- factor(rep(names(strata),strata), levels=names(strata))

  if (is.null(strata))
    CumHazard <- data.frame(hazard=H,time=sfit$time)
  else
    CumHazard <- data.frame(hazard=H,time=sfit$time,strata=strata)
 
  numTime <- length(CumHazard$time)
  levl <- NULL
  baselineNode <- NULL
  maxTim <- 0
 
  if(isStrata)
  {
    baseTable <- xmlNode("BaseCumHazardTables")
    for(i in 1:length(levels(CumHazard$strata)))
    {
      name <- levels(CumHazard$strata)[i]
      for(j in 1:length(CumHazard$hazard))
      {
        if(CumHazard$strata[j]==name)
        {
          levl <- gsub(strataVar,"",name)
          levl <- gsub("=","",levl)
          maxTim <- CumHazard$time[j]
        }
      }
      stratumNode <- xmlNode("BaselineStratum",attrs=c(value=levl,maxTime=maxTim))
      for(j in 1:length(CumHazard$hazard))
      {
	if(CumHazard$strata[j]==name)
	{
	 if(CumHazard$hazard[j] != Inf)
	 {
	  baseCell <- xmlNode("BaselineCell",attrs=c(time=CumHazard$time[j],cumHazard=CumHazard$hazard[j]))
	  stratumNode <- append.XMLNode(stratumNode,baseCell)
	 } else
	 {
	  baseCell <- xmlNode("BaselineCell",attrs=c(time=CumHazard$time[j],cumHazard=1.0E+10))
	  stratumNode <- append.XMLNode(stratumNode,baseCell)
	  warnNode <- xmlCommentNode("Cumulative Hazards approach Infinity after this time")
	  stratumNode <- append.XMLNode(stratumNode,warnNode)
	  break
	 } 
	}
      }
      baseTable <- append.XMLNode(baseTable,stratumNode) 
    }
  } else
  {
    baseTable <- xmlNode("BaseCumHazardTables",attrs=c(maxTime=CumHazard$time[numTime]))
    for(i in 1:length(CumHazard$time))
    {
     if(CumHazard$hazard[i] != Inf)
     {
       baseCell <- xmlNode("BaselineCell",attrs=c(time=CumHazard$time[i],
                                               cumHazard=CumHazard$hazard[i]))
       baseTable <- append.XMLNode(baseTable,baseCell)
     } else
     {
      baseCell <- xmlNode("BaselineCell",attrs=c(time=CumHazard$time[i],cumHazard=1.0E+10))
      baseTable <- append.XMLNode(baseTable,baseCell)
      warnNode <- xmlCommentNode("Cumulative Hazards approach Infinity after this time")
      baseTable <- append.XMLNode(baseTable,warnNode)
      break
     }
    }
  }
  the.model<- append.XMLNode(the.model,baseTable) 
 
  # Add to the top level structure.
  
  pmml <- append.XMLNode(pmml, the.model)
  
  return(pmml)
}
