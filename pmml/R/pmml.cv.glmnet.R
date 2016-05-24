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
# 
# Author: Tridivesh Jena
#
# Implemented: 080201 by Tridivesh Jena (info@zementis.com) to add the
# capability to export glmnet models.

pmml.cv.glmnet <- function(model,
                    model.name="Elasticnet_Model",
                    app.name="Rattle/PMML",
                    description="Generalized Linear Regression Model",
                    copyright=NULL,
                    transforms=NULL,
		    unknownValue=NULL,
		     dataset=NULL,
		     s=NULL,
                    ...)
{
  if (! inherits(model, "cv.glmnet")) stop("Not a legitimate cross-validated glmnet object")

  # Collect the required information.
  
  # the fitted model   
  fitmodel <- model$glmnet.fit 

  # the lambda sequence and the lambda resulting in the minimum covariance

  lambda <- model$lambda
  if(!is.null(s))
  {
    minlambda <- s
  } else
  {
    minlambda <- model$lambda.1se
  }

  # array index where lambda less or equals minlambda

  precision <- 0.000000001
  exact <- TRUE
  index <- 1
  index1 <- 1
  index2 <- 1
  if(lambda[1] < minlambda)
  {
    index <- 1
  } else if(lambda[length(lambda)] > minlambda)
  {
    index <- length(lambda)
  } else
  { 
    for(i in 1:length(lambda))
    {
      if(abs(lambda[i] - minlambda) < precision)
      {
        index <- i
        break
      }
      if(lambda[i] <= minlambda)
      {
        index1 <- i
        index0 <- i-1
	exact <- FALSE
        break 
      }
    }
  }

  # model distribution
  
  name <- attributes(model$name)$names
  if(grepl("deviance",name))
  {
    type <- "poisson"
  } else if(grepl("mse",name))
  {
    type <- "gaussian"
  } else
  {
    stop("Only poisson and gaussian family types supported")
  }

  # get regression information
   beta <- fitmodel$beta
   varnames <- attributes(beta)$Dimnames[[1]]

##new 
   if(!is.null(transforms))
   {
    for(i in 1:length(varnames))
    {
     varnames[i] <- row.names(transforms$fieldData)[i] 
    }
   }

  if(exact)
  {
    # intercept 
    intercept <- fitmodel$a0[index]

    # the coefficients
    coeffs <- beta[,index]
  } else
  {
    # intercept
    v0 <- fitmodel$a0[index0] 
    v1 <- fitmodel$a0[index1]
    l <- minlambda
    l0 <- lambda[index0]
    l1 <- lambda[index1]
    intercept <- (v1*(l-l0) - v0*(l-l1))/(l1-l0)

    # the coefficients
    v0 <- beta[,index0]
    v1 <- beta[,index1]
    l <- minlambda
    l0 <- lambda[index0]
    l1 <- lambda[index1]
    coeffs <- (v1*(l-l0) - v0*(l-l1))/(l1-l0)
  }

  #TODO: get alpha,df,weights? ; support offset, 
  # manual how to numerify categorical variables(NormDiscrete?)
  # if s between given s's...linearly interpolate
 
  class <- NULL  
  field <- NULL
  field$name <- c("predictedScore", varnames)

# get field class from the modelling dataset, if provided
# else assume it is numeric
  for(i in 1:length(field$name))
  {
    if(is.null(dataset))
    {
      class <- c(class,"numeric") 
    } else 
    {
      class <- c(class,class(dataset[,i]))
    }
  }

  field$class <- class
  names(field$class) <- field$name

  number.of.fields <- length(field$name)
  target <- "predictedScore" 

  # 110113 For the following grep, with certain Japanes characters we
  # see the string including a "[" in some encoding and causes the
  # grep to fail. We can't do the usual Encoding<- "UTF-8" trick since
  # the characters already look like UTF-8. But using enc2utf8 works -
  # does it hurt doing it always, or just when we have Japanese? Needs
  # testing.

  field$name <- enc2utf8(field$name)
  
# not checked; not used in present implementation of regression models
# kept for future use if multinomial models are implemented
  if(FALSE){
  ylevels <- FALSE
  numfac <- 0
  for (i in 1:number.of.fields)
  {
    # If target variable is binomial, get its categories
    # otherwise stop
    if (field$class[[field$name[i]]] == "factor")
    {
      numfac <- numfac + 1
      if (field$name[i] == target)
      {
        if(length(levels(model$data[[field$name[i]]])) != 2)
          stop("binomial family with more than two target categories is not
                currently supported by PMML")
        ylevels <- TRUE
        field$levels[[field$name[i]]] <- levels(model$data[[field$name[i]]])
      } else
      {
        field$levels[[field$name[i]]] <- model$xlevels[[field$name[i]]]

        #Tridi 2/16/12: remove any 'as.factor' from field names
        if (length(grep("^as.factor\\(", field$name[i])))
        {
          field$name[i] <- sub("^as.factor\\((.*)\\)", "\\1", field$name[i])
          names(field$class)[i] <- sub("^as.factor\\((.*)\\)", "\\1", names(field$class)[i])
          names(field$levels)[numfac] <- sub("^as.factor\\((.*)\\)", "\\1", names(field$levels)[numfac])
        }
      }
    }
  }}

  # PMML

  pmml <- .pmmlRootNode("4.2")

  # PMML -> Header

  pmml <- append.XMLNode(pmml, .pmmlHeader(description, copyright, app.name))

  # PMML -> DataDictionary

  pmml <- append.XMLNode(pmml, .pmmlDataDictionary(field, transformed=transforms))

# following not used or checked. Not used for present implementation of regression models
# kept for possible future use if multinomial general regression models are implemented

 if(FALSE){ 
# determine the distribution and link function to add. quasi distributions cannot be
#  listed. Certain link functions are not supported and an error must be thrown. Certain
# link function can be recast as a power link and must be constructed separately
  # PMML -> GeneralRegressionModel
  add <- FALSE
  addl <- FALSE
  addl2 <- FALSE
  if (model$call[[1]] == "glm"){
    model.type <- model$family$family
    model.link <- model$family$link
  }
  else
    model.type <- "unknown"

# Only binary categorical cases can be handled. For multinomial cases, glm assumes the first
#  category as one and all the rest together as one. The output then is the probability of the 
#  first category NOT be true. This case is not implemented.
  categ <- FALSE
  if(ylevels)
  {
   if(model.type == "binomial")
   {
     categ <- TRUE
     add <- TRUE
   }
  }

  if(model.type == "binomial"){
    add <- TRUE
  }
  if(model.type == "Gamma") {
    model.type <- "gamma"
    add <- TRUE
  }
  if(model.type == "inverse.gaussian") {
    model.type <- "igauss"
    add <- TRUE
  }
  if(model.type == "gaussian") {
    model.type <- "normal"
    add <- TRUE
  }
  if(model.type == "poisson") { 
    add <- TRUE
  } 

  if(model.link == "cloglog") {
    addl <- TRUE
  } else
   if(model.link == "identity") {
    addl <- TRUE
  } else
  if(model.link == "log") {
    addl <- TRUE
  } else
  if(model.link == "logit") {
    addl <- TRUE
  } else
  if(model.link == "probit") {
    addl <- TRUE
  } else
  if(model.link == "inverse") {
    addl <- TRUE
    addl2 <- TRUE
    d <- "-1"
  } else
  if(model.link == "sqrt") {
    addl <- TRUE
    addl2 <- TRUE
    d <- "0.5"
  } else {
    stop("link function currently not supported by PMML")
  }

  if(categ)
  {
     the.model <- xmlNode("GeneralRegressionModel",
                        attrs=c(modelName=model.name,
                          modelType="generalizedLinear",
                          functionName="classification",
                          algorithmName="glm",
                          distribution=model.type,
                          linkFunction=model.link))

  } else if(add && addl && addl2)
  {
    the.model <- xmlNode("GeneralRegressionModel",
                         attrs=c(modelName=model.name,
                           modelType="generalizedLinear",
                           functionName="regression",
                           algorithmName="glm",
                           distribution=model.type,
                           linkFunction="power",
                            linkParameter=d))
  } else if(add && addl && !addl2)
  {
    the.model <- xmlNode("GeneralRegressionModel",
                         attrs=c(modelName=model.name,
                           modelType="generalizedLinear",
                           functionName="regression",
                           algorithmName="glm",
                           distribution=model.type,
                           linkFunction=model.link))
  } else if(!add && addl && addl2)
  {
    the.model <- xmlNode("GeneralRegressionModel",
                         attrs=c(modelName=model.name,
                           modelType="generalizedLinear",
                           functionName="regression",
                           algorithmName="glm",
                           linkFunction="power",
                            linkParameter=d))
  } else if(!add && addl && !addl2)
  {
    the.model <- xmlNode("GeneralRegressionModel",
                         attrs=c(modelName=model.name,
                           modelType="generalizedLinear",
                           functionName="regression",
                           algorithmName="glm",
                           linkFunction=model.link))
  } else 
    stop("model type not supported")

 }

 the.model <- xmlNode("GeneralRegressionModel",
                         attrs=c(modelName=model.name,
                           modelType="generalLinear",
                           algorithmName="glmnet",
                           functionName="regression"))


  extensionNode <- xmlNode("Extension",attrs=c(name="lambda",value=minlambda))
  the.model <- append.XMLNode(the.model,extensionNode)

  # PMML -> RegressionModel -> MiningSchema

  the.model <- append.XMLNode(the.model, .pmmlMiningSchema(field, target, transforms, unknownValue=unknownValue))

  outn <- xmlNode("Output")
  outpn <- xmlNode("OutputField",attrs=c(name="predictedValue",feature="predictedValue"))
  outn <- append.XMLNode(outn,outpn)
  if(type=="poisson")
  {
    outpn <- xmlNode("OutputField",attrs=c(name="predictedMean",feature="transformedValue"))
    applyNode <- xmlNode("Apply",attrs=c("function"="exp"))
    fldNode <- xmlNode("FieldRef",attrs=c(field="predictedValue"))
    applyNode <- append.XMLNode(applyNode,fldNode)
    outpn <- append.XMLNode(outpn,applyNode)
    outn <- append.XMLNode(outn,outpn) 
  }

  the.model <- append.XMLNode(the.model, outn)

  #--------------------------------------------
  # PMML -> Model -> LocalTransforms

  # test of Zementis xform functions
  if(!is.null(transforms))
  {
    the.model <- append.XMLNode(the.model, .pmmlLocalTransformations(field, transforms, NULL))
  }

  plNode <- xmlNode("ParameterList")
  pnode <- xmlNode("Parameter",attrs=c(name="p0",label="Intercept"))
  plNode <- append.XMLNode(plNode,pnode)

  for(i in 1:length(coeffs))
  {
   pnode <- xmlNode("Parameter",attrs=c(name=paste("p",i,sep=""),label=varnames[i]))
   plNode <- append.XMLNode(plNode,pnode)
  }

  the.model <- append.XMLNode(the.model,plNode)

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
# numerical terms
  for(i in 1:length(coeffs)){
    ppcell <- xmlNode("PPCell",attrs=c(value="1",predictorName=varnames[i],
                      parameterName=paste("p",i,sep="")))
    ppm <- append.XMLNode(ppm,ppcell)
  }

  the.model <- append.XMLNode(the.model,ppm)

  pmNode <- xmlNode("ParamMatrix")
  pcNode <- xmlNode("PCell",attrs=c(parameterName="p0", df="1",beta=intercept[[1]]))
  pmNode <- append.XMLNode(pmNode,pcNode)
  for(i in 1:length(coeffs))
  {
   if((!is.na(coeffs[i])) && (coeffs[i] != 0) )
   {
     pcNode <- xmlNode("PCell",attrs=c(parameterName=paste("p",i,sep=""), df="1",
					beta=coeffs[[i]]))
     pmNode <- append.XMLNode(pmNode,pcNode)
   }
  } 

  the.model <- append.XMLNode(the.model,pmNode)

  # Add to the top level structure.
  
  pmml <- append.XMLNode(pmml, the.model)
  
  return(pmml)
}
