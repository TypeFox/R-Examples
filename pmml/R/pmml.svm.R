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
# Author: Wen Lin
# Date: Aug 2013
#-------------------------------------------------------------------------------------
pmml.svm <- function(model,
                        model.name="LIBSVM_Model",
                        app.name="R-PMML",
                        description="Support Vector Machine Model",
                        copyright=NULL,
                        transforms=NULL,
			unknownValue=NULL,
                        ...)
{
  if (! inherits(model, "svm")) stop("Not a legitimate svm object")

  #---------------------------------------------------
  # Collect the required fields information.
  # field: the list of all fields
  
  field <- NULL
  field$name <- as.character(attr(model$terms, "variables"))[-1]
  field$class <- attr(model$terms, "dataClasses")
  
  target <- as.character(attr(model$terms, "variables"))[-1][1]
  
  #-----------------------------------------------------
  # Determine functionName
  
  functionName <- "classification"
  if(field$class[[1]]=="numeric") {
     functionName <- "regression"
  }

  #----------------------------------------------------------
  # PMML
  
  pmml <- .pmmlRootNode("4.2")
  
  #----------------------------------------------------------
  # PMML -> Header

  pmml <- append.XMLNode(pmml, .pmmlHeader(description, copyright, app.name))

  #-----------------------------------------------------------
  # PMML -> DataDictionary

  pmml <- append.XMLNode(pmml, .pmmlDataDictionary(field, transformed=transforms))

  #------------------------------------------------------------
  # PMML -> SupportVectorMachineModel

  if(model$nclasses > 2) {

       xmlModel <- xmlNode("SupportVectorMachineModel",
                      attrs=c(modelName=model.name,
                        functionName=functionName, # Required
                        algorithmName="LIBSVM", classificationMethod="OneAgainstOne" ))                             
  } else {
       xmlModel <- xmlNode("SupportVectorMachineModel",
                      attrs=c(modelName=model.name,
                        functionName=functionName, # Required
                        algorithmName="LIBSVM" ))
  }
 

  #---------------------------------------------------------------
  # PMML -> SupportVectorMachineModel -> MiningSchema

  xmlModel <- append.XMLNode(xmlModel, .pmmlMiningSchema(field, target, transformed=transforms, unknownValue=unknownValue))

  #-----------------------------------------------------------------
  # PMML -> SupportVectorMachineModel -> Output

  xmlOutput <- NULL
   
  if(functionName == "regression") {

      xmlOutput <- xmlNode("Output")
      xmlOF_predicted <- xmlNode("OutputField", attrs=c(name="predictedValue",feature="predictedValue"))
      xmlOutput <- append.XMLNode(xmlOutput, xmlOF_predicted)
  
      xmlOF <- xmlNode("OutputField", attrs=c(name="svm_predict_function",feature="transformedValue"))
     
      ym<-model$y.scale$`scaled:center`[[1]]
      ys<-model$y.scale$`scaled:scale`[[1]]
      
      xmlApply <- xmlNode("Apply", attrs=c('function'="*"))
      xmlFR <-  xmlNode("FieldRef", attrs=c(field="predictedValue"))
      xmlConst <- xmlNode("Constant", -1*ys) 
      
      xmlApply <- append.XMLNode(xmlApply, xmlFR)
      xmlApply <- append.XMLNode(xmlApply, xmlConst)
      
      xmlApply_sum <- xmlNode("Apply", attrs=c('function'="+"))
      xmlConst_sum <- xmlNode("Constant", ym)
      xmlApply_sum <- append.XMLNode(xmlApply_sum, xmlApply)
      xmlApply_sum <- append.XMLNode(xmlApply_sum, xmlConst_sum )
      
      xmlOF <- append.XMLNode(xmlOF, xmlApply_sum)
      
      xmlOutput <- append.XMLNode(xmlOutput, xmlOF)
      
  } else {
      xmlOutput <- .pmmlOutput(field, target)
  }
  
  xmlModel <- append.XMLNode(xmlModel, xmlOutput)

  #------------------------------------------------------------------
  # PMML -> SupportVectorMachineModel -> LocalTransformations

  xmlLT <- NULL
  if(!is.null(transforms))
  {
     xmlLT <- .pmmlLocalTransformations(field, transforms)
  } else {
     xmlLT <- xmlNode("LocalTransformations")
  }
 
  
  if(is.null(model$x.scale) == FALSE) {
     #NormContinuous xform
     
     num.inputs <- length(model$x.scale[[1]])
     inputNames <- names(model$x.scale[[1]])

     inputDFNames <- array(NA, dim=num.inputs)
     for(i in 1:num.inputs) {
    
       dfName <- paste("algorithm_derived_nc_",inputNames[i],sep="")
       inputDFNames[i] <- dfName
       xmlDF <- xmlNode("DerivedField", attrs=c(name=dfName,dataType="double",optype="continuous"))
       xmlNC <- xmlNode("NormContinuous", attrs=c(field=inputNames[i]))
    
       m<-model$x.scale$`scaled:center`[[i]]
       s<-model$x.scale$`scaled:scale`[[i]]

       xmlLN1 <- xmlNode("LinearNorm", attrs=c(orig=0,norm=-m/s))
       xmlLN2 <- xmlNode("LinearNorm", attrs=c(orig=m,norm=0))
       
       if(m<0) {
          xmlNC <- append.XMLNode(xmlNC, xmlLN2)
          xmlNC <- append.XMLNode(xmlNC, xmlLN1)
       } else {
          xmlNC <- append.XMLNode(xmlNC, xmlLN1)
          xmlNC <- append.XMLNode(xmlNC, xmlLN2)
       }
       
       xmlDF <- append.XMLNode(xmlDF, xmlNC)
       xmlLT <- append.XMLNode(xmlLT, xmlDF)
     }
  }
  
  allVectorAttrName <- attr(model$SV,"dimnames")[2]
  num.vector.attr <- length(allVectorAttrName[[1]])
  
  #--------------------------------
  # vfNames : the array stores the names of all vector fields
  vfNames <- array(NA, dim=length(allVectorAttrName[[1]]))
  vfIndex <- 1
  
  for(i in 2:length(field$name)) {

     inputName <- NULL
     inputName <- field$name[[i]]
     
     if(field$class[[i]]=="numeric") {
        if(is.null(model$x.scale)) {
           vfNames[vfIndex] <- inputName
        } else {
           vfNames[vfIndex] <- paste("algorithm_derived_nc_",inputName,sep="")
        }
        
        vfIndex <- vfIndex + 1
        next
     }

     for(j in 1:num.vector.attr) {
        vectorAttr <- allVectorAttrName[[1]][j]
        
        if(grepl(inputName,vectorAttr) == TRUE) {

           ndValue <- NULL
           ndValue <- gsub(inputName,"",vectorAttr)
           
           if(grepl("_",ndValue) == TRUE) {
              # TO: check the first char instead of contain
              # also, check if ndValue is skipped from the previous one
              next
           }
           
           dfName <- NULL
           dfName <- paste("algorithm_derived_nd_",inputName,"_",ndValue,sep="")

           vfNames[vfIndex] <- dfName
           vfIndex <- vfIndex + 1
     
           xmlDF <- xmlNode("DerivedField", attrs=c(name=dfName,dataType="double",optype="continuous"))
           if (grepl("\\.", ndValue) == TRUE) {
              xmlND <- xmlNode("NormDiscrete", attrs=c(field=inputName,value=ndValue))
              xmlWarning <- newXMLCommentNode(" R Warning: The character '-' in the original data might have been replaced by the '.' character. Check the desired scoring data for consistency")
              xmlND <- append.XMLNode(xmlND, xmlWarning)
           } else {
              xmlND <- xmlNode("NormDiscrete", attrs=c(field=inputName,value=ndValue))
           }
           
           xmlDF <- append.XMLNode(xmlDF, xmlND)
           xmlLT <- append.XMLNode(xmlLT, xmlDF)
        }
     }
  }
 
  xmlModel <- append.XMLNode(xmlModel, xmlLT)
   
  #------------------------------------------------------------------
  # Kernel
  xmlKernel <- NULL
  if (model$kernel == 0) {
     xmlKernel <- xmlNode("LinearKernelType", attrs=c(description="Linear kernel type"))
  } else if (model$kernel == 1) {
     xmlKernel <- xmlNode("PolynomialKernelType", attrs=c(gamma=model$gamma,coef0=model$coef0,degree=model$degree,
     description="Polynomial kernel type"))
  } else if (model$kernel == 3) {
     xmlKernel <- xmlNode("SigmoidKernelType", attrs=c(gamma=model$gamma,coef0=model$coef0,
     description="Sigmoid kernel type"))
  } else {
     xmlKernel <- xmlNode("RadialBasisKernelType", attrs=c(gamma=model$gamma,description="Radial basis kernel type"))
  }
 
  xmlModel <- append.XMLNode(xmlModel, xmlKernel)
  #------------------------------------------------------------------
  # Vector Instances

  vectorSize <- length(model$SV[1,])
  
  xmlVD <- xmlNode("VectorDictionary",attrs=c(numberOfVectors=model$tot.nSV))
 
  xmlVF <- xmlNode("VectorFields",attrs=c(numberOfFields=vectorSize))
  
  for(i in 1:num.vector.attr) {
    xmlFR <- xmlNode("FieldRef", attrs=c(field=vfNames[i]))
    xmlVF <- append.XMLNode(xmlVF, xmlFR)
  }
  
  xmlVD <- append.XMLNode(xmlVD, xmlVF)
  for(i in 1:model$tot.nSV) {
    xmlVI <- xmlNode("VectorInstance", attrs=c(id=i))
    xmlRealSA <- xmlNode("REAL-SparseArray", attrs=c(n=num.vector.attr))
    
    indices <- NULL
    realEntries <- NULL
    
    for(j in 1:num.vector.attr) {
       indices <- paste(indices,j)
       realEntries <- paste(realEntries,model$SV[i,][[j]])
    }

    xmlIndices <- xmlNode("Indices", indices)
    xmlRealEntries <- xmlNode("REAL-Entries", realEntries)
    
    xmlRealSA <- append.XMLNode(xmlRealSA, xmlIndices)
    xmlRealSA <- append.XMLNode(xmlRealSA, xmlRealEntries)
    xmlVI <- append.XMLNode(xmlVI, xmlRealSA)
    xmlVD <- append.XMLNode(xmlVD, xmlVI)
  }
  
  xmlModel <- append.XMLNode(xmlModel, xmlVD)
    
  #------------------------------------------------------------------
  # Support Vector Machines
  
  startVectorIndex <- array(NA, dim=model$nclasses)
  startVectorIndex[1] <- 1
  
  for(i in 2:model$nclasses) {
      startVectorIndex[i] <- startVectorIndex[i-1]+model$nSV[i-1]
  }

  svmCount <- 0
  
  if(functionName == "classification") {
  
    for(i in 1:(model$nclasses-1)) {
       for(j in (i+1):model$nclasses) {
     
          svmCount <- svmCount+1
         
          xmlSVM = xmlNode("SupportVectorMachine", 
                        attrs=c(targetCategory=model$levels[i], alternateTargetCategory=model$levels[j]))

          si <- startVectorIndex[i]
          sj <- startVectorIndex[j]
          ci <- model$nSV[i]
          cj <- model$nSV[j]
         
          coef1Array <- model$coefs[,j-1]
          coef2Array <- model$coefs[,i]

          xmlSVs <- xmlNode("SupportVectors",attrs=c(numberOfAttributes=num.vector.attr, numberOfSupportVectors=ci+cj))
          xmlCFs <- xmlNode("Coefficients", attrs=c(absoluteValue=model$rho[svmCount], numberOfCoefficients=ci+cj))
          
          for(k in 0:(ci-1)) {
             xmlSV <- xmlNode("SupportVector", attrs=c(vectorId=si+k))
             xmlCF <- xmlNode("Coefficient",  attrs=c(value=-1*coef1Array[si+k]))
              
             xmlSVs <- append.XMLNode(xmlSVs, xmlSV)
             xmlCFs <- append.XMLNode(xmlCFs, xmlCF)
          }
          
          for(k in 0:(cj-1)) {
             xmlSV <- xmlNode("SupportVector", attrs=c(vectorId=sj+k))
             xmlCF <- xmlNode("Coefficient",  attrs=c(value=-1*coef2Array[sj+k]))
              
             xmlSVs <- append.XMLNode(xmlSVs, xmlSV)
             xmlCFs <- append.XMLNode(xmlCFs, xmlCF)
          }
     
          xmlSVM <- append.XMLNode(xmlSVM, xmlSVs)
          xmlSVM <- append.XMLNode(xmlSVM, xmlCFs)
          
          xmlModel <- append.XMLNode(xmlModel, xmlSVM) 
      }
    } 
  } else {

       xmlSVM = xmlNode("SupportVectorMachine")
      
       xmlSVs <- xmlNode("SupportVectors",
                   attrs=c(numberOfAttributes=num.vector.attr, numberOfSupportVectors=model$tot.nSV))
       xmlCFs <- xmlNode("Coefficients", 
                   attrs=c(absoluteValue=model$rho[1], numberOfCoefficients=model$tot.nSV))
       
       for(i in 1:model$tot.nSV) {
             xmlSV <- xmlNode("SupportVector", attrs=c(vectorId=i))
             xmlCF <- xmlNode("Coefficient",  attrs=c(value=-1*model$coefs[i]))
              
             xmlSVs <- append.XMLNode(xmlSVs, xmlSV)
             xmlCFs <- append.XMLNode(xmlCFs, xmlCF)
       }
       
       xmlSVM <- append.XMLNode(xmlSVM, xmlSVs)
       xmlSVM <- append.XMLNode(xmlSVM, xmlCFs)
          
       xmlModel <- append.XMLNode(xmlModel, xmlSVM) 
  }
  
  #-------------------------------------------------------------------
  pmml <- append.XMLNode(pmml, xmlModel)
  
  return(pmml)
}

