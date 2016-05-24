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
# Date: Sep 2013
#------------------------------------------------------------------------------------

pmml.ada <- function(model,
                       model.name="AdaBoost_Model",
                       app.name="R-PMML",
                       description="AdaBoost Model",
                       copyright=NULL,
                       transforms=NULL,
		       unknownValue=NULL,
                        ...)
{
  if (! inherits(model, "ada")) stop("Not a legitimate ada object")
    
   #---------------------------------------------------
   # Collect the required fields information.
 
   noTargetFields <- NULL
   noTargetFields$name <- as.character(attr(model$terms, "variables"))[-1]
   noTargetFields$class <- attr(model$terms, "dataClasses")
  
   noTargetFields$name[1] <- "ZementisHiddenTargetField"
   noTargetFields$class[1] <- "numeric" 
   names(noTargetFields$class)[1] <- "ZementisHiddenTargetField"
   
   target <- as.character(attr(model$terms, "variables"))[-1][1]
   target <- .removeAsFactor(target) 
   
   #----------------------------------------------------------
   # PMML
  
   pmml <- .pmmlRootNode("4.2")
  
   #----------------------------------------------------------
   # PMML -> Header

   pmml <- append.XMLNode(pmml, .pmmlHeader(description, copyright, app.name))

   #-----------------------------------------------------------
   # PMML -> DataDictionary

   pmml <- append.XMLNode(pmml, .pmmlDataDictionary(noTargetFields, transformed=transforms))
  
   #------------------------------------------------------------- 
   # PMML -> MiningModel
   
   #Wen: Even though ADA is a binary classifier, we set the functionName as regression here.
   #     This is because the output category needs to be obtained through post-processing
   #     the raw numeric Adapa output, the weighted sum of each small tree.
   
   functionName <- "regression"
   xmlModel <- xmlNode("MiningModel",
                      attrs=c(modelName=model.name,
                        functionName=functionName, # Required
                        algorithmName="ADA_BOOST"
                        ))
   
   #---------------------------------------------------------------
   # PMML -> MiningModel -> MiningSchema

   xmlModel <- append.XMLNode(xmlModel, .pmmlMiningSchema(noTargetFields, target=NULL, transformed=transforms, unknownValue=unknownValue))
   
   #---------------------------------------------------------------
   # PMML -> MiningModel -> Output
   
   weightSum <- 0.0
   for(i in 1:length(model$model$tree)) {
     weightSum <- weightSum + model$model$alpha[i]
   }

   xmlOutput <- xmlNode("Output")
   
   #-----------------------
   xmlOF_raw <- xmlNode("OutputField", attrs=c(name="rawValue",feature="predictedValue",dataType="double",optype="continuous"))
   xmlOutput <- append.XMLNode(xmlOutput, xmlOF_raw)
  
   #-----------------------
   xmlOF_boost <- xmlNode("OutputField", attrs=c(name="boostValue",feature="transformedValue"))
   xmlApply <- xmlNode("Apply", attrs=c('function'="*"))
   xmlFR <-  xmlNode("FieldRef", attrs=c(field="rawValue"))
   xmlConst <- xmlNode("Constant", weightSum) 
      
   xmlApply <- append.XMLNode(xmlApply, xmlFR)
   xmlApply <- append.XMLNode(xmlApply, xmlConst)
   xmlOF_boost <- append.XMLNode(xmlOF_boost, xmlApply)
      
   xmlOutput <- append.XMLNode(xmlOutput, xmlOF_boost)
   
   #------------------------
   newName <- paste("Predicted_",target,sep="")
   
   xmlOF_cate <- xmlNode("OutputField", attrs=c(name=newName,feature="transformedValue"))
   xmlApply_if <- xmlNode("Apply", attrs=c('function'="if"))
   
   xmlApply_compare <- xmlNode("Apply", attrs=c('function'="lessThan"))
   xmlFR_boost <-  xmlNode("FieldRef", attrs=c(field="boostValue"))
   xmlConst_boost <- xmlNode("Constant", "0.0")
   xmlApply_compare <- append.XMLNode(xmlApply_compare, xmlFR_boost)
   xmlApply_compare <- append.XMLNode(xmlApply_compare, xmlConst_boost)
   
   xmlApply_if <- append.XMLNode(xmlApply_if, xmlApply_compare)
   xmlConst_cate1 <- xmlNode("Constant", levels(model$fit)[1])
   xmlConst_cate2 <- xmlNode("Constant", levels(model$fit)[2])
   
   xmlApply_if <- append.XMLNode(xmlApply_if, xmlConst_cate1)
   xmlApply_if <- append.XMLNode(xmlApply_if, xmlConst_cate2)
   
   xmlOF_cate <- append.XMLNode(xmlOF_cate, xmlApply_if)
   xmlOutput <- append.XMLNode(xmlOutput, xmlOF_cate)
   
   #---------------------------
   prob1 <- paste("Probability_",levels(model$fit)[1],sep="")
   xmlOF_prob1 <- xmlNode("OutputField", attrs=c(name=prob1,feature="transformedValue"))
   
   xmlApply_multiply1 <- xmlNode("Apply", attrs=c('function'="*"))
   xmlConst_multiply1 <- xmlNode("Constant", "2.0")
   xmlFR_multiply1 <-  xmlNode("FieldRef", attrs=c(field="boostValue"))
   xmlApply_multiply1 <- append.XMLNode(xmlApply_multiply1, xmlConst_multiply1)
   xmlApply_multiply1 <- append.XMLNode(xmlApply_multiply1, xmlFR_multiply1)
  
   xmlApply_exp1 <- xmlNode("Apply", attrs=c('function'="exp"))
   xmlApply_exp1 <- append.XMLNode(xmlApply_exp1, xmlApply_multiply1)
   
   xmlApply_sum1 <- xmlNode("Apply", attrs=c('function'="+"))
   xmlConst_sum1 <- xmlNode("Constant", "1.0")
   xmlApply_sum1 <- append.XMLNode(xmlApply_sum1, xmlConst_sum1)
   xmlApply_sum1 <- append.XMLNode(xmlApply_sum1, xmlApply_exp1)
   
   xmlApply_divide1 <- xmlNode("Apply", attrs=c('function'="/"))
   xmlConst_divide1 <- xmlNode("Constant", "1.0")
   xmlApply_divide1 <- append.XMLNode(xmlApply_divide1, xmlConst_divide1)
   xmlApply_divide1 <- append.XMLNode(xmlApply_divide1, xmlApply_sum1)
  
   xmlOF_prob1 <- append.XMLNode(xmlOF_prob1, xmlApply_divide1)
   xmlOutput <- append.XMLNode(xmlOutput, xmlOF_prob1) 
   
   #-----------------------------------------------------
   prob2 <- paste("Probability_",levels(model$fit)[2],sep="")
   xmlOF_prob2 <- xmlNode("OutputField", attrs=c(name=prob2,feature="transformedValue"))
   
   xmlApply_minus <- xmlNode("Apply", attrs=c('function'="-"))
   xmlConst_minus <- xmlNode("Constant", "1.0")
   xmlFR_minus <-  xmlNode("FieldRef", attrs=c(field=prob1))
   xmlApply_minus <- append.XMLNode(xmlApply_minus, xmlConst_minus)
   xmlApply_minus <- append.XMLNode(xmlApply_minus, xmlFR_minus)
     
   xmlOF_prob2 <- append.XMLNode(xmlOF_prob2, xmlApply_minus)
   xmlOutput <- append.XMLNode(xmlOutput, xmlOF_prob2) 
   #-----------------------------
   xmlModel <- append.XMLNode(xmlModel, xmlOutput) 
     
   #------------------------------------------------------------------
   # PMML -> MiningModel -> LocalTransformations

   xmlLT <- NULL
   if(!is.null(transforms))
   {
      xmlLT <- .pmmlLocalTransformations(noTargetFields, transforms)
   } else {
      xmlLT <- xmlNode("LocalTransformations")
   }
   
   xmlModel <- append.XMLNode(xmlModel, xmlLT)
   #---------------------------------------------------------------
   
   # PMML -> MiningModel -> Segmentation
    
   xmlSegmentation <- xmlNode("Segmentation", attrs=c(multipleModelMethod="weightedAverage" ))
   
   for(i in 1:length(model$model$tree)) {
   
      xmlSegment <- xmlNode("Segment", attrs=c(id=i,weight=model$model$alpha[i]))
      xmlTrue <- xmlNode("True")
      xmlSegment <- append.XMLNode(xmlSegment, xmlTrue)
      
      #------------------------------------------------------
      # Wen: It is a bit tricky here: each segment tree is a classification model, return either -1 or 1.
      # However the MiningModel (ada algorithm) is doing weighted sum of the returning -1/1
      # So the functionName of MiningModel should be set as regression; the functionName of each
      # segment tree then has to be regression too (they have to be the same as the parent MiningModel
      # per schema request). Even though the tree nodes are written through classification means.
  
      xmlSegmentTree <- xmlNode("TreeModel", attrs=c(functionName=functionName, noTrueChildStrategy="returnLastPrediction"))

      # noTargetFields contains fields directly from the ADA fitting formula,
      # which are the fields will be used in each segment Tree
      
      xmlSegmentTree <- append.XMLNode(xmlSegmentTree, .pmmlMiningSchema(noTargetFields, target=NULL, transformed=NULL, unknownValue=unknownValue))
      
      adaTree <- model$model$tree[[i]]
      xmlNode <- .buildRpartTreeNode(adaTree,"classification")
      
      xmlSegmentTree <- append.XMLNode(xmlSegmentTree, xmlNode)
      xmlSegment <- append.XMLNode(xmlSegment, xmlSegmentTree) 
      xmlSegmentation <- append.XMLNode(xmlSegmentation, xmlSegment) 
   }
    
   xmlModel <- append.XMLNode(xmlModel, xmlSegmentation) 
    
   #-------------------------------------------------------------------
   pmml <- append.XMLNode(pmml, xmlModel)
  
   return(pmml)
}
