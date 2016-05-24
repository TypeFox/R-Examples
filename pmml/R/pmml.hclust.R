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

pmml.hclust <- function(model,
                        model.name="HClust_Model",
                        app.name="Rattle/PMML",
                        description="Hierarchical cluster model",
                        copyright=NULL,
                        transforms=NULL,
			unknownValue=NULL,
                        centers,
                        ...)
{

  if (! inherits(model, "hclust")) stop("Not a legitimate hclust object")

  # Collect the required information.

  field <- NULL
  field$name <-  colnames(centers)
  number.of.fields <- length(field$name)
  field$class <- rep("numeric", number.of.fields) # All fields are numeric
  names(field$class) <- field$name

  field2 <- NULL
  field2$name[1] <- "ZementisClusterIDPlaceHolder"
  field2$class[1] <- "ID"
  names(field2$class)[1] <- "ZementisClusterIDPlaceHolder"
  for(i in 1:number.of.fields)
  {
   field2$name[i+1] <- field$name[i]
   field2$class[i+1] <- field$class[i]
   names(field2$class)[i+1] <- names(field$class[i])
  }

  orig.fields <- field$name

  # 090822 Mark any categoric transforms as inactive since they won't
  # have been used in the clustering (at least not until we automate
  # the conversion to indicator variables.

#  for (i in which(sapply(transforms, function(x) x$type) %in%
#                  .TRANSFORMS.TO.CATEGORIC))
#    transforms[[i]]$status <- "inactive"

  number.of.clusters <- nrow(centers)
  cluster.names <- 1:number.of.clusters

  #-------------------------------------------------------------------
  # PMML

  pmml <- .pmmlRootNode("4.2")

  #-------------------------------------------------------------------
  # PMML -> Header

  pmml <- append.XMLNode(pmml, .pmmlHeader(description, copyright, app.name))

  #-------------------------------------------------------------------
  # PMML -> DataDictionary

  pmml <- append.XMLNode(pmml, .pmmlDataDictionary(field2,transformed=transforms))

  # PMML -> ClusteringModel

  cl.model <- xmlNode("ClusteringModel",
                      attrs=c(modelName=model.name,
                        functionName="clustering", # Required
                        algorithmName="HClust",
                        modelClass="centerBased", # Required
                        numberOfClusters=number.of.clusters)) # Required

  # PMML -> ClusteringModel -> MiningSchema

  cl.model <- append.XMLNode(cl.model, .pmmlMiningSchema(field2,transformed=transforms,unknownValue=unknownValue))

  #----------------------------------------------------------------------
  # Outputs
  output <- xmlNode("Output")
  out <- xmlNode("OutputField",attrs=c(name="predictedValue", feature="predictedValue"))
  output <- append.XMLNode(output, out)
  cl.model <- append.XMLNode(cl.model, output)

  #----------------------------------------------------------------------
  # PMML -> ClusteringModel -> LocalTransformations -> DerivedField -> NormContiuous

  # test of Zementis xform functions
  if(!is.null(transforms))
  {
    the.model <- append.XMLNode(the.model, .pmmlLocalTransformations(field2, transforms))
  }
  
  #------------------------------------------------------------------
  # PMML -> ClusteringModel -> ComparisonMeasure
  
  cl.model <- append.XMLNode(cl.model,
                             append.XMLNode(xmlNode("ComparisonMeasure",
                                                    attrs=c(kind="distance")),
                                            xmlNode("squaredEuclidean")))

  #-------------------------------------------------------------------
  # PMML -> ClusteringField

  for (i in orig.fields)
  {
    cl.model <- append.xmlNode(cl.model, xmlNode("ClusteringField", attrs=c(field=i)))
  }
  
  # PMML -> ClusteringModel -> Cluster -> Array
  
  clusters <- list()
  for (i in 1:number.of.clusters)
  {
    cl.model <- append.XMLNode(cl.model,
                               xmlNode("Cluster",  attrs=c(name=cluster.names[i],size=model$size[i]),
                               xmlNode("Array", attrs=c(n=number.of.fields, type="real"),
                                               paste(centers[i,],collapse=" "))))
  }
  pmml <- append.XMLNode(pmml, cl.model)

  return(pmml)
}
