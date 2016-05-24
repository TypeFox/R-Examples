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

pmml.itemsets <- function(model,
                       model.name="arules_Model",
                       app.name="Rattle/PMML",
                       description="arules frequent itemsets model",
                       copyright=NULL,
                       transforms=NULL, ...)
{
  
  if (! inherits(model, "itemsets")) stop("Not a legitimate arules itemsets rules object")

  requireNamespace("arules", quietly=TRUE)
  
  ## PMML
  pmml <- .pmmlRootNode("4.0")

  ## PMML -> Header

  pmml <- append.XMLNode(pmml, .pmmlHeader(description, copyright, app.name))

  ## PMML -> DataDictionary
  data.dictionary <- xmlNode("DataDictionary", attrs=c(numberOfFields = 2L))
  data.dictionary <- append.xmlNode(data.dictionary, list(
      xmlNode("DataField", attrs=c(name="transaction",
              optype = "categorical", dataType = "string")),
      xmlNode("DataField", attrs=c(name="item",
              optype = "categorical", dataType = "string"))
  ))

  pmml <- append.XMLNode(pmml, data.dictionary)


  ## model
  quality <- quality(model)
  is <- items(model)

  association.model <- xmlNode("AssociationModel", 
      attrs=c(functionName="associationRules",
          ## fixme: this is currently a hack
          numberOfTransactions=attr(quality(model), "size.data"), 
          numberOfItems=length(arules::itemLabels(model)),
          minimumSupport=min(quality$support),     
          minimumConfidence=0L,
          numberOfItemsets=length(is),     
          numberOfRules=0L))

  ## mining schema
  mining.schema <- xmlNode("MiningSchema")
  mining.schema <- append.xmlNode(mining.schema, list(
      xmlNode("MiningField", attrs = c(name = "transaction",usageType="group")),
      xmlNode("MiningField", attrs = c(name = "item",usageType="active"))
  ))
  
  association.model <- append.xmlNode(association.model, mining.schema)

  ## items
  items <- list()
  il <- .markupSpecials(arules::itemLabels(model))
  for (i in 1:length(il)) 
  items[[i]] <- xmlNode("Item", attrs = list(id = i, value = il[i]))

  association.model <- append.xmlNode(association.model, items)

  ## itemsets
  itemsets <- list()
  sizes <- arules::size(is)
  isl <- arules::LIST(is, decode=FALSE)
  for (i in 1:length(isl)){
      itemsets[[i]] <- xmlNode("Itemset", attrs = list(id = i, 
              numberOfItems = sizes[i], support=quality$support[i])) 
      
      items <- list()
      if(sizes[i] >0)
      for (j in 1:sizes[i])  
      items[[j]] <- xmlNode("ItemRef", attrs = list(itemRef = isl[[i]][j]))

      itemsets[[i]] <- append.xmlNode(itemsets[[i]], items)  
  }
  
  association.model <- append.xmlNode(association.model, itemsets)
  
  ## no rules

  pmml <- append.XMLNode(pmml, association.model)

  return(pmml)
}
