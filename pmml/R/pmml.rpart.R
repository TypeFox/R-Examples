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
# rpart PMML exporter
#
# Original by Togaware
# Updated by Zementis 080605 to add ScoreDistribution, missingValueStrategy 
# Updated by Zementis 090705 to add Output element
# Updated by Togaware 091003 to add other attributes to Output element
#
# To Do
#   Add the optional ModelStats, ModelExplanation and ModelVerification elements

pmml.rpart <- function(model,
                       model.name="RPart_Model",
                       app.name="Rattle/PMML",
                       description="RPart Decision Tree Model",
                       copyright=NULL,
                       transforms=NULL,
		       unknownValue=NULL,
                       dataset=NULL,
                        ...)
{
  if (! inherits(model, "rpart")) stop("Not a legitimate rpart object")
  requireNamespace("rpart", quietly=TRUE)

  function.name <- "classification"
  if (model$method != "class") function.name <- "regression"
  
  # Collect the required information.

  # We list all variables, irrespective of whether they appear in the
  # final model. This seems to be the standard thing to do with
  # PMML. It also adds extra information - i.e., the model did not
  # need these extra variables!

  field <- NULL
  field$name <- as.character(attr(model$terms, "variables"))[-1]
  field$class <- attr(model$terms, "dataClasses")

  # 100829 Record the variable used for a weight, if there is one.

  weights <- model$call$weights
  if (! is.null(weights))
    weights <- gsub("^\\(|\\)$", "",
                    gsub("crs\\$dataset\\$|\\[.*\\]", "",
                         capture.output(print(weights))))
  
  # 081229 Our names and types get out of sync for multiple transforms
  # on the one variable. TODO How to fix? For now, we will assume a
  # single transform on each variable.

  ofield <- field

  # 090617 Ensure that the list of fields includes those necessary for
  # the transforms. By this stage the transforms should have removed
  # from it any that are not needed in the model.

  number.of.fields <- length(field$name)

  target <- field$name[1]
  
  for (i in 1:number.of.fields)
  {
    # 090829 Zementis Move the test for factor to inside the test for
    # a target.  This will guarantee that for a non-string target in a
    # classification tree model, the probability can still be output
    # for each category.

    # if (field$class[[field$name[i]]] == "factor")

    if (field$name[i] == target)
        field$levels[[field$name[i]]] <- attr(model, "ylevels")
      else
        if (field$class[[field$name[i]]] == "factor")
          field$levels[[field$name[i]]] <- attr(model,"xlevels")[[field$name[i]]]
  }

  # 090519 Identify those variables that are not used in the model. We
  # need to deal with a model that has surrogates. Use maxsurrogate to
  # distinguish. However, when using surrogates perhaps we also need
  # to identify the smaller list of inactive, since there may still be
  # inactive variables. This is not yet implemented.

  if (model$control$maxsurrogate == 0)
  {
    frame <- model$frame
    leaves <- frame$var == "<leaf>"
    used <- unique(frame$var[!leaves])
    inactive <- setdiff(setdiff(levels(used), used), "<leaf>")
  }
  else
    inactive <- NULL
  
  # PMML

  pmml <- .pmmlRootNode("4.2")

  # PMML -> Header

  pmml <- append.XMLNode(pmml, .pmmlHeader(description, copyright, app.name))

  # PMML -> DataDictionary

  pmml <- append.XMLNode(pmml, .pmmlDataDictionary(field, dataset, weights=weights,transformed=transforms))

  # PMML -> TreeModel

  the.model <- xmlNode("TreeModel", attrs=c(modelName=model.name,
                                       functionName=function.name,
                                       algorithmName="rpart",
                                       splitCharacteristic="binarySplit",
                                       missingValueStrategy="defaultChild",
				       noTrueChildStrategy="returnLastPrediction"))

  # PMML -> TreeModel -> MiningSchema
  
  the.model <- append.XMLNode(the.model, .pmmlMiningSchema(field, target, transformed=transforms,unknownValue=unknownValue))

  # PMML -> TreeModel -> Output
  
  the.model <- append.XMLNode(the.model, .pmmlOutput(field, target, switch(function.name,
                             classification="categorical", regression="continuous")))

  # PMML -> TreeModel -> LocalTransformations -> DerivedField -> NormContiuous

  # test of Zementis xform functions
  if(!is.null(transforms))
  {
    the.model <- append.XMLNode(the.model, .pmmlLocalTransformations(field, transforms, NULL))
  }
  
  # PMML -> TreeModel -> Node

  # Collect information to create nodes.
  
  xmlTreeNode <- .buildRpartTreeNode(model,function.name)
  the.model <- append.XMLNode(the.model, xmlTreeNode)

  # Add to the top level structure.

  pmml <- append.XMLNode(pmml, the.model)

  return(pmml)
}


############################################################################
# .buildRpartTreeNode
.buildRpartTreeNode <- function(model,function.name)
{
#  depth <- rpart:::tree.depth(as.numeric(row.names(model$frame)))
  numnode <- as.numeric(row.names(model$frame))
  treedepth <- floor(log(numnode, base=2) + 1e-7)
  depth <- treedepth - min(treedepth)

  count <- model$frame$n
  label <- labels(model, minlength=0, digits=7)
  fieldLabel <- label[1]
  operator <- ""
  value <- "" #list("")

  # Added by Zementis: Information to create nodes.
  ff <- model$frame
  cp <- 0
  id <- as.integer(row.names(ff))
  parent.id <- ifelse(id==1,1, floor(id/2))
  parent.cp <- ff$complexity[match(parent.id, id)]
  rows <- (1:length(id))[parent.cp > cp]
  parent_ii <- 1

  # Check the function.name.
  
  score <- attr(model, "ylevels")[model$frame$yval]
  if(function.name == "regression") score <- ff$yval[rows]

  # Get the information for the primary predicates

  if (length(label) > 1)
  {
    for (i in 2:length(label))
    {
      fieldLabel <-  c(fieldLabel, strsplit(label[i], '>|<|=')[[1]][1])
      op <- substr(label[i], nchar(fieldLabel[i])+1, nchar(fieldLabel[i])+2)
      if (op == ">=")
      {
        operator <- c(operator, "greaterOrEqual")
        value <- c(value, substr(label[i], nchar(fieldLabel[i])+3, nchar(label[i])))
      }
      else if (op == "< ")
      {
        operator <- c(operator, "lessThan")
        value <- c(value, substr(label[i], nchar(fieldLabel[i])+3, nchar(label[i])))
      }
      else if (substr(op, 1, 1) == "=")
      {
        operator <- c(operator, "isIn")
        value <- c(value, substr(label[i], nchar(fieldLabel[i])+2, nchar(label[i])))
      }
    }
    node <- .genBinaryTreeNodes(depth, id, count, score, fieldLabel, operator, value,
                               model, parent_ii, rows,"right")
  }
  else
  {
    node <- .genBinaryTreeNodes(depth, id, count, score, fieldLabel, operator, value,
                               model, parent_ii, rows,"right")
  }


}
############################################################################
# FUNCTION: .genBinaryTreeNodes
#
# Goal: create nodes for the tree (a recursive function)

.genBinaryTreeNodes <- function(depths, ids, counts, scores, fieldLabels,
                               ops, values, model, parent_ii, rows,position)
{
  depth <- depths[1]
  count <- counts[1]
  score <- scores[1]
  fieldLabel <- fieldLabels[1]
  op <- ops[1]
  value <- values[1]

  # Added by zementis

  ff <- model$frame
  id <- ids[1]
  ii <- rows[1]

  # Assign the default child for non-leaf nodes.
  
  if(length(ids) > 1)  # Non-leaf node
  {
    sons <- 2*id + c(0,1)
    sons.n <- ff$n[match(sons, ids)]

    childId = sons[2]
    if(sons.n[1] >= sons.n[2]) childId = sons[1]
 
    node <- xmlNode("Node", attrs=c(id=id, score=score, recordCount=count,
                              defaultChild=childId))
  }
  else  # Leaf node
  {
    node <- xmlNode("Node", attrs=c(id=id,score=score, recordCount=count))
  }

  # Create the predicate for the node

  if (fieldLabel == "root")
  {
    predicate <- xmlNode("True")
  }
  else if (ff$nsurrogate[parent_ii] >0)  # When the node has surrogate predicates
  {  
    predicate <-  xmlNode("CompoundPredicate",attrs=c(booleanOperator = "surrogate"))

    # Add the primary predicate
    
    predicate <- append.XMLNode(predicate,.getPrimaryPredicates(fieldLabel,op,value))

    # Add the surrogate predicates
    
    predicate <- .getSurrogatePredicates(predicate, model, parent_ii, position)
  }
  else # When the node does not have surrogate predicates
  {
     # Add the primary predicate
    
     predicate <- .getPrimaryPredicates(fieldLabel, op, value)
  }
  node <- append.XMLNode(node, predicate) 

  # Add score distribution for classification case.
  
  if(model$method == "class")
  {
    ylevel <- attr(model,'ylevels')
    node <- .getScoreDistributions(node, ylevel, ff, ii)
  }

  # The recursive function to create child nodes.
  
  if (length(depths) == 1)
  {
    left <- NULL
    right <- NULL
  }
  else
  {
    split.point <- which(depths[c(-1,-2)] == depths[2]) + 1 # Binary tree
    lb <- 2:split.point
    rb <- (split.point + 1):length(depths)
    left <- .genBinaryTreeNodes(depths[lb], ids[lb], counts[lb], scores[lb], fieldLabels[lb],
                               ops[lb], values[lb], model, ii, rows[lb], "left")
    right <- .genBinaryTreeNodes(depths[rb], ids[rb],counts[rb], scores[rb], fieldLabels[rb],
                                ops[rb], values[rb], model, ii, rows[rb], "right")  
  }
 
  if (!is.null(left))
  {
    node <- append.XMLNode(node, left)
    node <- append.XMLNode(node, right)
  }

  return(node)
}

#############################################################################
# FUNCTION: .getPrimaryPredicates
#
# Goal: add the primary predicate for the node

.getPrimaryPredicates <- function(field,op,value)
{
  if (op %in% c("greaterOrEqual", "lessThan"))
  {
    predicate <- xmlNode("SimplePredicate",
                         attrs=c(field=field, operator=op, value=value))
  }
  else if (op == "isIn")
  {
    predicate <- .getSimpleSetPredicate(field,op,value)
  }

  return(predicate)
}

##############################################################################
# function: .getSurrogatePredicates
#
# goal: add the surrogate predicates to take care of the missing value cases
#
# author: Zementis, Inc.
#
# date: June, 2008
##############################################################################
.getSurrogatePredicates <- function(predicate, model,i,position)
{
    ff <- model$frame
    is.leaf <- (ff$var=='<leaf>')
    index <- cumsum(c(1, ff$ncompete + ff$nsurrogate + 1*(!is.leaf)))

    # j: indices of the surrogate predicates in the splits (list) for the current node
    j <- seq(1 +index[i] + ff$ncompete[i], length.out=ff$nsurrogate[i])

    predicateNameList <- dimnames(model$splits)[[1]]
    predicateSignList <- model$splits[,2]
    predicateValueList <- model$splits[,4]

    # n: number of surrogate predicates in the current node
    n<- length(predicateNameList[j])
    currentNodePredicateNameList <- predicateNameList[j]
    currentNodeSignList <- predicateSignList[j]
    currentNodeValueList <- predicateValueList[j]
     
    # for simple set predicate
    digits=getOption("digits") 
    temp <- model$splits[,2]
    cuts <- vector(mode='character', length=nrow(model$splits))
    for (k in 1:length(cuts)) 
    {
            if (temp[k] == -1)
                cuts[k] <-paste("<", format(signif(model$splits[k,4], digits=digits)))
            else if (temp[k] ==1)
                cuts[k] <-paste("<", format(signif(model$splits[k,4], digits=digits)))
            else cuts[k]<- paste(c("L", "-", "R")[model$csplit[model$splits[k,4], 1:temp[k]]], collapse='', sep='')
    }
    currentNodeCutsList <- cuts[j]
    field <- NULL
    field$name <- as.character(attr(model$terms, "variables"))[-1]
    number.of.fields <- length(field$name)
    field$class <- attr(model$terms, "dataClasses")
    target <- field$name[1]

    for (i in 1:number.of.fields)
    {
        if (field$class[[field$name[i]]] == "factor") 
        {
           if (field$name[i] == target) 
           {
               field$levels[[field$name[i]]] <- attr(model, "ylevels")
           } else
           {
               field$levels[[field$name[i]]] <- attr(model,"xlevels")[[field$name[i]]]
           }
        }
    }

   # generate simple set predicate
   for (k in 1:n)
   {
         if(position == "left")
         {
             if(currentNodeSignList[[k]]==1)
             {
                op <- "greaterOrEqual"
             } else if (currentNodeSignList[[k]]== -1)
             {
                op <- "lessThan"
             } else {
                op <- "isIn"
             }
         } else if(position == "right")
         {
             if(currentNodeSignList[[k]]==1)
             {
                op <- "lessThan"
             } else if (currentNodeSignList[[k]]==-1)
             {
                op <- "greaterOrEqual"
             } else
             {
                op <- "isIn"
             }
         }

         if (op == "isIn" && position == "left")
         {
             # simple set predicate for a left node
             value1 <- strsplit(currentNodeCutsList[k],"")
             value <- NULL
             for(s in 1:length(value1[[1]])) 
             {
                 if(value1[[1]][s] == "L") 
                 {
                     value <-  paste(value, field$levels[[currentNodePredicateNameList[k]]][s],sep="")
                     value <- paste(value,",",sep="") 
                 }
             }
             predicate <-  append.XMLNode(predicate,.getSimpleSetPredicate(currentNodePredicateNameList[k],op,value))
         } else if( op =="isIn" && position == "right")
         {
             # simple set predicate for a right node
             value1 <- strsplit(currentNodeCutsList[k],"")
             value <- NULL
             for(s in 1:length(value1[[1]])) {
                if(value1[[1]][s] == "R") {
                  value <-  paste(value, field$levels[[currentNodePredicateNameList[k]]][s],sep="")
                  value <- paste(value,",",sep="")
                }
             }
             predicate <-  append.XMLNode(predicate,.getSimpleSetPredicate(currentNodePredicateNameList[k],op,value))
         } else
         {
             predicate<- append.XMLNode(predicate,xmlNode("SimplePredicate",attrs=c(field=currentNodePredicateNameList[k], operator=op,
                                  value=currentNodeValueList[[k]])))
         }
   }
    return(predicate)
}

########################################################################
# Function: .getSimpleSetPredicate
#
# Goal: refactor the original code, which creates the simple set predicate
# 
# Original by: Togaware
# Refactored by: Zementis, Inc.
# Refactored date: June, 2008

.getSimpleSetPredicate <- function(field, op, value)
{
  predicate <- xmlNode("SimpleSetPredicate", attrs=c(field=field,
                                               booleanOperator=op))
  # 090112 gjw We need to account for embedded commans in the
  # values. This comes up when we have R generated bins which will
  # have a list of valeus like "[402,602],(602,697],(697,763]" or
  # "(101,165],(229,292]". This is getting split incorrectly when
  # splitting on commas. I should get the values directly, using the
  # variable levels, but for now, look out for this special case and
  # deal with it. Another solution would be to change how the values
  # are represented in the binning, replacing the embedded "," with a
  # "-"

  value <- value[[1]]
  if (length(grep("^[[(][[:digit:]\\.]*", value))>0)
    value <- sub(']]', ']', paste(strsplit(value, "],")[[1]], "]", sep=""))
  else
    value <- strsplit(value, ",")[[1]]

  # 081019 gjw We need quotes around the values since there may be
  # embedded spaces. So we ensure they have quotes. In the PMML this
  # comes out as "&quot;", which when read back into R comes as a '"',
  # so perhaps that is okay? The SPSS generated PMML has actual
  # quotes. We may need to change this to ensure we get actual quotes
  # rather than the XML code for a quote.
  
  vals <- paste('"', value, '"', collapse=" ", sep="")

  predicate <- append.XMLNode(predicate, xmlNode("Array", vals,
                                                 attrs=c(n=length(value),
                                                   type="string")))

  return(predicate)
}

#######################################################################
# function: .getScoreDistributions
#
# goal: extract the probability for each category (only for classification case)
#
# author: Zementis, Inc.
#
# date: June, 2008
#######################################################################
.getScoreDistributions <- function(node, ylevel, ff, ii)
{
   # ii: the sequential oreder of the current node

   # n: number of classification categories
   n = length(ylevel)
   recordCountMap <- ff$yval2[ii,2:(1+n) ,drop=TRUE]
   confidenceMap <- ff$yval2[ii,(n+2):(2*n+1) ,drop=TRUE]

   for(i in 1:n) 
   {
     recordCount <- ifelse(is.na(recordCountMap[i]), 0, recordCountMap[i])
     confidence <- ifelse(is.na(confidenceMap[i]), 0, confidenceMap[i])
     scoreDist <- xmlNode("ScoreDistribution", attrs = c(value=ylevel[i], recordCount = recordCount, confidence = confidence))
     node <- append.XMLNode(node,scoreDist)
   }

   return(node)
}

