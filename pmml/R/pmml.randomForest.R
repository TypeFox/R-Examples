# PMML: Predictive Model Markup Language
#
# Copyright (c) 2009-2015, some parts by Togaware Pty Ltd and other by Zementis, Inc. 
#
# This file is part of the PMML package for R.
#
# The PMML package is free software: you can redistribute it and/or 
# modify it under the terms of the GNU General Public License as
# published by the Free Software Foundation, either version 2 of 
#
# The PMML package is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. Please see the
# GNU General Public License for details (http://www.gnu.org/licenses/).
######################################################################################
#
# Author: Tridivesh Jena
#
# Implemented: by Tridivesh Jena (info@zementis.com) to add the
# capability to export random forest models.

pmml.randomForest <- function(model,
                              model.name="randomForest_Model",
                              app.name="Rattle/PMML",
                              description="Random Forest Tree Model",
                              copyright=NULL,
			      transforms=NULL,
			      unknownValue=NULL,
                              ...)

{
   if (! inherits(model, "randomForest"))
    stop("Not a legitimate randomForest object")

   requireNamespace("randomForest",quietly=TRUE)

   a <- randomForest::getTree(model,2)

  # Tridivesh: Collect the required information. We list all variables,
  # irrespective of whether they appear in the final model. This seems 
  # to be the standard thing to do with PMML. It also adds extra information
  # - i.e., the model did not need these extra variables!
  #
  # For a randomForest formula as currently used in Rattle, the
  # target is, for example, as.factor(Adjusted). Here, I need to
  # remove the as.factor(...). I wonder if I need to identify a
  # transformation in the PMML.

  field <- NULL
  tr.vars <- attr(model$terms, "dataClasses")
  var.names0 <- gsub("as\\.factor\\(","",names(tr.vars))
  var.names <- gsub("\\)","",var.names0)
  field$name <- var.names
  number.of.fields <- length(field$name)
  target <- var.names[1]

  # The following is a bit sus and does not really get the corect type
  # of the as.factor modified fields!
  # Tridi 2/8/12: modified to get category names correctly
  field$class <- attr(model$terms, "dataClasses")
  names(field$class) <- var.names

#Better implementation of the field constructor below:

  cat <- list() 
  for (i in 1:number.of.fields)
  {
    if (field$class[[field$name[i]]] == "factor")
    {
      if (field$name[i] == target)
      {
        field$levels[[field$name[i]]] <- model$classes
      }
      else
      {
          cat <- c(cat,model$forest$xlevels[field$name[i]])
      }
    }
  }
#  field <- c(field,list("levels"=cat))

  # PMML

  pmml <- .pmmlRootNode("4.2")

  # PMML -> Header

  pmml <- append.XMLNode(pmml, .pmmlHeader(description, copyright, app.name))

  # PMML -> DataDictionary

  pmml <- append.XMLNode(pmml, .pmmlDataDictionary(field,transformed=transforms))

  mmodel <- xmlNode("MiningModel",attrs=c(modelName=model.name,algorithmName="randomForest",functionName=model$type))
  mmodel <- append.XMLNode(mmodel,.pmmlMiningSchema(field,target,transforms,unknownValue))

  # Tridi: Add output fields
  mmodel <- append.XMLNode(mmodel, .pmmlOutput(field, target))

  #Tridi: If interaction terms do exist, define a product in LocalTransformations and use
  # it as a model variable. This step is rare as randomForest seems to avoid multiplicative
  # terms.
  ltNode <- xmlNode("LocalTransformations")
  interact <- FALSE
  for(fld in 1:number.of.fields)
  {
    if(length(grep(":",field$name[fld])) == 1)
    {
     interact <- TRUE
     drvnode <- xmlNode("DerivedField",attrs=c(name=field$name[fld],optype="continuous",
                                                               dataType="double"))
     applyNode <- xmlNode("Apply",attrs=c("function"="*"))
     for(fac in 1:length(strsplit(field$name[fld],":")[[1]]))
     {
       fldNode <- xmlNode("FieldRef",attrs=c(field=strsplit(field$name[fld],":")[[1]][fac]))
       if(length(grep("as\\.factor\\(",fldNode)) == 1)
         fldNode <- gsub("as.factor\\((\\w*)\\)","\\1", fldNode, perl=TRUE)
       applyNode <- append.XMLNode(applyNode, fldNode)
     }
     drvnode <- append.XMLNode(drvnode, applyNode)
    }
    if(interact)
      ltNode <- append.XMLNode(ltNode, drvnode)
  }
  if(interact && is.null(transforms))
  {
    mmodel <- append.XMLNode(mmodel, ltNode)
  }

  # test of Zementis xform functions
  if(interact && !is.null(transforms))
  {
    ltNode <- .pmmlLocalTransformations(field, transforms, ltNode)
    mmodel <- append.XMLNode(mmodel, ltNode)
  }
  if(!interact && !is.null(transforms))
  {
    mmodel <- append.XMLNode(mmodel,.pmmlLocalTransformations(field, transforms, ltNode))
  }

  if(model$type == "regression") 
  {
    segmentation <- xmlNode("Segmentation",attrs=c(multipleModelMethod="average"))
  } 
  if(model$type == "classification") 
  {
    segmentation <- xmlNode("Segmentation",attrs=c(multipleModelMethod="majorityVote"))
  }

  numTrees <-model$ntree
  
  segments <- lapply(1:numTrees,function(x){.makeSegment(x,model,model.name,field,target,unknownValue)})
  segmentation2 <- append.XMLNode(segmentation, segments)
  rm(segmentation)
  rm(segments)

  mmodel2 <- append.XMLNode(mmodel,segmentation2)
  rm(mmodel)
  rm(segmentation2)

  pmml2 <- append.XMLNode(pmml, mmodel2)
  rm(pmml)

  return(pmml2)
}

   .makeSegment <- function(b,model,model.name,field,target,unknownValue=NULL)
   {
    print(paste("Now converting tree ",b," to PMML"))
  # PMML -> TreeModel -> Node
  # Tridi: Tree structure information here as produced by the getTree function of the 
  # randomForest package
    if(model$type == "regression") 
    {
      tree <- cbind(model$forest$leftDaughter[,b],
                    model$forest$rightDaughter[,b],
                    model$forest$bestvar[,b],
                    model$forest$xbestsplit[,b],
                    model$forest$nodepred[,b])[1:model$forest$ndbigtree[b],]
    } else 
    {
      tree <- cbind(model$forest$treemap[,,b],
                    model$forest$bestvar[,b],
                    model$forest$xbestsplit[,b],
                    model$forest$nodepred[,b])[1:model$forest$ndbigtree[b],]
    }
    nodeList <- list()
    rowId <- which(tree[,2] == max(tree[,2]))
    nodeList <- lapply(1:tree[rowId,2],function(x){.makeNode(x,model,tree,field)})
    nodeF <- .makeTree(nodeList,tree,rowId)

    rm(nodeList)

  # PMML -> TreeModel
    if(model$type == "regression") 
    {
          tree.model <- xmlNode("TreeModel",
                        attrs=c(modelName=model.name,
                          functionName="regression",
                          algorithmName="randomForest",
                          splitCharacteristic="binarySplit"))
    }
    if(model$type == "classification") 
    {
  	tree.model <- xmlNode("TreeModel",
                        attrs=c(modelName=model.name,
                          functionName="classification",
                          algorithmName="randomForest",
                          splitCharacteristic="binarySplit"))
    }

  # PMML -> TreeModel -> MiningSchema
  tree.model <- append.XMLNode(tree.model, .pmmlMiningSchema(field, target,unknownValue=unknownValue))


  # Add to the top level structure.
     segment <- xmlNode("Segment",attrs=c(id=b))
     tru <- xmlNode("True")
     segment <- append.XMLNode(segment, tru)

  # Add to the top level structure.
     tree.model <- append.XMLNode(tree.model, nodeF)
     segment <- append.XMLNode(segment, tree.model)
     return(segment)
  }

.makeTree <- function(nodeLi,tre,rId)
{
    while(rId!=0)
    {
    if(tre[rId,1] != 0)
    {
      nodeR<-nodeLi[[tre[rId,2]]]
      nodeL<-nodeLi[[tre[rId,1]]]
      nodeT<-nodeLi[[rId]]

      nodeT<-append.XMLNode(nodeT,nodeL)

      nodeT<-append.XMLNode(nodeT,nodeR)

      nodeLi[[rId]]<-nodeT

      rm(nodeR)
      rm(nodeL)
      rm(nodeT)
    }
    rId=rId-1
  }

  return(nodeLi[[1]])
}


.makeNode <- function(n, mod, tinf, fieldInfo)
{

  if(n==1)
  {
    return(append.XMLNode(xmlNode("Node",attrs=c(id=1)),xmlNode("True")))
  } 
  else 
  {
    side <- 2
    if(n/2 == floor(n/2))
    {
      side <- 1
    }
    score  <- NULL
    if(tinf[n,1] == 0)
    {
      if(mod$type == "regression"){
        score <- tinf[n,5]
      } else{
	score <- mod$classes[tinf[n,5]]
      }
    }
 
    if(is.null(score))
    {
      rfNode <- xmlNode("Node",attrs=c(id=n))
    } else
    {
      rfNode <- xmlNode("Node",attrs=c(id=n,score=score))
    }
 # After the node, add the split info in pmml
     # ------------------------------------------------------------------------- 
     rowid <- which(tinf[,side]==n)
     fname <- names(mod$forest$xlevels[tinf[rowid,3]])
     # is the field categorical
     logical <- FALSE
     numeric <- FALSE
     fieldClass <- fieldInfo$class[fname]
     if(fieldClass == "numeric")
	numeric <- TRUE
     if(fieldClass == "logical")
	logical <- TRUE

  # split if var is numeric
     if(numeric) 
     {
       if(side == 1)
       { 
         splitNode <- xmlNode("SimplePredicate",attrs=c(field=fname,operator="lessOrEqual",
                        value=format(tinf[rowid,4],digits=18)))
       } else
       {
         splitNode <- xmlNode("SimplePredicate",attrs=c(field=fname,operator="greaterThan",
                        value=format(tinf[rowid,4],digits=18)))
       }
     } else if(logical)
     {
       bool = ifelse(tinf[rowid,4] <= 0.5, FALSE, TRUE)
       splitNode <- xmlNode("SimplePredicate",attrs=c(field=fname,operator="equal",
                        value=bool))
     } else  
     {
       if(tinf[rowid,4] >= 0)
       {
  # split if var is categorical
         binary <- .sdecimal2binary(tinf[rowid,4])
         ssp <- xmlNode("SimpleSetPredicate",attrs=c(field=fname,booleanOperator="isIn"))
         num1 <- 0
         scat <- NULL
         holder <- array(0,dim=c(1,mod$forest$ncat[fname][[1]]))
         for(k in 1:length(binary))
         {
          holder[k] = binary[k]
         }

 # for each category allowed, if value is 1 (from the binary conversion) then go left
         options(useFancyQuotes = FALSE)
         for(k in 1:mod$forest$ncat[fname][[1]]) 
         {
           if(side == 1)
           {
             if(holder[k]==1)
             {
              num1 <- num1 + 1
              catname <- mod$forest$xlevels[fname][[1]][k]
              scat <- paste(scat," ",dQuote(catname))
             }
           } else
           {
             if(holder[k]==0)
             {
              num1 <- num1 + 1
              catname <- mod$forest$xlevels[fname][[1]][k]
              scat <- paste(scat," ",dQuote(catname))
             }
           }
         }
             
 # all the gsubs are to strip intermediate, leading and trailing spaces. 
         scat <- gsub("^[ ]*","",scat)
         ap <- xmlNode("Array",attrs=c(n=num1,type="string"),scat)
         ssp <- append.XMLNode(ssp,ap)
         splitNode <- ssp
       } 
     }
     rfNode <- append.XMLNode(rfNode,splitNode)
   }
  return(rfNode)
}



.getRFTreeNodes2 <- function(recursiveObject, model, side, tinf, rowfrom, rownext, fieldInfo)
{
  if(!((model$type == "regression") || (model$type == "classification")))
     print("Model type not supported")
  treeSkip <- FALSE

 # Keep going over nodes; if leaf node, add score, else split and keep going
  if((rowfrom == 1) && (rownext == 1)) 
  {
#handle trees with 1 node only
    if(is.null(dim(tinf)))
    {
      rfNode <- xmlNode("Node",attrs=c(id="1",score=tinf[5]))
      nodeB <- xmlNode("True")
      rfNode <- append.XMLNode(rfNode,nodeB)
      recursiveObject$internalNode <- rfNode
      return(recursiveObject) 
    } else
    {
    # Add top node at first loop only
      rfNode <- xmlNode("Node",attrs=c(id="1"))
      nodeB <- xmlNode("True")
      rfNode <- append.XMLNode(rfNode,nodeB)
    }
  } else 
  {
      fname <- attributes(model$forest$xlevels[tinf[rowfrom,3]])[[1]]
 # Treat left and right leafs separately as their information is stored in separate column in tree
   if(side==-1)
   {
     if(tinf[rownext,1] == 0) 
     {
 # The score for classification must be translated from a number to the category name
       if(model$type == "regression") 
       {
 # The score for regresion can just be read off.
          rfNode <- xmlNode("Node",attrs=c(id=tinf[rowfrom,1],score=tinf[rownext,5]))
       } else 
       {
          rfNode <- xmlNode("Node",attrs=c(id=tinf[rowfrom,1],score=model$classes[tinf[rownext,5]]))
       }
     } else
     {
       rfNode <- xmlNode("Node",attrs=c(id=tinf[rowfrom,1]))
     } 
 # After the node, add the split info in pmml
     # ------------------------------------------------------------------------- 
  # left side, regression model, terminal node 
     # is the field categorical
     logical <- FALSE
     numeric <- FALSE
     if(is.numeric(model$forest$xlevels[[tinf[rowfrom,3]]][1])) 
     {
       name = names(model$forest$xlevels[tinf[rowfrom,3]]) 
       if(fieldInfo$class[name] == "logical")
       {
	 logical <- TRUE 
       } else
       {
         numeric <- TRUE
       }
     } else 
     {
       numeric <- FALSE
     }
  # split if var is numeric
     if(numeric) 
     { 
       splitNode <- xmlNode("SimplePredicate",attrs=c(field=fname,operator="lessOrEqual",
			value=tinf[rowfrom,4]))
     } else if(logical)
     {
       bool = ifelse(tinf[rowfrom,4] <= 0.5, FALSE, TRUE)
       splitNode <- xmlNode("SimplePredicate",attrs=c(field=fname,operator="equal",
                        value=bool))

     } else  
     {
       if(tinf[rowfrom,4] >= 0)
       {
  # split if var is categorical
         binary <- .sdecimal2binary(tinf[rowfrom,4])
         ssp <- xmlNode("SimpleSetPredicate",attrs=c(field=fname,booleanOperator="isIn"))
         num1 <- 0
         scat <- NULL
         holder <- array(0,dim=c(1,model$forest$ncat[fname][[1]]))
         for(k in 1:length(binary))
         {
          holder[k] = binary[k]
         }

 # for each category allowed, if value is 1 (from the binary conversion) then go left
         options(useFancyQuotes = FALSE)
         for(k in 1:model$forest$ncat[fname][[1]]) 
         {
           if(holder[k]==1)
           {
            num1 <- num1 + 1
            catname <- model$forest$xlevels[fname][[1]][k]
            scat <- paste(scat," ",dQuote(catname))
           }
         }
             
 # all the gsubs are to strip intermediate, leading and trailing spaces. 
         scat <- gsub("^[ ]*","",scat)
         ap <- xmlNode("Array",attrs=c(n=num1,type="string"),scat)
         ssp <- append.XMLNode(ssp,ap)
         splitNode <- ssp
            
       } else
       {
         treeSkip <- TRUE
         sknode <- xmlNode("skip")
         recursiveObject$internalNode <- "skip" 
         return(recursiveObject)
       }
     }
        if(treeSkip == FALSE) 
        rfNode <- append.XMLNode(rfNode,splitNode)
   } else  
   {
       # ----------------------------------------------------------------------------
  # right side, regression, terminal node
 # repeat all over for right side 
        if(tinf[rownext,1] == 0) 
        {
          if(model$type == "regression") 
          {
 # The only difference is where to read off the node info from the tree structure 
            rfNode <- xmlNode("Node",attrs=c(id=tinf[rowfrom,2],score=tinf[rownext,5]))
          } else 
          {
            rfNode <- xmlNode("Node",attrs=c(id=tinf[rowfrom,2],score=model$classes[tinf[rownext,5]]))
         }
        }
        else
        {
          rfNode <- xmlNode("Node",attrs=c(id=tinf[rowfrom,2]))
        }

       # is the field categorical
       logical <- FALSE
       numeric <- FALSE 
       if(is.numeric(model$forest$xlevels[[tinf[rowfrom,3]]][1]))      
       {
         name = names(model$forest$xlevels[tinf[rowfrom,3]])
         if(fieldInfo$class[name] == "logical")
         {
           logical <- TRUE 
         } else
         {
           numeric <- TRUE
         }
       } else
       {
         numeric <- FALSE
       }

       if(numeric)
       {
         splitNode <- xmlNode("SimplePredicate",attrs=c(field=fname,operator="greaterThan",
                        value=tinf[rowfrom,4]))
       } else if(logical)
       {
         bool = ifelse(tinf[rowfrom,4] <= 0.5, TRUE, FALSE)
         splitNode <- xmlNode("SimplePredicate",attrs=c(field=fname,operator="equal",
                        value=bool))

       } else
       {
         if(tinf[rowfrom,4] >= 0)
         {  
  # split if var is categorical
           binary <- .sdecimal2binary(tinf[rowfrom,4])
           ssp <- xmlNode("SimpleSetPredicate",attrs=c(field=fname,booleanOperator="isIn"))
           num1 <- 0
           scat <- NULL
           holder <- array(0,dim=c(1,model$forest$ncat[fname][[1]]))
           options(useFancyQuotes = FALSE)
           for(k in 1:length(binary))
           {
            holder[k] = binary[k]
           }
           for(k in 1:model$forest$ncat[fname][[1]]) 
           {
             if(holder[k]==0)
             {
               num1 <-  num1 + 1
               catname <- as.character(unlist(model$forest$xlevels[fname]))[k]
               scat <- paste(scat," ",dQuote(catname))
             }
           }

           scat <- gsub("^[ ]*","",scat)
           ap <- xmlNode("Array",attrs=c(n=num1,type="string"),scat)
           ssp <- append.XMLNode(ssp,ap)
           splitNode <- ssp
         } else
         {
           treeSkip <- TRUE
           sknode <- xmlNode("skip")
           recursiveObject$internalNode <- "skip" 
           return(recursiveObject)
         }
       }
       if(treeSkip == FALSE)
       rfNode <- append.XMLNode(rfNode,splitNode)
     } 
    } 

  if(tinf[rownext,5] == -1)
  {
    terminalFlag <- TRUE
  } else {
    terminalFlag <- FALSE
  }

  if (terminalFlag == TRUE) 
  {
#    only the predicted value for this node is the output
  }
  if(terminalFlag == FALSE) 
  {

    recursiveObject$internalNode <- NULL
    recursiveObject <- .getRFTreeNodes2(recursiveObject,model,-1,tinf,rownext,tinf[rownext,1],fieldInfo)
    if(!is.null(recursiveObject$internalNode) && (recursiveObject$internalNode == "skip")[[1]])
    {
      return(recursiveObject)
    } else
    {
    rfNode <- append.XMLNode(rfNode, recursiveObject$internalNode)
    }



    recursiveObject$internalNode <- NULL
    recursiveObject <- .getRFTreeNodes2(recursiveObject,model,1,tinf,rownext,tinf[rownext,2],fieldInfo)
    if(!is.null(recursiveObject$internalNode) && (recursiveObject$internalNode == "skip")[[1]])
    {
      return(recursiveObject)
    } else
    {
    rfNode <- append.XMLNode(rfNode, recursiveObject$internalNode)
  }

  }

  if(!is.null(recursiveObject$internalNode) && (recursiveObject$internalNode == "skip")[[1]])
  {
    return(recursiveObject)
  } else
  {
  recursiveObject$internalNode <- rfNode
  return(recursiveObject)
  }
}
