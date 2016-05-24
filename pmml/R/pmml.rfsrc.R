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
#####################################################################################
#
# Author: Tridivesh Jena

pmml.rfsrc <- function(model,
                     model.name="rsf_Model",
                     app.name="Rattle/PMML",
                     description="Random Survival Forest Model",
                     copyright=NULL,
		     transforms=NULL,
		     unknownValue=NULL, ...)
{
  # Based on RANDOM SURVIVAL FOREST 2.0.0, Copyright 2006, Cleveland Clinic
  # Original by Hemant Ishwaran and Udaya B. Kogalur
  # Unified with the pmml package by Graham Williams
  
  # Tridi 1/15/12
  # although seems logical, the constructed rsf model in R doesnt seem to require this.
  # remove this rather than force modeller to cast rsf object into a rsf,forest object  
  #  if (sum(inherits(model, c("rsf", "forest"), TRUE) == c(1, 2)) != 2)
  #    stop("Not a legitimate (rsf, forest) object")

  requireNamespace("randomForestSRC", quietly=TRUE)

  # Collect the required information.

  field <- NULL

  field$name <- model$xvar.names
  if (is.null(field$name))
    stop("RSF predictorNames is NULL.  Please ensure the object is valid.")
  number.of.fields <- length(field$name)

  field$class <- rep("numeric", number.of.fields) # All fields are numeric? 
  names(field$class) <- field$name

  nativeArray <- model$forest$nativeArray
  if (is.null(nativeArray))
    stop("RSF nativeArray content is NULL. Please ensure object is valid.")
  if(length(nativeArray[nativeArray[,"mwcpSZ"]!=0,"mwcpSZ"]) != 0)
   stop("Categorical predictor variables not yet supported for randomForestSRC models.")
  
  numTrees <- length(as.vector(unique(nativeArray$treeID))) # Trees in forest
  
  timeInterest = model$time.interest
  if (is.null(timeInterest))
    stop("RSF timeInterest content is NULL. Please ensure object is valid.")

  formula = model$call
  if (is.null(formula))
    stop("RSF formula is NULL.  Please ensure the object is valid.")

  forestSeed = model$forest$seed
  if (is.null(forestSeed))
    stop("RSF forestSeed content is NULL.  Please ensure object is valid.")

  # PMML

  pmml <- .pmmlRootNode("4.2")

  # PMML -> Header

  pmml <- append.XMLNode(pmml, .pmmlHeader(description, copyright, app.name))

#  # PMML -> MiningBuildTask
#
#  buildNode <- xmlNode("MiningBuildTask")
#
#  # PMML -> MiningBuildTask -> Extension
# 
#  extensionNode <- xmlNode("Extension")
#
#  # PMML -> MiningBuildTask -> Extension -> TimesOfInterest
#  timeVar <- colnames(model$yvar)[1]
#  tnode <- append.XMLNode(extensionNode, xmlNode("TimeVariable", attrs=c(name=timeVar)))
#
##  extensionNode <- append.XMLNode(extensionNode, 
##                             xmlNode("Array", attrs=c(n=length(timeInterest), 
##                                     type="double"),paste(timeInterest, collapse="  \n  ")))
#  extensionNode <- append.XMLNode(extensionNode,
#                             xmlNode("Array", attrs=c(n=length(timeInterest),               
#                                     type="double"),paste(timeInterest,collapse=" ")))
# 
## Add into the PMML.
#
#  pmml <- append.XMLNode(pmml, append.XMLNode(buildNode, extensionNode))
  
  # PMML -> DataDictionary

  pmml <- append.XMLNode(pmml, .pmmlDataDictionarySurv(field, timeName=model$yvar.names[1], statusName=model$yvar.names[2], dataset=NULL,weights=NULL, transformed=transforms))
  
  # Create a dummy XML node object to insert into the recursive
  # output object.

  internalNode <- xmlNode("Null")
  
  # Define the variables for the offset and leaf count in the
  # recursive output object.

  offset <- leafCount <- 1
  
  # Create the recursive output object.  This would be unnecessary if
  # it was possible to declare global variables in a package.

  recursiveOutput <- list(internalNode = internalNode,
                          offset = offset, leafCount = leafCount)

  # <MiningModel>
  miningModelNode <- xmlNode("MiningModel", attrs=c(modelName="RrsfModel",functionName="regression")) 

# if (.supportTransformExport(transforms))
# {
#    field <- .unifyTransforms(field, transforms)
#    transforms <- .activateDependTransforms(transforms)
#  }

  # MiningModel -> MiningSchema
  miningModelNode <- append.XMLNode(miningModelNode, .pmmlMiningSchemaSurv(field, timeName=model$yvar.names[1], statusName=model$yvar.names[2], target=NULL, inactive=NULL, transformed=transforms,unknownValue=unknownValue))

  #Tridi: If interaction terms do exist, define a product in LocalTransformations and use
  # it as a model variable. This step is rare as randomForest seems to avoid multiplicative
  # terms.
  ltNode <- xmlNode("LocalTransformations")
  interact <- FALSE
  for(fld in 1:number.of.fields){
    if(length(grep(":",field$name[fld])) == 1){
     interact <- TRUE
     drvnode <- xmlNode("DerivedField",attrs=c(name=field$name[fld],optype="continuous",dataType="double"))
     applyNode <- xmlNode("Apply",attrs=c("function"="*"))
     for(fac in 1:length(strsplit(field$name[fld],":")[[1]])){
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
    mmodel <- append.XMLNode(mmodel, ltNode)

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

  # ensemble method
  segmentationNode <- xmlNode("Segmentation",
			attrs=c(multipleModelMethod="average"))

  # Survival times 
  extensionNode <- xmlNode("Extension")

  # PMML -> MiningModel -> Segmentation -> Extension -> TimesOfInterest
  timeVar <- colnames(model$yvar)[1]
  tnode <- append.XMLNode(extensionNode, xmlNode("TimeVariable", attrs=c(name=timeVar)))

  extensionNode <- append.XMLNode(extensionNode,
                             xmlNode("Array", attrs=c(n=length(timeInterest),
                                     type="double"),paste(timeInterest,collapse=" ")))

  # Add into the PMML.
  segmentationNode <- append.XMLNode(segmentationNode, extensionNode)
  
  # Loop through all trees in the forest and extract the data.

  for (b in 1:numTrees)
  {
    print(paste("Converting Tree",b," to PMML",sep=""))
    segmentNode <- xmlNode("Segment",attrs=c(id=b))
    predicateNode <- xmlNode("True")
    segmentNode <- append.XMLNode(segmentNode, predicateNode)

    treeName <- paste("Tree",b,sep="")
    treeModelNode <- xmlNode("TreeModel",
                             attrs=c(modelName=treeName, functionName="regression",
                               algorithmName="rsf",
                               splitCharacteristic="binarySplit"))

    # PMML --> TreeModel [b] -> MiningSchema
    
    treeModelNode <- append.XMLNode(treeModelNode, .pmmlMiningSchemaSurv(field, timeName=model$yvar.names[1], statusName=model$yvar.names[2], NULL, unknownValue=unknownValue))
    
    
    # Initialize the root node.  This differs from the rest of the
    # internal nodes in the PMML structure.

    treeRoot <- xmlNode("Node", attrs=c(score=0,id=0))
    treeRoot <- append.XMLNode(treeRoot, xmlNode("True"))
    
    rootParmID <- nativeArray$parmID[recursiveOutput$offset] 
    rootSpltPT <- nativeArray$contPT[recursiveOutput$offset] 
    recursiveOutput$offset <- recursiveOutput$offset + 1
    recursiveOutput$leafCount <- 1
    
    # Check that the current tree is not a stump (root node only with
    # no branches)

    if (rootParmID != 0)
    {
      # find scored values
      member <- model$membership[,b]
      inbag <- model$inbag[,b]
      y0 <- model$yvar
      numNodes <- max(member)
      numData <- nrow(y0)
      #initialize array to hold chf data
      numtimes <- length(model$time.interest)
      nodescores <- vector("list",numNodes)

      for(node in 1:numNodes)
      {
        scorematrix <- matrix(0,numtimes,2)
        scorematrix[,1] <- model$time.interest

        bag1 <- inbag[member==node]
        y1 <- y0[member==node,]

        #check if number of inbag info equals number of data points
        sz <- length(bag1)
        if(sz != nrow(y1))
        {
          stop("Inconsistent data; inbag values information not complete. Please contact Technical Support.")
        }

        nelsonAalen <- 0

        # Yinint = total population in the node
        Yinit <- sum(bag1)

        # get inbag only data in time order
        bag2 <- bag1[bag1 != 0]
        y2 <- y1[bag1 != 0,]
        ybag1 <- cbind(y2,bag2)
        y <- ybag1[order(ybag1[,1]),]
        colnames(y) <- c("time","status","numInbag")
     
        for(bdat in 1:nrow(y))
        {
          nelsonAalen = nelsonAalen + (y[bdat,"status"]*y[bdat,"numInbag"]/Yinit)
          tindex <- which(scorematrix[,1] == y[bdat,"time"])
          scorematrix[tindex,2] <- nelsonAalen
          if(y[bdat,"status"] == 1)
          {
            Yinit = Yinit - y[bdat,"numInbag"]
          }
        }

        #fill up the rest of the chf matrix
        for(sindex in 2:numtimes)
        {
          if(scorematrix[sindex,2] == 0)
          {
            scorematrix[sindex,2] <- scorematrix[sindex-1,2]
          }
        }

        nodescores[[node]] <- scorematrix
      }

      # The tree must be created in two phases.  First, the root left
      # daughter branches are created.  Second, the root right
      # daughter branches are created.  This is due to the root node
      # having a slightly different structure using the PMML protocol.
      # The root node actually has no split information.  The split
      # information is encoded into the daughter nodes.  Thus, instead
      # of making a check for the root node in the recursive routine,
      # we call the recursive routine twice.
      
      # Create the left daughter nodes.  Note that the object node
      # content is irrelevant as input.

      recursiveOutput$internalNode <- NULL
      recursiveOutput <- .rsfMakeTree(recursiveOutput, nativeArray,
                                     field$name, b, -1, rootParmID,
                                     rootSpltPT,model,nodescores)
      
      treeRoot <- append.XMLNode(treeRoot, recursiveOutput$internalNode)
     
      recursiveOutput$leafCount <- recursiveOutput$leafCount + 1
      
      # Create the right daughter nodes.  Note that the object node
      # content is irrelevant as input.

      recursiveOutput$internalNode <- NULL
      recursiveOutput <- .rsfMakeTree(recursiveOutput, nativeArray,
                                     field$name, b, +1, rootParmID,
                                     rootSpltPT,model,nodescores)
      
      treeRoot <- append.XMLNode(treeRoot, recursiveOutput$internalNode)
     
    }
    
    # Add the current tree to the PMML data structure.

    treeModelNode <- append.XMLNode(treeModelNode, treeRoot)
    segmentNode <- append.XMLNode(segmentNode, treeModelNode)
    segmentationNode <- append.XMLNode(segmentationNode, segmentNode)
  }
  miningModelNode <- append.XMLNode(miningModelNode, segmentationNode)
  
  pmml <- append.XMLNode(pmml, miningModelNode)  
  return (pmml)
}

.rsfMakeTree <- function(recursiveObject, nativeArray, predictorNames, b,
                      daughter, splitParameter, splitValue, model, nodescores)
{
  # Node information encoded in a PMML TreeModel follows a slightly
  # different protocol than that encoded in our RSF matrix
  # representation.  Since the RSF representation is linear in
  # nature, each record containing node information must encode the
  # split information, particularly the split parameter and split
  # point, in the record itself.  In contrast, the PMML TreeModel
  # indicates a split by the presence of daughters in the node.  The
  # split parameter and split point are encoded by a SimplePredicate
  # tag in the daughters.  In creating a PMML tree from an RSF tree,
  # the recursive algorithm requires a "look back" to the previous
  # record in the RSF tree to determine the split parameter and
  # value.  This is accomplished via the parameters passed by the
  # parent call to this routine.
  
  # Weak consistency check to ensure that the iteration matches the
  # treeID in the nativeArray record.

  if(b != nativeArray$treeID[recursiveObject$offset])
    stop("Invalid nativeArray input record (treeID) at ",
         recursiveObject$offset, ".  Please contact Technical Support.")

  # Read the current nativeArray record, and determine whether this is
  # a terminal node.

  fwdSplitParameter <- nativeArray$parmID[recursiveObject$offset]
  fwdSplitValue <- nativeArray$contPT[recursiveObject$offset]

  # Create the node that will be returned on this call.
  ident <- nativeArray$nodeID[recursiveObject$offset]
  if (fwdSplitParameter == 0)
  {
    scorelist <- nodescores[[ident]]
    extensionNode <- xmlNode("Extension")
    extensionNode <- append.XMLNode(extensionNode,
                             xmlNode("Array", attrs=c(n=length(scorelist[,2]),
                                     type="double"),paste(scorelist[,2],collapse=" ")))
    rsfNode <- xmlNode("Node",attrs=c(score="HazardFunction",id=ident))
    rsfNode <- append.xmlNode(rsfNode,extensionNode)
    terminalFlag <- TRUE
  }
  else if (fwdSplitParameter != 0)
  {
    rsfNode <- xmlNode("Node")
 
    terminalFlag <- FALSE
  }
  
  # Determine whether this the left of right daughter.

  if (daughter == -1)
    parseString <- "lessOrEqual"
  else if (daughter == +1)
    parseString <- "greaterThan"
  else
    # Ensure that the function call is coherent to aid in debugging.
    stop("Invalid parse direction encountered during recursion.",
         "Please contact Technical Support.")

  # Add the split information to this node via the look back.
  pName <- predictorNames[splitParameter]
  if(length(grep("as\\.factor\\(",predictorNames[splitParameter])) == 1)
    pName <- gsub("as.factor\\((\\w*)\\)","\\1", predictorNames[splitParameter], perl=TRUE)

  rsfNode <- append.XMLNode(rsfNode,
                            xmlNode("SimplePredicate",
                                  attrs=c(field=pName,
                                    operator=parseString, value=splitValue)))

  # Increment the offset, always.

  recursiveObject$offset <- recursiveObject$offset + 1

  # Parse left and then right, if this is not a terminal node.

  if (terminalFlag == FALSE)
  {
    # Parse left: Do not increment the leafCount.  Internally
    # increment the offset, always.  Note that the object node content
    # is irrelevant as input.

    recursiveObject$internalNode <- NULL
    recursiveObject <- .rsfMakeTree(recursiveObject, nativeArray,
                                   predictorNames, b, daughter = -1,
                                   fwdSplitParameter, fwdSplitValue, model,nodescores)
    
    rsfNode <- append.XMLNode(rsfNode, recursiveObject$internalNode)
    
    # Parse right: Increment the leafCount.  Internally increment the
    # offset, always.  Note that the object node content is irrelevant
    # as input.
    
    recursiveObject$leafCount <- recursiveObject$leafCount + 1
    recursiveObject$internalNode <- NULL
    recursiveObject <- .rsfMakeTree(recursiveObject, nativeArray,
                                   predictorNames, b, daughter = +1,
                                   fwdSplitParameter, fwdSplitValue, model,nodescores)
    
    rsfNode <- append.XMLNode(rsfNode, recursiveObject$internalNode)
    
  }
  
  # Modify the recursive object with the new internal node structure.

  recursiveObject$internalNode <- rsfNode
  
  return (recursiveObject)
  
}

