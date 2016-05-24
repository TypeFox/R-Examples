#----------------------------------------------------------------------------------------------
#
# cutreeDynamic
#
#----------------------------------------------------------------------------------------------
# A wrapper function for cutreeHybrid and cutreeDynamicTree.

cutreeDynamic = function(
      dendro, cutHeight = NULL, minClusterSize = 20, 
                       
      method = "hybrid", distM = NULL, 
      deepSplit = (ifelse(method=="hybrid", 1, FALSE)), 

      # Advanced options
      maxCoreScatter = NULL, minGap = NULL,
      maxAbsCoreScatter = NULL, minAbsGap = NULL,

      minSplitHeight = NULL, minAbsSplitHeight = NULL,

      # External (user-supplied) measure of branch split
      externalBranchSplitFnc = NULL, minExternalSplit = NULL,
      externalSplitOptions = list(),
      externalSplitFncNeedsDistance = NULL,
      assumeSimpleExternalSpecification = TRUE,


      # PAM stage options
      pamStage = TRUE, pamRespectsDendro = TRUE,
      useMedoids = FALSE, maxDistToLabel = NULL,
      maxPamDist = cutHeight,
      respectSmallClusters = TRUE,

      # Various options
      verbose = 2, indent = 0)
{

  #if (!is.null(labelUnlabeled))
  #{
  #  pamStage = labelUnlabeled;
  #  warning("The argument 'labelUnlabeled' is deprecated. Please use 'pamStage' instead.");
  #}

  if (!is.null(maxDistToLabel)) 
  {
    printFlush("cutreeDynamic: maxDistToLabel in deprecated. Please use maxPamDist instead");
    maxPamDist = maxDistToLabel;
  }

  if (class(dendro)!="hclust") stop("Argument dendro must have class hclust.");
  methods = c("hybrid", "tree");
  met = charmatch(method, methods);
  if ( (met==1) && (is.null(distM)) )
  {
    warning(paste('cutreeDynamic: method "hybrid" requires a valid dissimilarity matrix "distM".',
                  '\nDefaulting to method "tree".'));
    met = 2;
  }
  if (is.na(met))
  {
    stop(paste("Invalid method argument. Accepted values are (unique abbreviations of)", 
                paste(methods, collapse = ", ")));
  } else if (met==1)
  {

    # if (is.null(distM)) stop('distM must be given when using method "hybrid"');
    return(cutreeHybrid(dendro = dendro, distM = distM, 
                      cutHeight = cutHeight, 
                      minClusterSize = minClusterSize, 
                      deepSplit = deepSplit,

                      maxCoreScatter = maxCoreScatter, minGap = minGap,
                      maxAbsCoreScatter = maxAbsCoreScatter, minAbsGap = minAbsGap,

                      minSplitHeight = minSplitHeight, minAbsSplitHeight = minAbsSplitHeight,

                      # External (user-supplied) measure of branch split
                      externalBranchSplitFnc = externalBranchSplitFnc,
                      minExternalSplit = minExternalSplit,
                      externalSplitOptions = externalSplitOptions,
                      externalSplitFncNeedsDistance = externalSplitFncNeedsDistance,
                      assumeSimpleExternalSpecification = assumeSimpleExternalSpecification,

                      pamStage = pamStage, pamRespectsDendro = pamRespectsDendro,
                      useMedoids = useMedoids, 
                      maxPamDist = maxPamDist, 
                      respectSmallClusters = respectSmallClusters, 
                      verbose = verbose, indent = indent)$labels);
  } else
  {
    return(cutreeDynamicTree(dendro = dendro, maxTreeHeight = cutHeight, deepSplit = deepSplit,
                             minModuleSize = minClusterSize)); 
  }
}

#----------------------------------------------------------------------------------------------
#
# merge2Clusters
#
#----------------------------------------------------------------------------------------------
# Manually merge 2 clusters.

merge2Clusters= function(labels, mainClusterLabel, minorClusterLabel)
{
  labels2 = ifelse(as.character(labels)==minorClusterLabel, mainClusterLabel,
                        as.character(labels))
  if (class(labels)=="numeric") labels2 = as.numeric(labels2);
  if (class(labels)=="factor") labels2 = factor(labels2)
  
  labels2;
}


    
