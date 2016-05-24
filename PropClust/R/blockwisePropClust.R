# Block-wise propensity clustering.

# Basic input should be an adjacency matrix, a logical indicating whether L2 or Poisson updates should be
# used, maximum block size, and most likely some pre-clustering arguments.

# The function may run something like hierarchical clustering coupled with dynamic tree cut to determine
# initial clusters. Then determine pairwise cluster dissimilarities and cluster the clusters. Then merge the
# clusters into blocks of size not exceeding the given maximum size.

# The advantage of this is relative simplicity and the fact that I get the blocks and clustering within each
# block in one step. The disadvantage is that the pre-clustering may not work all that well.

# The function requires WGCNA but the requirement may be removed.


#============================================================================================
#
# Helper functions
#
#============================================================================================

.spaste = function(...) { paste(..., sep="") }

.checkAdjMat = function (adjMat, min = 0, max = 1) 
{
    dim = dim(adjMat)
    if (is.null(dim) || length(dim) != 2) 
        stop("adjacency is not two-dimensional")
    if (!is.numeric(adjMat)) 
        stop("adjacency is not numeric")
    if (dim[1] != dim[2]) 
        stop("adjacency is not square")
    if (max(abs(adjMat - t(adjMat)), na.rm = TRUE) > 1e-12) 
        stop("adjacency is not symmetric")
    if (min(adjMat, na.rm = TRUE) < min || max(adjMat, na.rm = TRUE) > max) 
        stop("some entries are not between", min, "and", max)
}



.translateUsingTable = function(x, translationTable)
{
  translationTable[ match(x, translationTable[, 1]), 2]
}

.translationTable = function(from, to)
{
  if (length(from)!= length(to)) stop("Length of 'from' and 'to' must be the same.");

  tab = as.matrix(table(from, to));
  if (ncol(tab)!=nrow(tab)) {
     printFlush("Error in .translationTable: table is not 1-to-1.")
     if (ncol(tab)<10 && nrow(tab) < 10) print(tab);
     stop();
  }

  nonZeros.row = rowSums( tab!=0 )  
  nonZeros.col = colSums( tab!=0 )

  if (any(c(nonZeros.row, nonZeros.col) != 1)) {  
     printFlush("Error in .translationTable: table is not 1-to-1.")
     if (ncol(tab)<10 && nrow(tab) < 10) print(tab);
     stop();
  }

  unique.from = sort(unique(from));
  tt = cbind(unique.from, to [ match(unique.from, from)] );
  colnames(tt) = c("from", "to");
  tt;
}

#=====================================================================================================
#
# .minWhichMin: min() and which.min() of columns in a matrix
#
#=====================================================================================================

# This is a wrapper around my C-level function. 
.minWhichMin = function(x)
{
  x = as.matrix(x);
  nc= ncol(x);
  nr = nrow(x);
  min = rep(0, nc);
  which = rep(0, nc);
  whichmin = .C("minWhichMin", as.double(x),
                as.integer(nr), as.integer(nc),
                as.double(min), as.double(which), DUP = FALSE, NAOK = TRUE);
  cbind( min = whichmin[[4]], which = whichmin[[5]] + 1);
}


#===================================================================================================
#
# .clusterDissimilarity
#
#===================================================================================================

.clusterDissimilarity = function(dissimilarity, labels, method = c("average", "complete", "single"),
                                 unassignedLabel)
{
  levels0 = sort(unique(labels))
  levels = levels0 [levels0 != unassignedLabel]
  nClusters = length(levels);
  clusterDiss = matrix(0, nClusters, nClusters);

  distFunctions = c("mean", "max", "min");
  useFnc = match.fun(distFunctions[ match (match.arg(method), c("average", "complete", "single")) ]);

  if (nClusters > 1) 
  {
    for (c1 in 1:(nClusters-1)) for (c2 in (c1+1):nClusters)
      clusterDiss[c1, c2] = clusterDiss[c2, c1] = 
           useFnc(dissimilarity[ labels==levels[c1], labels==levels[c2]], na.rm = TRUE);
  }

  clusterDiss;
}

#===================================================================================================
#
# .nodeClusterDissimilarity
#
#===================================================================================================

# To make it easier to use the result, the function returns clusters in rows and nodes in columns.

.nodeClusterDissimilarity = function(dissim, labels, levels,
                                     method = c("average", "complete", "single"), unassignedLabel)
{
  unassigned = labels == unassignedLabel
  method = match.arg(method);

  nClusters = length(levels);
  nUnassigned = sum(unassigned);

  ncd = matrix(0, nClusters, nUnassigned);

  for (c in 1:nClusters)
  {
    if (method == "average") {
      ncd[ c, ] = colMeans( dissim[labels==levels[c], unassigned ], na.rm = TRUE);
    } else if (method == "complete") {
      ncd[ c, ] = apply( dissim[labels==levels[c], unassigned ], 2, max, na.rm = TRUE);
    } else 
      ncd[ c, ] = apply( dissim[labels==levels[c], unassigned ], 2, min, na.rm = TRUE);
  }
  ncd;
}
      

#===================================================================================================
#
# .mergeClustersIntoBlocks
#
#===================================================================================================

# .mergeClustersIntoBlocks: merge clusters into blocks. This code is adapted from WGCNA function
# projectiveKMeans. Assumes the labels are integers with no gaps starting from 1.

.mergeClustersIntoBlocks = function(dissim, labels, maxBlockSize, 
                                    method =  c("average", "complete", "single"),
                                    unassignedLabel)
{

  distFunctions = c("mean", "max", "min");
  useFnc = match.fun(distFunctions[ match (match.arg(method), c("average", "complete", "single")) ]);

  clusterDiss = .clusterDissimilarity(dissim, labels, method = method, unassignedLabel = unassignedLabel);
  diag(clusterDiss) = NA;

  clusterSizes = table(labels[ labels != unassignedLabel ] );
  nClusters = length(clusterSizes);
  clusterNames = names(clusterSizes);
  if (is.numeric(labels)) clusterNames = as.numeric(as.character(clusterNames));

  small = (clusterSizes < maxBlockSize);
  done = FALSE;
  while (!done & (sum(small)>1) & nClusters > 1)
  {
    smallIndex = c(1:nClusters)[small]
    nSmall = sum(small);

    distOrder = order(as.vector(clusterDiss[smallIndex, smallIndex]))[
                                seq(from=2, to = nSmall * (nSmall-1), by=2)];
    i = 1; canMerge = FALSE;
    while (i <= length(distOrder) && (!canMerge))
    {
      col = as.integer( (distOrder[i]-1)/nSmall + 1);
      whichJ = smallIndex[col];
      whichI = smallIndex[distOrder[i] - (col-1) * nSmall];
      canMerge = sum(clusterSizes[c(whichI, whichJ)]) < maxBlockSize;
      i = i + 1;
    }
    if (canMerge)
    {
      labels[labels==clusterNames[whichJ]] = clusterNames[whichI];
      clusterSizes[whichI] = sum(clusterSizes[c(whichI, whichJ)]);
      nClusters = nClusters -1;
      clusterSizes = clusterSizes[-whichJ];
      clusterNames = clusterNames[-whichJ];
      clusterDiss = clusterDiss[ -whichJ, -whichJ];

      for (c in 1:nClusters) if (c!=whichI)
        clusterDiss[whichI, c] = clusterDiss[c, whichI] = 
              useFnc(dissim[ labels==clusterNames[c], labels==clusterNames[whichI]], na.rm = TRUE)

      small = (clusterSizes < maxBlockSize);
    } else done = TRUE;
  }

  labels;
}

#===================================================================================================
#
# .assignToNearestCluster
#
#===================================================================================================

# Caution: if all labels equal unassignedLabel it will return the labels unchanged.

.assignToNearestCluster = function(dissim, labels, method, unassignedLabel)
{
  levels0 = sort(unique(labels));
  levels = levels0 [levels0 != unassignedLabel];
  if (length(levels) > 0)
  {
    nodeClusterSimilarity = .nodeClusterDissimilarity(dissim, labels=labels, levels=levels, method=method, 
                                                      unassignedLabel=unassignedLabel);
    nearest = .minWhichMin(nodeClusterSimilarity);
    labels[ labels == unassignedLabel ] = nearest[, "which"];
  }
  labels;
}

 
#===================================================================================================
#
# main user level function propensityClustering
#
#===================================================================================================

propensityClustering = function(adjacency, 
                              decompositionType = c("CPBA", "Pure Propensity"),
                              objectiveFunction = c("Poisson", "L2norm"),
                              fastUpdates = TRUE,
                              blocks = NULL,
                              initialClusters = NULL,
                              nClusters = NULL,
                              maxBlockSize = if (fastUpdates) 5000 else 1000,
                              clustMethod = "average",
                              cutreeDynamicArgs = list(deepSplit = 2,
                                                       minClusterSize = 20, verbose = 0),
                              dropUnassigned = TRUE,
                              unassignedLabel = 0,
                              verbose = 2,
                              indent = 0
                             )

{
  spaces = indentSpaces(indent);


  .checkAdjMat(adjacency, min = 0, max = max(adjacency, na.rm = TRUE));
  objectiveFunction = match.arg(objectiveFunction);
  nAllNodes = nNodes = ncol(adjacency);
  useNodes = rep(TRUE, nNodes);

  decompositionType = match.arg(decompositionType);

  if (decompositionType=="Pure Propensity")
  {
    # Run propensity decomposition on a single cluster that contains all nodes.
    nClusters = 1;
    initialClusters = rep(1, nNodes);
    return( CPBADecomposition(adjacency, clustering = initialClusters, 
                              objectiveFunction = objectiveFunction,
                              nClusters = nClusters) );
  }

  if (!is.null(nClusters)) 
  {
    # If the user supplies nClusters, use internal initialization.
    useInternalInit = TRUE
    # Check that the supplied nClusters makes sense.
    if (!is.finite(nClusters)) stop("The number of clusters 'nClusters' must be finite.");
    if (nClusters < 2) stop("If given, the number of clusters 'nClusters' must be at least 2.");
    # The following is necessary so splitting into blocks does not leave out any objects.
    dropUnassigned = FALSE
  } else {
    useInternalInit = FALSE;
  }

  # After this step the number of clusters is always non-null and is zero if it originally was NULL.

  if (!is.null(initialClusters))
  {
    initialClusters = as.vector(initialClusters);
    if (length(initialClusters)!=nNodes)
     stop(.spaste("Length of 'initialClusters' must equal the number of nodes\n",
                  "  (i.e., number of rows or columns of 'adjacency')."));
    tree = NULL;
  } else {
    dissim = 1-adjacency;
    # Cluster
    tree = hclust(as.dist(dissim), method = clustMethod);
    # Cut the tree
    cutreeDynamicArgs$dendro = tree;
    cutreeDynamicArgs$distM = dissim;
    initialClusters = do.call(cutreeDynamic, cutreeDynamicArgs);
    unassignedLabel = 0;
  } 

  if (all(initialClusters==unassignedLabel))
    stop(.spaste("All initial cluster labels are 'unassigned'.\n",
                "  Please supply a non-trivial initial clustering\n",
                "  or change the initial clustering arguments so hierarchical clustering\n",
                "  with Dynamic Tree Cut return clusters."));


  if (!is.null(blocks))
  {
    blocks = as.vector(blocks);
    if (length(blocks)!=nNodes)
     stop(.spaste("Length of 'blocks' must equal the number of nodes\n",
                  "  (i.e., number of rows or columns of 'adjacency')."));
  }

  if (dropUnassigned) 
  {
     useNodes = initialClusters != unassignedLabel;
     adjacency = adjacency[useNodes, useNodes];
     initialClusters = initialClusters[useNodes];
     if (!is.null(blocks))
        blocks = blocks[useNodes];
     nNodes = ncol(adjacency);
  } else {
    initialClusters = .assignToNearestCluster(1-adjacency, labels = initialClusters, 
                                method = clustMethod, unassignedLabel = unassignedLabel);
  }

  # Note: past this point initialClusters cannot contain any unassigned labels.
 
  if (is.null(blocks))
  {
    if (verbose > 0) printFlush(.spaste(spaces, " ..determining blocks.."));
    blocks = .mergeClustersIntoBlocks(1-adjacency, labels = initialClusters, maxBlockSize = maxBlockSize, 
                                    method = clustMethod, unassignedLabel = unassignedLabel);
  }

  # This code assumes there are no 0 labels among blocks.
  blocks = as.numeric(as.factor(blocks));
  nBlocks = length(unique(blocks));
  
  blockNodes = list();
  blockClusters = list();
  blockLevels = sort(unique(blocks));
  # Split given initial clusters by block for use below.
  for (b in 1:nBlocks)
  { 
    blockNodes1 = c(1:nNodes)[ blocks== blockLevels[b] ];
    blockNodes[[b]] = blockNodes1;
    blockClusters[[b]] = initialClusters[ blockNodes1 ];
  }

  propensityClusters = initialClusters;
  # Run propensity clustering on each block separately
  propClusts = list();
  if (verbose > 0) 
   printFlush(.spaste(spaces, " ..running propensityClustering in each block with non-trivial clustering.."));
  for (b in 1:nBlocks)
  {
    blockNodes1 = blockNodes[[b]];
    if (length(unique(initialClusters[ blockNodes1 ])) > 1 | useInternalInit)
    {
      if (verbose > 1) printFlush(.spaste(spaces, "   ..running propensityClustering in block ", b));
      initClust = initialClusters[ blockNodes1 ]
      initClust.norm = as.numeric(factor(initClust));
      norm2orig = .translationTable(initClust.norm, initClust);
      if (useInternalInit) 
      {
        ncl1 = nClusters;
      } else {
        ncl1 = length(unique(initClust.norm));
      }
    
      pc1 = .propensityClustering.internal( adjacency[ blockNodes1, blockNodes1 ], 
                                  initialClusters  = initClust.norm,
                                  l2bool = objectiveFunction=="L2norm",
                                  nClusters = ncl1,
                                  initbool = useInternalInit, fastUpdates = fastUpdates);
      if (useInternalInit)
      {
         # Each block has nClusters clusters, so this calculation is easy
         pc1$Clustering = pc1$Clustering + (b-1) * nClusters;
      } else {
         pc1$Clustering = .translateUsingTable(pc1$Clustering, norm2orig);
      }
      propensityClusters [ blockNodes1 ] = pc1$Clustering;
      propClusts[[b]] = pc1; 
    } else {
      propClusts[[b]] = NA;
      propensityClusters [ blockNodes1 ] = initialClusters [ blockNodes1 ];
    }
  }

  # Run one final propensity decomposition on the full clustering.

  if (verbose > 0) printFlush(.spaste(spaces, " ..running final propensity decomposition.."));
  propensity = rep(0, nAllNodes);
  propensityClusters.norm = as.numeric(as.factor(propensityClusters));
  pd = CPBADecomposition(adjacency, propensityClusters.norm, 
                         objectiveFunction = objectiveFunction,
                         nClusters = NULL);
  propensity[useNodes] = pd$Propensity;

  if (verbose > 0) printFlush(.spaste(spaces, " ..done (propensityClustering)."));

  blocks.all = initialClusters.all = propClusters.all = rep(0, nAllNodes);
  propClusters.all[useNodes] = propensityClusters;
  blocks.all[useNodes] = blocks;
  initialClusters.all[useNodes] = initialClusters;

  # Return value
  c(list(Clustering = propClusters.all, Propensity = propensity, NodeWasConsidered = useNodes),
          pd, list(Blocks = blocks.all, InitialClusters = initialClusters.all,
          InitialTree = tree));

}

                              
