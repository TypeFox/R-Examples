
# Attempt to speed up the function a bit. There is some on-the-fly growing of vectors that can be removed,
# and perhaps some variables could be deleted from the branch structure.
#  - remove Clusters, LastMerge
#  - add counters keeping track o number of merges, number of singletons, number of basic clusters

# ClusterTrim is removed for now; it would have been a bit more complicated to make it work with the
# enhanced PAM stage. May be re-introduced later if anyone shows any interest in this.

# Tree cut

# 1.11-2
#    . Bug fix in interpretation of deepSplit fixed

# 1.11-1
#    . If no merge lies below the cut height, simply exit with all labels=0 instead of throwing an
#    error.
#    . If cutHeight is above the highest merge, set it to the highest merge.

# 1.11
#    . change the default cutHeight to 99% of dendogram height range
 
# 1.10-02:
#    . if distM isn't given, method defaults to "tree".

# 1.10:
#    . Bug fix: since distM us necessary in the first stage as well, the function now complains if
#    distM is not given

# 1.09:
#    . Bug fix: number of Unassigned = 1 is now handled correctly

# 1.08: 
#    . Changing the meaning of minGap, maxCoreScatter: both now relative (i.e., fractions). minGap will
#    be the minimum cluster gap expressed as a fraction of the range between a certain quantile of the
#    merging heights and the cutHeight; maxCoreScatter will be interpreted in the same way. Adding two
#    more parameters, namely minAbsGap, maxAbsCoreScatter that can be used to supply hitherto used
#    absolute values. If they are given, they override minGap and maxCoreScatter, respectively. 
#    . Change of names: clusterMinSize -> minClusterSize
#    . Default value of minClusterSize now 20
#    . Fixed stage 2 labeling for non-medoids: ClusterDiam{eter} now given is the maximum of average
#      distances of points to the rest of the cluster.
#    . Default for maxDistToLabel is now cutHeight even though it may mean that some objects above the
#      cutHeight may be labeled. There doesn't seem to be a clean and simple way to only label objects
#      whose merging distance to (the continuation of) the cluster is below cutHeight
#    . If cutHeight not given, will be set to max(dendro$height) for the DynamicTree as well.

# 1.07: Changing parameter names to be more intuitive and in accord with generally used terminology
#       Also changing internal function names by prepending a . to them
#    . All extra functions (EvaluateCLusters etc) removed.
#    . Rename the function from GetClusters to cutreeHybrid

# 1.06: The tree cut is exactly the same as 1.05, but the color handling is relegated to NetworkFunctions.
#       The ColorsFromLabels in NetworkFunctions use a slightly different input format; 
#           in particular, 0 is considered grey. This means, among other things, that 
# 1.05:
#   . The PAM stage is changed: instead of calculating distances to medoids, will calculate average
#     distances to clusters. This is intuitively less desirable than the medoids, but simulations
#     seem to indicate that large clusters are too greedy. The average linkage may help that a bit.
#     It is not quite clear that cluster trimming makes a lot of sense in this scenario, but I'll
#     keep it in (assigning elements based on average distance to a trimmed cluster is not quite the
#     same as an untrimmed cluster, even for the elements that were trimmed). 
#   . Improve trimming: instead of the lowest joining heights, keep all singletons up to the first joined
#     object (branch or singleton) whose merging height is above the threshold; the rest (content of all
#     branches merged higher than threshold) is trimmed.

# 1.04:
#   . Implement Bin's idea that small branches unassigned in Stage 1 should not be broken
#     up by Stage 2. Implemented as follows: Only keep those branches unbroken whose only reason for not
#     being a cluster is that they don't have enough objects.
#     While going through the merge tree: mark branches that are (1) merged into composite clusters, (2)
#     not clusters themselves because of failing the minimum size requirement only.
#   . Change boudaries of clusters. Instead of regarding everything up to the merge automatically part
#     the respective cluster, only elements whose joining heights are less than a cutoff given for
#     the cluster are considered part of the cluster automatically; the rest is assigned in pam-like
#     manner. 

# 1.03: 
#   . Changing the definition of the core from the most connected points to the first points added
#     to the cluster. Note that the order depends on how branches are merged, so that better be
#     correct as well. This makes the core stable against adding outliers to the cluster.

# 1.02.01: fixing a bug in the main function that was referencing nonexistent heights member of
#          dendro.
#    . Fixing a correctness issue: the core average distance is the average distance between points
#      within the core, not the average distance of points within the core to all points in the
#      cluster. 

# 1.02: 
#    . Changing core size: instead of minClusterSize + sqrt(Size - minClusterSize) it is now
#      CoreSize  + sqrt(Size - CoreSize), with CoreSize = as.integer(minClusterSize) + 1.


# 1.01.03: Fixing memory usage and a bug in which singletons were added twice (and more times).
#  In this version, to simplify things, Singletons are only kept for basic clusters. For composite
#  clusters they are NULL; CoreScatter should never be called for composite clusters as that can get
#  hugely time and memory intensive.
# Another bug concerning couting unassigned points fixed.

# -02: adding distance matrix to the parameters of GetClusters.
# -03: Merging GetClusters and AssignLabel together; changes in variable names to make code more
# readable.


# Progress indicator...

.initProgInd = function( leadStr = "..", trailStr = "", quiet = !interactive())
{
  oldStr = " ";
  cat(oldStr);
  progInd = list(oldStr = oldStr, leadStr = leadStr, trailStr = trailStr);
  class(progInd) = "progressIndicator";
  .updateProgInd(0, progInd, quiet);
}

.updateProgInd = function(newFrac, progInd, quiet = !interactive())
{
  if (class(progInd)!="progressIndicator")
    stop("Parameter progInd is not of class 'progressIndicator'. Use initProgInd() to initialize",
         "it prior to use.");

  newStr = paste(progInd$leadStr, as.integer(newFrac*100), "% ", progInd$trailStr, sep = "");
  if (newStr!=progInd$oldStr)
  {
    if (quiet)
    {
      progInd$oldStr = newStr;
    } else {
      cat(paste(rep("\b", nchar(progInd$oldStr)), collapse=""));
      cat(newStr);
      if (exists("flush.console")) flush.console();
      progInd$oldStr = newStr;
    }
  }
  progInd;
}


# The following are supporting function for GetClusters. 

.CoreSize = function(BranchSize, minClusterSize)
{
  BaseCoreSize = minClusterSize/2 + 1;
  if (BaseCoreSize < BranchSize)
  {
    CoreSize = as.integer(BaseCoreSize + sqrt(BranchSize - BaseCoreSize));
  } else CoreSize = BranchSize;
  CoreSize;
}
  
# This assumes the diagonal of the distance matrix
# is zero, BranchDist is a square matrix whose dimension is at least 2.

.CoreScatter = function(BranchDist, minClusterSize)
{
  nPoints = dim(BranchDist)[1];
  PointAverageDistances = colSums(BranchDist) / (nPoints-1);
  CoreSize = minClusterSize/2 + 1;
  if (CoreSize < nPoints)
  {
    EffCoreSize = as.integer(CoreSize + sqrt(nPoints - CoreSize));
    ord = order(PointAverageDistances);
    Core = ord[c(1:EffCoreSize)];
  } else {
    Core = c(1:nPoints);
    EffCoreSize = nPoints;
  }
  CoreAverageDistances = colSums(BranchDist[Core, Core]) / (EffCoreSize-1);
  mean(CoreAverageDistances);
}

.interpolate = function(data, index)
{
  i = round(index);
  n = length(data);
  if (i<1) return(data[1]);
  if (i>=n) return(data[n]);

  r = index-i;
  data[i] * (1-r) + data[i+1] * r;
}

.chunkSize = 100;

#-------------------------------------------------------------------------------------------
#
# cutreeHybrid
#
#-------------------------------------------------------------------------------------------
# Traverses a given clustering tree and detects branches whose size is at least minClusterSize, average
# singleton joining height is at most maxCoreScatter and split (attaching height minus average
# height) is at least minGap. If cutHeight is set, all clusters are cut at that height.

# stabilityDistM is assumed normalized to the maximum possible distance, i.e. range between 0 and 1 and
# distance 1 means the two genes were put in separate modules in all sampled clusterings. 

cutreeHybrid = function(
      # Input data: basic
      dendro, distM, 

      # Branch cut criteria and options
      cutHeight = NULL, minClusterSize = 20, deepSplit = 1,

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
      useMedoids = FALSE,
      maxPamDist = cutHeight,
      respectSmallClusters = TRUE, 

      # Various options
      verbose = 2, indent = 0)
{

  spaces = indentSpaces(indent);

  # if (verbose>0) printFlush(paste(spaces, "cutreeHybrid cluster detection starting..."));

  # No. of merges in the tree
  nMerge = length(dendro$height);
  if (nMerge < 1)
    stop("The given dendrogram is suspicious: number of merges is zero.");

  if (is.null(distM)) stop("distM must be non-NULL")

  if (is.null(dim(distM)))
    stop("distM must be a matrix.");

  if ( (dim(distM)[1] != nMerge+1) | (dim(distM)[2]!=nMerge+1) ) 
    stop("distM has incorrect dimensions.");

  if (pamRespectsDendro & !respectSmallClusters)
    printFlush(paste("cutreeHybrid Warning: parameters pamRespectsDendro (TRUE)", 
                     "and respectSmallClusters (FALSE) imply contradictory intent.\n",
                     "Although the code will work, please check you really",
                     "intented these settings for the two arguments."));

  if (any(diag(distM!=0))) diag(distM) = 0;

  refQuantile = 0.05;
  refMerge = round(nMerge * refQuantile);
  if (refMerge < 1) refMerge = 1;
  refHeight = dendro$height[refMerge]; 
  
  if (is.null(cutHeight))
  {
    cutHeight = 0.99 * (max(dendro$height) - refHeight) + refHeight;
    if (verbose>0) 
      printFlush(paste(spaces, 
            "..cutHeight not given, setting it to", signif(cutHeight,3), 
            " ===>  99% of the (truncated) height range in dendro."));
  } else {
    if (cutHeight > max(dendro$height)) cutHeight = max(dendro$height);
  }

  # If maxPamDist is not set, set it to cutHeight

  if (is.null(maxPamDist)) maxPamDist = cutHeight;

  nMergeBelowCut = sum(dendro$height <= cutHeight);
  if (nMergeBelowCut < minClusterSize) 
  {
    if (verbose>0) printFlush(paste(spaces, "cutHeight set too low: no merges below the cut."));
    return(list(labels = rep(0, times = nMerge+1)))
  }

  # Check external branch split function(s), if given.

  nExternalSplits = length(externalBranchSplitFnc);
  if (nExternalSplits>0)
  {
    if (length(minExternalSplit)<1) stop("'minExternalBranchSplit' must be given.");
    if (assumeSimpleExternalSpecification && nExternalSplits==1)
    {
      externalSplitOptions = list(externalSplitOptions);
    }
    externalBranchSplitFnc = lapply(externalBranchSplitFnc, match.fun);
    for (es in 1:nExternalSplits)
    {
      externalSplitOptions[[es]]$tree = dendro;
      if (length(externalSplitFncNeedsDistance)==0 || 
           externalSplitFncNeedsDistance[es]) externalSplitOptions[[es]]$dissimMat = distM;
    }
  } 

  MxBranches = nMergeBelowCut
  branch.isBasic = rep(TRUE, MxBranches);
  branch.isTopBasic = rep(TRUE, MxBranches);
  branch.failSize = rep(FALSE, MxBranches);
  branch.rootHeight = rep(NA, MxBranches);
  branch.size = rep(2, MxBranches);
  branch.nMerge = rep(1, MxBranches);
  branch.nSingletons = rep(2, MxBranches);
  branch.nBasicClusters = rep(0, MxBranches);
  branch.mergedInto = rep(0, MxBranches);
  branch.attachHeight = rep(NA, MxBranches);
  branch.singletons = vector(mode = "list", length = MxBranches);
  branch.basicClusters = vector(mode = "list", length = MxBranches);
  branch.mergingHeights = vector(mode = "list", length = MxBranches);
  branch.singletonHeights = vector(mode = "list", length = MxBranches);

  # Note: size equals the number of singletons on a basic branch, but not on a composite branch.
                  
  nBranches = 0;

  spyIndex = NULL;
  if (file.exists(".dynamicTreeCutSpyFile"))
  {
     spyIndex = read.table(".dynamicTreeCutSpyFile", header = FALSE);
     printFlush("Found 'spy file' with indices of objects to watch for.");
     spyIndex= as.numeric(spyIndex[, 1]);
     printFlush(paste(spyIndex, collapse = ", "));
  }

  # Default values for maxCoreScatter and minGap:

  defMCS = c(0.64, 0.73, 0.82, 0.91, 0.95); 	# Default max core scatter
  defMG = (1-defMCS)*3/4;			# Default minimum gap

  nSplitDefaults = length(defMCS);

  # Convert deep split to range 1..5
  if (is.logical(deepSplit)) deepSplit = as.integer(deepSplit)*(nSplitDefaults - 2);
  deepSplit = deepSplit + 1;
  if ((deepSplit<1) | (deepSplit>nSplitDefaults))
    stop(paste("Parameter deepSplit (value", deepSplit, 
               ") out of range: allowable range is 0 through", nSplitDefaults-1));

  # If not set, set the cluster gap and core scatter according to deepSplit.
  if (is.null(maxCoreScatter)) maxCoreScatter = .interpolate(defMCS, deepSplit)
  if (is.null(minGap)) minGap = .interpolate(defMG, deepSplit);

  # Convert (relative) minGap and maxCoreScatter to corresponding absolute quantities if the latter were
  # not given.

  if (is.null(maxAbsCoreScatter))
     maxAbsCoreScatter = refHeight + maxCoreScatter * (cutHeight - refHeight);
  if (is.null(minAbsGap))
     minAbsGap = minGap * (cutHeight - refHeight);

  # if minSplitHeight was not given, set it to 0
  if (is.null(minSplitHeight)) minSplitHeight = 0;
  # Convert (relative) minSplitHeight to corresponding minAbsSplitHeight if the latter was not given.
  if (is.null(minAbsSplitHeight))
     minAbsSplitHeight = refHeight + minSplitHeight * (cutHeight - refHeight);
    
  nPoints = nMerge+1;

  # For each merge, record the cluster that it belongs to
  IndMergeToBranch = rep(0, times = nMerge)

  # For each object that joins a composite branch, record the number of the branch
  onBranch = rep(0, nPoints);

  # The root
  RootBranch = 0;
  if (verbose>2) 
  {
     printFlush(paste(spaces, "..Going through the merge tree"));
     pind = .initProgInd();
  }

  mergeDiagnostics = data.frame(smI = rep(NA, nMerge), smSize = rep(NA, nMerge), 
                                smCrSc = rep(NA, nMerge), smGap = rep(NA, nMerge), 
                                lgI = rep(NA, nMerge), lgSize = rep(NA, nMerge), 
                                lgCrSc = rep(NA, nMerge), lgGap = rep(NA, nMerge),
                                merged = rep(NA, nMerge));

  if (nExternalSplits > 0)
  {
     externalMergeDiags = matrix(NA, nMerge, nExternalSplits);
     colnames(externalMergeDiags) = paste("externalBranchSplit", 1:nExternalSplits, sep = ".");
  }

  extender = rep(0, .chunkSize);

  for (merge in 1:nMerge) if (dendro$height[merge]<=cutHeight)
  {
    # are both merged objects sigletons?
    if (dendro$merge[merge,1]<0 & dendro$merge[merge,2]<0)
    {
      # Yes; start a new branch.
      nBranches = nBranches + 1;
      branch.isBasic[nBranches] = TRUE;
      branch.isTopBasic[nBranches] = TRUE;
      branch.singletons[[nBranches]] = c(-dendro$merge[merge,], extender);
      branch.basicClusters[[nBranches]] = extender;
      branch.mergingHeights[[nBranches]] = c(rep(dendro$height[merge], 2), extender);
      branch.singletonHeights[[nBranches]] = c(rep(dendro$height[merge], 2), extender);
      IndMergeToBranch[merge] = nBranches;
      RootBranch = nBranches;
    } else if (sign(dendro$merge[merge,1]) * sign(dendro$merge[merge,2]) <0)
    {
      # merge the sigleton into the branch
      clust = IndMergeToBranch[max(dendro$merge[merge,])];
      if (clust==0) stop("Internal error: a previous merge has no associated cluster. Sorry!");
      gene = -min(dendro$merge[merge,]);
      ns = branch.nSingletons[clust] + 1;
      nm = branch.nMerge[clust] + 1;
      if (branch.isBasic[clust]) 
      {
        if (ns>length(branch.singletons[[clust]]))
        {
          branch.singletons[[clust]] = c(branch.singletons[[clust]], extender);
          branch.singletonHeights[[clust]] = c(branch.singletonHeights[[clust]], extender)
        }
        branch.singletons[[clust]] [ns] = gene;
        branch.singletonHeights[[clust]] [ns] = dendro$height[merge]
      } else {
        onBranch[gene] = clust;
      }
      if (nm >= length(branch.mergingHeights[[clust]]))
         branch.mergingHeights[[clust]] = c(branch.mergingHeights[[clust]], extender)
      branch.mergingHeights[[clust]] [nm] = dendro$height[merge];
      branch.size[clust] = branch.size[clust] + 1;
      branch.nMerge[clust] = nm;
      branch.nSingletons[clust] = ns;
      IndMergeToBranch[merge] = clust;
      RootBranch = clust;
    } else {
      # attempt to merge two branches:
      clusts = IndMergeToBranch[dendro$merge[merge,]];
      sizes = branch.size[clusts];
      # Note: for 2 elements, rank and order are the same.
      rnk = rank(sizes, ties.method = "first");
      small = clusts[rnk[1]]; large = clusts[rnk[2]];
      sizes = sizes[rnk];
      branch1 = branch.singletons[[large]] [1:sizes[2]];
      branch2 = branch.singletons[[small]] [1:sizes[1]];
      spyMatch = FALSE;
      if (!is.null(spyIndex))
      {
        n1 = length(intersect(branch1, spyIndex));
        if ( (n1/length(branch1) > 0.99 && n1/length(spyIndex) > 0.99) )
        {
          printFlush(paste("Found spy match for branch 1 on merge", merge))
          spyMatch = TRUE
        }
        n2 = length(intersect(branch2, spyIndex));
        if ( (n2/length(branch1) > 0.99 && n2/length(spyIndex) > 0.99) )
        {
          printFlush(paste("Found spy match for branch 2 on merge", merge))
          spyMatch = TRUE
        }
      }

      if (branch.isBasic[small])
      {
         coresize = .CoreSize(branch.nSingletons[small], minClusterSize);
         Core = branch.singletons[[small]] [c(1:coresize)];
         # SmAveDist = mean(apply(distM[Core, Core], 2, sum)/(coresize-1)); 
         SmAveDist = mean(colSums(distM[Core, Core, drop = FALSE])/(coresize-1));        
      } else { SmAveDist = 0; }
       
      if (branch.isBasic[large])
      {
         coresize = .CoreSize(branch.nSingletons[large], minClusterSize);
         Core = branch.singletons[[large]] [c(1:coresize)];
         LgAveDist = mean(colSums(distM[Core, Core])/(coresize-1)); 
      } else { LgAveDist = 0; }

      mergeDiagnostics[merge, ] = c(small, branch.size[small], SmAveDist, 
                                     dendro$height[merge] - SmAveDist,
                                     large, branch.size[large], LgAveDist, 
                                     dendro$height[merge] - LgAveDist,
                                     NA);
      # We first check each cluster separately for being too small, too diffuse, or too shallow:
      SmallerScores = c(branch.isBasic[small], 
                        branch.size[small] < minClusterSize,
                        SmAveDist > maxAbsCoreScatter, 
                        dendro$height[merge] - SmAveDist < minAbsGap,
                        dendro$height[merge] < minAbsSplitHeight );

      
      if ( SmallerScores[1] * sum(SmallerScores[-1]) > 0 )
      {
        DoMerge = TRUE;
        SmallerFailSize = !(SmallerScores[3] | SmallerScores[4]);  # Smaller fails only due to size
      } else 
      {
        LargerScores = c(branch.isBasic[large], 
                         branch.size[large] < minClusterSize,
                         LgAveDist > maxAbsCoreScatter, 
                         dendro$height[merge] - LgAveDist < minAbsGap,
                         dendro$height[merge] < minAbsSplitHeight );
        if ( LargerScores[1] * sum(LargerScores[-1]) > 0 )
        { # Actually: the large one is the one to be merged
          DoMerge = TRUE;
          SmallerFailSize = !(LargerScores[3] | LargerScores[4]);  # cluster fails only due to size
          x = small; small = large; large = x;
          sizes = rev(sizes);
        } else {
          DoMerge = FALSE; # None of the two satisfies merging criteria
        }
      }

      if (DoMerge)
      {
        mergeDiagnostics$merged[merge] = 1;
      } #else
        #browser();

      # Still not merging? If user-supplied criterion is given, check whether the criterion is below the
      # specified threshold.
      
      if (!DoMerge && (nExternalSplits > 0) && branch.isBasic[small] && branch.isBasic[large])
      {
        if (verbose > 4) printFlush(paste0("Entering external split code on merge ", merge));
        branch1 = branch.singletons[[large]] [1:sizes[2]];
        branch2 = branch.singletons[[small]] [1:sizes[1]];
          
        if (verbose > 4 | spyMatch) printFlush(paste0("  ..branch lengths: ", sizes[1], ", ", sizes[2]))
        #if (any(is.na(branch1)) || any(branch1==0)) browser();
        #if (any(is.na(branch2)) || any(branch2==0)) browser();
        es = 0;
        while (es < nExternalSplits && !DoMerge)
        {
          es = es + 1;
          args = externalSplitOptions[[es]];
          args = c(args, list(branch1 = branch1, branch2 = branch2));
          extSplit = do.call(externalBranchSplitFnc[[es]], args);
          if (spyMatch)
            printFlush(" .. external criterion ", es, ": ", extSplit);
          DoMerge = extSplit < minExternalSplit[es];
          externalMergeDiags[merge, es] = extSplit;
          mergeDiagnostics$merged[merge] = if (DoMerge) 2 else 0;
        }
      }
        
      if (DoMerge)
      {
        # merge the small into the large cluster and close it.
        branch.failSize[[small]] = SmallerFailSize;
        branch.mergedInto[small] = large; 
        branch.attachHeight[small] = dendro$height[merge];
        branch.isTopBasic[small] = FALSE;
        nss = branch.nSingletons[small];
        nsl = branch.nSingletons[large];
        ns = nss + nsl;
        if (branch.isBasic[large]) 
        {
          nExt = ceiling(  (ns - length(branch.singletons[[large]]))/.chunkSize  );
          if (nExt > 0) 
          {
            if (verbose > 5) 
              printFlush(paste("Extending singletons for branch", large, "by", nExt, " extenders."));
            branch.singletons[[large]] = c(branch.singletons[[large]], rep(extender, nExt));
            branch.singletonHeights[[large]] = c(branch.singletonHeights[[large]], rep(extender, nExt));
          }
          branch.singletons[[large]] [(nsl+1):ns] = branch.singletons[[small]][1:nss];
          branch.singletonHeights[[large]] [(nsl+1):ns] = branch.singletonHeights[[small]][1:nss];
          branch.nSingletons[large] = ns;
        } else {
          if (!branch.isBasic[small])
            stop("Internal error: merging two composite clusters. Sorry!");
          onBranch[ branch.singletons[[small]] ] = large;
        }
        nm = branch.nMerge[large] + 1;
        if (nm > length(branch.mergingHeights[[large]]))
           branch.mergingHeights[[large]] = c(branch.mergingHeights[[large]], extender);
        branch.mergingHeights[[large]] [nm] = dendro$height[merge];
        branch.nMerge[large] = nm;
        branch.size[large] = branch.size[small] + branch.size[large];
        IndMergeToBranch[merge] = large;
        RootBranch = large;
      } else {
        # start or continue a composite cluster.

        # If large is basic and small is not basic, switch them.
        if (branch.isBasic[large] & !branch.isBasic[small])
        {
          x = large; large = small; small = x;
          sizes = rev(sizes);
        }

        # Note: if pamRespectsDendro, need to start a new composite cluster every time two branches merge,
        # otherwise will not have the necessary information.
        # Otherwise, if the large cluster is already composite, I can simply merge both clusters into 
        # one of the non-composite clusters.

        if (branch.isBasic[large] | (pamStage & pamRespectsDendro))
        {
          nBranches = nBranches + 1;
          branch.attachHeight[c(large, small)] = dendro$height[merge];
          branch.mergedInto[c(large, small)] = nBranches;
          if (branch.isBasic[small])
          {
            addBasicClusters = small; # add basic clusters
          } else 
            addBasicClusters = branch.basicClusters[[small]];
          if (branch.isBasic[large])
          {
            addBasicClusters = c(addBasicClusters, large);
          } else 
            addBasicClusters = c(addBasicClusters, branch.basicClusters[[large]]);
          # print(paste("  Starting a composite cluster with number", nBranches));
          branch.isBasic[nBranches] = FALSE;
          branch.isTopBasic[nBranches] = FALSE;
          branch.basicClusters[[nBranches]] = addBasicClusters;
          branch.mergingHeights[[nBranches]] = c(rep(dendro$height[merge], 2), extender);
          branch.nMerge[nBranches] = 2;
          branch.size[nBranches] = sum(sizes);
          branch.nBasicClusters[nBranches] = length(addBasicClusters);
          IndMergeToBranch[merge] = nBranches;
          RootBranch = nBranches;
        } else {
          # Add small branch to the large one 
          addBasicClusters = if (branch.isBasic[small]) small else branch.basicClusters[[small]];
          nbl = branch.nBasicClusters[large];
          nb = branch.nBasicClusters[large] + length(addBasicClusters);
          if (nb > length(branch.basicClusters[[large]]))
          {
            nExt = ceiling(  ( nb - length(branch.basicClusters[[large]]))/.chunkSize);
            branch.basicClusters[[large]] = c(branch.basicClusters[[large]], rep(extender, nExt));
          }
          branch.basicClusters[[large]] [(nbl+1):nb] = addBasicClusters;
          branch.nBasicClusters[large] = nb;   
          branch.size[large] = branch.size[large] + branch.size[small];
          nm = branch.nMerge[large] + 1;
          if (nm > length(branch.mergingHeights[[large]]))
            branch.mergingHeights[[large]] = c(branch.mergingHeights[[large]], extender); 
          branch.mergingHeights[[large]] [nm] = dendro$height[merge];
          branch.nMerge[large] = nm;
          branch.attachHeight[small] = dendro$height[merge];
          branch.mergedInto[small] = large; 
          IndMergeToBranch[merge] = large;
          RootBranch = large;
        }
      }
    }
    if (verbose > 2) pind = .updateProgInd(merge/nMerge, pind);
  }
  if (verbose > 2) 
  {
    pind = .updateProgInd(1, pind);
    printFlush("");
  }


  if (verbose>2) printFlush(paste(spaces, "..Going through detected branches and marking clusters.."));
  
  isCluster = rep(FALSE, times = nBranches);
  SmallLabels = rep(0, times = nPoints);
  for (clust in 1:nBranches)
  {
    if (is.na(branch.attachHeight[clust])) branch.attachHeight[clust] = cutHeight;
    if (branch.isTopBasic[clust])
    {
       coresize = .CoreSize(branch.nSingletons[clust], minClusterSize);
       Core = branch.singletons[[clust]] [c(1:coresize)];
       CoreScatter = mean(colSums(distM[Core, Core])/(coresize-1)); 
       isCluster[clust] = branch.isTopBasic[clust] & (branch.size[clust] >= minClusterSize) & 
                  (CoreScatter < maxAbsCoreScatter) &
                  (branch.attachHeight[clust] - CoreScatter > minAbsGap);
    } else { CoreScatter = 0; }
    if (branch.failSize[clust]) SmallLabels[branch.singletons[[clust]]] = clust;
  }
  if (!respectSmallClusters) SmallLabels = rep(0, times = nPoints);

  if (verbose>2) printFlush(paste(spaces, "..Assigning Tree Cut stage labels.."));

  Colors = rep(0, times = nPoints);
  coreLabels = rep(0, times = nPoints);
  clusterBranches = c(1:nBranches)[isCluster];
  branchLabels = rep(0, nBranches);
  color = 0;
  for (clust in clusterBranches)
  {
    color = color+1;
    Colors[branch.singletons[[clust]]] = color;
    SmallLabels[branch.singletons[[clust]]] = 0;
    coresize = .CoreSize(branch.nSingletons[clust], minClusterSize);
    Core = branch.singletons[[clust]] [c(1:coresize)];
    coreLabels[Core] = color;
    branchLabels[clust] = color;
  } 

  Labeled = c(1:nPoints)[Colors!=0];
  Unlabeled = c(1:nPoints)[Colors==0];
  nUnlabeled = length(Unlabeled);
  UnlabeledExist = (nUnlabeled>0);
  if (length(Labeled)>0)
  {
    LabelFac = factor(Colors[Labeled]);
    nProperLabels = nlevels(LabelFac);
  } else
    nProperLabels = 0;

  if (pamStage & UnlabeledExist & nProperLabels>0)
  {
     if (verbose>2) printFlush(paste(spaces, "..Assigning PAM stage labels.."));
     nPAMed = 0;
     # Assign some of the grey genes to the nearest module. Define nearest as the distance to the medoid,
     # that is the point in the cluster that has the lowest average distance to all other points in the
     # cluster. First get the medoids.
     if (useMedoids)
     {
        Medoids = rep(0, times = nProperLabels);
        ClusterRadii = rep(0, times = nProperLabels);	
        for (cluster in 1:nProperLabels)
        {
          InCluster = c(1:nPoints)[Colors==cluster];
          DistInCluster = distM[InCluster, InCluster];
          DistSums = colSums(DistInCluster);
          Medoids[cluster] = InCluster[which.min(DistSums)];
          ClusterRadii[cluster] = max(DistInCluster[, which.min(DistSums)])
        }
        # If small clusters are to be respected, assign those first based on medoid-medoid distances.
        if (respectSmallClusters)
        {
          FSmallLabels = factor(SmallLabels);
          SmallLabLevs = as.numeric(levels(FSmallLabels));
          nSmallClusters = nlevels(FSmallLabels) - (SmallLabLevs[1]==0);
          if (nSmallClusters>0) for (sclust in SmallLabLevs[SmallLabLevs!=0])
          {
            InCluster = c(1:nPoints)[SmallLabels==sclust];
            if (pamRespectsDendro)
            {
              onBr = unique(onBranch[InCluster]);
              if (length(onBr)>1)
                stop(paste("Internal error: objects in a small cluster are marked to belong",
                           "\nto several large branches:", paste(onBr, collapse = ", ")));
              if (onBr > 0)
              {
                 basicOnBranch = branch.basicClusters[[onBr]];
                 labelsOnBranch = branchLabels[basicOnBranch]
              } else {
                 labelsOnBranch = NULL;    
              }
            } else {
                labelsOnBranch = c(1:nProperLabels)
            }
            # printFlush(paste("SmallCluster", sclust, "has", length(InCluster), "elements."));
            DistInCluster = distM[InCluster, InCluster, drop = FALSE];
            if (length(labelsOnBranch) > 0)
            {
               if (length(InCluster)>1)
               {
                 DistSums = apply(DistInCluster, 2, sum);
                 smed = InCluster[which.min(DistSums)];
                 DistToMeds = distM[Medoids[labelsOnBranch], smed];
                 closest = which.min(DistToMeds);
                 DistToClosest = DistToMeds[closest];
                 closestLabel = labelsOnBranch[closest]
                 if ( (DistToClosest < ClusterRadii[closestLabel]) | (DistToClosest <  maxPamDist) )
                 {
                   Colors[InCluster] = closestLabel;
                   nPAMed = nPAMed + length(InCluster);
                 } else Colors[InCluster] = -1;  # This prevents individual points from being assigned later 
               }
            } else
              Colors[InCluster] = -1;
          }
        }
        # Assign leftover unlabeled objects to clusters with nearest medoids
        Unlabeled = c(1:nPoints)[Colors==0];
        if (length(Unlabeled>0)) for (obj in Unlabeled)
        {
          if (pamRespectsDendro)
          {
            onBr = onBranch[obj]
            if (onBr > 0)
            {
              basicOnBranch = branch.basicClusters[[onBr]];
              labelsOnBranch = branchLabels[basicOnBranch]
            } else {
              labelsOnBranch = NULL;
            }
          } else {
              labelsOnBranch = c(1:nProperLabels)
          }
          if (!is.null(labelsOnBranch))
          {
             UnassdToMedoidDist = distM[Medoids[labelsOnBranch], obj];
             nearest= which.min(UnassdToMedoidDist)
             NearestCenterDist = UnassdToMedoidDist[nearest];
             nearestMed = labelsOnBranch[nearest]
             if ( (NearestCenterDist < ClusterRadii[nearestMed]) |
                  (NearestCenterDist < maxPamDist))
             {
               Colors[obj] = nearestMed;
               nPAMed = nPAMed + 1;
             }
          }
        }
        UnlabeledExist = (sum(Colors==0)>0);
     } else # Instead of medoids, use average distances
     { # This is the default method, so I will try to tune it for speed a bit.
        ClusterDiam = rep(0, times = nProperLabels);	
        for (cluster in 1:nProperLabels)
        {
          InCluster = c(1:nPoints)[Colors==cluster];
          nInCluster = length(InCluster)
          DistInCluster = distM[InCluster, InCluster];
          if (nInCluster>1) {
             AveDistInClust = colSums(DistInCluster)/(nInCluster-1);
             ClusterDiam[cluster] = max(AveDistInClust);
          } else {
             ClusterDiam[cluster] = 0;
          }
        }
        # If small clusters are respected, assign them first based on average cluster-cluster distances.
        ColorsX = Colors;
        if (respectSmallClusters)
        {
          FSmallLabels = factor(SmallLabels);
          SmallLabLevs = as.numeric(levels(FSmallLabels));
          nSmallClusters = nlevels(FSmallLabels) - (SmallLabLevs[1]==0);
          if (nSmallClusters>0)
          {
            if (pamRespectsDendro)
            {
              for (sclust in SmallLabLevs[SmallLabLevs!=0])
              {
                InCluster = c(1:nPoints)[SmallLabels==sclust];
                onBr = unique(onBranch[InCluster]);
                if (length(onBr)>1)
                  stop(paste("Internal error: objects in a small cluster are marked to belong",
                             "\nto several large branches:", paste(onBr, collapse = ", ")));
                if (onBr > 0)
                {
                  basicOnBranch = branch.basicClusters[[onBr]];
                  labelsOnBranch = branchLabels[basicOnBranch]
                  useObjects = ColorsX %in% unique(labelsOnBranch)
                  DistSClustClust = distM[InCluster, useObjects, drop = FALSE];
                  MeanDist = colMeans(DistSClustClust);
                  useColorsFac = factor(ColorsX[useObjects])
                  MeanMeanDist = tapply(MeanDist, useColorsFac, mean);
                  nearest = which.min(MeanMeanDist);
                  NearestDist = MeanMeanDist[nearest];
                  nearestLabel = as.numeric(levels(useColorsFac)[nearest])
                  if ( ((NearestDist < ClusterDiam[nearestLabel]) | (NearestDist <  maxPamDist)) )
                  {
                    Colors[InCluster] = nearestLabel;
                    nPAMed = nPAMed + length(InCluster);
                  } else Colors[InCluster] = -1;  # This prevents individual points from being assigned later
                }
              }    
            } else {
              labelsOnBranch = c(1:nProperLabels)
              useObjects = c(1:nPoints)[ColorsX!=0];
              for (sclust in SmallLabLevs[SmallLabLevs!=0])
              {
                InCluster = c(1:nPoints)[SmallLabels==sclust];
                DistSClustClust = distM[InCluster, useObjects, drop = FALSE];
                MeanDist = colMeans(DistSClustClust);
                useColorsFac = factor(ColorsX[useObjects])
                MeanMeanDist = tapply(MeanDist, useColorsFac, mean);
                nearest = which.min(MeanMeanDist);
                NearestDist = MeanMeanDist[nearest];
                nearestLabel = as.numeric(levels(useColorsFac)[nearest])
                if ( ((NearestDist < ClusterDiam[nearestLabel]) | (NearestDist <  maxPamDist)) )
                {
                  Colors[InCluster] = nearestLabel;
                  nPAMed = nPAMed + length(InCluster);
                } else Colors[InCluster] = -1;  # This prevents individual points from being assigned later
              }
            }
          }
        }
        # Assign leftover unlabeled objects to clusters with nearest medoids
        Unlabeled = c(1:nPoints)[Colors==0];
        #ColorsX = Colors;
        if (length(Unlabeled)>0) 
        {
          if (pamRespectsDendro)
          {
            unlabOnBranch = Unlabeled[onBranch[Unlabeled] > 0];
            for (obj in unlabOnBranch)
            {
              onBr = onBranch[obj]
              basicOnBranch = branch.basicClusters[[onBr]];
              labelsOnBranch = branchLabels[basicOnBranch]
              useObjects = ColorsX %in% unique(labelsOnBranch);
              useColorsFac = factor(ColorsX[useObjects])
              UnassdToClustDist = tapply(distM[useObjects, obj], useColorsFac, mean);
              nearest = which.min(UnassdToClustDist);
              NearestClusterDist = UnassdToClustDist[nearest];
              nearestLabel = as.numeric(levels(useColorsFac)[nearest])
              if ((NearestClusterDist < ClusterDiam[nearestLabel]) |
                  (NearestClusterDist < maxPamDist) )
              {
                Colors[obj] = nearestLabel
                nPAMed = nPAMed + 1;
              }
            }
          } else {
            useObjects = c(1:nPoints)[ColorsX !=0];
            useColorsFac = factor(ColorsX[useObjects])
            nUseColors = nlevels(useColorsFac);
            UnassdToClustDist = apply(distM[useObjects, Unlabeled, drop = FALSE], 2, tapply, useColorsFac, mean);
            # Fix dimensions for the case when there's only one cluster
            dim(UnassdToClustDist) = c(nUseColors, length(Unlabeled));
            nearest = apply(UnassdToClustDist, 2, which.min);
            nearestDist = apply(UnassdToClustDist, 2, min);
            nearestLabel = as.numeric(levels(useColorsFac)[nearest])
            assign = (nearestDist < ClusterDiam[nearestLabel]) | (nearestDist < maxPamDist)
            Colors[Unlabeled[assign]] = nearestLabel[assign];
            nPAMed = nPAMed + sum(assign);
          }
        }
     }
     if (verbose>2) printFlush(paste(spaces, "....assigned", nPAMed, "objects to existing clusters."));
  }

  # Relabel labels such that 1 corresponds to the largest cluster etc.
  Colors[Colors<0] = 0;
  UnlabeledExist = (sum(Colors==0)>0);
  NumLabs = as.numeric(as.factor(Colors));
  Sizes = table(NumLabs);
  if (UnlabeledExist)
  {
    if (length(Sizes)>1)
    {
        SizeRank = c(1, rank(-Sizes[2:length(Sizes)], ties.method="first")+1);
    } else {
        SizeRank = 1;
    }
    OrdNumLabs = SizeRank[NumLabs];
  } else {
    SizeRank = rank(-Sizes[1:length(Sizes)], ties.method="first");
    OrdNumLabs = SizeRank[NumLabs];
  }
  ordCoreLabels = OrdNumLabs-UnlabeledExist;
  ordCoreLabels[coreLabels==0] = 0; 

  if (verbose>0) printFlush(paste(spaces, "..done."));

  list(labels = OrdNumLabs-UnlabeledExist,
       cores = ordCoreLabels,
       smallLabels = SmallLabels,
       onBranch = onBranch,
       mergeDiagnostics = if (nExternalSplits==0) mergeDiagnostics else 
                              cbind(mergeDiagnostics, externalMergeDiags),
       mergeCriteria = list(maxCoreScatter = maxCoreScatter, minGap = minGap, 
             maxAbsCoreScatter = maxAbsCoreScatter, minAbsGap = minAbsGap, 
             minExternalSplit = minExternalSplit),
       branches  = list(nBranches = nBranches, # Branches = Branches, 
                        IndMergeToBranch = IndMergeToBranch,
                        RootBranch = RootBranch, isCluster = isCluster, 
                        nPoints = nMerge+1));
} 


