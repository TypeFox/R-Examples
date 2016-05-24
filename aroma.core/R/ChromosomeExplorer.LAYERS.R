setMethodS3("writeAxesLayers", "ChromosomeExplorer", function(this, ...) {
  model <- getModel(this);

  path <- getPath(this);
  path <- filePath(getParent(path), "axes,chrLayer");
  path <- Arguments$getWritablePath(path);
  plotAxesLayers(model, path=path, imageFormat="png", transparent=TRUE, ...);

  invisible(path);
}, private=TRUE)


setMethodS3("writeGridHorizontalLayers", "ChromosomeExplorer", function(this, ...) {
  model <- getModel(this);

  path <- getPath(this);
  path <- filePath(getParent(path), "gridH,chrLayer");
  path <- Arguments$getWritablePath(path);
  plotGridHorizontalLayers(model, path=path, imageFormat="png", transparent=TRUE, ...);

  invisible(path);
}, private=TRUE)



setMethodS3("writeCytobandLayers", "ChromosomeExplorer", function(this, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get the model
  model <- getModel(this);

  path <- getPath(this);
  path <- filePath(getParent(path), "cytoband,chrLayer");
  path <- Arguments$getWritablePath(path);
  plotCytobandLayers(model, path=path, imageFormat="png", transparent=TRUE, ...);

  invisible(path);
}, private=TRUE)


setMethodS3("getSampleLayerName", "Explorer", function(this, name, class="sampleLayer", ...) {
  layer <- c(getSampleLayerPrefix(this), name, class);
  layer <- paste(layer, collapse=",");
  layer;
}, private=TRUE)


setMethodS3("writeRawCopyNumberLayers", "ChromosomeExplorer", function(this, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get the model
  model <- getModel(this);

  path <- getPath(this);
  layer <- getSampleLayerName(this, "rawCNs");
  path <- filePath(getParent(path), layer);

  path <- Arguments$getWritablePath(path);
  plotRawCopyNumbers(model, path=path, imageFormat="png", transparent=TRUE, ...);

  invisible(path);
}, private=TRUE)


setMethodS3("writeCopyNumberRegionLayers", "ChromosomeExplorer", function(this, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get the model
  model <- getModel(this);

  # Not supported?
  if (!inherits(model, "CopyNumberSegmentationModel"))
    return(NULL);

  path <- getPath(this);
  tags <- getAsteriskTags(model, collapse=",");
  path <- filePath(getParent(path), sprintf("%s,sampleLayer", tags));

  path <- Arguments$getWritablePath(path);
  plotCopyNumberRegionLayers(model, path=path, imageFormat="png", transparent=TRUE, ...);

  invisible(path);
}, private=TRUE)



##############################################################################
# HISTORY:
# 2009-12-30
# o CLEAN UP: Now ChromosomeExplorer is unaware of the data set tuple in
#   the model.  Instead all queries goes to the model.
# o CLEAN UP: Removed getSetTuple() from ChromosomeExplorer.
# 2009-12-25
# o Added getSampleLabels() to ChromosomeExplorer.
# o Now getNames() of ChromosomeExplorer is based on the getNames()
#   output, and no longer on getArrays().
# o BUG FIX: translateFullNames() of ChromosomeExplorer would translate the
#   full names, but return the original ones.
# 2009-12-22
# o BUG FIX: getFullNames() of ChromosomeExplorer reported: Error in
#   UseMethod("translateFullNames") : no applicable method for
#   'translateFullNames' applied to an object of class "character".
# 2009-05-19
# o Now testing for file permissions for creat-/writ-/updating files/dirs.
# 2008-12-13
# o Added argument 'zooms' to the constructor of ChromosomeExplorer.
#   Added methods get- and setZooms().
# 2008-06-05
# o Made updateSamplesFile(), writeAxesLayers(), writeGridHorizontalLayers(),
#   writeCytobandLayers(), writeRegions(), setup(), process() parallel safe.
# 2008-05-08
# o Now one can pass argument 'aliased=TRUE' to process() which causes the
#   ChromosomeExplorer and coupled CopyNumberSegmentationModel to return
#   tags that inferred from aliased full names.
# o Now updateSamplesFile() of ChromosomeExplorer passes arguments '...'
#   to getFullNames().
# o Now getFullNames() of ChromosomeExplorer passes arguments '...' along.
# 2007-10-17
# o Now the  accepts CopyNumberChromosomalModel:s.ChromosomeExplorer
# o Added support for specifying the output version (at least for now).
# 2007-10-10
# o Added support for layers.  This should be backward compatible.
# 2007-10-09
# o Now process() writes cytoband/*.png layer images.
# o Added writeCytobandLayers() etc.
# 2007-09-30
# o Now the main HTML file for ChromosomeExplorer is ChromosomeExplorer.html,
#   which is in analogue to how the ArrayExplorer works.  This HTML now
#   loads to a different Javascript file with a different name so that already
#   existing index.html ChromosomeExplorer files will still work.
# 2007-09-04
# o Now ChromosomeExplorer recognizes CopyNumberSegmentationModel:s.
# 2007-05-08
# o Added default zoom levels to updateSamplesFile() for ChromosomeExplorer.
#   This is applies the first time process() is called.
# 2007-03-19
# o Now ChromosomeExplorer extends Explorer.
# 2007-02-06
# o Now templates are in reports/templates/ and includes in reports/includes/.
# o Updated the path to <rootPath>/<dataSetName>/<tags>/<chipType>/<set>/.
# 2007-01-17
# o Now all 'arrays' arguments can contain array names.
# o Added getArrays() and setArrays() in order to focus on a subset of the
#   arrays in the model.
# 2007-01-16
# o Added getAlias() and setAlias(), and updated getName() accordingly.
#   This makes it easy to change the name of output set for subsets of
#   arrays.
# 2007-01-15
# o Added some more Rdoc comments.
# 2007-01-10
# o BUG FIX: setup() would only add index.html if includes/ were missing.
# 2007-01-08
# o Added display().
# 2007-01-07
# o TODO: The region filter for writeRegions() is hardwired for now.
# o Update to include regions.xls file on the CE output page too.
# o Now all HTML, CSS, and Javascript files are created too.
# 2007-01-04
# o Created.
##############################################################################
