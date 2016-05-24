library("aroma.affymetrix");
verbose <- Arguments$getVerbose(-8, timestamp=TRUE);

cdf <- AffymetrixCdfFile$byChipType("Mapping50K_Hind240");
print(cdf);

unitTypes <- getUnitTypes(cdf);
print(table(unitTypes));

units <- which(unitTypes == 2);
str(units);

# Group units by the number of groups and cells per group
sets <- groupUnitsByDimension(cdf, units=units, verbose=verbose);
str(sets);

# Apply cell indices to sets containing "rectangular" unit
cellSets <- lapply(sets$nestedSets, FUN=function(set) {
  nbrOfGroups <- set$nbrOfGroups;
  set$sets <- lapply(set$sets, FUN=function(subset) {
    nbrOfCellsPerGroup <- subset$nbrOfCellsPerGroup;
    nbrOfCellsPerGroup <- unique(nbrOfCellsPerGroup);
    if (length(nbrOfCellsPerGroup) == 1) {
      nbrOfUnits <- length(subset$units);
      dim <- c(nbrOfCellsPerGroup, nbrOfGroups, nbrOfUnits);
      cells <- getCellIndices(cdf, units=subset$units, 
                                   unlist=TRUE, useNames=FALSE);
      dim(cells) <- dim;
      subset$cells <- cells;
    }
    subset;
  });
  set;
});

str(cellSets);
