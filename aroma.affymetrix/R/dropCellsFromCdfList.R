# @author "HB"
dropCellsFromCdfList <- function(cdfList, maxNbrOfCells, cellsToExclude, ..., verbose=0) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'cdfList':
  if (!is.list(cdfList)) {
    stop("Argument 'cdfList' is not a list: ", mode(cdfList));
  }

  # Argument 'maxNbrOfCells':
  maxNbrOfCells <- as.integer(maxNbrOfCells);
  if (maxNbrOfCells < 1) {
    stop("Argument 'maxNbrOfCells' is non-positive: ", maxNbrOfCells);
  }

  # Argument 'cellsToExclude':
  cellsToExclude <- as.integer(cellsToExclude);
  if (any(cellsToExclude < 1 | cellsToExclude > maxNbrOfCells)) {
    stop(sprintf("Argument 'cellsToExclude' is out of range [1,%d].", maxNbrOfCells));
  }


  nbrOfUnits <- length(cdfList);

  if(verbose >= 1) {
    cat("Excluding cells from CDF list structure...\n");
    cat("Number of units: ", nbrOfUnits, "\n", sep="");
    cat("Number of cells to exclude: ", length(cellsToExclude), "\n", sep="");
  }

  toExclude <- logical(maxNbrOfCells);
  toExclude[cellsToExclude] <- TRUE;

  # Sanity check (cannot drop all cells)
  stopifnot(!all(toExclude));

  # Nothing todo?
  if (!any(toExclude)) {
    if(verbose >= 1) {
      cat("No cells to be exclude. Skipping.\n");
      cat("Excluding cells from CDF list structure...done\n");
    }
    return(cdfList);
  }

  if(verbose >= 1) {
    cat("Number of units left: ");
  }

  # Number of cells dropped from CDF
  nDrop <- 0L;

  # Total number of cells
  nCells <- 0L;

  fields <- c("x", "y", "indices", "pbase", "tbase", "atom", "indexpos");
  for (jj in seq_len(nbrOfUnits)) {
    if(verbose >= 1 && jj %% 10e3 == 1) {
      cat(nbrOfUnits-jj+1L, ", ", sep="");
    }

    cdfUnit <- cdfList[[jj]];

    # Total number of cells in unit
    nbrOfAtoms <- 0L;
    nbrOfCells <- 0L;

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Update the groups
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    groups <- cdfUnit[["groups"]];

    # Number of cells dropped for this unit
    nDropUU <- 0L;
    nCellsUU <- 0L;

    nbrOfGroups <- length(groups);
    for (gg in seq_len(nbrOfGroups)) {
      group <- groups[[gg]];

      indices <- group$indices;
      nCellsGG <- length(indices);
      nCellsUU <- nCellsUU + nCellsGG;
      toExcludeGG <- toExclude[indices];
      drop <- which(toExcludeGG);
##print(list(indices=indices, drop=drop, toExcludeGG=toExcludeGG));
      nDropGG <- length(drop);

      # Sanity check (assuming not all cells are dropped)
      ## stopifnot(nDropGG < nCellsGG);

      nbrOfCellsGG <- nCellsGG - nDropGG;


      if (nDropGG > 0) {
        keep <- which(!toExcludeGG);
##print(list(unit=jj, group=gg, keep=keep));
        for (field in fields) {
          value <- group[[field]];
          value <- value[keep];
          group[[field]] <- value;
        }

        # The 'atom' and 'indexpos' fields must be "complete" (I think)
        group$atom <- seq_len(nbrOfCellsGG);
        group$indexpos <- seq_len(nbrOfCellsGG);
#print(group);

        nDropUU <- nDropUU + nDropGG;
      }

      # Hmm... (is the following true)
      nbrOfAtomsGG <- nbrOfCellsGG;

      # Update the group count
      group$natoms <- nbrOfAtomsGG;
      group$ncells <- nbrOfCellsGG;

      nbrOfAtoms <- nbrOfAtoms + nbrOfAtomsGG;
      nbrOfCells <- nbrOfCells + nbrOfCellsGG;

      groups[[gg]] <- group;
    } # for (gg ...)

    cdfUnit[["groups"]] <- groups;

    nDrop <- nDrop + nDropUU;
    nCells <- nCells + nCellsUU;

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Update the unit counts
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Total number of cells in unit
    cdfUnit$natoms <- nbrOfAtoms;
    cdfUnit$ncells <- nbrOfCells;

    cdfList[[jj]] <- cdfUnit;
  } # for (jj ...)

  if(verbose >= 1) {
    cat("0.\n");
    cat("Total number of cells: ", nCells, "\n");
    cat("Number of cells dropped: ", nDrop, "\n");
    cat("Excluding cells from CDF list structure...done\n");
  }

  cdfList;
} # dropCellsFromCdfList()


############################################################################
# HISTORY:
# 2011-09-17
# o Note that dropCellsFromCdfList() may create empty unit groups.
# o The dropCellsFromCdfList() method was designed such that it can
#   be moved to the affxparser package as-is.
# o Note that the method does not assume that a cell belong to unique unit.
# o Now also adjusting the counts.
# 2011-09-15
# o Added dropCellsFromCdfList() for excluding cells from a CDF list
#   structure by their cell indices.
############################################################################
