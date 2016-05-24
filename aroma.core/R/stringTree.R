#     # Find the longest common suffix of CEL file names.
#     names <- getNames(this);
#     names <- strsplit(names, split="");
#     names <- lapply(names, FUN=rev);
#     names <- lapply(names, FUN=paste, collapse="");
#     names <- unlist(names, use.names=FALSE);
#     suffix <- stringTree(names, maxDepth=1);
#     suffix <- suffix[1];
#     suffix <- names(suffix);
#     if (identical(suffix, "EOS")) {
#       suffix <- "";
#     } else {
#       suffix <- strsplit(suffix, split="")[[1]];
#       suffix <- rev(suffix);
#       suffix <- paste(suffix, collapse="");
#     }
#
setMethodS3("stringTree", "character", function(strs, maxDepth=-1, countRemaining=TRUE, ...) {
  if (length(strs) == 0)
    return(list(EOS=0));

  if (all(nchar(strs) == 0))
    return(list(EOS=length(strs)));

  if (maxDepth == 0) {
    if (countRemaining) {
      return(list(remaining=length(strs)));
    } else {
      return(list(remaining=strs));
    }
  }

  # Extract longest common prefix and remaining "heads"
  prefix <- c();
  pos <- 1;
  while (TRUE) {
    heads <- substring(strs, first=pos, last=pos);
    uHeads <- sort(unique(heads));

    # Number of unique heads
    nbrOfUHeads <- length(uHeads);

    # If some strings are completed, then stop.
    emptyPos <- which(nchar(uHeads) == 0);
    hasCompleteStrings <- (length(emptyPos) > 0);
    if (hasCompleteStrings) {
      uHeads <- uHeads[-emptyPos];
      nbrOfUHeads <- nbrOfUHeads - 1;
      break;
    }

    pos <- pos + 1;

    # String headers differ now?
    if (nbrOfUHeads > 1) {
      break;
    }

    # Update most common prefix
    prefix <- c(prefix, uHeads);
  } # while(TRUE)


  if (hasCompleteStrings) {
    # All strings share same prefix (may be empty)
    prefix <- paste(prefix, collapse="");

    # Are all strings identical?
    if (nbrOfUHeads == 0) {
      children <- list(list(EOS=length(strs)));
      names(children) <- prefix;
      return(children);
    }

    # No, some are completed and some are not.
    isCompleted <- (nchar(heads) == 0);
    nbrOfCompleted <- sum(isCompleted);

    children <- vector("list", nbrOfUHeads+1);
    names(children) <- c("EOS", uHeads);
    children$EOS <- nbrOfCompleted;

    # Remaining strings
    heads <- heads[!isCompleted];
    strs <- strs[!isCompleted];
    strs <- substring(strs, first=pos);

    for (uHead in uHeads) {
      idx <- (heads == uHead);
      value <- stringTree(strs[idx], maxDepth=maxDepth-1);
      if (length(value) == 1)
        value <- value[[1]];
      children[[uHead]] <- value;
    }

    subtree <- list(children);
    names(subtree) <- prefix;

    return(subtree);
  }

  # Here we know that there are no complete strings and that all
  # heads now share the same prefix.
  prefix <- paste(prefix, collapse="");

  # All unique heads
  uCommonHeads <- paste(prefix, uHeads, sep="");

  # Allocate subtree
  children <- vector("list", nbrOfUHeads+1);
  names(children) <- c("EOS", uCommonHeads)
  children$EOS <- as.integer(0);

  # Unique remain string tails
  strs <- substring(strs, first=pos);

  for (kk in seq_along(uHeads)) {
    uHead <- uHeads[kk];
    uCommonHead <- uCommonHeads[kk];
    idx <- (heads == uHead);
    children[[uCommonHead]] <- stringTree(strs[idx], maxDepth=maxDepth-1);
  }

  children;
}, private=TRUE)


############################################################################
# HISTORY:
# 2006-11-30
# o Not working. Broken.
# 2006-??-??
# o Created.
############################################################################
