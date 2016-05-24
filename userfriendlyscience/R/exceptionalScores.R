exceptionalScores <- function(dat, items=NULL,
                              exception=.025, totalOnly=TRUE, append=TRUE,
                              both=TRUE, silent=FALSE, suffix = "_isExceptional",
                              totalVarName = "exceptionalScores") {
  
  if (is.data.frame(dat)) {
    if (is.null(items)) {
      items <- names(dat);
      if (!silent) {
        cat("No items specified: extracting all variable names in dataframe.\n");
      }
    }
    exceptionalScores <- dat[, items];
  } else {
    ### Vector provided; store in dataframe.
    exceptionalScores <- data.frame(dat);
    names(exceptionalScores) <- deparse(substitute(dat));
  }
  
  originalCols <- ncol(exceptionalScores);
  exceptionalScores <- data.frame(exceptionalScores[, unlist(lapply(exceptionalScores, is.numeric))]);
  if ((originalCols > ncol(exceptionalScores) & !silent)) {
    cat0("Note: ", originalCols - ncol(exceptionalScores), " variables ",
         "were not numeric and will not be checked for exceptional values.\n");
  }
  
  namesToUse <- paste0(colnames(exceptionalScores), suffix);
  
  exceptionalScores <- apply(exceptionalScores, 2,
                             exceptionalScore, prob = exception, both=both, silent=silent);
  
  colnames(exceptionalScores) <- namesToUse;
  
  if (totalOnly) {
    totalTrues <- rowSums(exceptionalScores, na.rm=TRUE);
    if (append) {
      dat[, totalVarName] <- totalTrues;
      return(dat);
    } else {
      return(totalTrues);
    }    
  } else {
    if (append) {
      return(data.frame(dat,
                        exceptionalScores));
    } else {
      return(exceptionalScores);
    }
  }
  
}
