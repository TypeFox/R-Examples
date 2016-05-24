calcLRTs <- function(testPointList, cleanResu="") {
  if( ! is.null(testPointList) & length(testPointList)>0) {
    message.redef("Computing likelihood ratio tests...")
    LRTlist <- list()
    for(testPoint in testPointList) { ## for each successive sublist of testPointList
      ##Nb can be a correct testPoint variable absent from fittedNames ?
      notValidParams <- (names(testPoint) %w/o% c(blackbox.getOption("ParameterNames"), "latt2Ns2", "Nratio", "NactNfounderratio", "NfounderNancratio", "condS2")) ## names(LRTfixedvals) must be standardNames
      LRTfixedvals <- testPoint[(names(testPoint) %w/o% notValidParams)] ## removeparameters not considered in the analyses but that user would have input in TestPoint setting
      if (length(notValidParams)>0) {
        mess <- paste("   this parameter (", notValidParams, ") and its associated value wil be ignored.")
        message.redef("Some parameters in a TestPoint is not a parameter of the model at hand, ")
        message.redef(mess)
      }
      toBeRemoved <- intersect(names(testPoint), blackbox.getOption("constantNames")) ## names(LRTfixedvals) must be standardNames
      LRTfixedvals <- testPoint[(names(LRTfixedvals) %w/o% toBeRemoved)] ## removes (canonical) values fixed in kriging but that user would have input in TestPoint setting
      if (length(LRTfixedvals)==0) {
        message.redef("A testPoint has been ignored")
        message.redef("   because none of its component parameters has been fitted.")
        ##next
      } else {
        usernames <- sapply(names(LRTfixedvals), userunit, format="ASCII") ## may be nicer than those actually entered by user
        readablename <- paste("LRT", usernames, LRTfixedvals, sep="_", collapse="_") ## before log transfo !;
        for (st in names(LRTfixedvals)) if(islogscale(st)) {LRTfixedvals[st] <- log(as.numeric(LRTfixedvals[st]))}
        tmp <- LRTfn(LRTfixedvals, cleanResu=cleanResu)
        if (!is.na(tmp$LRT)) { ## storing only successful computations for later use
          tmplist <- list(tmp) ## a list with a list as single element
          names(tmplist) <- readablename ## giving a name, say 'LRTtwoNm_2', to the single element
          ## stacking sublists. Usage: '...LRTlist$LRTtwoNm$proftp$par', ...LRTlist[[1]]$LRT
          LRTlist[names(tmplist)] <- tmplist
        }
      }
    }
    .blackbox.data$options$LRTlist[names(LRTlist)] <- LRTlist
  } ##
  return(blackbox.getOption("LRTlist"))
}
