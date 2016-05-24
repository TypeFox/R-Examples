`DiscretizeImpreciseModel` <-
function(ListOfFunctions,UpperPrevisions,x,omega,s){

  RangeEpsilon <- 0.0005      # This is the RangeEpsilon, which corresponds to the accuracy of the discrete model

  # The values of the functions at the additional nodes omega:
    omega.nodevalues <- fevaluation(omega,ListOfFunctions, useApply = FALSE)
    omega.nodevalues <-  omega.nodevalues[,!duplicated(omega.nodevalues,MARGIN=2)] # a simple but very
                                                                                   # efficient reduction of
                                                                                   # the number of nodes
  # The values of the functions at the observations:
    x.nodevalues <- fevaluation(x,ListOfFunctions, useApply = FALSE)

  # All present nodevalues as a matrix:
    nodevalues <- cbind(x.nodevalues,omega.nodevalues)

  # Standardizing the functions f
    f.max <- apply(nodevalues,1,max)
    f.min <- apply(nodevalues,1,min)
    x.nodevalues <- (x.nodevalues-f.min)/(f.max-f.min)
    omega.nodevalues <- (omega.nodevalues-f.min)/(f.max-f.min)

  # Standardizing the upper previsions
    UpperPrevisions <- (UpperPrevisions-f.min)/(f.max-f.min)+RangeEpsilon

  # Discretizing the range of the functions f
    x.nodevalues <- ceiling(x.nodevalues/RangeEpsilon)*RangeEpsilon
    omega.nodevalues <- ceiling(omega.nodevalues/RangeEpsilon)*RangeEpsilon

    DiscretePrevision <- list(x.nodevalues,omega.nodevalues,UpperPrevisions)

  DiscretePrevision
}

