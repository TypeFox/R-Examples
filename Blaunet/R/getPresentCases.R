getPresentCases <-
function(blauObj, presentCases){
    blauObj$ids <- blauObj$ids[presentCases, , drop=FALSE]
    blauObj$memberships <- blauObj$memberships[presentCases, , drop=FALSE]
    blauObj$dimensions <- blauObj$dimensions[presentCases, , drop=FALSE]
    blauObj$weights <- blauObj$weights[presentCases, , drop=FALSE]

    if (!is.null(blauObj$isInNiche)) {
      blauObj$isInNiche <- blauObj$isInNiche[presentCases, , drop=FALSE]
    }
    
    if (!is.null(blauObj$primaryMembership)) {
      blauObj$primaryMembership <- blauObj$primaryMembership[presentCases, , drop=FALSE]
    }

    return(blauObj)

    #we don't cut down connections, we keep the full network object
    #we don't cut down top/lowbounds because they're computed in a manner specified with the input function and aren't the same size
    #we don't cut down 'final' elements such as nodalLocal because they've been computed already according to user specification
}
