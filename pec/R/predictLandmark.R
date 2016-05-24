# methods for dynamic 
# --------------------------------------------------------------------

predictLandmark <- function(object,newdata,times,landmark,cause,...){
  UseMethod("predictLandmark",object)
}

predictLandmark.jointPenal <- function(object,newdata,times,landmark,cause,...){
    ## xxx
    ## stopifnot(
}
