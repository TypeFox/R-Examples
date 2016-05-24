#' Create bins from a set of binned distances and a set of cutpoints.
#'
#' This is an internal routine and shouldn't be necessary in normal analyses.
#'
#' @param data \code{data.frame} with at least the column \code{distance}.
#' @param cutpoints vector of cutpoints for the bins
#'
#' @return data \code{data} with two extra columns \code{distbegin} and
#'        \code{distend}.
#'
#' @author David L. Miller
#' @export
create.bins <- function(data,cutpoints){

  # don't do anything if there are NAs in the distance column
  if(any(is.na(data$distance))){
    stop("Some distances are NA, can't create bins.")
  }

  # lazy typist
  cp <- cutpoints

  # remove distances outside bins
  in.cp.ind <- data$distance>=cp[1] & data$distance<=cp[length(cp)]
  if(!all(in.cp.ind)){
    warning("Some distances were outside bins and have been removed.")
  }
  data <- data[in.cp.ind,]

  # pull out the distances (removing the NAs for now)
  na.ind <- is.na(data$distance)
  d <- data$distance[!na.ind]

  distbegin<-rep(NA,length(d))
  distend<-rep(NA,length(d))

  for(i in 1:(length(cp)-1)){
    # which elements of d lie between cutpoints i and i+1
    ind <- which(d>=cp[i] & d<cp[i+1])

    distbegin[ind] <- cp[i]
    distend[ind]   <- cp[i+1]
  }
  # last cutpoint, include those observations AT the truncation point
  ind <- which(d>=cp[i] & d<=cp[i+1])

  distbegin[ind] <- cp[i]
  distend[ind]   <- cp[i+1]


  # handle NA distances, that we need to preserve
  distbegin.na <- rep(NA,length(data$distance))
  distend.na <- rep(NA,length(data$distance))
  distbegin.na[!na.ind] <- distbegin
  distend.na[!na.ind] <- distend

  # put all that together and make a data.frame
  data <- cbind(data,
                distbegin=distbegin.na,
                distend=distend.na)
  data <- data.frame(data)

  return(data)
}
