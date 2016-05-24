##' A method of the generic function by for objects of class \'trackdata\'
##' 
##' A given function 'FUN' is applied to the data corresponding to each segment
##' of data.
##' 
##' It is the same as trapply but with the extension to subsume calculation to
##' groups of segments. Note, if you do not want to apply the function fun to a
##' special group of segments, use \link{trapply} instead.
##' 
##' @aliases by.trackdata by
##' @param data a track data object
##' @param INDICES a list of segment indices, like a label vector
##' @param FUN a function that is applied to each segment
##' @param \dots arguments of the function fun
##' @param simplify simplify = TRUE , output is a matrix; simplify = FALSE a
##' list is returned
##' @return list or vector
##' @author Jonathan Harrington
##' @seealso \code{\link{trapply}}, \code{\link{by}}, \code{\link{trackdata}}
##' \code{\link{dapply}} \code{\link{smooth}} \code{\link{apply}}
##' @keywords methods
##' @examples
##' 
##'   data(demo.vowels)
##'   data(demo.vowels.fm)
##' 
##' 
##'    #mean F1 subsumed for each vowel
##'    ################################
##'    lab = label(demo.vowels)
##'    by(demo.vowels.fm[,1], lab ,sapply,mean,simplify=FALSE)
##' 
##' 
##'    #mean F1 subsumed for segment onsets mids and offsets
##'    ##############################################
##'    data = demo.vowels.fm
##'    llabs = NULL
##'    for (ind in 1:dim(data$ftime)[1]) {
##'      seglabs = rep("mid",data$index[ind,2]-data$index[ind,1]+1)
##'      seglabs[1] = "on"
##'      seglabs[length(seglabs)] = "off"
##'      llabs = as.vector(c(llabs , seglabs))
##'    }
##' 
##'    by(demo.vowels.fm[,1], llabs , sapply, mean , simplify=FALSE)
##' 
##'    #mean F1 subsumed for segment onsets mids and offsets subsumed for each vowel
##'    #####################################################################
##'    by(demo.vowels.fm[,1], list(lab = lab, llabs = llabs) , sapply, mean , simplify=FALSE)
##' 
##' 
##' 
##' @export
`by.trackdata` <- function (data, INDICES = NULL, FUN, ..., simplify = FALSE) 
{
  #there might be a problem with data$data therefore data is replaced by abitrary datam
  # data=m.fdat.int[,1]; FUN=targtime; simplify=TRUE;INDICES = NULL
  
  datam = data
  
  orgindices = INDICES
  fun = FUN
  arr.indices = function(orgindices, datam) {
    indices = NULL
    if (is.null(orgindices) || is.vector(orgindices)) {
      n = 1:nrow(datam$index)
      if (is.null(orgindices)) {
        indices = rep(n, datam$index[, 2] - datam$index[, 
                                                        1] + 1)
      }
      if (length(orgindices) == dim(datam$ftime)[1]) {
        indices = NULL
        for (ind in 1:dim(datam$ftime)[1]) {
          indices = c(indices, rep(orgindices[ind], datam$index[ind, 
                                                                2] - datam$index[ind, 1] + 1))
        }
      }
      if (length(orgindices) == dim(datam$data)[1]) {
        indices = orgindices
      }
      return(indices)
    } else {
      warning("Can not arrange INDICES!")
      return(orgindices)
    }
  }
  indices = NULL
  if (is.list(orgindices)) {
    for (var in 1:length(names(orgindices))) {
      orgindices[[var]] = arr.indices(orgindices[[var]],datam)
    }
    indices = orgindices
  } else {
    indices = arr.indices(orgindices,datam)
  }
  result <- o <- by(I(datam$data), indices, fun, ...)
  if (simplify) {
    if (is.null(attributes(summary(o))$dim)) 
      result <- c(unlist(o))
    else {
      result <- NULL
      for (j in 1:length(o)) {
        result <- rbind(result, o[[j]])
      }
    }
  }
  result
}
