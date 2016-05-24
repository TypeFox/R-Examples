#' Create a model frame for ddf fitting
#'
#' Creates a model.frame for distance detection function fitting. It includes
#' some pre-specified and computed variables with those included in the model
#' specified by user (formula)
#'
#' The following fields are always included: detected, observer, binned, and
#' optionally distance (unless null), timesdetected (if present in data). If
#' the distance data were binned, include distbegin and distend point fields.
#' If the integration width varies also include int.begin and int.end and
#' include an offset field for an iterative glm, if used.  Beyond these fields
#' only fields used in the model formula are included.
#'
#' @param xmat dataframe for ddf
#' @param scale.formula user specified formula for scale of distance detection
#'   function
#' @param meta.data user-specified meta.data (see \code{\link{ddf}}
#' @param shape.formula user specified formula for shape parameter of distance
#'   detection function
#' @return model frame for analysis
#' @note Internal function and not called by user
#' @author Jeff Laake
#' @keywords utility
create.model.frame <- function(xmat,scale.formula,meta.data,shape.formula=NULL){
  # Create data frame with variables used in the formula and object #
  # (Object id); this code fix: jll 18-Nov-04; a different approach was used
  # to allow for use of as.factor(x) etc in formula
  if(!is.null(shape.formula)){
    varlist <- unique(c(all.vars(scale.formula),all.vars(shape.formula)))
  }else{
    varlist <- all.vars(scale.formula)
  }

  if(any(! varlist %in% names(xmat))){
     stop("The following variables used in the formula are not in the data: ",
          varlist[!varlist%in%names(xmat)])
  }

  data <- cbind(object = xmat$object,
                xmat[,varlist,drop=FALSE])
  colnames(data) <- c("object",varlist)
  data <- as.data.frame(data)

  # fix: ljt 21-Sep-05; row.names of xmat need to be copied over to data,
  #     otherwise some later stuff in e.g., ddf.ds doesn't work
  row.names(data) <- row.names(xmat)

  # Include other fields such as detected, observer, binned, and optionally
  # distance (unless null)
  if(is.null(data$distance)){
    data <- cbind(data,
                  distance=xmat$distance,
                  binned=xmat$binned)
  }else{
    data <- cbind(data,
                  detected=xmat$detected,
                  binned=xmat$binned,
                  observer=xmat$observer)
  }
  if(!is.null(xmat$observer)){
    data$observer <- xmat$observer
  }

  if(!is.null(xmat$detected)){
    data$detected <- xmat$detected
  }

  # Also include # times the object was detected if it is in data
  if(!is.null(xmat$timesdetected)){
    data <- cbind(data,
                  timesdetected=xmat$timesdetected)
  }

  # If the distance data were binned, include begin and end point fields
  if(meta.data$binned){
    data <- cbind(data,
                  distbegin=xmat$distbegin,
                  distend=xmat$distend)
  }

  # If the integration width varies (int.begin/int.end are in data frame)
  # also include those fields
  if(!is.null(xmat$int.begin)&!is.null(xmat$int.end)){
    data <- cbind(data,
                  int.begin=xmat$int.begin,
                  int.end=xmat$int.end)
  }

  # If data frame contains an offset for glm also include it.
  if(!is.null(xmat$offsetvalue)){
    data <- cbind(data,
                  offsetvalue=xmat$offsetvalue)
  }

  # return completed dataframe
  return(data)
}
