#' Process data for fitting distance sampling detection function
#'
#' Sets up dataframe and does some basic error checking. Adds needed fields to
#' dataframe and to \code{meta.data}.
#'
#' The function does a number of error checking tasks, creating fields and
#' adding to \code{meta.data} including:
#'
#' 1) If \code{check=TRUE}, check to make sure the record structure is okay for
#' mrds data. The number of primary records (observer=1) must equal the number
#' of secondary records (observer=2). Also, a field in the dataframe is created
#' \code{timesseen} which counts the number of times an object was detected
#' 0,1,2; if \code{timesseen=0} then the record is tossed from the analysis.
#' Also if there are differences in the data (distance, size, covariates) for
#' observer 1 and 2 a warning is issued that the analysis may fail.  The code
#' assumes these values are the same for both observers.
#'
#' 2) Based on the presence of fields \code{distbegin} and \code{distend}, a
#' determination is made of whether the data analysis should be based on binned
#' distances and a field \code{binned} is created, which is \code{TRUE} if the
#' distance for the observation is binned.  By assigning for each observation
#' this allows an analysis of a mixture of binned and unbinned distances.
#'
#' 4) Data are restricted such that distances are not greater than \code{width}
#' and not less than \code{left} if those values are specified in
#' \code{meta.data}.  If they are not specified then \code{left} defaults to 0
#' and \code{width} defaults to the largest distance measurement.
#'
#' 5) Determine if an integration range (\code{int.begin} and \code{int.end}
#' has been specified for the observations.  If it has, add the structure to
#' \code{meta.data}.  The integration range is typically used for aerial
#' surveys in which the altitude varies such that the strip width (left to
#' width) changes with a change in altitude.
#'
#' 6) Fields defined as factors are cleaned up such that any unused levels are
#' eliminated.
#'
#' 7) If the restrictions placed on the data, eliminated all of the data, the
#' function stops with an error message
#'
#' @param data dataframe object
#' @param meta.data meta.data options; see \code{\link{ddf}} for a description
#' @param check if \code{TRUE} check data for errors in the mrds structure; for
#'   \code{method="ds" check=FALSE}
#' @return \item{xmat}{processed \code{data.frame} with added fields}
#'   \item{meta.data}{meta.data list}
#' @author Jeff Laake
#' @keywords utility
process.data <- function(data,meta.data=list(),check=TRUE){

  set.default.width=function(data,meta.data){
  # set.default.width - sets default transect width when none was specified
  #  Arguments:
  #  data      - dataframe
  #  meta.data - meta.data list
  # Values:  width of transect
    if(meta.data$binned){
      width <- max(c(data$distend,data$distance),na.rm=TRUE)
    }else{
      width <- max(data$distance)
    }
    return(width)
  }

  # assign dataframe to data

  # Check to make sure the record structure is ok. Number of primary
  # records = number of secondary
  if(check){
    if(length(data$detected[data$observer==1]) !=
        length(data$detected[data$observer==2])){
      stop("number of records for primary observer not equal to number for secondary observer")
    }
  }

  # Create field which counts the number of times an object was detected 0,1,2
  if(check){
    timesdetected <- data$detected[data$observer==1] +
                     data$detected[data$observer==2]
    data$timesdetected <- rep(0,dim(data)[1])
    data$timesdetected[data$observer==1] <- timesdetected
    data$timesdetected[data$observer==2] <- timesdetected

    # If any 00 (not detected by either observer), stop and issue error message
    if(any(data$timesdetected==0)){
      stop("following objects were never detected:",
            paste(data$object[data$observer==1&data$timesdetected==0],
                  collapse=","),"\n")
    }
  }

  # also check for mrds that the data fields have the same value for both 
  # observers for example same distance, size etc.  This is only a warning as 
  #some fields may be validly different
  #  if(any(apply(data[data$observer==1,names(data)!="observer"&names(data)!="detected"],1,paste,collapse="")!=
  #    apply(data[data$observer==2,names(data)!="observer"&names(data)!="detected"],1,paste,collapse="")))
  #    warning("If analysis fails it may be due to difference in data between observer 1 and 2;\n fields such as distance, size and covariates should be the same")

  # Determine if data are binned by presence of distbegin and distend fields
  if(is.null(data$distend)|is.null(data$distbegin)){
    binned <- FALSE
  }else{
    if(all(is.null(data$distend))|all(is.null(data$distbegin))){
      binned <- FALSE
    }else{
      if(any(is.null(data$distend) & !is.null(data$distbegin)) |
         any(is.null(data$distbegin)&!is.null(data$distend))){
        stop("mismatched distance intervals - one or more endpoints are missing")
      }else{
        binned <- TRUE
      }
    }
  }

  if(meta.data$binned & !binned){
    stop("binned set to TRUE in meta.data but distbegin and distend fields are missing")
  }

  if(!meta.data$binned & binned){
    warning("data contain distbegin and distend fields but binned=FALSE. Analyzing as not binned",immediate.=TRUE)
    binned <- FALSE
  }

  meta.data$binned <- binned

  if(meta.data$binned & is.null(meta.data$breaks)){
    stop("breaks must be set in meta.data for binned data")
  }

  # Fill in distance field for binned observations and create logical variable
  data$binned <- rep(FALSE,dim(data)[1])
  if(binned){
    meta.data$binned <- TRUE
    data$distance[!is.na(data$distend)]<-(data$distbegin[!is.na(data$distend)]+
                                          data$distend[!is.na(data$distend)])/2
    data$binned[!is.na(data$distbegin)] <- TRUE
  }

  # Restrict data to width interval 
  # If no width set, use largest measured distance as width
  if(is.na(meta.data$width)){
    width <- set.default.width(data,meta.data)
    meta.data$width <- width
    xmat <- data
    warning("no truncation distance specified; using largest observed distance",immediate.=TRUE)
  }else{
    # change: jll 2 June 05; ref to width changed to meta.data$width
    # This piece of code makes sure that the set width is as large as the
    # largest bin end point for binned data.
    if(meta.data$binned){
      if(any(data$binned & data$distend > meta.data$width)){
        stop("width must exceed largest interval end point")
      }else{
        xmat <- data[data$binned |
                     (!data$binned&data$distance<=meta.data$width),]
      }
    }else{
      xmat <- data[data$distance <= meta.data$width,]
    }
  }

  # Determine if integration range has been specified
  if(is.null(xmat$int.begin)|is.null(xmat$int.end)){
    if(any(is.na(meta.data$int.range))){
      meta.data$int.range <- c(meta.data$left,meta.data$width)
    }
  }else{
      meta.data$int.range <- rbind(c(meta.data$left,meta.data$width),
                                   cbind(xmat$int.begin,xmat$int.end))
  }

  # If left >0 perform left truncation by restricting values
  if(meta.data$left >0){
    if(binned){
      if(any(data$binned&data$distbegin < meta.data$left)){
        stop("left truncation must be smaller than the smallest interval begin point")
      }else{
        xmat <- data[data$binned|(!data$binned&data$distance>=meta.data$left),]
      }
    }else{
      xmat <- xmat[xmat$distance>=meta.data$left,]
    }
  }

  # Clean up factor levels
  b <- dim(xmat)[2]
  for(i in 1:b){
    if(is.factor(xmat[,i])){
      xmat[,i] <- factor (xmat[,i])
    }
  }

  # If the exclusion eliminated all of the data, stop with error message
  if(dim(xmat)[1]==0){
    stop("no data to analyze")
  }

  return(list(xmat=xmat,meta.data=meta.data))
}
