#' create an object of class 'dataObj' required for further processing
#' 
#' create an standardized input object of class 'dataObj' containing
#' information on weights, household ids, household sizes, person ids and
#' optionally strata. Outputs of this function are typically used in
#' \code{\link{simStructure}}.
#' 
#' @name specifyInput
#' @param data a \code{data.frame} or \code{data.table} featuring sample data.
#' @param hhid character vector of length 1 specifying variable containing
#' household ids within slot \code{data}.
#' @param hhsize character vector of length 1 specifying variable containing
#' household sizes within slot \code{data}. If NULL, household sizes are
#' automatically calculated.
#' @param pid character vector of length 1 specifying variable containing
#' person ids within slot \code{data}. If NULL, person ids are automatically
#' calculated.
#' @param weight character vector of length 1 specifying variable holding
#' sampling weights within slot \code{data}.
#' @param strata character vector of length 1 specifing variable name within
#' slot \code{data} of variable holding information on strata, e.g. regions or
#' NULL if such variable does not exist.
#' @param population TRUE/FALSE vector of length 1 specifing if the data object is a sample or a population
#' NULL if such variable does not exist.
#' @export
#' @author Bernhard Meindl
#' @keywords method
#' @examples
#' data(eusilcS)
#' inp <- specifyInput(data=eusilcS, hhid="db030", weight="rb050", strata="db040")
#' class(inp)
#' inp
specifyInput <- function(data, hhid, hhsize=NULL, pid=NULL, weight=NULL, strata=NULL, population=FALSE) {
  if ( !class(hhid)=="character" | length(hhid) != 1 | is.na(match(hhid, colnames(data)))) {
    stop("hhid must be a character defining the variable holding household ids and must be of length 1!\n")
  }
  ## check if weight var is defined well for sample data
  if(!population){
    if ( !class(weight)=="character" | length(weight) != 1 | is.na(match(weight, colnames(data)))) {
      stop("weight must be a character defining the variable holding sampling weights and must be of length 1!\n")
    }
  }else{ # for population data no weight car is accepted
    if(!is.null(weight))
      stop("weight must not be defined for a population")
  }
  if ( !is.null(strata) ) {
    if ( !class(strata)=="character" | length(strata) != 1 | is.na(match(strata, colnames(data)))) {
      stop("strata must be a character defining the variable holding information on stratas and must be of length 1!\n")
    }
  }

  data <- as.data.table(data)
  setkeyv(data, hhid)

  if ( !is.null(hhsize) ) {
    if ( !class(hhsize)=="character" | length(hhsize) != 1 | is.na(match(hhsize, colnames(data)))) {
      stop("strata must be a character defining the variable holding information on stratas and must be of length 1!\n")
    }
  } else {
    hhsize <- "hhsize"
    sizes <- data[,.N,by=key(data)]
    setnames(sizes, c(key(data), hhsize))
    data <- data[sizes]
  }
  if ( !is.null(pid) ) {
    if ( !class(pid)=="character" | length(pid) != 1 | is.na(match(pid, colnames(data)))) {
      stop("strata must be a character defining the variable holding information on stratas and must be of length 1!\n")
    }
  } else {
    pid <- "pid"
    sizes <- data[,.N,by=key(data)]
    data$pid <- paste(data[[hhid]], ".",unlist(sapply(sizes[["N"]], function(x) { seq(1, x) })), sep="")
  }
  classes <- sapply(data,class)
  for(i in seq_along(classes)){
    if(classes[[i]][1]=="labelled"){ 
      class(data[[i]]) <- classes[[i]][-1]
    }
  }
  invisible(new("dataObj", data=data, hhid=hhid, hhsize=hhsize, pid=pid, weight=weight, strata=strata, ispopulation=population))
}

