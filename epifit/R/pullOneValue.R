##' Pull one set of values from variables included in data which are thought to include essentially the same information.
##'
##' When some part of the data is missing, the missing information may be recovered from another source of information. This function scans all variables which are thought to include essentially the same information, and pulls value from a variable which is not missing.
##' @param data a data.frame in which variables are included
##' @param varlist a character vector of variable list which is assumed to contain the same information.
##' @param check a logical value specifying whether to scan all variables in varlist for incompatible values or not.
##' @param choose a character value specifying which to choose among incompatible values. To choose \sQuote{first} or \sQuote{last} depends on the order of the data. \sQuote{lower} and \sQuote{upper} cannot be used for character variables. This option works only when check option is enabled (set to TRUE).
##' @return a vector including one set of values obtained from varlist.
##' @examples
##' dat <- data.frame(a1=c(NA,2,3), a2=c(1,NA,2), a3=c(1,2,NA), b=c(10,11,20))
##' dat
##' pullOneValue(dat, c("a1", "a2", "a3"))
##' @export
pullOneValue <- function(data=NULL, varlist=c(""), check=TRUE, choose=c("first", "last", "lower", "upper")){

  res <- NULL
  
  if(is.null(data) || !is.data.frame(data))
    stop("data must be data.frame")
  
  allvar <- colnames(data)
  
  for(i in 1:length(varlist)){
    if(!varlist[i] %in% allvar)
      stop("at least one variable in varlist is not included in data")
  }
  
  idx <- GetParamPosition(varlist, allvar)
  
  if(length(idx) < 2)
    return(data[,idx])
  
  n <- nrow(data)

  # check variable mode
  if(is.numeric(data[,idx[1]])){
    res <- as.numeric(rep(NA, n))
    for(i in 2:length(idx)){
      if(!is.numeric(data[,idx[i]]))
        stop("mode of variables are different")
    }
  } else if(is.character(data[,idx[1]])){
    res <- as.character(rep(NA, n))
    for(i in 2:length(idx)){
      if(!is.character(data[,idx[i]]))
        stop("mode of variables are different")
    }
  } else if(is.factor(data[,idx[1]])){
    res <- as.factor(rep(NA, n))
    levels(res) <- levels(data[, idx[1]])
    for(i in 2:length(idx)){
      if(!is.factor(data[,idx[i]]))
        stop("mode of variables are different")
    }
    warning("factor is not fully supported")
  } else {
    stop("unsupported variable type")
  }

  for(i in 1:n){
    for(j in idx){
      
      if(!is.na(data[i,j]) && !is.nan(data[i,j])){
        if(is.na(res[i])){
          res[i] <- data[i,j]
          if(!check)
            next
        } else {
          if(res[i] != data[i,j]){
            cat("Imcompatible data exists in row:", i,"\n")
            for(j in idx){
              cat(allvar[j], " = ", data[i,j], "\n")
            }
            if(choose=="last"){
              res[i] <- data[i,j]
            } else if(choose=="lower"){
              if(res[i] > data[i,j])
                res[i] <- data[i,j]
            } else if(choose=="upper"){
              if(res[i] < data[i,j])
                res[i] <- data[i,j]
            } else if(choose!="first"){
              stop("unknown option is specified in choose argument")
            }
          }
        }
      }
    }
  }
  return(res)
}
