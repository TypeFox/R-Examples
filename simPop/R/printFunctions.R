setMethod("show", "dataObj", function(object){
  if ( is.null(object@ispopulation) ) {
    cat("this object does not contain any data yet!\n")
  } else {
    if(object@ispopulation){
      dname <- "population"
    } else {
      dname <- "survey sample"
    }
    cat("\n -------------- \n")
    cat(dname, "of size", nrow(object@data), "x", ncol(object@data), "\n")
    cat("\n Selected important variables: \n")
    cat("\n household ID:", object@hhid)
    cat("\n personal ID:", object@pid)
    cat("\n variable household size:", object@hhsize)
    cat("\n sampling weight:", object@weight)
    cat("\n strata:", object@strata)
    cat("\n -------------- \n")
    cat("\n")
  }
})

setMethod("show", "simPopObj", function(object){
  if ( is.null(popObj(object)) ) {
    cat("synthetic population has not been generated!\n")
  } else {
    dname <- "synthetic population"
    cat("\n")
    cat("-------------- \n")
    dd <- dim(popData(object))
    cat(dname, " of size \n", dd[1], "x", dd[2], "\n")
    cat("\n")
    if(sampleData(object)[,all(object@sample@weight==1),with=TRUE]){
      cat("build from a population of size \n")
    }else{
      cat("build from a sample of size \n")  
    }  
    
    dd <- dim(sampleData(object))
    cat(dd[1],"x",dd[2])
    cat("\n")
    cat("-------------- \n\n")
    cat("variables in the population:\n")
    cat(paste(colnames(popData(object)), collapse=","))
    cat("\n")
  }
})
