plot.HMDdata <-
function(x, ...){
  xx <- x
  at <- attributes(xx)
  type <- strsplit(at$`country-data-sex`,
                   split="-")[[1]][2]
  if(is.null(dim(xx))){
    ## one-dimensional graph
    x <- as.numeric(at$names)
    if(type=="Rates"){
      logrates <- log(as.vector(xx))
      plot(x, logrates, ...)
    }
    if(type=="Population"){
      population <- as.vector(xx)
      plot(x, population, ...)
    }
    if(type=="Deaths"){
      deaths <- as.vector(xx)
      plot(x, deaths, ...)
    }
    if(type=="Exposures"){
      exposures <- as.vector(xx)
      plot(x, exposures, ...)
    }
  }else{
    x1 <- as.numeric(at$dimnames[[1]])
    x2 <- as.numeric(at$dimnames[[2]])
    listMSF <- list(ages=x1, years=x2)
    gridMSF <- expand.grid(listMSF)
    if(type=="Rates"){
      gridMSF$logrates <- log(as.vector(xx))
    }
    if(type=="Population"){
      gridMSF$population <- as.vector(xx)
    }
    if(type=="Deaths"){
      gridMSF$deaths <- as.vector(xx)
    }
    if(type=="Exposures"){
      gridMSF$exposures <- as.vector(xx)
    }
    levelplot(gridMSF[,3] ~ years * ages, gridMSF, ...)
  }
}
