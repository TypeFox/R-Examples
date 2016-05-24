convert.report <- function(rep,inmodel=NULL,outmodel="MAgPIE",full=FALSE,as.list=TRUE) {
  .convert <- function(input,inmodel=NULL,outmodel="MAgPIE",full=FALSE) {
    map  <- magclassdata$map
    if(!(outmodel %in% names(map)))            stop("No existing transformation rules for output model \"",outmodel,"\"!",call.=FALSE)
    if(!(inmodel %in% names(map[[outmodel]]))) stop("No existing transformation rules for input model \"",inmodel,"\" in combination with output model \"",outmodel,"\"!",call.=FALSE)
    if(outmodel %in% names(input))             stop("Input already contains data for model \"",outmodel,"\"",call.=FALSE)
    map <- map[[outmodel]][[inmodel]]
    mag <- input[[inmodel]]
    if("GLO" %in% getRegions(mag)) map$GLO <- "GLO"
    outmag <- mag[rep(1,length(map)),,unlist(magclassdata$trans)[unlist(magclassdata$trans) %in% getNames(mag)]]
    outmag[,,] <- NA
    dimnames(outmag)[[1]] <- names(map)   
    for(reg in names(map)) {
      #sum
      elem <- getNames(mag)[getNames(mag) %in% magclassdata$trans$sum]
      if(length(elem)>0) if(!is.na(map[[reg]])) suppressWarnings(outmag[reg,,elem] <- colSums(mag[map[[reg]],,elem]))
      #mean
      elem <- getNames(mag)[getNames(mag) %in% magclassdata$trans$mean]
      if(length(elem)>0) if(!is.na(map[[reg]])) suppressWarnings(outmag[reg,,elem] <- colMeans(mag[map[[reg]],,elem]))    
    }
    out <- list();
    if(full) out[[inmodel]]=mag;
    out[[outmodel]]=outmag;
    return(out)
  }
  if(is.character(rep)) rep <- read.report(rep)
  if(is.null(inmodel)) {
    if(length(names(rep[[1]]))==1) inmodel <- names(rep[[1]])
    else stop("Not clear which model should be used as input!")
  }
  if(inmodel!=outmodel) {
    rep <- lapply(rep,.convert,inmodel,outmodel,full)  
  } 
  if(!as.list) {
    for(scenario in names(rep)) {
      for(model in names(rep[[scenario]])) {
        getNames(rep[[scenario]][[model]]) <- paste(scenario,model,getNames(rep[[scenario]][[model]]),sep=".")
      }
    }
    rep <- mbind(unlist(rep,recursive=FALSE))
    names(dimnames(rep))[3] <- "scenario.model.value"
  }
  return(rep)  
}