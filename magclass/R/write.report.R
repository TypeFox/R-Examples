write.report <- function(x,file=NULL,model="MAgPIE",scenario="default",unit=NA,ndigit=4,append=FALSE,skipempty=TRUE) {
  if(is.list(x)) {
    if(is.list(x[[1]])) {
      for(scenario in names(x)){
        for(model in names(x[[scenario]])) {
          write.report(x[[scenario]][[model]],file=file,model=model,scenario=scenario,unit=unit,ndigit=ndigit,append=append)  
          append <- TRUE
        }
      }  
    } else {
      stop("Wrong format. x must be either a list of lists or a MAgPIE object! Only single list found!")
    }
  } else {

    if(!is.magpie(x)) stop("Input is not a MAgPIE object!")
    dimnames(x)[[1]] <- sub("^GLO(\\.{0,1}[0-9]*)$","World\\1",dimnames(x)[[1]])
    
    # If data was read in by read.report there is an attribute $dimnames$scenario.model.variable.
    # This will be used below to structure the output of write.report exactly like the input.
    # That means: write.report automatically recognizes models and scenarios and 
    # does not put it into the data-dimension
    if (scenario == "default" & length(names(attr(x,"dimnames"))) >2) {
      if (names(attr(x,"dimnames"))[[3]] == "scenario.model.variable") {
        scenario <- fulldim(x)[[2]]$scenario
        model <- fulldim(x)[[2]]$model
      }
    }

    unitdef<-unit
    ii<-1
    for (mod in model) {
      for (scen in scenario) {
        if (length(fulldim(x)[[2]]) == 5 ) {
          if (length(strsplit(getNames(x)[1],split="\\.")[[1]]) > 1 
              & scen %in% unlist(lapply(strsplit(getNames(x),split="\\."),'[[',1))
              & mod %in% unlist(lapply(strsplit(getNames(x),split="\\."),'[[',2))) {
          xtemp<-x[,,paste(scen,mod,sep=".")] 
          } else {
          xtemp <- x
          }
        } else {
          xtemp<-x           
        }
        ndata <- ndata(xtemp)
        nregions <- nregions(xtemp)
        nyears <- nyears(xtemp)
        regions <- getRegions(xtemp)
        years <- gsub(" ",0,prettyNum(getYears(xtemp,as.integer=TRUE),width=4))
        if(length(unit)==1) {
          nelem_with_brackets <- length(grep("\\(*\\)$",getNames(xtemp)))
          if(nelem_with_brackets==dim(xtemp)[3]) {
            tmp <- getNames(xtemp)
            dimnames(xtemp)[[3]] <- sub(" ?\\(.*\\)$","",tmp)
            unit <- sub("^[^\\(]*\\((.*)\\)$","\\1",tmp)    
          } else {
            if(nelem_with_brackets > 0) warning("Some but not all variable entries provide information in brackets which might be a unit information. To have it detected as unit all entries must provide this information!")
            unit <- rep(unit,ndata)
          }
        }
        
        output <- matrix(NA,nregions*ndata,5+nyears)
        colnames(output) <- c("Model","Scenario","Region","Variable","Unit",years)
        output[,"Model"] <- mod
        output[,"Scenario"] <- scen
        output[,"Region"] <- rep(regions,ndata)
        
        for(i in 1:ndata){
          if (length(fulldim(x)[[2]]) == 5){
            if (length(strsplit(getNames(xtemp)[i],split="\\.")[[1]]) > 1 
                & strsplit(getNames(xtemp)[i],split="\\.")[[1]][[1]]==scen 
                & strsplit(getNames(xtemp)[i],split="\\.")[[1]][[2]]==mod) {
              output[(i-1)*nregions + 1:nregions,"Variable"] <- strsplit(getNames(xtemp)[i],split="\\.")[[1]][[3]]
            } else {
              output[(i-1)*nregions + 1:nregions,"Variable"] <- gsub(".","|",getNames(xtemp)[i],fixed=TRUE)
            }
          } else {
            output[(i-1)*nregions + 1:nregions,"Variable"] <- gsub(".","|",getNames(xtemp)[i],fixed=TRUE)          
          }
          output[(i-1)*nregions + 1:nregions,"Unit"] <- unit[i]
          output[(i-1)*nregions + 1:nregions,5+1:nyears] <- round(xtemp[,,i],ndigit)
        }
        
        if(skipempty) {
          toskip <- which(rowSums(!is.na(output[,5+(1:nyears),drop=FALSE]))==0)
          if(length(toskip)>0) output <- output[-toskip,,drop=FALSE]
        }
        output[is.na(output)] <- "N/A"
        output[which(output=="NaN")]<-"N/A"
        if(is.null(file)) {
          print(output)
        } else {
          if(!file.exists(file)) append <- FALSE
          if(ii > 1) append <-TRUE
          write.table(output,file,quote=FALSE,sep=";",row.names=FALSE,col.names=!append,append=append,eol=";\n")
          ii<-ii+1
        }
        unit<-unitdef
      }
    }

    

  }
}
