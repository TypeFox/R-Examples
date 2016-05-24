
read.report <- function(file,as.list=TRUE) {
  
  .trim <- function(a) return(gsub("(^ +)|( +$)", "",as.character(a)))
  
  .return_magpie <- function(tmp,scenario,model) {
    regions <- unique(as.character(tmp$Region))
    names(regions) <- regions
    years <- sub("X","y",grep("^X[0-9]{4}$",dimnames(tmp)[[2]],value=TRUE))
    names <- unique(paste(tmp$Variable, "#SPLITHERE# (",tmp$Unit,")",sep =""))
    units <- sub("^.*#SPLITHERE# \\((.*)\\)$","\\1",names)
    names(names) <- sub("#SPLITHERE#", "", names)
    names <- sub("#SPLITHERE#","",names)
    #delete dots if they are aparrently not used as dimension separator
    ndots <- nchar(gsub("[^\\.]*","",names))
    if(any(ndots!=ndots[1])) names <- gsub("\\.","",names)
    mag <- new.magpie(sub("ZZZZZZGLO","GLO",(sort(sub("GLO","ZZZZZZGLO",regions)))),years,names)
    yearelems <- grep("^X[0-9]{4}$",dimnames(tmp)[[2]])
    regions[order(sub("GLO","ZZZZZZGLO",regions))] <- dimnames(mag)[[1]]
    mag <- as.array(mag)
    coord <- cbind(regions[tmp$Region],rep(years,each=dim(tmp)[1]),names[paste(tmp$Variable, " (",tmp$Unit,")",sep ="")])
    if(dim(coord)[1]>length(mag)) {
      duplicates <- duplicated(coord)
      warning("Duplicate entries found for model \"",model,"\" and scenario \"",scenario,"\" and only the last entry found in the data will be used (duplicate entries: ",paste(apply(rbind(NULL,unique(coord[duplicates,c(1,3)])),1,paste,collapse="|"),collapse=", "),")!")    
    }

    mag[coord] <- suppressWarnings(as.numeric(as.vector(as.matrix(tmp[,yearelems]))))
    names(dimnames(mag)) <- c("region","year","variable")
return(as.magpie(mag,spatial=1,temporal=2))  
  }

  .readmif <- function(file) {
    default_header <- c("Model","Scenario","Region","Variable","Unit","X2005",
                        "X2010","X2020","X2030","X2040","X2050","X2060","X2070",
                        "X2080","X2090","X2100")
    #determine seperator
    s <- read.table(file,sep=";",header=FALSE,nrows=1) 
    if (all(names(s) == "V1")) sep <- "," else sep <- ";"
    #recognize header
    s <- read.table(file,sep=sep,header=FALSE,nrows=1) 
    header <- (.trim(s[,1]) == "Model" | .trim(s[,1]) == "MODEL")
    #read in raw data
    raw <- read.table(file,sep=sep,header=header,stringsAsFactors=FALSE,na.strings = "N/A",fileEncoding = "UTF8")
    ugly_format <-  all(is.na(raw[,dim(raw)[2]]))
    if(ugly_format) raw <- raw[,-dim(raw)[2]]
    
    if("number of items read is not a multiple of the number of columns" %in% names(warnings())) {
      stop("Inconsistent input data! At least one line is incomplete!")
    }
    
    #rename from uppercase to lowercase
    if (header & .trim(s[,1]) == "MODEL") {
      names(raw)[1:5] <- default_header[1:5]
    }
    
    if(!header) {
      if(dim(raw)[2]==length(default_header)) dimnames(raw)[[2]] <- default_header
      else stop("Cannot read report. No header given and report has not the standard size!")   
    }
    
    output <- list()
    raw$Scenario <- .trim(raw$Scenario)
    raw$Model    <- .trim(raw$Model) 
    raw$Region   <- .trim(raw$Region)
    raw$Unit     <- .trim(raw$Unit)
    raw$Variable <- .trim(raw$Variable)
    
    raw$Region <- sub("R5\\.2","",raw$Region)
    raw$Region <- sub("World|glob","GLO",raw$Region)
    models <- unique(raw$Model)
    scenarios <- unique(raw$Scenario)
    for(scenario in scenarios) {
      output[[scenario]] <- list()
      for(model in models) {  
        if (nrow(raw[raw$Model==model & raw$Scenario==scenario,]) > 0) {
          output[[scenario]][[model]] <- .return_magpie(raw[raw$Model==model & raw$Scenario==scenario,],scenario,model)
          if(!as.list) getNames(output[[scenario]][[model]]) <- paste(scenario,model,getNames(output[[scenario]][[model]]),sep=".")          
        }
      }
    }
    return(output)
  }
  
  #expand wildcards
  file_name_unexpanded <- file  
  file <- Sys.glob(file)
  if(length(file)>1) {
    output <- NULL
    for(f in file) {
      output <- c(output,.readmif(f))
    }
  } else if(length(file)==0) {
    stop("File ",file_name_unexpanded," could not be found!")
  } else {
    output <- .readmif(file)
  }
    
  if(!as.list) {
     output <- mbind(unlist(output,recursive=FALSE))
     names(dimnames(output))[3] <- "scenario.model.variable"
  }
  return(output)
}













