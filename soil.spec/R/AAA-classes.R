#' Class and methods for spectra data
#'
#' @author Tomislav Hengl \email{tom.hengl@wur.nl} and Andrew Sila \email{asila@cgiar.org}
#' Note: similar classes are also available via the 'inspectr' package [http://r-forge.r-project.org/projects/inspectr/]


#################### CLASS DEFINITIONS #########################

## class for pure measurements
setClass("Spectra", representation(samples = "data.frame", ab = "data.frame"), validity = function(object) {
   ## check the column names:
  knames <- c("SAMPLEID", "MID", "DateTime","LWN","DateTime","Material","Zero.Filing","Resolution")
   if(!any(names(object@samples) %in% knames)){
      return(paste("Expecting column names:", paste(knames, collapse=", ")))
   }
   if(!any(class(object@samples$DateTime) %in% "POSIXct")){
      return("'DateTime' of class 'POSIXct' required")
   }
   cnames <- c("SAMPLEID", "wavenumber", "value")
   if(!any(names(object@ab) %in% cnames)&ncol(object@ab)==3){
      return(paste("Expecting column names:", paste(cnames, collapse=", ")))
   }
   if(!any(names(object@ab)==cnames[1])&!ncol(object@ab)==3){
      return(paste("Expecting sample number column:", cnames[1]))
   }
   ## check the range of ir values
   ab.r <- sapply(c("MIR", "NIRS", "NIRP", "VISNIR1", "VISNIR2", "VISNIR3"), get, envir=spec.opts)
   if(any(names(object@ab)=="wavenumber")){  ## long format
     if(length(object@ab$wavenumber)>0&sum(!is.na(object@ab$wavenumber))>0){
        if(any(object@ab$wavenumber <= min(ab.r[1,]) | object@ab$wavenumber >= max(ab.r[2,]))){
          return(paste("Wavenumber not in the", paste(min(ab.r[1,]), max(ab.r[2,]), collapse=" - "), "IR range"))
        }
     }
   } else { ## for the wide format
     cn <- names(object@ab)[-which(names(object@ab)=="SAMPLEID")]
     wn <- as.numeric(sapply(cn, function(x){gsub("[^0-9.]", "", x)}))
     if(any(wn <= min(ab.r[1,]) | wn >= max(ab.r[2,]))){
       return(paste("Wavenumber not in the", paste(min(ab.r[1,]), max(ab.r[2,]), collapse=" - "), "IR range"))
     }
   }
   ## check if all IDs in the samples match other tables:
   if(length(object@ab$SAMPLEID)>0&sum(!is.na(object@ab$SAMPLEID))>0){
      if(any(!(levels(as.factor(paste(object@samples$SAMPLEID))) %in% levels(as.factor(paste(object@ab$SAMPLEID)))))){
        return("'SAMPLEID' sample serial numbers in the samples table and data tables do not match")
      }
   }
})


## Complete class for soil spectroscopy measurements:
setClass("SpectraPoints", representation(metadata = "data.frame", data = "Spectra", sp = "SpatialPoints"), validity = function(object) {
   ## check if the object has standard names:
   mnames <- get("mdnames", envir=spec.opts)
   if(!any(names(object@metadata) %in% mnames)){
      return(paste("Expecting column names:", paste(mnames, collapse=", ")))
   }
   ## check if the metadata i.d.'s corresponds:
   if(any(!(levels(as.factor(paste(object@data@samples$MID))) %in% levels(as.factor(paste(object@metadata$MID)))))){
       return("'MID' metadata IDs in the samples table and metadata tables do not match")
   }
})

## output class for spectral models
setClass("SpectraModel", representation(variable = "character", Space = "matrix", model = "ANY"), validity = function(object) {
   if(length(object@variable)>1){
     return("Single variable name expected")
   }
   if(!all(attr(object@model$fitted.values, "dimnames")[[1]] %in% attr(object@Space, "dimnames")[[1]])){
     return("Sample ID's between the 'Space' and 'model' slots do not match.")
   }
})


## end of script;