#' Generate objects of class Spectra or SpectraPoints
#'
#' @author Tomislav Hengl \email{tom.hengl@wur.nl} and Andrew Sila \email{asila@cgiar.org}

Spectra <- function(samples, ab, idcol="SAMPLEID", replace.prefix="m"){
    ## if missing absorbance generate empty data frames:
    if(missing(ab)){
       ab = data.frame(SAMPLEID=NA, wavenumber=NA, value=NA)
    }
    ## subset the spectra to valid range:
    ab.r <- sapply(c("MIR", "NIRS", "NIRP", "VISNIR1", "VISNIR2", "VISNIR3"), get, envir=spec.opts)
    ## long format:
    if(!any(names(ab) %in% c(idcol, "wavenumber", "value"))&ncol(ab)==3){
       sel <- ab$wavenumber>min(ab.r[1,]) & ab$wavenumber<max(ab.r[2,])
       if(sum(sel)<nrow(ab)){
          warning("Spectral absorbances outside the range. See '?spec.env' for more info.")
          ab <- ab[sel,]
       }
    }
    ## replace prefix:
    names(ab)[-which(names(ab)==idcol)] <- sapply(names(ab)[-which(names(ab)==idcol)], function(x){gsub(replace.prefix, "X", x)})
    ## rename the id column:
    names(ab)[which(names(ab)==idcol)] <- "SAMPLEID"
    out <- new("Spectra", samples=samples, ab=ab)
    return(out)
}

SpectraPoints <- function(metadata, Spectra, sp, silent=TRUE){

    ## check if all IDs in the samples table match spatial ID's:
    if(silent==FALSE){
      ssn <- attr(coordinates(sp), "dimnames")[[1]]
      ssn.m <- !(levels(as.factor(paste(Spectra@samples$SAMPLEID))) %in% levels(as.factor(ssn)))
      if(any(ssn.m)){
        warning(paste(sum(ssn.m), "sample(s) missing coordinates."), immediate.=TRUE)
      }
    }

    ## if missing metadata generate it:
    if(missing(metadata)){
       mnames <- get("mdnames", envir=spec.opts)
       metadata <- data.frame(as.list(rep(NA, length(mnames))))
       names(metadata) <- mnames
       metadata$MID <- levels(Spectra@samples$MID)
    }
    ## subset the spectra to only points with coords?
    out <- new("SpectraPoints", metadata=metadata, data=Spectra, sp=sp)
    return(out)
}

## print a summary of "SpectraPoints" object:
setMethod("summary", signature(object = "SpectraPoints"), function(object){
  cat("Object of class Spectra\n")
  ## print metadata which is not NA:
  for(i in 1:nrow(object@metadata)){
    if(!is.na(object@metadata[,i])){
      txt <- names(object@metadata)[i]
      cat(" ", txt, rep(" ", 17-nchar(txt)), ": ", paste(object@metadata[,i], collapse=", "), "\n", sep="")
    }
  }
  ## print summary of ab table:
  cat(" Number of bands  :", ncol(object@data@ab)-1, "\n")
  ab <- slot(slot(object, "data"), "ab")
  ## select region of interest:
  idcol <- names(object@data@samples)[1] 
  ab.id <- ab[idcol]
  ab <- ab[,-which(names(ab)==idcol)]
  ## extract wavenumbers from column names:
  wavenumber <- as.numeric(sapply(names(ab), function(x){gsub("[^0-9.]", "", x)}))
  cat(" Start wavenumber :", min(wavenumber, na.rm=TRUE), "\n")
  cat(" End wavenumber   :", max(wavenumber, na.rm=TRUE), "\n")
  ## print summary of sp slot:
  summary(object@sp)
})

## end of script;