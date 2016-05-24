#' Rotate spectral absorbances from wide to long format and back
#'
#' @author Andrew Sila \email{asila@cgiar.org} and Tomislav Hengl \email{tom.hengl@wur.nl}

## rotate from wide to lon format:
reshape.SpectraPoints <- function(sampled, idcol="SAMPLEID", prefix="X"){
  if(!class(sampled)=="SpectraPoints"){
    stop("Object of class 'SpectraPoints' expected")
  }
  ## guess the format:
  if(!(all(c(idcol, "wavenumber", "value") %in% names(sampled@data@ab))&ncol(sampled@data@ab)==3)){
    xd = as.list(sampled@data@ab[,-which(names(sampled@data@ab)==idcol)])
    ab = lapply(names(xd), function(x){data.frame(SAMPLEID=sampled@data@ab[,idcol], wavenumber=as.numeric(strsplit(x, "m")[[1]][2]), value=xd[[x]])})
    ab = do.call(rbind, ab)
  } else {
    ab <- stats::reshape(sampled@data@ab, timevar="wavenumber", idvar=idcol, direction="wide")
    if(prefix==""){
      names(ab)[2:ncol(ab)] <- as.numeric(sapply(names(ab)[2:ncol(ab)], function(x){strsplit(x, "value.")[[1]][2]}))
    } else {
      names(ab)[2:ncol(ab)] <- sapply(names(ab)[2:ncol(ab)], function(x){paste(prefix, strsplit(x, "value.")[[1]][2], sep="")})
    }
  }
  sampled@data@ab <- ab
  return(sampled)
}

## end of script;