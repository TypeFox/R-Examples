# Internal functions
UseDistMatrix <- function(matdist, knownpts, unknownpts){
  i <- factor(row.names(knownpts), levels = row.names(knownpts))
  j <- factor(row.names(unknownpts), levels = row.names(unknownpts))
  matdist <- matdist[levels(i), levels(j)]
  return(round(matdist, digits = 8))
}

ComputeInteractDensity <- function(matdist, typefct, beta, span)
{
  if(typefct == "pareto") {
    alpha  <- (2 ^ (1 / beta) - 1) / span
    matDens <- (1 + alpha * matdist) ^ (-beta)
  } else if(typefct == "exponential") {
    alpha  <- log(2) / span ^ beta
    matDens <- exp(- alpha * matdist ^ beta)
  } else {
    stop("Please choose a valid interaction function argument (typefct)")
  }
  matDens <- round(matDens, digits = 8)
  return(matDens)
}

ComputeOpportunity <- function(knownpts, matdens, varname = varname)
{
  matOpport <- knownpts@data[, varname] * matdens
  return(round(matOpport, digits = 8))
}

ComputePotentials <- function(unknownpts, matopport)
{
  unknownpts@data$OUTPUT <- apply(matopport, 2, sum, na.rm = TRUE)
  return(unknownpts)
}

ComputeReilly <- function(unknownpts, matopport)
{
  unknownpts@data$OUTPUT <- row.names(matopport)[apply(matopport, 2, which.max)]
  return(unknownpts)
}

ComputeHuff <- function(unknownpts, matopport)
{
  sumCol <- colSums(x = matopport, na.rm = TRUE)
  matOpportPct <- 100 * t(t(matopport) / sumCol)
  matOpportPct[is.na(matOpportPct) | is.infinite(matOpportPct)] <- 0
  unknownpts@data$OUTPUT <- apply(matOpportPct, 2, max, na.rm = TRUE)
  return(unknownpts)
}

TestSp <- function(x){
  if (substr(class(x),1,7) != "Spatial"){
    stop(paste("Your input (",quote(x),") is not a spatial object.", sep=""),
         call. = F)
  }
  if (is.na(x@proj4string)){
    stop(
      paste(
        "Your input (", quote(x),
        ") does not have a valid coordinate reference system.", sep=""),
      call. = F)
  }
}
