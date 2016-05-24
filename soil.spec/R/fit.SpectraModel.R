#' Fit SpectraModel using calibration data
#'
#' @author Andrew Sila \email{asila@cgiar.org} and Tomislav Hengl \email{tom.hengl@wur.nl}


.getTrans <- function(sampled, idcol="SAMPLEID", CO2.band = get("CO2.band", envir=spec.opts), ab.r=get("MIR", envir=spec.opts), order=1, gap=21){  ## works only with wide format
  ab <- slot(slot(sampled, "data"), "ab")
  ## select region of interest:
  ab.id <- ab[idcol]
  ab <- ab[,-which(names(ab)==idcol)]
  ## remove columns with missing values:
  ab <- ab[,sapply(ab, function(x){ifelse(any(is.na(x)), FALSE, TRUE)})]
  ## extract wavenumbers from column names:
  wavenumber <- as.numeric(sapply(names(ab), function(x){gsub("[^0-9.]", "", x)}))
  names(ab) <- wavenumber
  ## remove values in the CO2.band:
  sel <- wavenumber>=ab.r[1]&wavenumber<=ab.r[2]&(wavenumber<CO2.band[1]|wavenumber>CO2.band[2]) 
  ab <- cbind(ab.id, ab[,sel])
  ## remove noise from SpectraPoints:
  raw <- as.matrix(ab[,-which(names(ab)==idcol)])
  row.names(raw) <- ab[,idcol]
  der <- data.frame(signif(trans(raw, tr="derivative", order=order, gap=gap)$trans, 3))
  return(der)
}

if(!isGeneric("fit")){
  setGeneric("fit", function(formulaString, sampled, reference, ...){standardGeneric("fit")})
}

## fit SpectraModel using calibration ('reference') data
fit.SpectraModel <- function(formulaString, sampled, reference, idcol = "SAMPLEID", ab.r = get("MIR", envir=spec.opts), CO2.band = get("CO2.band", envir=spec.opts), ncomp, ncomp.max=15, repl=5, segment.type="interleaved", prefix="X", ...){

  if(nrow(reference)>200){  warning("Fitting 'SpectraModel' for large data set can be time consuming", call. = FALSE, immediate. = TRUE) }
  ## check if the soil column is available in the reference table:
  covs <- all.vars(formulaString)[-1]
  tvar <- all.vars(formulaString)[1]
  if(!any(names(reference)==tvar)){
     stop(paste("Target variable", tvar, "not found in the 'reference' object"))
  }
  ## if in long format, reshape to wide format:
  if(all(names(sampled@data@ab) %in% c("SAMPLEID", "wavenumber", "value"))){
     sampled <- reshape.SpectraPoints(sampled, prefix=prefix)
  }
  der <- .getTrans(sampled=sampled, idcol=idcol, CO2.band=CO2.band, ab.r=ab.r)
  ## merge cleaned Spectral values and reference (wet chemistry):
  calibset <- merge(reference, cbind(SAMPLEID=row.names(der), der), by="SAMPLEID")
  response <- eval(formulaString[[2]], calibset)
  ## prepare a data.frame for mvr_dcv function:
  if(any(names(attributes(response))=="na.action")){
    mis.row <- attr(response, "na.action") 
    calibset <- calibset[-mis.row,]
  }  
  calibset.df <- data.frame(response, train=rep(1, nrow(calibset)))
  ## TH: "X" added by default to column name!
  ## TH: subset to existing columns
  ## (just in case if some columns are ommitted due to missing values or similar)
  sel <- match(names(calibset), covs)
  calibset.df$IR <- as.matrix(calibset[,covs[sel[!is.na(sel)]]])
  row.names(calibset.df$IR) <- calibset[,idcol]
  ## Model selection...
  if(missing(ncomp)){
    if(!nrow(calibset)>0){
      stop("Empty calibration data found. See '?soil.spec::fit' for more info.")
    }
    require(chemometrics)
    message("Performing double cross-validation using 'mvr_dcv'...")
    try( res.pls <- chemometrics::mvr_dcv(formula=response~IR, data=calibset.df, ncomp=ncomp.max, method="svdpc", repl=repl, selstrat="relchange", segment0.type=segment.type, na.action=na.omit, segment.type=segment.type) )
    if(class(.Last.value)[1]=="try-error"){
      warning("Double cross-validation unsuccesfull")
      ncomp <- 5
    } else {
      ## The optimal number of PCs, minimum should be at least 3!
      ncomp <- ifelse(res.pls$afinal<3, 3, res.pls$afinal)
    }
  }
  ## Fit model:
  m <- plsr(formula=response~IR, data=calibset.df, ncomp=ncomp, validation="CV", ...)
  ## Copy model coefficients and SAMPLEID's for fitted.values:
  out.m <- list(fitted.values=m$fitted.values[,,ncomp])
  attr(out.m$fitted.values, "names") = attr(m$model$IR, "dimnames")[[1]]
  out.m$RMSEP <- signif(pls::RMSEP(m)$val[,,ncomp], 3)
  out.m$opt.coef <- coef(m, comps = 1:ncomp)
  out.m$Xmeans <- m$Xmeans
  out.m$Ymeans <- m$Ymeans 
  ## save to an object:
  message("Deriving principal components for 'Space' slot...")
  pca.der <- prcomp(as.formula(paste("~", paste(names(der), collapse="+"))), der)
  out.m$pca <- pca.der[c("sdev","rotation","center","scale","call")]
  out.m$pca$rotation <- out.m$pca$rotation[,1:3]
  class(out.m$pca) <- "prcomp"
  out <- new("SpectraModel", variable=tvar, Space=pca.der$x[,1:3], model=out.m)
  return(out)
}
                              
setMethod("fit", signature(formulaString = "formula", sampled = "SpectraPoints", reference = "data.frame"), fit.SpectraModel)

## list of models i.e. multiple soil properties:
setMethod("fit", signature(formulaString = "list", sampled = "SpectraPoints", reference = "data.frame"), function(formulaString, sampled, reference, ...){
   out <- NULL
   for(i in 1:length(formulaString)){
     out[[i]] <- fit.SpectraModel(formulaString[[i]], sampled, reference, ...)
   }
   return(out)
})

## end of script;