#' Predict SpectraModel using new data
#'
#' @author Tomislav Hengl \email{tom.hengl@wur.nl} and Andrew Sila \email{asila@cgiar.org}

if(!isGeneric("predict")){
  setGeneric("predict", function(object, ...){standardGeneric("predict")})
}

## Predict soil properties:
predict.SpectraPoints <- function(object, idcol = "SAMPLEID", model, variable=model@variable, output, validate = FALSE, model.class = "mvr", confidence.band = TRUE, prob. = .90, signif.digits = 3, st.wavenumbers=wavenumbers, instr.range = c("ten-mir", "alp-mir", "mpa-nir")[1], ...){
   
   ## remove duplicate bands:
   object@data@ab <- object@data@ab[,!duplicated(names(object@data@ab))]
   ## derive first der:
   suppressWarnings( IR <- .getTrans(object) )
   ## rename the column names using standard band names:
   if(missing(st.wavenumbers)){
     st.wavenumbers <- wavenumbers[wavenumbers$TYPE==instr.range,]
   }
   w.n <- as.numeric(sapply(names(IR), function(x){gsub("[^0-9.]", "", x)}))
   rn <- rank(w.n, ties.method="first")
   w.s <- cut(w.n[rn], breaks=c(st.wavenumbers$LOWER, st.wavenumbers$UPPER[nrow(st.wavenumbers)]), labels=st.wavenumbers$BAND, include.lowest=TRUE)
   w.s <- as.numeric(sapply(w.s, function(x){gsub("[^0-9.]", "", x)}))
   names(IR)[rn] <- paste0("X", ifelse(is.na(w.s), round(w.n[rn], 1), w.s))
   object@data@ab <- cbind(object@data@ab[idcol], IR)
   
   ## check that all bands in newdata required by the model are available:
   if(!all(attr(model@model$opt.coef, "dimnames")[[1]] %in% names(object@data@ab)[-1])){		
     stop("Missing column names (wavenumbers) in the 'object'. See '?predict.SpectraPoints' for more info.")
   } 
   ## Check overlap in feature space:
   if(validate==TRUE){
     outliers <- validate(object, model, silent=FALSE)
   }
   if(model.class=="mvr"){
     ## selection of bands from the regression model:
     bands <- match(attr(model@model$opt.coef, "dimnames")[[1]], names(object@data@ab)[-1])
     IR <- as.matrix(IR[,bands]) ## TH: this sorts IR matrix so the coefficient names match column names
     ## attach unique row names (it is fine to have duplicate samples but they should have different unique name) otherwise a warning is issued!
     row.names(IR) <- make.unique(paste(object@data@ab[,idcol]))
     B <- rowSums(model@model$opt.coef, dims = 2)
     B0 <- model@model$Ymeans - model@model$Xmeans %*% B
     pred <- IR %*% B + rep(B0, each = dim(IR)[1])  
     attr(pred, "dimnames")[[2]] <- variable
     ## return values as "data.frame"
     if(!(model@variable=="PHIHOX"|model@variable=="SNDLDF")){
       if(confidence.band == TRUE){
          t.v <- qt(prob.+(1-prob.)/2, df=10000)
          lower <- exp(pred - t.v * model@model$RMSEP[2])
          upper <- exp(pred + t.v * model@model$RMSEP[2])
       }
       pred <- exp(pred)
     } else{
       if(confidence.band == TRUE){
          t.v <- qt(prob.+(1-prob.)/2, df=10000)
          lower <- pred - t.v * model@model$RMSEP[2]
          upper <- pred + t.v * model@model$RMSEP[2]
       }
     }
     if(confidence.band == TRUE){
        attr(lower, "dimnames")[[2]] <- paste(variable, "_lower", sep="")
        attr(upper, "dimnames")[[2]] <- paste(variable, "_upper", sep="")
        out <- data.frame(signif(pred,signif.digits), signif(lower,signif.digits), signif(upper,signif.digits))
     } else{
        out <- data.frame(signif(pred,signif.digits))
     }
   }
   if(!missing(output)){
     warning("Not available at the moment")
   }
   ## return object of a SoilProfileCollection class:
   return(out)
      
}

setMethod("predict", signature(object = "SpectraPoints"), predict.SpectraPoints)

setMethod("predict", signature(object = "data.frame"), function(object, idcol = "SAMPLEID", ...){
   ## check if the table is usable
   if(ncol(object)<3|sum(names(object)==idcol)==0){
     stop("Missing 'idcol' or absorbances")
   }
   ## coordinates missing so generate "0"s:
   sp <- SpatialPoints(data.frame(lat=0,lon=0), proj4string=CRS(as.character(NA)))
   ## prepare 'samples' table
   samples <- cbind(object[idcol], MID="NA", DateTime=Sys.time())
   ## convert to "SpectraPoints"
   object.sp <- SpectraPoints(Spectra=Spectra(samples, object), sp=sp)
   out <- predict.SpectraPoints(object.sp, ...)
   return(out)
})
