#' Validate that the new data matches referent model
#'
#' @author Tomislav Hengl \email{tom.hengl@wur.nl} and Andrew Sila \email{asila@cgiar.org}

if(!isGeneric("validate")){
  setGeneric("validate", function(obj, ...){standardGeneric("validate")})
}

validate.SpectraPoints <- function(obj, model, silent=FALSE){
     if(!class(model)=="SpectraModel"){
       stop("Model of class 'SpectraModel' expected")
     }
     IR <- as.matrix(.getTrans(obj))
     message("Analyzing overlap of absorbance values in feature space (using PC1-3)...")
     IR.df <- data.frame(IR)
     Space1 <- model@Space ## referent space
     x <- which(! names(IR.df) %in% attr(model@model$pca$rotation, "dimnames")[[1]])
     if(length(x)>0){ stop(paste("Values missing for bands:", paste(attr(model@model$pca$rotation, "dimnames")[[1]][x], collapse=", ", sep=""))) }
     Space0 <- predict(model@model$pca, IR.df) ## new points
     s.r <- range(c(Space1[,1],Space0[,1]))
     ## 3D point pattern analysis:
     require(spatstat)
     Space1 <- spatstat::pp3(Space1[,1], Space1[,2], Space1[,3], spatstat::box3(c(s.r[1]-.1, s.r[2]+.1)))
     Space0 <- spatstat::pp3(Space0[,1], Space0[,2], Space0[,3], spatstat::box3(c(s.r[1]-.1, s.r[2]+.1)))
     ## distance between two 3D point patterns:
     d. <- spatstat::nncross.pp3(Space0, Space1)
     nd. <- spatstat::nndist.pp3(Space1)
     if(any(d.$dist > 2*quantile(nd., .95))){
       outliers <- which(d.$dist > 2*quantile(nd., .95))
       if(silent==TRUE){
         warning(paste("Possible outliers in feature space detected (PC1-3):", paste(outliers, collapse=", ")))
       }
       return(outliers)
     } else {
        message("... no outliers detected.")
     }
}

setMethod("validate", signature(obj = "SpectraPoints"), validate.SpectraPoints)

# end of script;