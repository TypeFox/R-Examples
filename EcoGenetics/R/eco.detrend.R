#' Detrending spatial data with polynomial interpolation
#' 
#' @description  This program performs a Trend Surface Analysis
#' (Borcard et al. 2011, Legendre and Legendre 2012, Lichstein et al 2002) 
#' for the data Z and the given coordinates, projected or in decimal degrees format,
#' in which case will be projected with \code{\link[SoDA]{geoXY}}. 
#' 
#' @param Z Data frame, matrix or vector with dependent variables.
#' @param XY Data frame, matrix or vector with projected coordinates (X, XY or XYZ).
#' For longitude-latitude data in decimal degrees format, use the option latlon = TRUE.
#' @param degree Polynomial degree.
#' @param center Should the data be centered? Default TRUE
#' @param scale Should the data be scaled? Default FALSE
#' @param raw Use raw and not orthogonal polynomials? Default FALSE
#' @param latlon Are the coordinates in decimal degrees format? Defalut FALSE. If TRUE,
#' the coordinates must be in a data.frame/matrix with the longitude in the first
#' column and latitude in the second. The position is projected onto a plane in
#' meters with the function \code{\link[SoDA]{geoXY}}.
#' 
#' @return An object of class "eco.detrend" with the following slots:
#' @return > POLY.DEG polynomial degree used in the analysis
#' @return > RES detrended data
#' @return > XY projected coordinates
#' @return > MODEL models selected with the Akaike criterion
#' @return > ANALYSIS object of class "eco.mlm" with 
#' the regression results for each variable 
#' 
#'
#' \strong{ACCESS TO THE SLOTS}
#' The content of the slots can be accessed 
#' with the corresponding accessors, using
#' the generic notation of EcoGenetics 
#' (<ecoslot.> + <name of the slot> + <name of the object>).
#' See help("EcoGenetics accessors") and the Examples
#' section below.
#' 
#' @examples
#' \dontrun{
#' 
#' data(eco2)
#' 
#' # original data
#' data1 <- matrix(eco2[["P"]][,1], 30, 30)
#' image(data1)
#' 
#' # original data + trend
#' data2 <- matrix(eco2[["P"]][,2], 30, 30)
#' image(data2)
#' 
#' # data detrending
#' data2.det <- eco.detrend(Z = eco2[["P"]][,2], XY =  eco2[["XY"]], degree =  1)
#' 
#' 
#' #-----------------------
#' # ACCESSORS USE EXAMPLE
#' #-----------------------
#' 
#' # the slots are accesed with the generic format 
#' # (ecoslot. + name of the slot + name of the object). 
#' # See help("EcoGenetics accessors")
#' 
#' data2.det <- ecoslot.RES(data2.det)       # detrended data in slot RES
#' 
#' data2.det <- matrix(data2.det[,1], 30, 30)
#' image(data2.det)
#' 
#' 
#'}
#'
#' 
#' @references 
#' 
#' Borcard D., F. Gillet, and P. Legendre. 2011. Numerical ecology with R. 
#' Springer Science & Business Media.
#' 
#' Legendre P., and L. Legendre. 2012. Numerical ecology. Third English edition.
#' Elsevier Science, Amsterdam, Netherlands.
#' 
#' Lichstein J., T. Simons, S. Shriner, and K. Franzreb. 2002. 
#' Spatial autocorrelation and autoregressive models in ecology. 
#' Ecological monographs, 72: 445-463.
#' 
#' @export

setGeneric("eco.detrend", 
           function(Z, XY, degree, center = TRUE, 
                    scale = FALSE, raw = FALSE, 
                    latlon = FALSE) {
	
  if(latlon) {
    XY <- SoDA::geoXY(XY[,2], XY[,1], unit=1)
  }
  
  XY.scale <- scale(XY, center = center, scale = scale)
	
	polinomio <- data.frame(poly(as.matrix(XY.scale), degree = degree, raw = raw))
	
	
	capture.output(surf.ls <- eco.lmtree(Z, polinomio, analysis = "mlm"))
		
	
	residuos <- surf.ls@RESIDUALS
	formulas  <- lapply(1:length(surf.ls@MLM), function(i) surf.ls@MLM[[i]]$call)
	names(formulas) <- colnames(Z)
	
	res <- new("eco.detrend")
	res@POLY.DEG <- degree
	res@RES <- as.data.frame(residuos)
	res@XY <- XY
	res@MODEL <- formulas
	res@ANALYSIS <- surf.ls
	
	res
	
})

