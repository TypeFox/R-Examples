"MetaTable" <- 
function (x)
{
 
 rb <- rbar (x)
 vr <- varr (x)
 ve <- vare (x)
 pv <- pvse (x)[1]
 lclhet <- CIrb(x,,F)[1]
 uclhet <- CIrb(x,,F)[2]
 lclhom <- CIrb(x)[1]
 uclhom <- CIrb(x)[2]
 rho <- rhoCA(x)
 vrho <- varRCA(x)
 pva <- pvaaa(x)
 clcl <- CredIntRho(x, level=.8)[[1]]
 cucl <- CredIntRho(x, level=.8)[[2]]
mat <- data.frame(rbar = rb, Variance.rbar = vr, VarianceSamplingError = ve,
	PercentDueError = pv, HET95LCL = lclhet, HET95UCL = uclhet, HOM95LCL = lclhom,
	HOM95UCL = uclhom, RHO = rho, VarianceRho = vrho, PercentDueErrorCorrect = pva,
	CredInt80LCL = clcl, CredInt80UCL = cucl)
 return(mat)
 }
 
 
 
 