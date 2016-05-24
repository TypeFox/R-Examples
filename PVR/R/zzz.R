.onLoad <-function(libname, pkgname)
{
	ver <- read.dcf(file.path(libname, pkgname, "DESCRIPTION"), "Version")
	ver <- as.character(ver)	
	
	packageStartupMessage("PVR ", ver, " loaded\n")
}

setOldClass("phylo")

setClass("PVR", representation(
				Eigen = "list", phyDist = "matrix", phylo = "phylo", Selection = "list", PVR = "list", VarPart = "list")
)

setClass("PSR", representation(
				PSRarea = "data.frame", PSR = "data.frame", Expect.area.values = "list", nullPSR = "matrix", BrownianPSR = "matrix"),
			contains = "PVR" 
)

setMethod("plot",
		signature(x = "PSR"),
		function (x, y, expected, ...) 
		{
			plot(slot(x, "PSR"), 
					xlab = "cumulative eigenvalues(%)", 
					ylab = "r squared", ylim = c(0,1), 
					xlim = c(0,1), bty = "L", expected = segments(0, 0, 1, 1), pch = 16, ...)
			
#			coords = data.frame(x = c(0,slot(x, "PSR")[,1]), y = c(0,slot(x, "PSR")[,2]))
#			polygon(coords)
		}
)

setMethod("show",
		signature(object = "PVR"),
		function (object) 
		{
			if(!is.null(object@PVR$R2)){
				
				cat("\n", "Selected vectors", "\n", sep = "")
				cat("Used selection method: ", object@Selection$Method, "\n", sep = "")
				cat("R2: ",  object@PVR$R2, "\n", sep = "")
				cat("Residuals: ",  "\n", object@PVR$Residuals, "\n", sep = "")	
			} else{
				
				cat("\n", "Unselected vectors", "\n", sep = "")
				print(object@Eigen)
			}					
		}
)

setMethod("show",
		signature(object = "PSR"),
		function (object) 
		{
			cat("\n", "Phylogenetic Signal Representation (PSR) curve", "\n", sep = "")
			cat("\t", "PSR area: ", object@PSRarea$PSR.area, "\n", sep = "")
			cat("\t", "Probability of PSR area been equal to the null (random) expectancy: ", object@PSRarea$null.p, "\n", sep = "")
			cat("\t", "Probability of PSR area been equal to the Brownian expectancy**: ", object@PSRarea$Brownian.p, "\n", sep = "")
			cat("\t", "Iterations: ", object@PSRarea$iterations, "\n", sep = "")
			cat("\t", "Trait type: ", attr(object, "trait.type"), "\n", sep = "")
			cat("\n","\t", "**: ", "If trait is binary than the null expectancy is equal to the Brownian expectancy", "\n", sep = "")
		}
)
