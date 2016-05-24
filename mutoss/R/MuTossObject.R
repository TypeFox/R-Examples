setClass("ErrorControl",
		representation  = representation(
				type           = "character", # Type of error rate controlled for (FWER, FWER.weak, FDR, FDX, gFWER, perComparison (?))
				alpha          = "numeric",   # Error rate of the procedure
				k              = "numeric",   # Additional parameter for  generalized FWE control
				q              = "numeric"    # Additional paramter for FDX control
		)
)

##TODO: add slots for "model": formula, link, family
##TODO: write some header, it is the main object that comes out after four weeks of hard work!!

setClass("Mutoss",
		representation    = representation(
				data			= "ANY",            # Raw data used in model
				model			= "ANY",            # link function,error family and design 
				description     = "character",
				statistic       = "numeric",      	# For Z, T or F statistics (maybe different slots?)
				hypotheses		= "ANY",
				hypNames        = "character",    	# Identifiers for the hypotheses tested
				criticalValues	= "numeric",	    # Procedure-specific critical values
				pValues         = "numeric",      	# Raw p-values. Either imported or calculated with data and model 
				adjPValues      = "numeric",      	# Procedure-specific adjusted p-values
				errorControl    = "ErrorControl", 	# Details of the multiplicity control procedure used.				  
				rejected        = "logical",      	# Logical vector of the output of a procedure at a given error rate
				qValues         = "numeric",      	# Storey's estimates of the supremum of the pFDR
				locFDR          = "numeric",      	# Efron's local fdr estimates (by which method?)
				pi0             = "numeric",      	# Estimate of the proportion of null hypotheses (by which method?)
				confIntervals   = "matrix",       	# Confidence intervals for selected parameters (of which kind? selected how?)
				commandHistory  = "character"
		)
)
