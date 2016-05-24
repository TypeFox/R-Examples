meanvar.plot <- function(x, ...)
UseMethod("meanvar.plot")

################################################################################
# PLOT.MEANVAR: Constructs mean/variance plots for an mvabund objects          #
# (multivariate abundance data)                                                #
################################################################################
meanvar.plot.mvabund <- function(x, 
  				n.vars = NULL, 
				var.subset=NULL, 
				subset=NULL, 
				table=FALSE, 
				... ) {

	#Pass all of the arguments to the default plot functions.
	MeanVarTable <- default.meanvar.plot.mvabund(x=x, 
		  				n.vars=n.vars, 
						var.subset=var.subset, 
						subset=subset, 
						table=table, 
						...)

	#Return the table if requested from the default plot call.
	if (table) { return(MeanVarTable) }

}
######## END MEANVAR.MVABUND FUNCTION #########

################################################################################
# PLOT.MEANVAR: Constructs mean/variance plots for an mvabund objects          #
# (multivariate abundance data)                                                #
# Separate means and variances are calculated for different samples as defined #
# by the formula, groups are chosen by factors in the independent data         #
# if there are no factors an overall Mean-Variance Plot is drawn               #
################################################################################
meanvar.plot.mvformula <- function(  x,  
				n.vars = NULL,
				var.subset=NULL, 
				subset=NULL, 
				table=FALSE, 
  				overall.main=NULL, 
				overlay=TRUE, 
				... ) {

	#Pass all of the arguments to the default plot functions.
	MeanVarTable <- default.meanvar.plot.mvformula(x=x,  
				n.vars = n.vars,
				var.subset=var.subset, 
				subset=subset, 
				table=table, 
  				overall.main=overall.main, 
				overlay=overlay, 
				... )

	#Return the table if requested from the default plot call.
	if (table) return(MeanVarTable)

}
###### END THE MEANVAR FORMULA PLOT ######

setGeneric("meanvar.plot")
setMethod("meanvar.plot", "mvabund", meanvar.plot.mvabund)
setMethod("meanvar.plot", "matrix", meanvar.plot.mvabund)
setMethod("meanvar.plot", "data.frame", meanvar.plot.mvabund)
setMethod( "meanvar.plot", "mvformula", meanvar.plot.mvformula )
setMethod( "meanvar.plot", "formula", meanvar.plot.mvformula )

