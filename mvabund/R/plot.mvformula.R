################################################################################
## PLOT.mvformula: Plot functions for mvabund objects                          #
## (multivariate abundance data)                                               #
################################################################################

plot.mvformula <- function(	x,
				y=NA,
				type="p", 
				var.subset=NA, 
				n.vars=if(any(is.na(list(var.subset)))) 12 else length(var.subset),
				xvar.select=TRUE,
				xvar.subset=NA, 
				n.xvars=NA,
        transformation="log", 
				... ) 
{

	default.plot.mvformula(x,
				y=NA,
				type=type, 
				n.vars=n.vars,
				var.subset=var.subset,
				xvar.select=xvar.select,
				n.xvars=n.xvars,
				xvar.subset=xvar.subset,  
                                transformation.formula=transformation,
				...) 
 
}

####### END DEFINITION OF FUNCTION ########
setMethod("plot", "mvformula", plot.mvformula)
setMethod("plot", signature("formula","logical"), plot.mvformula)