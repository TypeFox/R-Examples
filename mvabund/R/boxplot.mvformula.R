################################################################################
# BOXPLOT.mvformula: Display the abundances of several species and their       #
# relation to explanatory variables in boxplots                                #
################################################################################
boxplot.mvformula <- function (x, 
				n.vars=12, 
				overall.main="", 
				var.subset=NA,
				... )
{

	#pass all arguments to default function call.
	default.boxplot.mvformula(x=x, 
				n.vars=n.vars, 
				overall.main=overall.main, 
				var.subset=var.subset,
				...)

}

