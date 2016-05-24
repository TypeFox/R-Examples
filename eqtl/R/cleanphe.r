#####################################################################
#
# cleanphe.R
#
# copyright (c) 2008-3, Ahmid A Khalili
# 
# last modified Jul, 2008
# first written Mar, 2008
# Licensed under the GNU General Public License version 2 (June, 1991)
#
# Part of the R/eqtl package
# Contains: cleanphe 
#
######################################################################

######################################################################
#
# cleanphe: Remove undesired phenotypes and LOD results from cross
#               and scanone object respectively
#
######################################################################

`cleanphe` <-
function(x,string='Buffer')
{
	require(qtl)

	if ( ! any(class(x) == 'cross') & ! any(class(x) == 'scanone') )
		stop("Input should have class \"cross\" or class \"scanone\".")
	if ( ! is.character(string) & !is.vector(string))
		stop("Expecting a characte  r vector for string")

	if ( class(x)[1] == 'scanone') {
		coord <- grep(string,names(x))
		cat("Drop ",length(coord),"lodcolumn\n");
		if(!length(coord)) return(x);
		nf <- x[,-coord]
	} else {
		coord <- grep(string,names(x$pheno))
		cat("Drop ",length(coord),"phenotypes\n");
		if(!length(coord)) return(x);
		nf <- x
		nf$pheno <- x$pheno[,-coord]
  }
	attributes(nf)$class <- class(x)
	try(return(nf),silent=FALSE)
}

