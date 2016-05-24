#####################################################################
#
# gpt.R
#
# copyright (c) 2008-3, Ahmid A Khalili
# 
# last modified Jul, 2008
# first written Mar, 2008
# Licensed under the GNU General Public License version 2 (June, 1991)
# 
# Karl W Broman wrote the scanone function part of R/qtl package
# licensed under the GNU General Public License version 2 (June, 1991)
#
# Part of the R/eqtl package
# Contains: gpt
#
######################################################################

######################################################################
#
# gpt: Calculate the global permutation threshold.
#
######################################################################

`gpt` <-
function(cross,n_etrait=100,n_perm=1000)
{

	require(qtl)

	if ( class(cross)[2] != 'cross')
		stop("Input should have class \"cross\".")
	if ( !is.vector(n_etrait))
		stop("Expecting a numeric VECTOR of length 1 for n_etrait")
	if ( !is.numeric(n_etrait))
		stop("Expecting a NUMERIC vector of length 1 for n_etrait")
	if ( length(n_etrait) !=1)
		stop("Expecting a numeric vector of length 1! for n_etrait")
	if ( !is.vector(n_perm))
		stop("Expecting a numeric VECTOR of length 1 for n_perm")
	if ( !is.numeric(n_perm))
		stop("Expecting a NUMERIC vector of length 1 for n_perm")
	if ( length(n_perm) !=1)
		stop("Expecting a numeric vector of length 1! for n_perm")
	if ( n_etrait > nind(cross))
		stop("Too many individuals: n_etrait > nind(cross)")

	cat(	'WARNINGS Global Permutation Threshold is computing: this calcul is long and may run for few days\n',
		'Parameters:',n_perm,'permutations on',n_etrait,'individuals\n'	)

	nphe <- nphe(cross)
	permutation <- suppressWarnings(scanone(	cross,
				 model='normal', method='em', n.perm=n_perm,
				 pheno.col=sample(1:nphe,n_etrait,replace=FALSE)
                        ))

	return(permutation)
}

