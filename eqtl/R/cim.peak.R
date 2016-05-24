#####################################################################
#
# cim.peak.R
#
# copyright (c) 2008-3, Ahmid A Khalili
# 
# last modified Jul, 2008
# first written Mar, 2008
# Licensed under the GNU General Public License version 2 (June, 1991)
# 
# Karl W Broman wrote the functions scanone, find.flanking, pull.geno
# parts of R/qtl package
# licensed under the GNU General Public License version 2 (June, 1991)
#
# Part of the R/eqtl package
# Contains: cim.peak
#
######################################################################

######################################################################
#
# cim.peak: genome scan, single genome scan with additive covariate
#               previously defined as QTL.
#
######################################################################

`cim.peak` <-
function( cross, peak )
{

	require(qtl)

	if (length(class(cross)) < 2 || class(cross)[2] != "cross")
		stop("Input should have class \"cross\".")
	if ( ! all(attr(peak,'class',exact=TRUE) %in% c('peak','list')) )
		stop("Input should have class \"peak\".")

	scanone <- get( attr(peak,'scanone',exact=TRUE) )

	if ( ! all(levels(scanone$chr) == names(cross$geno)) )
		stop("Arguments peak or cross misspecified: attributes(peak)$scanone should describe a scanone object which has been performed on cross object")
	if ( ! all(scanone$pos == pseudo.map(cross)) )
		stop("Arguments peak or cross misspecified: attributes(peak)$scanone should describe a scanone object which has been performed on cross object")

	cimlist <- list(un=NA)
	markerdata <- pull.geno(cross)

	for (i in seq(names(peak)) ){

		pheno.col <- grep(paste(names(peak[i])),names(cross$phe))

		if (length(pheno.col) <= 0)
			stop('peak trait should have phenotype in cross')

		covar <- NA

		for (y in seq(length(peak[[i]])) ){
			if (!is.na(peak[[i]][y])){
				pos <- peak[[i]][[y]]$peak.cM
				covar <- c(covar,paste(find.flanking(cross,pos=pos,chr=names(peak[[i]][y]))$close))
			}
		}

		if (length(covar) > 1){
			covar <- covar[-1]

			# INFORMATIVE SCREEN MESSAGE
			# cat('trait:',names(peak[i]),'\tcovar:',covar,'\n')

			covar <- markerdata[,paste(covar)]
			cim <- suppressWarnings(scanone(cross,pheno.col=pheno.col,addcovar=covar))
			attributes(cim)$names <- c('chr','pos',names(peak[i]))
			cimlist <- c(cimlist,list(cim))
		}
			# INFORMATIVE SCREEN MESSAGE
			# else cat('trait:',names(peak[i]),'\tno QTL\n')
	}

	if (length(cimlist) > 1){
		cim <- cimlist[[2]]
		if (length(cimlist)>2){
			for (i in 3:length(cimlist)){
				n <- c(names(cim),names(cimlist[[i]])[3])
				cim <- suppressWarnings(qtl:::c.scanone(cim,cimlist[[i]]))
				attributes(cim)$names <- n
			}
		}
		attributes(cim)$class <- c('scanone','data.frame')
		return(cim)
	} else stop('No QTL define in peak object. No additive covariates to compute composite interval mapping.')
}

