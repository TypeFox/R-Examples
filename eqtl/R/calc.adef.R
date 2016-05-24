#####################################################################
#
# calc.adef.R
#
# copyright (c) 2008-3, Ahmid A Khalili
# 
# last modified Mar, 2012
# first written Mar, 2008
# Licensed under the GNU General Public License version 2 (June, 1991)
# 
# Karl W Broman wrote the effect.plot function part of R/qtl package
# licensed under the GNU General Public License version 2 (June, 1991)
#
# Part of the R/eqtl package
# Contains: calc.adef
#
######################################################################

######################################################################
#
# calc.adef: calculate additive effect at each QTL marker, meaning
#               the phenotype for each genotypic group
#
######################################################################

`calc.adef` <-
function( cross, scanone, peak, round, ...) 
{

	require(qtl)

	if ( length(class(cross)) < 2 || class(cross)[2] != "cross")
		stop("Input should have class \"cross\".")
	if ( !all(attr(scanone,'class',exact=TRUE) %in% c('scanone','data.frame')) )
		stop("Input should have class \"scanone\"")
	if ( !all(attr(peak,'class',exact=TRUE) %in% c('peak','list')) )
		stop("Input should have class \"peak\".")
	if ( attr(peak,'scanone',exact=TRUE) != deparse(substitute(scanone)) )
		stop("Arguments scanone or peak misspecified: peak should have attribute scanone identical to the name of scanone object.\n",
			"peak should be computed from scanone.")
	if ( !all(levels(scanone$chr) == names(cross$geno)) )
		stop("Arguments scanone or cross misspecified: scanone should have been performed on cross object")
	if ( !all(scanone$pos == pseudo.map(cross)) )
		stop("Arguments scanone or cross misspecified: scanone should have been performed on cross object")
	if ( ! 'draws' %in% names(cross$geno[[1]]) ) stop('First running sim.genoprob')
	if ( ! 'prob' %in% names(cross$geno[[1]]) ) stop('First running calc.genoprob')

	bool <- attributes(cross$geno[[1]]$draws) %in% attributes(cross$geno[[1]]$prob)

        # next line added by Karl Broman
        bool[2] <- all(colnames(cross$geno[[1]]$draws) == colnames(cross$geno[[1]]$prob))

	if ( ! all(bool[-1]) )
		stop("cross object error: calc.genoprob and sim.genoprob should be computed with the same parameters ( step, map.function, etc. )")
	if ( "additive.effect" %in% attr(peak,'features',exact=TRUE) ){
		cat("WARNING: The additive effect is already computed. It will be replace by the new one.\n")
		peak <- drop.peakfeat(peak,'additive.effect')
	}

	if ( !missing(round) && round < 0 && !is.numeric(round) )
		stop("Argument round misspecified: round should be an integer >= 0")

	res <- list(un=NA)

	for ( i in seq(length(peak)) ){

		trait <- names(peak[i])
		lodcolumn <- grep(trait,names(scanone))
		lodcolumn <- as.numeric(lodcolumn)-2

		if (length(lodcolumn) > 1) {
			lodcolumn=lodcolumn[1]
			cat("double id:",trait,"\n")
		}

		if ( length(lodcolumn) < 1 ){
			stop(trait," is not a trait found\n")
		}

		# INFORMATIVE SCREEN MESSAGE
		# cat("trait ",trait,"\n")

		resbytrait <- list(un=NA)

		for (y in 1:length(peak[[i]])){

			bool=TRUE

			if(!is.na(peak[[i]][y])){

				# debug hk
				#echo(names(peak[[i]]))

				pic <- peak[[i]][[y]]$mname.peak
				ad<-1;m<-1;pos<-1

				for (j in seq(length(pic))) {

					marker <-as.character(pic[j])
					z <- grep(paste('^',marker,'$',sep=''),row.names(scanone))

					# INFORMATIVE SCREEN MESSAGE
					# cat("add. effect computed :",marker," chr.",names(peak[[i]][y]),"\n")

					adef<-effectplot( cross, pheno.col=lodcolumn, mname1=marker, draw=FALSE )
					pos <- c(pos,scanone$pos[z])
					m <- c(m,marker)
					ad <- c(ad,adef$Means[[2]]-adef$Means[[1]])
				}

				if ( ! missing(round) ) ad <- round(ad,round)
				ae<-data.frame('additive.effect'=ad[-1])
				ae <- cbind(peak[[i]][[y]],ae)
				ae<-list(ae)
			}else {
				ae <- NA
			}
			attributes(ae)$names<-names(peak[[i]][y])
			resbytrait<-c(resbytrait,ae)
		}
		resbytrait <- list(resbytrait[-1])
		attributes(resbytrait)$names <- trait
		res <- c(res,resbytrait)
	}

	res <- res[-1]
	attributes(res)$class <- c('peak',paste(class(res)[1]))
	attributes(res)$features <-  c(attributes(peak)$features,'additive.effect')
	attributes(res)$scanone <- attributes(peak)$scanone
	attributes(res)$lod.th <- attributes(peak)$lod.th
	attributes(res)$si <- attributes(peak)$si
	attributes(res)$window <- attributes(peak)$window

	return(res)
}

