#####################################################################
#
# cover.peak.R
#
# copyright (c) 2008-3, Ahmid A Khalili
# 
# last modified Jul, 2008
# first written Mar, 2008
# Licensed under the GNU General Public License version 2 (June, 1991)
# 
# Part of the R/eqtl package
# Contains: cover.peak 
#
######################################################################

######################################################################
#
# cover.peak: List QTLs within a genetical region. 
#
######################################################################

`cover.peak` <-
function( peak, pos, chr, pre=0 )
{

	if ( ! all(attr(peak,'class',exact=TRUE) %in% c('peak','list')) )
		stop("Input should have class \"peak\".")

	scanone <- get( attr(peak,'scanone',exact=TRUE) )

	if ( ! is.vector(chr) ) stop("Argument chr misspecified: expecting a vector")
	if ( ! all(chr %in% levels(scanone$chr)) )
		stop("Argument chr misspecified: could not find all chr \"",chr,"\" in scanone$chr")
	if ( length(chr) > length(levels(scanone$chr)) )
		stop("Argument chr misspecified: chr could not be longer than the number of chr. in scanone",length(levels(scanone$chr)) )
	if ( !is.vector(pos) || !is.numeric(pos) || length(pos)>1) stop("pos should be a numeric vector of length 1\n")
	if ( !is.vector(pre) || !is.numeric(pre) || length(pre)>1) stop("pre should be a numeric vector of length 1\n")
	if ( pre > pos ) stop("pre should be <= to pos")

	for ( i in seq(names(peak)) ){
		for (j in 1:length(peak[[i]])) {
       			if ( is.data.frame(peak[[i]][[j]])){
                		feat <- names(peak[[i]][[j]])
                		break
			}
       		}
		if (is.data.frame(peak[[i]][[j]])) break
	}

	cat("peak features:\t",feat,"\n")

	data.qtl <- data.frame(t(rep(NA,length(feat))))
	attributes(data.qtl)$names <- feat

	trait <- NA; chrom <- NA;
	for ( i in seq(names(peak)) ){
		if ( !any(chr %in% attributes(factor(names(peak[[i]])))$levels) )
			stop("chr should be a factor level of names(peak[[i]]) \n")
		for ( y in chr ){
			if ( is.data.frame(peak[[i]][[paste(y)]]) ) {
				for ( z in seq(peak[[i]][[paste(y)]]$mname.peak)){
					if ( peak[[i]][[paste(y)]]$sup.cM[z] >= pos+pre && peak[[i]][[paste(y)]]$inf.cM[z] <= pos-pre ) {
						trait <- c(trait,names(peak[i]))
						chrom <- c(chrom,names(peak[[i]][paste(y)]))
						data.qtl <- rbind(data.qtl,peak[[i]][[paste(y)]][z,])
					} else next 
				}
			} else next
		}
	}

	out <- cbind(trait=trait[-1],chr=chrom[-1],data.qtl[-1,])
	attributes(out)$class <- c('peak.array','data.frame')

	return(out)
}

