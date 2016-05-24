#####################################################################
#
# Rsq.2.array.R
#
# copyright (c) 2008-3, Ahmid A Khalili
# 
# last modified Jul, 2008
# first written Mar, 2008
# Licensed under the GNU General Public License version 2 (June, 1991)
# 
# Part of the R/eqtl package
# Contains: Rsq.2.array
#
######################################################################

######################################################################
#
# Rsq.2.array: Add R square result to QTL result array
#               of class peak.array.
#
######################################################################

`Rsq.2.array` <-
function(rsq,peak.array)
{

	if ( ! all(attr(rsq,'class',exact=TRUE) %in% c('rsq','data.frame')) )
		stop("rsq should have class \"rsq\" \"data.frame\" .")

	if ( !any(names(rsq) %in% c('qtl','rsq','pF')) )
		stop("rsq should be compute by the function calc.Rsq and therefore should have columns: 'qtl','rsq','pF'")
	if ( ! all(attr(peak.array,'class',exact=TRUE) %in% c("peak.array","data.frame")) )
		stop('peak.array should have class \"peak.array\" \"data.frame\". ')

	interaction <- grep(':',rsq$qtl)
	single.Rsq <- rsq$rsq[-interaction]
	Rsq.pF <- rsq$pF[-interaction]

	cat("length of the R square object:",length(rsq$rsq),"\n")
	cat("number of individual qtl R square:",length(single.Rsq),"\n")
	cat("number of detected qtl:",length(peak.array$mname.peak[!is.na(peak.array$mname.peak)]),"\n")

	Rpf <- NA
	Rsq <- NA
	y <- 1
	for(i in 1:nrow(peak.array)){
		if( is.na(peak.array$mname.peak[i])){
			Rpf <- c(Rpf,NA)
			Rsq <- c(Rsq,NA)
			next
		} else {
			Rpf <- c(Rpf,paste(Rsq.pF[y]))
			Rsq <- c(Rsq,paste(single.Rsq[y]))
			y<-y+1
		}
	}

	Rpf <- Rpf[-1]
	Rsq <- Rsq[-1]
	
	array <- cbind(peak.array,Rsq,Rpf)
        attributes(array)$class <- c('peak.array','data.frame')
	return(array)
}

