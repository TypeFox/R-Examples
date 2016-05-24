#####################################################################
#
# plotRsq.R
#
# copyright (c) 2008-3, Ahmid A Khalili
# 
# last modified Jul, 2008
# first written Mar, 2008
# Licensed under the GNU General Public License version 2 (June, 1991)
# 
# Part of the R/eqtl package
# Contains: plotRsq 
#
######################################################################

######################################################################
#
# plotRsq: histograms of R square values distribution
#
######################################################################

`plotRsq` <-
function(rsq, par=c(2,2), ...)
{

if ( ! any(class(rsq) == "rsq"))
	stop("rsq should have class \"rsq\" \"data.frame\" .")
if ( !is.vector(par) & !is.integer(par) & length(par) != 2 )
	stop("par is the mfrow parameter of the par() function. par should be a vector of 2 integers.")


interaction <- grep(':',rsq$qtl)

if ( ! length(interaction) > 0 ) par(mfrow=c(1,1))
else par(mfrow=par)

hist(   as.numeric(as.vector(rsq$rsq)),
        nclass=100,
        xaxp=c(0,1,10),
        main='',
        ylab='Number of eQTLs and significant eQTL-eQTL interactions',
        xlab='R square values',
	...
)

if ( ! length(interaction) > 0 ) return('No qtl interaction detected')

hist(   as.numeric(as.vector(rsq$rsq[-interaction])),
        nclass=100,
        xaxp=c(0,1,10),
        main='',
        ylab='Number of eQTLs',
        xlab='R square values',
	...
)

hist(   as.numeric(as.vector(rsq$rsq[interaction])),
        nclass=20,
        xaxp=c(0,1,100),
        main='',
        ylab='Number of eQTL-eQTL interactions',
        xlab='R square values',
	...
)

}

