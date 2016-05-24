setClassUnion('characterORNULL', c('character', 'NULL'))
setClassUnion('listORNULL', c('list', 'NULL'))
setClassUnion('numericORNULL', c('numeric', 'NULL'))


#' @name mf-class
#' @title Class mf
#' @usage mf$new(nboot, alpha, seed, compare, rng)
#' @description Parent class for package MF data objects.
#' @docType class
#' @section Fields: 
#' \itemize{
#' \item{\code{nboot: } }{numeric value specifying number of samples}
#' \item{\code{alpha: }}{numeric value specifying complement of confidence interval}
#' \item{\code{seed: }}{vector of integers specifying seed for pseudo-random number generator used}
#' \item{\code{compare: }}{vector of character strings naming groups compared}
#' \item{\code{rng: }}{character string naming type of random number generator}
#' }
#' @keywords documentation
#' @family mf
#' @author Marie Vendettuoli \email{marie.c.vendettuoli@@aphis.usda.gov}
mf <- setRefClass('mf', 
		fields = list(nboot = 'numeric', alpha = 'numeric', seed = 
		'integer', compare = 'character', rng = 'character'))


#' @name mfboot-class
#' @title Class mfboot
#' @usage mfboot$new(nboot, alpha, seed, compare, rng, sample, stat, stuff)
#' @description class for data objects produced by MFBoot, contains class mf with 
#' the two additional fields \emph{stat} and \emph{stuff}.
#' @docType class
#' @section Fields:
#' \itemize{
#' \item{\code{nboot: } }{numeric value specifying number of samples}
#' \item{\code{alpha: }}{numeric value specifying complement of confidence interval}
#' \item{\code{seed: }}{vector of integers specifying seed for pseudo-random number generator used}
#' \item{\code{compare: }}{vector of character strings naming groups compared}
#' \item{\code{rng: }}{character string naming type of random number generator}
#' \item{\code{sample: }}{ what is this?}
#' \item{\code{stat:} }{matrix of estimates}
#' }
#' @section Contains:
#' \code{\link{mf-class}}
#' @keywords documentation
#' @family mf
#' @seealso \code{\link{MFBoot}}
#' @author Marie Vendettuoli \email{marie.c.vendettuoli@@aphis.usda.gov}	
mfboot <- setRefClass('mfboot', contains = 'mf', fields = list(stat = 'matrix',
	sample = 'numericORNULL'))

#' @name mfhlboot-class
#' @title Class mfhlboot
#' @usage mfhlboot$new(nboot, alpha, seed, compare, rng, sample, MFstat, HLstat,
#' QDIFstat, QXstat, QYstat)
#' @description class for data objects produced by HLBoot, contains class mf with
#' addtional fields \emph{MFstat, HLstat, QDIFstat, QXstat, QYstat}
#' @docType class
#' @section Fields:
#' \itemize{
#' \item{\code{nboot: } }{numeric value specifying number of samples}
#' \item{\code{alpha: }}{numeric value specifying complement of confidence interval}
#' \item{\code{seed: }}{vector of integers specifying seed for pseudo-random number generator used}
#' \item{\code{compare: }}{vector of character strings naming groups compared}
#' \item{\code{rng: }}{character string naming type of random number generator}
#' \item{\code{sample: }}{ what is this?}
#' \item{\code{MFstat}}{matrix with columns \emph{observed, median, lower, upper} for 
#' Equal Tailed and Highest Density estimates of mitigated fraction (MF)}
#' \item{\code{HLstat}}{matrix with columns \emph{observed, median, lower, upper} for 
#' Equal Tailed and Highest Density estimates of Hodge-Lehmann estimator (HL)}
#' \item{\code{QDIFstat}}{matrix with columns \emph{observed, median, lower, upper} for 
#' estimates of Quartile Differences}
#' \item{\code{QXstat}}{matrix with columns \emph{observed, median, lower, upper} for 
#' quartiles of treatments}
#' \item{\code{QYstat}}{matrix with columns \emph{observed, median, lower, upper} for 
#' quartiles of responses}
#' \item{\code{sample: }}{ what is this?}
#' }
#' @section Contains:
#' \code{\link{mf-class}}
#' @keywords documentation
#' @family mf
#' @seealso \code{\link{HLBoot}}
#' @author Marie Vendettuoli \email{marie.c.vendettuoli@@aphis.usda.gov}	
mfhlboot <- setRefClass('mfhlboot', contains = 'mf', fields = list(MFstat = 'matrix',
	HLstat = 'matrix', QDIFstat = 'matrix', QXstat = 'matrix', QYstat = 'matrix', 
	sample = 'listORNULL'))
	
#' @name mfmp-class
#' @title Class mfmp
#' @usage mfmp$new(ci, x, what, alpha, tdist, df)
#' @description Class mfmp is created from output of funtion MFmp
#' @docType class
#' @section Fields:
#' \itemize{
#' \item{\code{ci:} }{numeric vector of point and interval estimates}
#' \item{\code{x: } }{numeric vector of length three holding data}
#' \item{\code{what: }}{text string describing interval type}
#' \item{\code{alpha: }}{numeric value specifying complement of confidence interval}
#' \item{\code{tdist: }}{Logical indicating if t distribution(TRUE) or gaussian (FALSE)}
#' \item{\code{df: }}{numeric value indicating degrees freedom}
#' }
#' @keywords documentation
#' @family mfmp
#' @author Marie Vendettuoli \email{marie.c.vendettuoli@@aphis.usda.gov}
#' @seealso \code{\link{MFmp}}
mfmp <- setRefClass('mfmp', 
		fields = list(ci = 'numeric', x = 'numeric', what = 'character', alpha = 
		'numeric', tdist = 'logical', df = 'numeric'))

#' @name mfbootcluster-class
#' @title Class mfbootcluster
#' @usage mfbootcluster$new(nboot, alpha, seed, compare, rng, stat, what, excludedClusters,
#' call, sample)
#' @description Class mfbootcluster is created from output of function MFClusBoot
#' @docType class
#' @section Fields:
#' \itemize{
#' \item{\code{nboot: } }{numeric value specifying number of samples}
#' \item{\code{alpha: }}{numeric value specifying complement of confidence interval}
#' \item{\code{seed: }}{vector of integers specifying seed for pseudo-random number generator used}
#' \item{\code{compare: }}{vector of character strings naming groups compared}
#' \item{\code{rng: }}{character string naming type of random number generator}
#' \item{\code{stat: }}{matrix matrix with columns \emph{observed, median, lower, upper} for 
#' estimates}
#' \item{\code{what: }}{character vector naming what was resampled: \emph{clusters}, \emph{units}, \emph{both}}
#' \item{\code{excludedClusters: }}{character vector naming clusters excluded because of missing treatment(s)}
#' \item{\code{call: }}{the call to \code{MFClusBoot}}
#' \item{\code{sample: }}{what is this?}
#' }
#' @section Contains:
#' \code{\link{mf-class}}
#' @keywords documentation
#' @family mf
#' @seealso \code{\link{MFClusBoot}}
#' @author Marie Vendettuoli \email{marie.c.vendettuoli@@aphis.usda.gov}
mfbootcluster <- setRefClass('mfbootcluster', contains = 'mf', fields = list(
	stat = 'matrix', what = 'character', excludedClusters = 'character', call = 'call',
	sample = 'numericORNULL'))
	
#' @name mfcluster-class
#' @title Class mfcluster
#' @usage mfcluster$new(All, bycluster, excludedClusters, call, compare)
#' @description Class mfcluster is created from output of function MFClus
#' @docType class
#' @section Fields:
#' \itemize{
#' \item{\code{All: }}{vector with elements:
#' \itemize{
#' \item{\emph{w }}{Wilcoxon statistic}
#' \item{\emph{u }}{Mann-Whitney statistic}
#' \item{\emph{r }}{mean ridit}
#' \item{\emph{n1 }}{size of group 1}
#' \item{\emph{n2 }}{size of group 2}
#' \item{\emph{mf }}{mitigated fraction}
#' }}
#' \item{\code{byCluster: }}{As for All, by clusters}
#' \item{\code{excludedClusters: }}{character vector naming clusters excluded because of missing treatment}
#' \item{\code{call: }}{the call to \code{MFClus}}
#' \item{\code{compare: }}{character vector naming groups compared}
#' }
#' @keywords documentation
#' @family mfcluster
#' @seealso \code{\link{MFClus}}
#' @author Marie Vendettuoli \email{marie.c.vendettuoli@@aphis.usda.gov}
mfcluster <- setRefClass('mfcluster', fields = list(All = 'data.frame', byCluster = 
	'matrix', excludedClusters = 'characterORNULL', call = 'call', compare = 'character'))
	
#' @name mfcomponents-class
#' @title Class mfcomponents
#' @usage mfcomponents$new(mf, x, y, subj, compare)
#' @description Class mfcomponens is created from output of function MFSubj
#' @docType class
#' @section Fields:
#' \itemize{
#' \item{\code{mf: }}{numeric estimator for mitigated fraction}
#' \item{\code{x: }}{numeric vector containing responses of group 1}
#' \item{\code{y: }}{numeric vector containing responses of group 2}
#' \item{\code{subj: }}{matrix where \code{mf.j} are the subject components}
#' \item{\code{compare: }}{character vector naming groups being compared}
#' }
#' @keywords documentation
#' @family mfcomponents
#' @seealso \code{\link{MFSubj}}
#' @author Marie Vendettuoli \email{marie.c.vendettuoli@@aphis.usda.gov}
mfcomponents <- setRefClass('mfcomponents', fields = list(mf = 'numeric', x = 'numeric',
	y = 'numeric', subj = 'matrix', compare = 'character'))