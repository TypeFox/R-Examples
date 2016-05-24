##' Analysis of Multivariate Events 
##' 
##' Implementation of various statistical models for multivariate
##' event history data. Including multivariate cumulative incidence models,
##' and bivariate random effects probit models (Liability models) 
##'
##' @name mets-package
##' @docType package
##' @author Klaus K. Holst and Thomas Scheike
##' @useDynLib mets
##' @import stats lava splines timereg Rcpp
##' @importFrom survival Surv is.Surv
##' @keywords package
##' @examples
##' 
##' ## To appear
##' 
NULL

##' np data set
##'
##' @name np
##' @docType data
##' @keywords data
##' @source Simulated data
NULL

##' Migraine data
##'
##' @name migr
##' @docType data
##' @keywords data
NULL

##' Dermal ridges data (families)
##'
##' Data on dermal ridge counts in left and right hand in (nuclear) families
##' @name dermalridges
##' @docType data
##' @keywords data
##' @format Data on 50 families with ridge counts in left and right
##' hand for moter, father and each child. Family id in 'family' and
##' gender and child number in 'sex' and 'child'.
##' @source Sarah B. Holt (1952). Genetics of dermal ridges: bilateral
##' asymmetry in finger ridge-counts.  Annals of Eugenics 17 (1),
##' pp.211--231. DOI: 10.1111/j.1469-1809.1952.tb02513.x
##' @examples
##' data(dermalridges)
##' fast.reshape(dermalridges,id="family",varying=c("child.left","child.right","sex"))
NULL

##' Dermal ridges data (monozygotic twins)
##'
##' Data on dermal ridge counts in left and right hand in (nuclear) families
##' @name dermalridgesMZ
##' @docType data
##' @keywords data
##' @format Data on dermal ridge counts (left and right hand) in 18
##' monozygotic twin pairs.
##' @source Sarah B. Holt (1952). Genetics of dermal ridges: bilateral
##' asymmetry in finger ridge-counts.  Annals of Eugenics 17 (1),
##' pp.211--231. DOI: 10.1111/j.1469-1809.1952.tb02513.x
##' @examples
##' data(dermalridgesMZ)
##' fast.reshape(dermalridgesMZ,id="id",varying=c("left","right"))
NULL

##' Menarche data set
##'
##' @name mena
##' @docType data
##' @keywords data
##' @source Simulated data
NULL

##' Multivariate Cumulative Incidence Function example data set
##'
##' @name multcif
##' @docType data
##' @keywords data
##' @source Simulated data
NULL

##' Stutter data set
##'
##' Based on nation-wide questionnaire answers from 33,317 Danish twins
##' @format
##' tvparnr: twin-pair id
##' zyg: zygosity, MZ:=mz, DZ(same sex):=dz, DZ(opposite sex):=os
##' stutter: stutter status (yes/no)
##' age: age
##' nr: number within twin-pair
##' @name twinstut
##' @docType data
##' @keywords data
NULL

##' BMI data set
##'
##' @format
##' Self-reported BMI-values on 11,411 subjects
##' 
##' tvparnr: twin id
##' bmi: BMI (m/kg^2)
##' age: Age
##' gender: (male/female)
##' zyg: zygosity, MZ:=mz, DZ(same sex):=dz, DZ(opposite sex):=os
##' @name twinbmi
##' @docType data
##' @keywords data
NULL

##' Prostate data set
##'
##' @name prt
##' @docType data
##' @keywords data
##' @source Simulated data
NULL

##' For internal use
##' 
##' @title For internal use
##' @name npc
##' @rdname internal
##' @author Klaus K. Holst
##' @keywords utilities
##' @export
##' @aliases plotcr npc nonparcuminc simnordic corsim.prostate
##' alpha2kendall alpha2spear coefmat piecewise.twostage surv.boxarea
##' cluster.index familycluster.index faster.reshape piecewise.data
##' simBinPlack simBinFam simBinFam2 simSurvFam corsim.prostate.random
##' simnordic.random simCox sim pmvn pbvn
##' loglikMVN scoreMVN grouptable jumptimes folds
NULL


