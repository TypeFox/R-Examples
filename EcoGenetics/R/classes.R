#################################
### ECOGENETICS CLASSES
#################################

#------------------------------------------------------------------------------#
#' eco.correlog-class
#' 
#' @name eco.correlog-class 
#' @slot OUT analysis output
#' @slot IN analysis input data
#' @slot BEAKS breaks
#' @slot CARDINAL number of elements in each class
#' @slot NAMES variables names
#' @slot METHOD analysis method 
#' @slot DISTMETHOD method used in the construction of breaks
#' @slot TEST test method used (bootstrap, permutation)
#' @slot NSIM number of simulations
#' @slot PADJUST P-values adjust method for permutation tests
#' @aliases eco.correlog-class
#' @author Leandro Roser \email{leandroroser@@ege.fcen.uba.ar} 
#' @keywords internal

setClass("eco.correlog", 
         
         representation(OUT = "list",
                        IN = "list",
                        BREAKS = "numeric", 
                        CARDINAL = "numeric", 
                        NAMES = "character",
                        METHOD = "character", 
                        DISTMETHOD = "character",
                        TEST = "character",
                        NSIM = "numeric",
                        PADJUST = "character")
         )


#------------------------------------------------------------------------------#
#' eco.variogram class
#' @name eco.variogram-class 
#' @aliases eco.variogram-class
#' @author Leandro Roser \email{leandroroser@@ege.fcen.uba.ar} 
#' @keywords internal

setClass("eco.variogram", contains = "data.frame")


#------------------------------------------------------------------------------##
#' eco.gsa class
#' @slot METHOD method used in the analysis 
#' @slot OBS observed value when a single variable is tested
#' @slot EXP expected value when a single variable is tested
#' @slot PVAL P-value when a single variable is tested
#' @slot ALTER alternative hypotesis when a single variable is tested
#' @slot NSIM number of simulations
#' @slot MULTI table with observed and expected values, P-values and alternative
#' hypoteses when multiple variables are tested
#' @aliases eco.gsa-class 
#' @author Leandro Roser \email{leandroroser@@ege.fcen.uba.ar} 
#' @keywords internal

setClass("eco.gsa",  
         representation(METHOD = "character",
                        OBS = "numeric",
                        EXP = "numeric",
                        PVAL = "numeric",
                        ADJUST = "character",
                        ALTER = "character",
                        NSIM ="numeric",
                        MULTI = "list")
         )


#------------------------------------------------------------------------------#
#' eco.lsa class
#' @slot OUT results
#' @slot METHOD method used in the analysis 
#' @slot TEST test method used (bootstrap, permutation)
#' @slot NSIM number of simulations
#' @slot PADJUST P-values adjust method for permutation tests
#' @slot COND conditional randomization (logical)
#' @slot XY input coordinates  
#' @aliases eco.lsa-class
#' @author Leandro Roser \email{leandroroser@@ege.fcen.uba.ar}
#' @keywords internal

setClass("eco.lsa",  
         representation(OUT = "list",
                        METHOD = "character",
                        TEST = "character",
                        NSIM ="numeric",
                        PADJ = "character",
                        COND = "logical",
                        XY = "data.frame")
         )


#------------------------------------------------------------------------------#
#' eco.weight class
#' @slot METHOD weights construction method
#' @slot PAR parameters used for the construction of breaks
#' @slot PAR.VAL values of the parameters used for the construction of breaks
#' @slot ROW.SD row standardization (logical)
#' @slot SELF data self-included (logical)
#' @slot W weights list
#' @slot XY input coordinates
#' @slot NONZERO number non-zero connections
#' @slot NONZEROIND number of individuals
#' with non-zero connections (as percentage)
#' @slot AVERAGE average number of connection (as percentage)
#' @aliases eco.weight-class 
#' @author Leandro Roser \email{leandroroser@@ege.fcen.uba.ar}
#' @keywords internal

setClass("eco.weight",  
         
         representation(W = "matrix",
                        XY = "data.frame",
                        METHOD = "character",
                        PAR = "character",
                        PAR.VAL = "numeric",
                        ROW.SD = "logical",
                        SELF ="logical",
                        NONZERO = "numeric",
                        NONZEROIND = "numeric",
                        AVG = "numeric")
         )


#------------------------------------------------------------------------------#
#' eco.lagweight class
#' @slot PAR parameters used for the construction of breaks
#' @slot PAR.VAL values of the parameters used for the construction of breaks
#' @slot ROW.SD row standardization (logical)
#' @slot SELF data self-included (logical)
#' @slot W weights list
#' @slot XY input coordinates
#' @slot MEAN mean class distances
#' @slot LOGMEAN mean of the class distances logarithm
#' @slot CARDINAL number of elements in each class
#' @slot BREAKS breaks
#' @slot METHOD breaks construction method
#' @aliases eco.lagweight-class
#' @author Leandro Roser \email{leandroroser@@ege.fcen.uba.ar}
#' @keywords internal

setClass("eco.lagweight",  
         
         representation(W = "list",
                        XY = "data.frame",
                        PAR = "character",
                        PAR.VAL = "numeric",
                        ROW.SD = "logical",
                        SELF ="logical",
                        CUMMUL = "logical",
                        MEAN = "numeric",
                        LOGMEAN = "numeric",
                        CARDINAL = "numeric",
                        BREAKS = "numeric",
                        METHOD = "character")
         )


#------------------------------------------------------------------------------#
#' eco.mlm-class
#' @name eco.mlm-class
#' @keywords internal
#' @slot MLM mlm results
#' @slot SUMMARY.MLM summary of the results
#' @slot ANOVA.MLM anovas for the results
#' @slot DF1 data frame
#' @slot DF2 data frame
#' @author Leandro Roser \email{leandroroser@@ege.fcen.uba.ar} 
#' @aliases eco.mlm-class

setClass("eco.mlm",
         
         representation( MLM = "list",
                         SUMMARY.MLM = "list", 
                         ANOVA.MLM = "list",
                         PREDICTED = "data.frame",
                         RESIDUALS = "data.frame",
                         DF1 = "dataframeORmatrix",
                         DF2 = "dataframeORmatrix")
         )


#------------------------------------------------------------------------------#
#' eco.mctree-class 
#' @name eco.mctree-class 
#' @keywords internal
#' @slot TREES trees obtained
#' @slot PREDICTIONS predictions of the analysis
#' @slot FREQUENCIES frequencies of individuals per class in nodes
#' @slot DF1 data frame
#' @slot DF2 data frame
#' @author Leandro Roser \email{leandroroser@@ege.fcen.uba.ar}
#' @aliases eco.mctree-class


setClass("eco.mctree",
         
         representation( TREES = "list",
                         CLASSPREDICT = "list", 
                         FREQUENCIES = "list",
                         PREDICTED = "data.frame",
                         RESIDUALS = "data.frame",
                         DF1 = "dataframeORmatrix",
                         DF2 = "dataframeORmatrix") 
         )


#------------------------------------------------------------------------------#
#' eco.detrend class
#' @slot POLY.DEG polynomial degree used in the analysis
#' @slot RES detrended data
#' @slot XY projected coordinates
#' @slot MODEL models selected with the Akaike criterion
#' @slot ANALYSIS object of class "eco.mlm" with 
#' the regression results for each variable 
#' @aliases eco.detrend-class
#' @author Leandro Roser \email{leandroroser@@ege.fcen.uba.ar}
#' @keywords internal

setClass("eco.detrend",  
         
         representation(POLY.DEG = "numeric",
                        RES = "data.frame",
                        XY = "data.frame",
                        MODEL = "list",
                        ANALYSIS ="eco.mlm")
         )

#------------------------------------------------------------------------------#
#' int.multiplot class
#' @author Leandro Roser \email{leandroroser@@ege.fcen.uba.ar}
#' @keywords internal

setClass("int.multiplot")

#------------------------------------------------------------------------------#
#' eco.IBD class
#' @author Leandro Roser \email{leandroroser@@ege.fcen.uba.ar}
#' @keywords internal

setClass("eco.IBD", representation(SP = "list"), contains = "eco.correlog")

#------------------------------------------------------------------------------#
