#' @name blocksdesign-package
#' @title Blocks design package
#' @aliases blocksdesign
#' @docType package
#' 
#' @description The \code{blocksdesign} package provides functionality for the construction of nested or crossed block designs for unstructured 
#' treatment sets with arbitrary levels of replication and arbitrary depth of nesting.
#'  
#' @details
#' 
#' Block designs group experimental units into homogeneous blocks to provide maximum precision for treatment comparisons. 
#' The most basic type of block design is the complete randomized blocks design where each block contains one or more complete
#' sets of treatments. Complete randomized blocks must contain numbers of treatment plots proportional to the treatment replication 
#' so the maximum possible number of complete randomized blocks is the highest common factor (hcf) of the replication numbers.
#' 
#' Complete randomized block designs give fully efficient 'within-blocks' treatment estimates and are usually the best choice for small experiments.
#' For larger experiments, however, the average variability within blocks may be too large for good precision of estimation and then it may 
#' be desirable to sub-divide the large complete blocks into smaller incomplete blocks. Sometimes it is advantageous to use a double 
#' blocking system in which one set of blocks, usually called row blocks, are crossed with a second set of blocks, usually called column blocks. Double
#' blocking systems can be valuable for controlling trends in two dimensions simultaneously.
#' 
#' The analysis of incomplete block designs is complex 
#' and requires the combination of treatment information from different strata using appropriately weighted treatment estimates. Before the 
#' advent of modern electronic computing, the computational complexity of the analysis meant that only a single level of nesting was feasible for 
#' a practicable incomplete block designs. Nowadays, however, the availability of modern computers and modern software such as the \code{lme4} mixed model   
#' package (Bates et al 2014) have largely eliminated such restrictions and the analysis of nested block designs with any feasible depth of
#' nesting is now routine. 
#' 
#' The \code{blocksdesign} package provides for the construction of general block designs where treatments can have any arbitrary
#' levels of replication and simple or crossed blocks can be nested repeatedly to any feasible depth of nesting. The advantage of block designs with more 
#' than a single level of nesting is that random variability can be controlled by a range of block sizes that provide for the control of 
#' positional effects or correlated plot effects at different scales of variability. The design algorithm optimizes the nested blocks of each
#' stratum hierarchically with each successive set of nested blocks optimized within the blocks of the preceding stratum. Block sizes within 
#' each stratum are as equal as possible and never differ by more than a single plot.
#'
#' The main design function is \code{\link[blocksdesign]{blocks}} which is used to generate the actual required design. The output from 
#' \code{blocks} includes a data frame of the block and treatment factors for each plot, a set of plans showing the allocation of 
#' treatments to block factors in the bottom stratum of the design and a  table showing the achieved D- and A-efficiency factors
#' for each blocks stratum together with A-efficiency upper bounds, where available. 
#'  
#' The subsidiary function \code{\link[blocksdesign]{upper_bounds}} estimates A-efficiency upper bounds for regular block designs with equally replicated 
#' treatments and equal block sizes. 
#'
#' Further discussion of designs with repeatedly nested strata can be found in the package vignette at: vignette("blocksdesign")
#' 
#' @seealso \url{http://www.expdesigns.co.uk}
#' 
#' @references
#' 
#' Bates, D., Maechler, M., Bolker, B. and Walker, S. (2014). lme4: Linear mixed-effects models using Eigen and S4. 
#' R package version 1.1-6. http://CRAN.R-project.org/package=lme4
#' 
NULL


