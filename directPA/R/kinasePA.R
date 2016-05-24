#' Direction Analysis for Kinases
#'
#' This is a wrapper for runing directPA for kinase perturbation analysis (kinasePA)
#' 
#' @usage kinasePA(Tc, direction, annotation, minSize=5, substrate.method="OSP", 
#' kinase.method="Stouffer", visualize=TRUE, ...)
#' 
#' @param Tc a numeric matrix. The columns are phosphorylation sites and the columns are treatments vs 
#' control statistics.
#' @param direction the direction to be tested for enrichment. Either specified as a degree for 
#' two-dimensional analysis or as contrast (in a triplet) for three-dimensional analysis.
#' @param annotation a list with names correspond to kinases and elements correspond to substrates belong 
#' to each kinase, respectively.
#' @param minSize the size of annotation groups to be considered for calculating enrichment. Groups 
#' that are smaller than the minSize will be removed from the analysis.
#' @param substrate.method the method to be used for integrating statistics across treatments for each 
#' substrate (phosphorylation site). Available methods are Stouffer, OSP, Fisher, and maxP. 
#' Default method is OSP.
#' @param kinase.method the method to be used for integrating statistics of all phosphorylation
#' sites that belongs to a kinase. Available methods are Stouffer, OSP, Fisher, and maxP. 
#' Default method is Stouffer. 
#' @param visualize whether to visualize the plot.
#' @param ... other visualization parameters to pass on. 
#' @return a list that contains directional p-values for each substrate and directional enrichment for 
#' each kinase.
#' @export
#' @examples
#' 
#' # load the phosphoproteomics dataset
#' data(HEK)
#' 
#' # load the kinase-substrate annoations
#' data(PhosphoSite)
#' 
#' # direction pathway analysis in 2-dimensional space. Implemented as rotating by degree 
#' # (1) test combined effect of Torin1 and Rapamycin vs insul both on "down-regulation"
#' # (180 degree to original direction)
#' kPA <- kinasePA(Tc=HEK, direction=pi, annotation=PhosphoSite.mouse)
#' kPA$kinase[order(unlist(kPA$kinase[,1])),][1:20,]
#' # rank substrates on the direciton of interest
#' sort(kPA$substrate.pvalues)[1:20]
#' 
#' # (2) test combined effect of Torin1 and Rapamycin vs insul on "no change and down-regulation"
#' # (135 degree to the original direction) 
#' kPA <- kinasePA(Tc=HEK, direction=pi*3/4, annotation=PhosphoSite.mouse)
#' kPA$kinase[order(unlist(kPA$kinase[,1])),][1:20,]
#' 
#' # (3) test combined effect of Torin1 and Rapamycin vs insul on "down-regulation and no change"
#' # (225 degree to the original direction) 
#' kPA <- kinasePA(Tc=HEK, direction=pi*5/4, annotation=PhosphoSite.mouse)
#' kPA$kinase[order(unlist(kPA$kinase[,1])),][1:20,]
#' 
kinasePA <- function(Tc, direction, annotation, minSize=5, substrate.method="OSP", kinase.method="Stouffer", visualize=TRUE, ...){
   dPA <- directPA(Tc, direction, annotation, minSize, gene.method=substrate.method, path.method=kinase.method, visualize, ...)
   kPA <- list()
   kPA$substrate.pvalues <- dPA$gene.pvalues 
   kPA$kinase <- dPA$gst
   return(kPA)
}
