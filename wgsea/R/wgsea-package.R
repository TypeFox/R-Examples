

##' Gene set enrichment analysis using Wilcoxon rank tests
##' 
##' Gene set enrichment analysis (GSEA) is typically based on tests derived
##' from the Kolmogorov-Smirnov, which is underpowered and a need for simpler
##' methods has been identified.  The wgsea package contains functions for
##' conducting GSEA using a Wilcoxon test to test for differences in the
##' distribution of p values between SNPs within the gene set under test and a
##' control set of SNPs.
##' 
##' \tabular{ll}{ Package: \tab wgsea\cr Type: \tab Package\cr Version: \tab
##' 1.0\cr Date: \tab 2012-04-18\cr License: \tab GPL\cr }
##' 
##' See the vignette for further details.
##' 
##' @name wgsea-package
##' @aliases wgsea-package wgsea
##' @docType package
##' @author Chris Wallace <chris.wallace@@cimr.cam.ac.uk>
##' @references Irizarry, R. A.; Wang, C.; Zhou, Y. & Speed, T. P. Gene set
##' enrichment analysis made simple. Stat Methods Med Res 2009, 18, 565-575
##' 
##' Heinig, M.; Petretto, E.; Wallace, C.; et al. A trans-acting locus
##' regulates an anti-viral expression network and type 1 diabetes risk. Nature
##' 2010, 467, 460-464
##' @keywords package
##' @examples
##' 
##' vignette(package="wgsea")
##' 
NULL



