## $Id: generics.R 4896 2014-02-21 18:07:14Z bpikouni $

setGeneric("fit",
           function(data, type, ...) {
             standardGeneric("fit")
           })

setGeneric("pointGraph",
           function(data, ...) {
             standardGeneric("pointGraph")
           })

setGeneric("profileGraph",
           function(data, ...) {
             standardGeneric("profileGraph")
           })

setGeneric("diffGraph",
           function(data,...) {
             standardGeneric("diffGraph")
           })

setGeneric("descriptiveTable",
           function(data, display="print", ...) {
             standardGeneric("descriptiveTable")
           })

setGeneric("correlationTable",
           function(data, display="print", ...) {
             standardGeneric("correlationTable")
           })

setGeneric("globalTest",
           function(fit, display="print", ...) {
             standardGeneric("globalTest")
           })

setGeneric("comparisonsTable",
           function(fit, 
                    type="pairwisereflect", 
                    alpha=0.05, addpct=FALSE, display="print", ...) {
             standardGeneric("comparisonsTable")
           })

setGeneric("grpSummaryTable",
           function(fit, mcadjust=FALSE, alpha=0.05, display="print", ...) {
             standardGeneric("grpSummaryTable")
           })


setGeneric("varianceTable",
           function(fit, display="print", ...) {
             standardGeneric("varianceTable")
           })

setGeneric("errorBarGraph",
           function(fit, mcadjust=FALSE,
             alpha=0.05, cgtheme=TRUE, device="single", ...) {
             standardGeneric("errorBarGraph")
           })

setGeneric("comparisonsGraph",
           function(compstable, cgtheme=TRUE, device="single",
                    wraplength=20, cex.comps=0.7, ...) {
             standardGeneric("comparisonsGraph")
           })

setGeneric("varianceGraph",
           function(fit, trend=NULL, cgtheme=TRUE,  device="single", ...) {
             standardGeneric("varianceGraph")
           })

setGeneric("qqGraph",
           function(fit, line=TRUE, cgtheme=TRUE,  device="single", ...) {
             standardGeneric("qqGraph")
           })

setGeneric("downweightedTable",
           function(fit, cutoffwt, display="print",  ...) {
             standardGeneric("downweightedTable")
           })

setGeneric("kmGraph",
           function(data, cgtheme=TRUE, distfcn="survival",
                    ylab=NULL, title=NULL, ...) {
             standardGeneric("kmGraph")
           })

setGeneric("samplesizeTable",
           function(fit, direction, mmdvec, power=0.80, alpha=0.05,
                    nmax=1000, display="print", ...) {
             standardGeneric("samplesizeTable")
           })


setGeneric("samplesizeGraph",
           function(sstable, Nscale = "log", mmdscale = "log",
                    ...) {
             standardGeneric("samplesizeGraph")
           })


setGeneric("showObj",
           function(object) {
             standardGeneric("showObj")
           })
                                           
##############################

if(FALSE) {

setGeneric("marginalProfile",
           function(data,...) {
             standardGeneric("marginalProfile")
           })

setGeneric("assessSplineFit",
           function(fit, ...) {
             standardGeneric("assessSplineFit")
           })

}
		   
