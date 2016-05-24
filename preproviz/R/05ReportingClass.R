
## USES: Analysis

#' A ReportClass is an class containing ggplot2 objects created based on an AnalysisClass object.  
#'
#' @slot barplot
#' @slot heatmap
#' @slot multidimensionalscaling
#' @slot variableclusters
#' @slot outliers
#' @slot varimp
#' @slot lofsum
#' @export 

setClass("ReportClass", representation(barplot="list", heatmap="list", multidimensionalscaling="list",  variableclusters="list", outliers="list", varimp="list", lofsum="list"))

### CHARACTERIZATION PLOTS

#' initializeReportClass is a constructor function for initializing a ReportClass object. 
#' @param object (AnaysisClass)
#' @export

initializeReportClass <- function(object){

# Creating an environment in order to input object name to ggplot object outside data
  
genv <- new.env()    
name <- object@objectname
genv$name <- name
genv$object <- object

## Bar plots
g_bar <- ggplot2::ggplot(getlongformatconstructeddata(object), aes (value), environment=genv) 
g_bar <- g_bar + geom_bar() + facet_wrap(~Var2, scales="free") + theme_bw() + ggtitle(name)

# TEST THIS: theme_set(theme_bw())
# g_scattermatrix <- GGally::ggpairs(getminmaxconstructeddata(object), diag = list(continuous = "density"), lower=list(continuous = "smooth"), upper="blank")

## Heat map
g_heatmap <- ggplot2::ggplot(getlongformatminmaxconstructeddata(object), aes(y=Var1,x=Var2), environment=genv) 
g_heatmap <- g_heatmap + geom_tile(aes(fill=value)) + scale_fill_gradient(low="white", high="black") + theme_bw() + ggtitle(name) + coord_flip()

### CLUSTERING PLOTS

## Multidimensional scaling

g_scatter <- ggplot2::ggplot(getcmdsdata(object), aes(x=X1, y=X2), environment=genv) + geom_point(shape=1)
g_scatter <- g_scatter + theme_bw() + ggtitle(name)

## Hieararchical clustering ## MINOR ISSUE: NOTE ANALYSIS DONE HERE AND NOT IN ANALYSIS CLASS

g_dendro <- suppressWarnings(ggdendro::ggdendrogram(ClustOfVar::hclustvar(getminmaxconstructeddata(object)), rotate = FALSE, size = 2) + labs(title=genv$name))
g_dendro <- g_dendro + coord_flip() 

## LOF Scores

g_outlier <- ggplot2::ggplot(getlofscores(object), aes(x=object.lofscores), environment=genv) 
g_outlier <- g_outlier + geom_density()+ theme_bw() + ggtitle(name)

## Variable importance

g_varimp <- ggplot2::ggplot(getvariableimportancedata(object), aes(y=MeanDecreaseAccuracy, x=features), environment=genv) 
g_varimp <- g_varimp + geom_bar(stat="identity") + coord_flip() + theme_bw() + ggtitle(name)

##

## Lofsum

g_lofsum <- ggplot2::ggplot(getlofsumdata(object), aes(x=seq, y=variable), environment=genv) + geom_tile(aes(fill=value)) 
g_lofsum <- g_lofsum + scale_fill_gradient(low="white", high="black") + theme_bw() + ggtitle(name) + ylab("") + theme(axis.title.y = element_blank())

ReportClass <- new("ReportClass", barplot=list(g_bar), heatmap=list(g_heatmap), multidimensionalscaling=list(g_scatter), variableclusters=list(g_dendro), outliers=list(g_outlier), varimp=list(g_varimp), lofsum=list(g_lofsum))
return(ReportClass)
}

## METHODS

#' plotCMDS  
#'
#' plotCMDS is a generic function for plotting classical multidimensional scaling of constructed features
#' @param object (ReportClass or RunClass)
#' @rdname plotCMDS
#' @export

setGeneric("plotCMDS", function(object) {
  standardGeneric("plotCMDS")
})

#' plotCMDS ReportClass 
#' @describeIn plotCMDS

setMethod("plotCMDS", signature(object = "ReportClass"), function(object) {
  object@multidimensionalscaling}
)

#' plotCMDS RunClass 
#' @describeIn plotCMDS

setMethod("plotCMDS", signature(object = "RunClass"), function(object) {
  listcmds <- lapply(object@reports, function(x) slot(x, "multidimensionalscaling"))
  listcmds <- lapply(listcmds, `[[`, 1)
  do.call(gridExtra::grid.arrange,  listcmds)
  }
)

###

#' plotBAR  
#'
#' plotBAR is a generic function for for plotting barplots of constructed features
#' @param object (ReportClass or RunClass)
#' @rdname plotBAR
#' @export

setGeneric("plotBAR", function(object) {
  standardGeneric("plotBAR")
})

#' plotBAR ReportClass 
#' @describeIn plotCMDS

setMethod("plotBAR", signature(object = "ReportClass"), function(object) {
  object@barplot}
)

#' plotBAR RunClass 
#' @describeIn plotCMDS

setMethod("plotBAR", signature(object = "RunClass"), function(object) {
  listcmds <- lapply(object@reports, function(x) slot(x, "barplot"))
  listcmds <- lapply(listcmds, `[[`, 1)
  do.call(gridExtra::grid.arrange,  listcmds)
}
)

## OUTLIERS

#' plotOUTLIERS 
#'
#' plotOUTLIERS is a generic function for plotting density of LOF scores of constructed features
#' @param object (ReportClass or RunClass)
#' @rdname plotOUTLIERS
#' @export

setGeneric("plotOUTLIERS", function(object) {
  standardGeneric("plotOUTLIERS")
})

#' plotOUTLIERS ReportClass 
#' @describeIn plotOUTLIERS

setMethod("plotOUTLIERS", signature(object = "ReportClass"), function(object) {
  object@outliers}
)

#' plotOUTLIERS RunClass 
#' @describeIn plotOUTLIERS

setMethod("plotOUTLIERS", signature(object = "RunClass"), function(object) {
  listcmds <- lapply(object@reports, function(x) slot(x, "outliers"))
  listcmds <- lapply(listcmds, `[[`, 1)
  do.call(gridExtra::grid.arrange,  listcmds)
}
)

##



#' plotVARCLUST 
#'
#' plotVARCLUST is a generic function for plotting variable clusters of constructed features
#' @param object (ReportClass or RunClass)
#' @rdname plotVARCLUST
#' @export

setGeneric("plotVARCLUST", function(object) {
  standardGeneric("plotVARCLUST")
})

#' plotVARCLUST ReportClass 
#' @describeIn plotVARCLUST

setMethod("plotVARCLUST", signature(object = "ReportClass"), function(object) {
  object@variableclusters}
)

#' plotVARCLUST RunClass 
#' @describeIn plotVARCLUST

setMethod("plotVARCLUST", signature(object = "RunClass"), function(object) {
  listcmds <- lapply(object@reports, function(x) slot(x, "variableclusters"))
  listcmds <- lapply(listcmds, `[[`, 1)
  do.call(gridExtra::grid.arrange,  listcmds)
}
)

### HEATMAP

#' plotHEATMAP 
#'
#' plotHEATMAP is a generic function for plotting heatmap of constructed features
#' @param object (ReportClass or RunClass)
#' @rdname plotHEATMAP
#' @export

setGeneric("plotHEATMAP", function(object) {
  standardGeneric("plotHEATMAP")
})

#' plotHEATMAP ReportClass 
#' @describeIn plotHEATMAP

setMethod("plotHEATMAP", signature(object = "ReportClass"), function(object) {
  object@heatmap}
)

#' plotHEATMAP RunClass 
#' @describeIn plotHEATMAP

setMethod("plotHEATMAP", signature(object = "RunClass"), function(object) {
  listcmds <- lapply(object@reports, function(x) slot(x, "heatmap"))
  listcmds <- lapply(listcmds, `[[`, 1)
  do.call(gridExtra::grid.arrange,  listcmds)
}
)

## VARIMP

#' plotVARIMP 
#'
#' plotHEATMAP is a generic function for variable importance (predicting the class labels in original data) of constructed features
#' @param object (ReportClass or RunClass)
#' @rdname plotVARIMP
#' @export

setGeneric("plotVARIMP", function(object) {
  standardGeneric("plotVARIMP")
})

#' plotVARIMP ReportClass 
#' @describeIn plotVARIMP

setMethod("plotVARIMP", signature(object = "ReportClass"), function(object) {
  object@varimp}
)

#' plotVARIMP RunClass 
#' @describeIn plotVARIMP

setMethod("plotVARIMP", signature(object = "RunClass"), function(object) {
  listcmds <- lapply(object@reports, function(x) slot(x, "varimp"))
  listcmds <- lapply(listcmds, `[[`, 1)
  do.call(gridExtra::grid.arrange,  listcmds)
}
)

## LOFSUM

#' plotLOFSUM 
#'
#' plotLOFSUM is a generic function for lof sum of constructed features
#' @param object (ReportClass or RunClass)
#' @rdname plotLOFSUM
#' @export

setGeneric("plotLOFSUM", function(object) {
  standardGeneric("plotLOFSUM")
})

#' plotLOFSUM ReportClass 
#' @describeIn plotLOFSUM

setMethod("plotLOFSUM", signature(object = "ReportClass"), function(object) {
  object@lofsum}
)

#' plotLOFSUM RunClass 
#' @describeIn plotLOFSUM

setMethod("plotLOFSUM", signature(object = "RunClass"), function(object) {
  listcmds <- lapply(object@reports, function(x) slot(x, "lofsum"))
  listcmds <- lapply(listcmds, `[[`, 1)
  do.call(gridExtra::grid.arrange,  listcmds)
}
)

