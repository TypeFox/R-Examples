#' Class CrossoverDesign
#' 
#' A S4 class for Crossover designs: CrossoverDesign
#' 
#' 
#' @name CrossoverDesign-class
#' @aliases CrossoverDesign-class CrossoverDesign show,CrossoverDesign-method
#' print,CrossoverDesign-method
#' @docType class
#' @section Slots: \describe{ \item{list("design")}{ Matrix specifying the design. Rows represent periods and columns the subjects.}
#' \item{list("s")}{ Number of sequences. }
#' \item{list("p")}{ Number of periods. }
#' \item{list("v")}{ Number of treatments. }
#' \item{list("model")}{ A numeric specifying the model the design was searched
#' for or -1 if unknown. }
#' \item{list("description")}{ Optional
#' description of design or reference. }
#' \item{list("attr")}{ List with attributes. }
#' \item{list("misc")}{ List with miscellaneous stuff - not used yet. } }
#' @author Kornelius Rohmeyer \email{rohmeyer@@small-projects.de}
#' @keywords graphs
#' @examples
#' 
#' 
#' design <- t(rbind(c(1,1,2,2),
#'                   c(2,2,1,1),
#'                   c(1,1,2,2),
#'                   c(2,2,1,1),
#'                   c(1,2,2,1),
#'                   c(2,1,1,2)))
#'                    
#' new("CrossoverDesign", design)
#' 
#' 
setClass("CrossoverDesign",	
		representation(design="matrix", 
				s="numeric", p="numeric", v="numeric", 
        model="numeric",
        description="character",
				attr="list", 
				misc="list"),
		validity=function(object) validDesign(object))

setMethod("initialize", "CrossoverDesign",
		function(.Object, design, v, model, description="", attr=list(), misc=list()) {			
			if (missing(design)) {			
				stop("Please specify missing design.")
			}						
			.Object@design <- design
      if (missing(v)) {
        v <- length(levels(as.factor(design)))
      }
			.Object@v <- v
      .Object@s <- dim(design)[2]
			.Object@p <- dim(design)[1]
      if (!missing(model)) {
        .Object@model <- model
      }
			.Object@description <- description
			.Object@attr <- attr
			.Object@misc <- misc      
			validObject(.Object)
			return(.Object)
		})

validDesign <- function(object) {
	if (max(length(object@s),length(object@p),length(object@v))>1) return(FALSE)
	return(TRUE)
}

setMethod("print", "CrossoverDesign",
          function(x, ...) {
            callNextMethod(x, ...)
          })

setMethod("show", "CrossoverDesign",
          function(object) {
            # callNextMethod(x, ...)            
            cat(paste(object@description, " (s=", object@s, ", p=",  object@p, ", v=",  object@v ,")\n", sep=""))
            print(object@design)
            cat("\n")
            cat("Av.eff.trt.pair.adj: ", design.efficiency(object@design)$av.eff.trt.pair.adj, "\n")
            if (length(object@model)>0) {
              print(general.carryover(object@design, model=object@model))
            }
          })

#' Class CrossoverSearchResult
#' 
#' A S4 class for the search result for Crossover designs:
#' CrossoverSearchResult
#' 
#' 
#' @name CrossoverSearchResult-class
#' @aliases CrossoverSearchResult-class CrossoverSearchResult
#' show,CrossoverSearchResult-method print,CrossoverSearchResult-method
#' @docType class
#' @section Slots: \describe{ \item{list("design")}{An object of class
#' \code{CrossoverDesign} describing the best design that was
#' found.}
#' \item{list("startDesigns")}{A list of start
#' designs to search from.}
#' \item{list("model")}{A numeric specifying the model the design was searched
#' for or -1 if unknown.}
#' \item{list("eff")}{List, Progress of the
#' algorithm. TODO: Explain further.}
#' \item{list("search")}{List, TODO}
#' \item{list("time")}{Named numeric with the time in seconds the
#' algorithm was searching.}
#' \item{list("misc")}{List - in the moment not used.}}
#' @author Kornelius Rohmeyer \email{rohmeyer@@small-projects.de}
#' @examples
#' 
#' # n=c(100,10) is very small, but it's just an example and should not take much time
#' x <- searchCrossOverDesign(s=9, p=5, v=4, model=4, n=c(100,10))
#' print(x)
#' 
setClass("CrossoverSearchResult",		
		representation(design="CrossoverDesign",
        startDesigns="list",
		    model="numeric",
				eff="list",
				search="list",
				time="numeric",
        misc="list")
)

setMethod("print", "CrossoverSearchResult",
		function(x, ...) {
			callNextMethod(x, ...)
		})

setMethod("show", "CrossoverSearchResult",
		function(object) {
			# callNextMethod(x, ...)
			cat("Crossover Search Result\n")			
			cat("\nDesign:\n")
			print(object@design)
			
		})

#' Plots information about the search algorithm and its process.
#' 
#' Plots information about the search algorithm and its process.
#' 
#' The x-axis corresponds to the consecutive simulation runs and the y-axis to
#' the design criterion \var{E} that depending on the model is either a
#' weighted average of efficiency factors or standardized pairwise variances
#' and described in detail in the vignette of this package.  Also see the
#' vignette for a few examples and a discussion what can be derived from this
#' plots.
#' @name plot
#' @aliases plot,CrossoverSearchResult,missing-method plot
#' @param x Result from searchCrossOverDesign.
#' @param y Missing.
#' @param type Type of plot. Number 1 is
#' more colorful, but number 2 perhaps a bit easier to understand.
#' @param show.jumps If \code{TRUE} vertical lines will show where the specified
#' jumps occured.
#' @return Returns a ggplot object of the plot.
#' @author Kornelius Rohmeyer \email{rohmeyer@@small-projects.de}
#' @examples
#' 
#' \dontrun{
#' x <- searchCrossOverDesign(s=9, p=5, v=4, model=4)
#' plot(x)
#' }
#' 
#' x <- searchCrossOverDesign(s=9, p=5, v=4, model=4, n=c(50,10), jumps=c(10, 10))
#' plot(x, show.jumps=TRUE)
#' plot(x, type=2)
#' 
#' @export
#' @docType methods
#' @rdname plot-methods
setMethod("plot", c(x="CrossoverSearchResult", y="missing") ,
          function(x, y, type=1, show.jumps=FALSE) { #function(x, ) {
            eff <- unlist(x@eff)
            run <- as.factor(rep(1:length(x@eff), each=length(x@eff[[1]])))
            n <- 1:(length(x@eff[[1]])*length(x@eff))
            n2 <- rep(1:length(x@eff[[1]]), times=length(x@eff))
            d <- data.frame(eff=eff, run=run, n=n, n2=n2)
            d <- d[complete.cases(d),]
            n <- x@search$n
            jumps <- x@search$jumps
            if (type==1) {
              plot <- ggplot(d, aes(x=n, y=eff, colour=run)) + geom_point()
              plot <- plot + geom_line(aes(x=n, y=eff, group=run, colour=run))              
              if (show.jumps) plot <- plot + geom_vline(xintercept = 1:((n[1]*n[2])/jumps[2])*jumps[2], colour="grey")              
            } else {
              plot <-ggplot(d, aes(x=n2, y=eff)) + geom_point(colour="#444499", size=1) + facet_wrap( ~ run)
              plot <- plot + geom_line(aes(x=n2, y=eff, group=run, colour=run))
            }            
            plot <- plot + geom_abline(intercept = max(d$eff), slope = 0, colour="#525252")
            plot <- plot + xlab("Simulation run") + ylab("E")
            plot <- plot + theme(axis.text.x = element_text(angle = 90)) + theme(legend.position="none")
            return(plot)
          }
)

setGeneric("getDesign", function(object, ...) standardGeneric("getDesign"))

#' Extract Design from a CrossoverSearchResult
#' 
#' Extract Design from a CrossoverSearchResult
#' 
#' @name getDesign
#' @aliases getDesign,CrossoverSearchResult-method getDesign getDesign,matrix-method getDesign,character-method
#' @param object A searchCrossOverDesign object from which the design should be extracted.
#' @param ... Possible parameters for subclasses (not yet used).
#' @return Returns a numeric matrix representing the crossover design.
#' Rows represent periods, columns represent sequences.
#' @author Kornelius Rohmeyer \email{rohmeyer@@small-projects.de}
#' @examples
#' 
#' # n=c(100,10) is very small, but it's just an example and should not take much time
#' x <- searchCrossOverDesign(s=9, p=5, v=4, model=4, n=c(100,10))
#' getDesign(x)
#' 
#' getDesign("williams4t")
#' 
#' @export
#' @docType methods
#' @rdname getDesign-methods
setMethod("getDesign", c("CrossoverSearchResult"),
          function(object, ...) {			
              return(object@design@design)
          })

setMethod("getDesign", c("matrix"),
          function(object, ...) {  		
            return(object)
          })

setMethod("getDesign", c("character"),
          function(object, ...) {    	
            design <- get(object, envir=Crossover.env)            
            rownames(design) <- paste("p", 1:dim(design)[1], sep="")
            colnames(design) <- paste("s", 1:dim(design)[2], sep="")
            return(design)
          })
