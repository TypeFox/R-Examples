#' Plot MCMC histogram
#'
#' @usage \method{plot}{admb_hist}(x,type=c("lattice","ggplot"),dtype=c("hist","density"),pars,...)
#' @param x plotting data
#' @param type only "lattice" at present 
#' @param dtype either "hist" or "density"
#' @param pars passed to rhist
#' @param \dots additional parameters for compatibility
#' @return plot object 
#' @export
#' @importFrom lattice xyplot

## don't know how to structure this properly.
##  I would like to have plot.admb_hist plot
##  a graph (as a side effect) and invisibly return a
##  data frame restructured for plotting,
##  but that causes problems when incorporating lattice and ggplot
##  plots into Sweave documents: one would normally say
##  print(plotfun(...)), but that doesn't work when plotfun()
##  doesn't actually return the plot (I think)
## I broke the restructuring function out so it can be
## used separately, but at the moment am still 
## returning the restructured data.  A possible workaround for
## ggplot() plots:
##  plotfun(...); print(last_plot())

plot.admb_hist <- function(x,type=c("lattice","ggplot"),
		dtype=c("hist","density"),
		pars,
		...) { ## dots for generic compat
	type <- match.arg(type)
	dtype <- match.arg(dtype)
	if (dtype=="hist") {
		xx <- rhist(x$hists,pars)
		if (type=="ggplot") {
			stop("ggplot disabled to avoid dependency")
			## if (!require(ggplot2)) stop("must install ggplot2 package")
			## X1 <- ""; X2 <- ""  ## hack to circumvent NOTE in R CMD check
			## vplot <- ggplot2::ggplot(xx,aes(x=X1,y=X2))+
			##   geom_step()+
			##     ## geom_bar(stat="identity",fill="darkgray")+
			##     facet_wrap(~param,scales="free")+
			##       labs(y="Frequency",x="")
		} else if (type=="lattice") {
			##barchart(X2~X1|param,data=xx,horiz=FALSE,
			vplot <- xyplot(X2~X1|param,type="s",data=xx,
					scales=list(x=list(relation="free",tick.number=3),
							y=list(relation="free")),
					as.table=TRUE)
		}
	}
	vplot
	## invisible(xx)
}

rhist <- function(x,pars) {
	nbars <-sapply(x,nrow)
	xx <- data.frame(do.call(rbind,x),
			param=rep(names(x),nbars))
	if (!missing(pars)) {
		if (is.numeric(pars))
			pars <- levels(xx$param)[pars]
		xx <- subset(xx,as.character(xx$param) %in% pars)
	}
	xx
}
