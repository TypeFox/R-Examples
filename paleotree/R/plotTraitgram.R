#' Plot a Traitgram for Continuous Traits
#' 
#' \code{plotTraitgram} plots a traitgram showing the evolution of a continuous trait.
#' If node values are not given (i.e. the data is empirical data collected from tips,
#' rather than simulated data), maximum-likelihood ancestral trait estimation is used
#' to calculate node values. (Ackerly, 2009) given a tree and a set of continuous trait
#' values.
#' 

#' @details By default, this function will use \code{\link{ace}} from the library \code{ape} to
#' reconstruct ancestral traits and confidence intervals using the PIC method, if internal
#' node values (i.e. ancestral node values) are not given.
#' 
#' As with many functions in the paleotree library, absolute time is always
#' decreasing, i.e. the present day is zero.
#' 

#' @param trait A vector of continuous trait values. If the length of \code{trait}
#' is equal to the number of tips and number of nodes, than it is presumed internal
#' node values are given. If the length of \code{trait} is equal to the number of tips
#' of \code{tree} then these will be treated as tip values and ancestral trait reconstruction
#' will be used to reconstruct the missing ancestral values. If \code{trait} is not
#' named, or if internal node values are given, then values will be presumed to be in
#' the order of tip/node numbering in \code{tree$edge}, which for tips is the same as
#' the ordering in \code{tree$tip.label}.

#' @param tree A \code{phylo} object.

#' @param main Main title of traitgram plot.

#' @param conf.int If \code{TRUE} (the default), confidence intervals are plotted.

#' @param lwd The line width used for branches in the figure.

#' @return Return no value, just plot the traitgram.

#' @note 
#' One should probably never do ancestral trait estimation without
#' looking at the confidence intervals, as these reconstructed estimates tend
#' to be very uncertain.

#' @author David W. Bapst

#' @seealso \code{\link{ace}}
#' 
#' Also see the functions \code{traitgram} in the library picante and
#' \code{phenogram} in the library phytools.

#' @references 
#' Ackerly, D. 2009 Conservatism and diversification of plant
#' functional traits: Evolutionary rates versus phylogenetic signal.
#' \emph{Proceedings of the National Academy of Sciences} \bold{106}(Supplement
#' 2):19699--19706.

#' @examples
#' 
#' set.seed(444)
#' tree <- rtree(10)
#' trait <- rTraitCont(tree)
#' 
#' #first, traitgram without conf intervals
#' plotTraitgram(trait,tree,conf.int=FALSE)
#' 
#' #now, with
#' plotTraitgram(trait,tree)
#' #not much confidence, eh?
#' 
#' # plotting simulated data
#'     # with values for ancestral nodes as input
#' trait <- rTraitCont(tree, ancestor=TRUE)
#' plotTraitgram(tree=tree,trait=trait)
#'


#' @export plotTraitgram
plotTraitgram<-function(trait,tree,main="",conf.int=TRUE,lwd=1.5){
	# traitgram plotted using ML ASR from geiger (or ace() from ape if ci=TRUE)
		# or traitgram plotted with ancestral values
	# checks
	if(!inherits(tree,"phylo")){
		stop("tree must be of class 'phylo'")
		}
	# get root time
	if(is.null(tree$root.time)){
		tree$root.time<-max(node.depth.edgelength(tree)[1:Ntip(tree)])
		}
	# get times
	times<-tree$root.time-node.depth.edgelength(tree)
	# are the node values already included?
	if(length(trait)==Ntip(tree)){
		#sort trait, if not sorted already
		if(is.null(names(trait))){
			message("No names for trait data, assuming in same order as tree$tip.label")
		}else{
			trait<-trait[tree$tip.label]
			}
		noNodes<-TRUE
		message("Plotting ancestral reconstructions for nodes.")
	}else{
		noNodes<-FALSE
		if(length(trait)!=(Ntip(tree)+Nnode(tree))){
			stop("length(trait) not equal to number of tips, nor number of tips+nodes")
			}
		message("Plotting internal node values included in trait data")
		message("assuming trait data in same order as tip/node numbering")
		}
	# do ancestral reconstruction if necessary
	if(noNodes){
		asr<-ace(trait,tree,method="pic")
		tr1<-c(trait,asr$ace)
	}else{
		tr1<-trait
		conf.int<-FALSE
		}
	# calculate x lims
	if(conf.int){		
		ci<-asr$CI95
		xlims<-c(min(c(tr1,ci))-0.1,max(c(tr1,ci))+0.1)
	}else{
		xlims<-c(min(tr1)-0.1,max(tr1)+0.1)
		}
	# initial plot
	plot(1,1,type="n",
		xlim=xlims,
		ylim=c(max(times),min(times)),
		xlab="Trait Values",
		ylab="Time (Before Present)",
		main=main)
	# plot branches	
	for(i in 1:nrow(tree$edge)){
		anc<-tree$edge[i,1]
		desc<-tree$edge[i,2]
		lines(c(tr1[anc],tr1[desc]),
			c(times[anc],times[desc]),
			lwd=lwd)
		}
	# plot confidence intrvals		
	if(conf.int){	
		for(i in 1:Nnode(tree)){
			lines(c(ci[i,1],ci[i,2]),
				c(times[i+Ntip(tree)],times[i+Ntip(tree)]),
				lwd=lwd)
			}
		}	
	# done
	}
