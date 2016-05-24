#' Transform a Topology into a Set of Constraint Commands for MrBayes
#'
#' Takes a phylogeny in the form of an object of class \code{phylo} and
#' outputs a set of topological constraints for \emph{MrBayes} as a set of character
#' strings, either printed in the R console or in a named text file, which
#' can be used as commands in the \emph{MrBayes} block of a NEXUS file for use with 
#' (you guessed it!) \emph{MrBayes}.

#' @details
#' \code{partial = TRUE} may be useful if the reason for using
#' \code{createMrBayesConstraints} is to constrain a topology containing
#' some of the taxa in an analysis, while allowing other taxa to freely
#' vary. For example, Slater (2013) constrained an analysis so extant
#' taxon relationships were held constant, using a molecular-based topology, 
#' while allowing fossil taxa to freely vary relative to their morphological
#' character data. 

#' @param tree An object of class \code{phylo}.

#' @param partial If \code{TRUE} (the default), then constraints will
#' be defined as partial constraints with respect to the rest of the taxa
#' in the input \code{tree}. If \code{FALSE}, constraints will be
#' defined as hard clade membership constraints (i.e. no additional
#' taxa will be allowed to belong to the nodes present in the observed tree.
#' Depending on your analysis, \code{partial = TRUE} may require additionally
#' defining an outgroup.

#' @param file Filename (possibly with path) as a character string
#' to a file which will be overwritten with the output constraint lines.
#' If not null, not constraint lines are output to the console.

#' @return
#' If argument \code{file} is \code{NULL}, then the constrain commands
#' are ouput as a series of character strings.

#' @author
#' David W. Bapst, with some inspiration from Graham Slater.

#' @references
#' Slater, G. J. 2013. Phylogenetic evidence for a shift in the mode of mammalian
#' body size evolution at the Cretaceous-Palaeogene boundary.
#' \emph{Methods in Ecology and Evolution} 4(8):734-744.

#' @examples
#' set.seed(444)
#' tree<-rtree(10)
#' createMrBayesConstraints(tree)
#' createMrBayesConstraints(tree,partial=FALSE)
#' 
#' \dontrun{
#' 
#' createMrBayesConstraints(tree,file="topoConstraints.txt")
#' 
#' }

#' @name createMrBayesConstraints
#' @rdname createMrBayesConstraints
#' @export
createMrBayesConstraints<-function(tree,partial=TRUE,file=NULL){
	#checks
	if(!inherits(tree,"phylo")){
		stop("tree must be of class 'phylo'")
		}
	########################################################
	#get the splits
	#
	splits<-prop.part(tree)	#prop.part
	#remove split with all taxa
	splits<-splits[sapply(splits,length)!=Ntip(tree)]
	#
	# seperate into 'includes' and 'excludes'
	includes<-lapply(splits,function(x) tree$tip.label[x])
	excludes<-lapply(splits,function(x) tree$tip.label[-x])
	#
	#turn into character strings
	includes<-sapply(includes,function(x) paste0(x,collapse=" "))
	excludes<-sapply(excludes,function(x) paste0(x,collapse=" "))
	#
	########################################################
	if(partial){	# if they should be formulated as partial constraints
		# convert into partial constraint commands
		constraints<-sapply(1:length(includes),function(i)
			paste0("constraint node",i," partial = ",
			includes[i]," : ",excludes[i],";"))
	}else{
		# convert into monophyletic constraint commands
		constraints<-sapply(1:length(includes),function(i)
			paste0("constraint node",i," = ",includes[i],";"))
		}
	#
	#########################################################
	# create prset line
	constList<-paste0("node",1:length(splits),collapse=",")
	prsetLine<-paste0("prset topologypr = constraints(",
		constList,");")
	#
	############################################################
	#
	# create final text block
	finalText<-c(constraints,"",prsetLine)
	#
	if(!is.null(file)){
		write(finalText,file)
	}else{
		return(finalText)
		}
	}

