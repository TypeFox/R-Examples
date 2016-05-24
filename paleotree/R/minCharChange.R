#' Estimating the Minimum Number of Character Transitions Using Maximum Parsimony
#'
#' \code{minCharChange} is a function which takes a cladogram and a discrete trait and finds the
#' solutions of inferred character states for ancestral nodes that minimizes the number of
#' character state transitions (either gains or losses/reversals) for a given topology and a set of
#' discrete character data. \code{minCharChange} relies on \code{ancPropStateMat}, which is a wrapper
#' for \code{phangorn}'s function \code{ancestral.pars}.

#' @param trait A vector of trait values for a discrete character, preferably named with taxon names
#' identical to the tip labels on the input tree.

#' @param tree A cladogram of type 'phylo'. Any branch lengths are ignored.

#' @param randomMax The maximum number of cladograms examined when searching a large number of solutions
#' consistent with the reconstructed ancestral states from \code{ancestral.pars} with the minimum number
#' of character state transitions. If the number of potential solutions is less than randomMax, then
#' solutions are exhaustively searched.

#' @param maxParsimony If maxParsimony is TRUE (the default) then only solutions with the smallest
#' number of total transitions examined will be returned. Note that since solutions are stochastically
#' 'guessed' at, and the number of possible solutions may not be exhaustively searched, there may have
#' been solutions not examined with a lower number of transitions even if \code{maxParsimony = TRUE}. Regardless,
#' one may want to do \code{maxParsimony = FALSE} if one is interested in whether there are solutions with a
#' smaller number of gains or losses and thus wants to return all solutions.

#' @param printMinResult If TRUE (the default), a summary of the results is printed to the terminal. The
#' information in this summary may be more detailed if the results of the analysis are simpler (i.e. 
#' fewer unique solutions).

#' @param orderedChar If TRUE (not the default), then the character will be reconstructed with a cost (step)
#' matrix of a linear, ordered character. This is not applicable if \code{type = "ACCTRAN"}, as cost matrices cannot
#' be used with ACCTRAN in \code{ancestral.pars}, and an error will be returned if \code{orderedChar = TRUE} but
#' a cost matrix is given, as the only reason to use orderedChar is to produce a cost matrix automatically.

#' @param type The parsimony algorithm applied by \code{ancestral.pars}, which can apply one of two:
#' "MPR" (the default) is a relatively fast algorithm developed by Hamazawa et al. (1995) and Narushima
#' and Hanazawa (1997), which relies on reconstructing the states at each internal node by re-rooting at
#' that node.  "ACCTRAN", the 'accelerated transitions' algorithm (Swofford and Maddison, 1987), favors
#' character reversal over independent gains when there is ambiguity. The "ACCTRAN" option in
#' ancestral.pars avoids repeated rerooting of the tree to search for a smaller set of maximum-parsimony
#' solutions that satisfy the ACCTRAN algorithm, but does so by assigning edge weights.
#' As of phangorn v1.99-12, both of these algorithms apply
#' the Sankoff parsimony algorithm, which allows multifurcations (polytomies).
 
#' @param cost A matrix of the cost (i.e. number of steps) necessary to change between states of the input
#' character trait. If NULL (the
#' default), the character is assumed to be unordered with equal cost to change from any state to another.
#' Cost matrices only impact the "MPR" algorithm; if a cost matrix is given but \code{type = "ACCTRAN"}, an error
#' is issued.

#' @param ambiguity A vector of values which indicate ambiguous (i.e. missing or unknown) character state codings
#' in supplied \code{trait} data. Taxa coded ambiguously as treated as being equally likely to be any state coding.
#' By default, \code{NA} values and "?" symbols are treated as ambiguous character codings, in agreement with behavior 
#' of functions in packages \code{phangorn} and \code{Claddis}. 
#' This argument is designed to mirror an hidden argument with an identical name in function \code{phyDat} in package {phangorn}.

#' @param dropAmbiguity A logical. If \code{TRUE} (which is not the default), all taxa with ambiguous codings as defined
#' by argument \code{ambiguity} will be dropped prior to ancestral nodes being inferred. This may result in too few taxa.

#' @param polySymbol A single symbol which separates alternative states for polymorphic codings; the default symbol is
#' \code{"&"}, following the output by \code{Claddis}'s \code{ReadMorphNexus} function, where polymorphic taxa are indicated
#' by default with a string with state labels separated by an \code{"&"} symbol.
#' For example, a taxon coded as polymorphic for states 1 or 2, would be indicated by the string "1&2".
#' \code{polySymbol} is used to break up these strings and automatically construct a fitting \code{contrast} table
#' for use with this data, including for ambiguous character state codings.
	
#' @param contrast A matrix of type integer with cells of 0 and 1, where each row is labelled with a string value
#' used for indicating character states in \code{trait}, and each column is labelled with the formal state label to
#' be used for assign taxa to particular character states. A value of 1 indicates that the respective coding string for
#' that row should be interpreted as reflecting the character state listed for that column. A coding could reflect
#' multiple states (such as might occur when taxa are polymorphic for some morphological character), so the sums of
#' rows and columns can sum to more than 1. 
#' If \code{contrast} is not \code{NULL} (the default), the arguments will nullify 
#' This argument is designed to mirror an hidden argument with an identical name in function \code{phyDat} in package {phangorn}.
#' This structure is based on the \code{\link{contrasts}} tables used for statistical evaluation of factors.
#' See the \code{phangorn} vignette "Special features of phangorn" for more details
#' on its implementation within \code{phangorn} including an example.
#' See examples below for the construction of an example contrast matrix for character data with polymorphisms, 
#' coded as character data output by \code{Claddis}'s \code{ReadMorphNexus} function, where polymorphic taxa are indicated
#' with a string with state labels separated by an \code{"&"} symbol.

#' @param returnContrast If TRUE, the contrast table used by \code{ancestral.pars} will be output instead for
#' user evaluation that polymorphic symbols and ambiguous states are being parsed correctly.

#' @details
#' The wrapper function \code{ancPropStateMat} simply automates the application of functions
#' \code{ancestral.pars} and \code{phyDat} from \code{phangorn}, along with several additional checks
#' and code to present the result as a matrix, rather than a specialized list. 
#' 
#' Note that although the default \code{orderedChar} argument assumes that multistate characters are unordered,
#' the results of character change will always be reported as gains and losses relative to the numbering of the
#' states in the output \code{transitionSumChanges}, exactly as if they had been ordered. In the case
#' where the character is actually ordered, this may be
#' considered a conservative approach, as using a parsimony algorithm for unordered character states allows fewer
#' gains or losses to be counted on branches where multiple gains and losses are reported. If the character is
#' presumably unordered \emph{and multistate}, however, then the gains and losses division
#' is \emph{arbitrary nonsense} and should be combined to
#' to obtain the total number of character changes.

#' @return
#' By default, \code{ancPropStateMat} returns a matrix, with rows corresponding to the ID numbers of tips and nodes in
#' \code{$edge}, and columns corresponding to character states, with the value representing the proportional
#' weight of that node being that state under the algorithm used (known tip values are always 1).
#' If argument \code{returnContrast} is \code{TRUE} then \code{ancPropStateMat} will instead return the final
#' contrast table used by \code{phyDat} for interpreting character state strings.
#'
#' \code{minCharChange} invisibly returns a list containing the following elements, several of which are printed
#' by default to the console, as controlled by argument \code{printMinResult}:
#'
#' \describe{

#'  \item{\code{message}}{Describes the performance of \code{minCharChange} at searching for a minimum solution.}

#' \item{\code{sumTransitions}}{A vector recording the total number of necessary transitions (sum total of gains
#' and losses/reversal) for each solution; effectively the parsimony cost of each solution.}

#' \item{\code{minTransitions}}{A symmetrical matrix with number of rows and columns equal to the number of
#' character states, with values in each cell indicating the minimum number of transitions from one ancestral
#' state (i.e. the rows) to a descendant state (i.e. the columns), taken across the set of kept solutions
#' (dependent on which are kept as decided by argument \code{maxParsimony}).  Generally guaranteed not to
#' add up to the number of edges contained within the input tree, and thus may not represent any realistic
#' evolutionary scenario but does represent a conservative approach for asking 'what is the smallest possible
#' number of transitions from 0 to 1' or 'smallest possible number of transitions from 1 to 0', independently
#' of each other.}

#' \item{\code{solutionArray}}{A three-dimensional array, where for each solution, we have a matrix with edges
#' for rows and two columns indicating the ancestral and child nodes of that edge, with values indicating the
#' states inferred for those nodes in a particular solution.}

#' \item{\code{transitionArray}}{A labelled three-dimensional array where for each solution we have a symmetrical
#' matrix with number of rows and columns equal to the number of character states, with values in each cell
#' indicating the total number of transitions from one ancestral state (i.e. the rows) to a descendant state
#' (i.e. the columns).}

#' \item{\code{transitionSumChanges}}{Which is a three column matrix with a row for every solution, with the
#' values in the three columns measuring the number of edges (branches) inferred to respectively have gains,
#' no change or losses (i.e. reversals), as calculated relative to the order of character states.}

#' }  

#' @aliases ancPropStateMat

#' @seealso
#' The functions described here are effectively wrapers of \code{phangorn}'s function
#' \code{ancestral.pars}.

#' @author David W. Bapst

#' @references
#' Hanazawa, M., H. Narushima, and N. Minaka. 1995. Generating most parsimonious reconstructions on
#' a tree: A generalization of the Farris-Swofford-Maddison method. Discrete Applied Mathematics
#' 56(2-3):245-265.
#' 
#' Narushima, H., and M. Hanazawa. 1997. A more efficient algorithm for MPR problems in phylogeny.
#' Discrete Applied Mathematics 80(2-3):231-238.
#'
#' Schliep, K. P. 2011. phangorn: phylogenetic analysis in R. \emph{Bioinformatics} 27(4):592-593.
#'
#' Swofford, D. L., and W. P. Maddison. 1987. Reconstructing ancestral character states under
#' Wagner parsimony. Mathematical Biosciences 87(2):199-229.


#' @examples
#' # let's write a quick & dirty ancestral trait plotting function
#' 
#' quickAncPlotter<-function(tree,ancData,cex){
#'	ancCol<-(1:ncol(ancData))+1
#' 	plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="upwards")
#' 	tiplabels(pch=16,pie=ancData[(1:Ntip(tree)),],cex=cex,piecol=ancCol,
#'		col=0)
#' 	nodelabels(pie=ancData[-(1:Ntip(tree)),],cex=cex,piecol=ancCol)	
#' 	}
#'
#' # example with retiolitid graptolite data
#' 
#' data(retiolitinae)
#' 
#' #unordered, MPR
#' ancMPR<-ancPropStateMat(retioTree, trait=retioChar[,2], type="MPR")
#' #unordered, ACCTRAN
#' ancACCTRAN<-ancPropStateMat(retioTree, trait=retioChar[,2], type="ACCTRAN")
#' 
#' #let's compare MPR versus ACCTRAN results
#' layout(1:2)
#' quickAncPlotter(retioTree,ancMPR,cex=0.5)
#' text(x=4,y=5,"type='MPR'",cex=1.5)
#' quickAncPlotter(retioTree,ancACCTRAN,cex=0.5)
#' text(x=5,y=5,"type='ACCTRAN'",cex=1.5)
#' 
#' minCharChange(retioTree,trait=retioChar[,2],type="MPR")
#' minCharChange(retioTree,trait=retioChar[,2],type="ACCTRAN")
#' 
#' # with simulated data
#' 
#' set.seed(444)
#' tree<-rtree(50)
#' #simulate under a likelihood model
#' char<-rTraitDisc(tree,k=3,rate=0.7)
#' tree$edge.length<-NULL
#' tree<-ladderize(tree)
#' 
#' #unordered, MPR
#' ancMPR<-ancPropStateMat(tree, trait=char, type="MPR")
#' #unordered, ACCTRAN
#' ancACCTRAN<-ancPropStateMat(tree, trait=char, type="ACCTRAN")
#' #ordered, MPR
#' ancMPRord<-ancPropStateMat(tree, trait=char, orderedChar=TRUE, type="MPR")
#' 
#' #let's compare MPR versus ACCTRAN results
#' layout(1:2)
#' quickAncPlotter(tree,ancMPR,cex=0.3)
#' text(x=8,y=15,"type='MPR'",cex=1.5)
#' quickAncPlotter(tree,ancACCTRAN,cex=0.3)
#' text(x=9,y=15,"type='ACCTRAN'",cex=1.5)
#' #MPR has much more uncertainty in node estimates
#' 	#but that doesn't mean ACCTRAN is preferable
#'
#' #let's compare unordered versus ordered under MPR
#' layout(1:2)
#' quickAncPlotter(tree,ancMPR,cex=0.3)
#' text(x=8,y=15,"unordered char\nMPR",cex=1.5)
#' quickAncPlotter(tree,ancMPRord,cex=0.3)
#' text(x=9,y=15,"ordered char\nMPR",cex=1.5)
#' layout(1)
#' 
#' \dontrun{
#' # what ancPropStateMat automates (with lots of checks):
#'
#' require(phangorn)
#' char1<-matrix(char,,1)
#' rownames(char1)<-names(char)
#' #translate into something for phangorn to read
#' char1<-phyDat(char1,type="USER",levels=sort(unique(char1)))
#' x<-ancestral.pars(tree,char1,type="MPR")
#' y<-ancestral.pars(tree,char1,type="ACCTRAN")
#' }
#' 
#' #estimating minimum number of transitions with MPR 
#' minCharChange(tree,trait=char,type="MPR")
#'
#' #and now with ACCTRAN
#' minCharChange(tree,trait=char,type="ACCTRAN")
#'
#' #POLYMORPHISM IN CHARACTER DATA
#' 
#' 
#' # example trait data with a polymorphic taxon
#'      # separated with '&' symbol
#' # similar to polymorphic data output by ReadMorphNexus from package Claddis
#' charPoly<-as.character(c(1,2,NA,0,0,1,"1&2",2,0,NA,0,2,1,1,"1&2"))
#' #simulate a tree with 16 taxa
#' set.seed(444)
#' tree<-rtree(15)
#' tree$edge.length<-NULL
#' tree<-ladderize(tree)
#' names(charPoly)<-tree$tip.label
#' charPoly
#' 
#' # need a contrast matrix that takes this into account
#'     #can build row by row, by hand
#' 
#' #first, build contrast matrix for basic states
#' contrast012<-rbind(c(1,0,0),c(0,1,0),c(0,0,1))
#' colnames(contrast012)<-rownames(contrast012)<-0:2
#' contrast012
#' 
#' #add polymorphic state and NA ambiguity as new rows
#' contrastPoly<-c(0,1,1)
#' contrastNA<-c(1,1,1)
#' contrastNew<-rbind(contrast012,'1&2'=contrastPoly,contrastNA)
#' rownames(contrastNew)[5]<-NA
#' 
#' #let's look at contrast
#' contrastNew
#' 
#' # now try this contrast table we've assembled
#'     # default: unordered, MPR
#' ancPoly<-ancPropStateMat(tree, trait=charPoly, contrast=contrastNew)
#' 
#' # but...!
#' # we can also do it automatically, 
#'     # by default, states with '&' are automatically treated
#'     # as polymorphic character codings by ancPropStateMat
#' ancPolyAuto<-ancPropStateMat(tree, trait=charPoly, polySymbol="&")
#'
#' # but does this match what the table we constructed?
#' ancPropStateMat(tree, trait=charPoly,
#' 		polySymbol="&", returnContrast=TRUE)
#' 
#' # compare to contrastNew above!
#' # only difference should be the default ambiguous
#' 	# character '?' is added to the table
#' 
#' #compare reconstructions
#' layout(1:2)
#' quickAncPlotter(tree,ancPoly,cex=0.5)
#' text(x=3.5,y=1.2,"manually-constructed\ncontrast",cex=1.3)
#' quickAncPlotter(tree,ancPolyAuto,cex=0.5)
#' text(x=3.5,y=1.2,"auto-constructed\ncontrast",cex=1.3)
#' layout(1)
#' 
#' #look pretty similar!
#' 
#' #i.e. the default polySymbol="&", but could be a different symbol
#'      #such as "," or "\"... it can only be *one* symbol, though
#' 
#' # all of this machinery should function just fine in minCharChange
#'		# again, by default polySymbol="&" (included anyway here for kicks)
#' minCharChange(tree, trait=charPoly, polySymbol="&")
#' 



#' @name minCharChange
#' @rdname minCharChange
#' @export
minCharChange<-function(trait, tree, randomMax=10000, maxParsimony=TRUE, orderedChar=FALSE,
		type="MPR", cost=NULL, printMinResult=TRUE,  ambiguity= c(NA, "?"),
		dropAmbiguity=FALSE, polySymbol="&", contrast=NULL){
	#randomMax=100;maxParsimony=TRUE;printMinResult=TRUE;type="MPR";cost=NULL
	#print result gives back a reasonable 
	ancMat<-ancPropStateMat(trait, tree, orderedChar=orderedChar, type=type, cost=cost)
	#num of potential solutions
	taxSol<-apply(ancMat,1,function(x) sum(x>0))	#taxSol = solution length of each taxon
	nSol<-prod(taxSol)
	#supposedly charN (my trait vector to be sampled) can be character, its fine
	charN<-colnames(ancMat)
	if(nSol>randomMax){
		solMat<-t(apply(ancMat,1,function(x) sample(charN[x>0],randomMax,replace=T)))
	}else{	
		#exhaustive search needed
		#first, build matrix of non-changing taxa
		noChange<-which	(taxSol==1)
		solMat<-matrix(sapply(noChange,function(x) charN[ancMat[x,]>0]),,1)
		rownames(solMat)<-noChange
		if(nSol>1){
			for(i in 2:max(taxSol)){
				changers<-which(taxSol==i)
				for(j in changers){
					solMat2<-lapply(charN[ancMat[j,]>0],function(x) rbind(solMat,x))
					solMat1<-solMat2[[1]]
					for(k in 2:length(solMat2)){solMat1<-cbind(solMat1,solMat2[[k]])}
					colnames(solMat1)<-NULL
					rownames(solMat1)<-c(rownames(solMat),j)
					solMat<-solMat1
					}
				}
			}
		solMat<-solMat[order(as.numeric(rownames(solMat))),,drop=FALSE]
		}
	#are all solMats unique? (yes, if TRUE)
	solUnq<-all(!sapply(1:ncol(solMat),function(x) 
		any(sapply((1:ncol(solMat))[-x],function(y) identical(solMat[,x],solMat[,y])))))
	#do I need to stop if not all solutions are unique???
	if(!solUnq){
		if(nSol>randomMax){
			#if random, then okay, I guess you might have non unique solutions
			solDup<-c(FALSE,sapply(2:ncol(solMat),function(x) 
				any(sapply((1:ncol(solMat))[1:(x-1)],function(y) identical(solMat[,x],solMat[,y])))))
			solMat<-solMat[,!solDup]
		}else{
			#if not random, then stop cause something is wrong!
			stop("Not all solutions are unique, as calculated, despite random permutations not used. Please investigate or contact Dave Bapst.")
			}
		}
	#edgeSol is a 3D array, where for each solution, we have a matrix with edges for
		# rows and two columns indicating the ancestral node of that edge
		# and the child node of that edge, with values indicating the states
		# inferred for those nodes in a particular solution
	edgeSol<-array(,dim=c(Nedge(tree),2,ncol(solMat)))
	for(i in 1:ncol(solMat)){
		xSol<-solMat[,i]
		#rearrange as an edge matrix of transitions
		edgeSol[,,i]<-cbind(xSol[sapply(tree$edge[,1],function(x) which(x==names(xSol)))],
			xSol[sapply(tree$edge[,2],function(x) which(x==names(xSol)))])
		}
	#tranMat is a 3D array where for each solution we have a symmetrical matrix
		#equal to the number of character states, with values indicating the total
		#number of transitions from one ancestral state (given as the rows) to
		#a descendant state (given as columns)
	tranMat<-array(,dim=c(length(charN),length(charN),ncol(solMat)))
	rownames(tranMat)<-paste("anc.",colnames(ancMat),sep="")
	colnames(tranMat)<-paste("desc.",colnames(ancMat),sep="")
	#sumTran is the parsimony cost: number of gains+losses
	sumTran<-numeric()	
	for(i in 1:ncol(solMat)){
		edgeTran<-edgeSol[,,i]
		#turn into transition matrix
		tranMat1<-t(sapply(charN,function(x) sapply(charN,function(y) 
			sum(edgeTran[,1]==x & edgeTran[,2]==y))))
		tranMat[,,i]<-tranMat1
		diag(tranMat1)<-0
		sumTran[i]<-sum(tranMat1)
		#rows are the ancestor state, columns are the desc state
		}
	#are all tranMats unique? generally not
	#	unqTran<-sapply(1:length(tranMat),function(x) 
	#		any(sapply((1:length(tranMat))[-x],function(y) identical(tranMat[,,x],tranMat[,,y]))))	
	#hist(sumTran)
	if(nSol==1){
		maxPars<-1
		tranMat<-tranMat[,,1,drop=FALSE]
		edgeSol<-edgeSol[,,1,drop=FALSE]
	}else{
		maxPars<-which(sumTran==min(sumTran))
		if(maxParsimony){
			#select only most parsimonious solutions
			solMat<-solMat[,maxPars,drop=FALSE]
			tranMat<-tranMat[,,maxPars,drop=FALSE]
			sumTran<-sumTran[maxPars]
			}
		}
	#get # of gains and # of losses and # of no-change for each transition matrix
	tranSumChange<-t(sapply(lapply(1:dim(tranMat)[3],function(y) tranMat[,,y]),function(x) 
		c(sum(x[upper.tri(x)]),sum(diag(x)),sum(x[lower.tri(x)]))))
	colnames(tranSumChange)<-c("Gains","NoChanges","Losses")
	#get the minimum solution
	minTran<-apply(tranMat,c(1,2),min)
	#
	funcMess<-c(paste0(nSol," potential solutions under ",type,", ",length(maxPars)," most parsimonious solutions found"),
		ifelse(nSol>randomMax,"Solutions sampled stochastically","Solutions exhaustively checked"))
	if(printMinResult){
		if(length(maxPars)<6){
			print(list(message=funcMess,sumTransitions=sumTran,
				transitionArray=tranMat,minTransitions=minTran))
		}else{
			print(list(message=funcMess,sumTransitions=sumTran,minTransitions=minTran))
			}
		}
	return(invisible(list(message=funcMess,sumTransitions=sumTran,minTransitions=minTran,
		solutionArray=edgeSol,transitionArray=tranMat,transitionSumChanges=tranSumChange))) #
	}



#' @rdname minCharChange
#' @export
ancPropStateMat<-function(trait, tree, orderedChar=FALSE, type="MPR", cost=NULL, ambiguity= c(NA, "?"),
	dropAmbiguity=FALSE, polySymbol="&", contrast=NULL, returnContrast=FALSE){
	#wrapper for phangorn's ancestral.pars that returns a fully labeled matrix indicating
		#the relative frequency of a node being reconstructed under a given state
	#require(phangorn)
	#convert trait to a character vector
	saveNames<-names(trait)
	trait<-as.character(trait)
	names(trait)<-saveNames
	#check trait
	if(!is.vector(trait) | !is.character(trait)){
		stop("trait must be vector of state data for a single character, that can be coerced to type 'character'")}
	#check orderedChar
	if(!is.logical(orderedChar)){
		stop("orderedChar must be a logical class element")
		}
	if(length(orderedChar)!=1){
		stop("orderedChar must be a single logical element")
		}
	#return error if cost is not null and type=ACCTRAN
	if(type=="ACCTRAN" & !is.null(cost)){
		stop("cost matrix is inapplicable if ACCTRAN algorithm is used")}
	#return error if cost is not null and orderedChar=TRUE
	if(orderedChar & !is.null(cost)){
		stop("Cannot treat character as ordered; cost matrix inapplicable under ACCTRAN")}
	#check polySymbol
	if(length(polySymbol)!=1){stop("polySymbol must be length 1, multiple (or zero) polySymbols not allowed")}
	#check names
	if(is.null(names(trait))){
		if(Ntip(tree)!=length(trait)){
			stop("names(trait) missing and length(trait) isn't same as number of tips in tree!")
		}else{
			names(trait)<-tree$tip.label
			message("names(trait) missing \n","trait values will be assigned to taxa exactly as in tree$tip.label")
			}
		}
	if(dropAmbiguity){
		#if dropAmbiguity, drop all taxa with ambiguity
		isAmbig<-sapply(trait,function(x) any(sapply(ambiguity,identical,x)))
		whichAmbig<-names(trait)[isAmbig]
		#make message
		message(paste0("dropping following taxa with ambiguious codings: ",
			whichAmbig,collapse="\n    "))
		#drop from trait
		trait<-trait[!isAmbig]
		#drop from tree
		tree<-drop.tip(phy=tree,tip=whichAmbig)
		if(Ntip(tree)<2){stop("Too few non-ambiguous taxa remain on the tree to continue analysis")}
		}
	#if contrast isn't null, ignore polySymbol and ambiguity
	if(!is.null(contrast)){
		message("contrast supplied, thus ignoring states indicated in arguments ambituity and polySymbol")
	}else{
		#if contrast isn't null..
		#if anything is polymorphic, need to build contrast table
		anyPoly<-any(sapply(unique(c(trait)),grepl,pattern=polySymbol))
		if(anyPoly){
			#Need to build a contrasts matrix
			contrast<-buildContrastPoly(trait=trait,
					polySymbol=polySymbol, ambiguity=ambiguity)
			}
		}
	#now if contrast exists, get trueStates from it, otherwise figure them out relative to ambiguity
	if(is.null(contrast)){
		unqState<-unique(c(trait))
		#identify unique non-ambiguous states
		isAmbigState<-sapply(unqState,function(x) any(sapply(ambiguity,identical,x)))
		trueStates <- sort(unqState[!isAmbigState], na.last = TRUE)
	}else{
		#get trueStates from contrast
		trueStates<-colnames(contrast)
		}
	#basic data structure setup
	char1<-matrix(trait,,1)
	rownames(char1)<-names(trait)
	#translate into something for phangorn to read
	char1<-phyDat(char1,type="USER",levels=trueStates,
		ambiguity=ambiguity,contrast=contrast,compress=FALSE)
	#if ordered
	if(orderedChar){
		if(!is.null(cost)){stop("Do not give cost matrix if you set argument cost = TRUE")}
		#if orderedChar 
		nStates<-length(trueStates)
		cost<-matrix(,nStates,nStates)
		for(i in 1:nStates){for(j in 1:nStates){
			cost[i,j]<-abs(i-j)
			}}
		colnames(cost)<-rownames(cost)<-trueStates
		}
	#get anc states
	anc1<-ancestral.pars(tree,char1,type=type,cost=cost)
	#check to make sure trait data isn't empty
	if(length(anc1[[1]])<1){
		stop("Ancestral reconstruction returned by ancestral.pars is empty, check arguments involving state codings")}
	#pull final contrast table
	contrastTable <- attr(anc1, "contrast")
	dimnames(contrastTable) <- list(attr(anc1, "allLevels"), attr(anc1, "levels"))
	#turn into a col-per-state matrix with each row a node or tip, numbered as in edge
	anc2<-matrix(unlist(anc1),,length(attr(anc1, "levels")),byrow=T)
	#based on conversation with Klaus on 04-17-15
		#will treat output as if it was always ordered exactly as tips and nodes
		#are numbered in $edge; should be as basic as numbering 1:nrow
	rownames(anc2)<-1:nrow(anc2)
	#does that make sense for the tree
	if(nrow(anc2) != (Nnode(tree)+Ntip(tree))){
		stop("ancestral state matrix has wrong number of rows??")}
	#and now name the columns by the levels
	colnames(anc2)<-attributes(anc1)$levels
	if(returnContrast){
		result<-contrastTable
	}else{
		result<-anc2
		}
	return(result)
	}
	
#hidden internal function
buildContrastPoly<-function(trait, polySymbol="&", ambiguity=c(NA, "?")){
	#test if any states with polySymbol
	unqState<-unique(c(trait))
	containPoly<-sapply(unqState,grepl,pattern=polySymbol)
	#identify unique non-poly, non-ambiguous states
	isAmbig<-sapply(unqState,function(x) any(sapply(ambiguity,identical,x)))
	trueStates <- sort(unqState[!containPoly & !isAmbig], na.last = TRUE)
	#WAIT also need to include all partial states from poly codings 06-04-15
		# a state might not be listed 'alone'
		# might only be listed for polymorphic taxon!
	#unique polymorphic codings
	polyState<-sort(unqState[containPoly])
	#break the strings up
	polyStateBr<-strsplit(polyState,polySymbol)	#so polySymbol needs to be length=1
	#get unique true states listed by poly codings
	uniquePolyTrueStates<-unique(unlist(polyStateBr))
	trueStates<-sort(unique(c(trueStates,uniquePolyTrueStates)))
	#now get number of trueStates
	nTrue<-length(trueStates)	
	#
	## now create contrast table
	#create contrastTrue
	contrastTrue<-matrix(0,nTrue,nTrue)
	diag(contrastTrue)<-1
	colnames(contrastTrue)<-trueStates
	rowNamesGo<-trueStates
	#now ambiguous characters
	if(length(ambiguity)>0){
		#now create rows for ambiguous characters
		contrastAmbig<-matrix(1,length(ambiguity),nTrue)
		rowNamesGo<-c(rowNamesGo,ambiguity)
	}else{
		contrastAmbig<-NULL
		}
	#now polymorphic characters
		#break polymorphic states down into contrast rows
	polyMatch<-sapply(polyStateBr,function(x) match(x, trueStates))
	#number of unique polymorphic codings
	nPoly<-length(polyState)
	contrastPoly<-matrix(0,nPoly,nTrue)
	contrastPolyRaw<-numeric(nTrue) #get a bunch of zeroes
	for(i in 1:nPoly){
		newContrast<-contrastPolyRaw
		newContrast[polyMatch]<-1		#fill in with 1s
		contrastPoly[i,]<-newContrast
		}
	rowNamesGo<-c(rowNamesGo, polyState)
	#now combine	
	contrast<-rbind(contrastTrue,contrastAmbig,contrastPoly)
	rownames(contrast)<-rowNamesGo
	#return contrast table
	return(contrast)
	}
	
