#' Example Species Abundances Tables
#'
#' A totally fictional example of species abundance data, for testing functions
#' that require a site-by-taxon table of community ecology data.

#' @name kanto
#' @rdname kanto
#' @aliases kanto

#' @details
#' A classic dataset of ecological data collected by Satoshi and Okido, consisting of
#' individual counts for 54 terrestrial faunal and floral species,
#' fron 23 sites across the mainland Kanto region.
#'
#' Different ontogenetic stages were compounded and recorded by the common name for the
#' first ontogenetic stage, with some inconsistency for species whose earliest stage have
#' only been recently recognized. When separate names are commonly applied to sexual
#' dimorphic forms, these were also combined and a single common name was used.
#'
#' \emph{Note: This data is a totally made-up, satirical homage to
#' a well-known video game series (thus constituting fair-use).}

#' @format 
#' A table of type integer, representing terrestrial fauna and flora abundance counts.

#' @source 
#' Pokemon And All Respective Names are Trademark and Copyright of Nintendo 1996-2015.

#' @keywords datasets

#' @docType data

#' @examples
#'
#' data(kanto)
#'
#' #visualize site abundances as barplots
#' barplotAbund<-function(x){
#' 	x<-x[,colSums(x)>0]
#' 	layout(1:(nrow(x)+1))
#' 	xpar<-par(mar=c(0,7,2,0))
#' 	for(i in 1:(nrow(x)-1)){
#' 		barplot(x[i,],ylab=rownames(x)[i],
#' 			names.arg="")
#' 		}
#' 	barplot(x[nrow(x),],
#' 		ylab=rownames(x)[nrow(x)],las=3)
#' 	par(xpar)
#' 	layout(1)
#' 	mtext("Abundances",side=2,line=3,adj=0.8)
#' 	}
#' 
#' #first five sites
#' kanto5<-kanto[1:5,]
#' barplotAbund(kanto5)
#' 
#' #get pairwise Spearman rho coefficients
#' rhoCoeff<-pairwiseSpearmanRho(kanto,dropAbsent="bothAbsent")
#' 
#' #what are the nearest-neighbor rhos (largest rho correlations)?
#' diag(rhoCoeff)<-NA
#' rhoNearest<-apply(rhoCoeff,1,max,na.rm=TRUE)
#' rhoNearest
#' 
#' # We can see the power plant sample is extremely different from the rest
#' 
#' # measure evenness: Hurlbert's PIE
#' 
#' kantoPIE<-HurlbertPIE(kanto)
#' 
#' # compare to dominance (relative abundance of most abundant taxon)
#' 
#' dominance<-apply(kanto,1,function(x) max(x)/sum(x) )
#' 
#' plot(kantoPIE,dominance)
#' 
#' # relatively strong relationship!
#' 
#' 
#' \dontrun{
#'
#' #get bray-curtis distances
#' library(vegan)
#' bcDist <- vegdist(kanto,method="bray")
#' 
#' #do a PCO on the bray-curtis distances
#' pcoRes <- pcoa(bcDist,correction="lingoes")
#' scores <- pcoRes$vectors
#' #plot the PCO
#' plot(scores,type="n")
#' text(labels=rownames(kanto),scores[,1],scores[,2],cex=0.5)
#' 
#' #the way the power plant the the pokemon tower converge
#' 	# is very suspicious: may be distortion due to a long gradient
#' 
#' #do a DCA instead with vegan's decorana
#' dcaRes<-decorana(kanto)
#' #plot using native vegan functions
#' 	#will show species scores in red
#' plot(dcaRes,cex=0.5)
#' #kind of messy
#' 
#' #show just the sites scores
#' plot(dcaRes,cex=0.5,display="sites")
#' 
#' #show just the species scores
#' plot(dcaRes,cex=0.5,display="species")
#' 
#' #well, that's pretty cool
#' 
#' 
#' #get the nearest neighbor for each site based on pair-wise rho coefficients
#' rhoNeighbor<-apply(rhoCoeff,1,function(x)
#' 	rownames(kanto)[tail(order(x,na.last=NA),1)])
#' 
#' #let's plot the nearest neighbor connections with igraph
#' NNtable<-cbind(rownames(kanto),rhoNeighbor)
#' 
#' # now plot with igraph
#' library(igraph)
#' NNlist <- graph.data.frame(NNtable)
#' plot(NNlist)
#' 
#' #arrows point at the nearest neighbor of each sample
#' 	# based on maximum Spearman rho correlation
#'
#' #testing for differences between groups of sites
#' 
#' #is there a difference between routes and non-routes
#' groups<-rep(0,nrow(kanto))
#' groups[grep(rownames(kanto),pattern="Route")]<-1
#' 
#' #anosim (in vegan)
#' 	#are distances within groups smaller than distances between?
#' #we could also use adonis from vegan instead 
#' library(vegan)
#' 
#' anosim(dat=kanto,grouping=groups)
#' adonis(kanto~factor(groups))
#' #both are very significant
#' 
#' #alternative: using multivariate GLMs in mvabund
#' 
#' library(mvabund)
#' 
#' ft <- manyglm(formula=kanto~factor(groups))
#' anova(ft)
#' #also highly significant!
#' 
#'
#' }
#' 


#'
NULL