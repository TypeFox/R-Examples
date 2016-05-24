#' Principal Coordinates of Phylogenetic Structure
#' 
#' Function to generate Principal Coordinates of Phylogenetic Structure (PCPS).
#' 
#' The function obtains a matrix containing phylogeny-weighted species composition 
#' (\code{\link{matrix.p}}) and is submitted to principal coordinates analysis (PCoA). 
#' This method generates the principal coordinates of phylogenetic structure 
#' (PCPS) (Duarte, 2011).
#'
#' The function scores.pcps re-scales the correlation values for \code{\link{biplot}} 
#' graphics. The function plot.pcps draws a simple biplot and represent clades as 
#' "spider" graphs (see \code{\link{ordispider}}).
#' 
#' @encoding UTF-8
#' @import SYNCSA
#' @importFrom stats as.formula cor fitted gaussian glm lm quantile summary.lm
#' @importFrom graphics plot points segments text
#' @importFrom vegan ordispider wcmdscale ordilabel vegdist scores
#' @aliases pcps print.pcps summary.pcps scores.pcps plot.pcps
#' @param comm Community data, with species as columns and sampling units as rows. 
#' This matrix can contain either presence/absence or abundance data.
#' @param dist.spp Matrix containing phylogenetic distances between species.
#' @param method Dissimilarity index, as accepted by \code{\link{vegdist}} (Default dist="bray").
#' @param squareroot Logical argument (TRUE or FALSE) to specify if use square root of 
#' dissimilarity index (Default squareroot = TRUE).
#' @param object An object of class pcps.
#' @param x An object of class pcps.
#' @param choices Axes for re-scaling. Choices must have length equal to two (Default choices = c(1,2)).
#' @param display Display text or points for the sampling units.
#' @param groups Factor giving the groups (Clades) for each species.
#' @param showlabel Label the groups by their names in the centroid of the object.
#' @param ... Other parameters for the respective functions.
#' @return \item{P}{Phylogeny-weighted species composition matrix.} \item{values}{The eigenvalues, 
#' relative eigenvalues and cumulative relative eigenvalues.} \item{vectors}{The principal coordinates
#' of phylogenetic structure (PCPS).} \item{correlations}{Correlations between a PCPS axis and 
#' phylogenetically weighted species abundances or frequencies.}
#' @note \strong{IMPORTANT}: The sequence species show up in the community data matrix MUST be the 
#' same as they show up in the phylogenetic distance matrix. See \code{\link{organize.syncsa}}.
#' @author Vanderlei Julio Debastiani <vanderleidebastiani@@yahoo.com.br>
#' @seealso \code{\link{matrix.p}}, \code{\link{wcmdscale}} 
#' @references Duarte, L.S. (2011). Phylogenetic habitat filtering influences forest nucleation 
#' in grasslands. Oikos, 120, 208:215.
#' @keywords PCPS
#' @examples
#' data(flona)
#' res<-pcps(flona$community, flona$phylo)
#' res
#' summary(res)
#' scores(res)
#' plot(res, display = "text", groups = c(rep("Clade-A", 30), rep("Clade-B", 29)))
#'
#' @export
pcps<-function(comm,dist.spp,method="bray",squareroot=TRUE){
	dis<-dist.spp
	P<-matrix.p(comm,dis,notification=FALSE)$matrix.P
	P.dist<-vegan::vegdist(P,method=method)
	if(squareroot==TRUE){
		P.dist<-sqrt(P.dist)
	}
	ordi.P<-vegan::wcmdscale(P.dist,eig=TRUE)
	vectors<-ordi.P$points
	colnames(vectors)<-NULL
	colnames(vectors)<-colnames(vectors,do.NULL=FALSE,prefix="pcps.")
	values<-ordi.P$eig[which((ordi.P$eig>=0)==TRUE)]
	if(length(unique(ordi.P$eig<0))>1){
		warning("Warning: Negative eigenvalues are present in the decomposition result, but only positive eigenvalues were considered",call.=FALSE)
	}
	relative<-values/sum(values)
	cumulative<-as.vector(rep(NA,length(values)))
	for (i in 1:length(values)){
		cumulative[i]<-sum((values/sum(values))[1:i])
	}
	Values<-cbind(values,relative,cumulative)
	colnames(Values)=c("Eigenvalues","Relative_eig","Cumul_eig")
	rownames(Values)=1:length(values)
	n.col<-dim(P)[2]
	n.axis<-dim(vectors)[2]
	correlations<-matrix(NA,nrow=n.col,ncol=n.axis)
	for (i in 1:n.col){
		for (j in 1:n.axis){
			correlations[i,j]<-cor(P[,i],vectors[,j])
		}
	}		
	colnames(correlations)<-colnames(vectors)
	rownames(correlations)<-colnames(P,do.NULL=FALSE,prefix="Uni.")
	Res<-list(call= match.call(), P=P, values=Values, vectors=vectors, correlations=correlations)
	class(Res) <- "pcps"
	return(Res)
}