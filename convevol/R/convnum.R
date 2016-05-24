#'Quantify convergence by the number of convergent events
#'
#'This program takes in a set of taxa that are already suspected to be convergent in a particular area of  morphospace.  It then counts the number of times that a lineage has invaded that region of morphospace.  
#'
#'@param phyl The phylogeny of interest in phylo format
#'@param phendata Phenotypic data for all tips
#'@param convtips A list consisting of the names of all convergent taxa
#'@param ellipse Optional.  An ellipse defining the region of interest, into which groups may or may not converge.
#'@param plot Whether or not to plot a phylomorphospace with lineages that cross into the region of interest highlighted as red arrows.  Default=TRUE
#'@param plotellipse Optional.  The ellipse defining the region of interest in the first two dimensions.  
#'
#'@details This function will construct an ellipse around all convergent taxa.  Then it will reconstruct ancestral states throughout the phylogeny, and use those to determine how many lineages have crossed into this ellipse from the outside.  
#'
#'@return The number of lineages that have crossed into the region of trait space occupied by the convergent taxa.  
#'
#'@import ape geiger MASS phytools
#'
#'@importFrom cluster ellipsoidhull
#'
#'@export
#'
#'@references Maechler, M., Rousseeuw, P., Struyf, A., Hubert, M., Hornik, K.(2013).
#'cluster: Cluster Analysis Basics and Extensions. R package version 1.14.4.
#'
#'Paradis, E., J. Claude, and K. Strimmer (2004) APE: Analyses of phylogenetics
#'and evolution in R langauge. Bioinformatics, 20, 289-290.
#'
#'Revell, L. J. (2012) phytools: An R package for phylogenetic comparative 
#'biology (and other things). Methods Ecol. Evol. 3 217-223.
#'
#'@examples
#'
#'phyl<-rtree(10)
#'phendata<-fastBM(phyl,nsim=2)
#'convtips<-c("t1","t2","t3")
#'answer<-convnum(phyl,phendata,convtips,plot=TRUE,ellipse=NULL,plotellipse=NULL)

convnum <-
function(phyl,phendata,convtips,plot=TRUE,ellipse=NULL,plotellipse=NULL)

{

#Error checking

if (class(phyl) != "phylo") 
	stop("your tree must be class 'phylo.'")

if (nrow(phendata) != length(phyl$tip)) 
	stop("your data matrix must have the same number of rows as tips in the tree.")
    
if (is.list(convtips)==TRUE) {convtips<-unlist(convtips)} 

if (length(convtips)<=ncol(phendata)) 
	stop("You must have fewer variables than putatively convergent taxa")


#The function.  First, defining the area of morphospace by an ellipse, if the ellipse option is not used..

phendata<-as.matrix(phendata)

phendata<-phendata[phyl$tip.label,]

convtaxa<-phendata[convtips ,]

if (is.null(ellipse)) {convell<-ellipsoidhull(convtaxa,0.0001)}
	else {convell<-ellipse}

if (is.null(plotellipse)) {plotell<-ellipsoidhull(convtaxa[,1:2],0.0001)}
	else {plotell<-plotellipse}

allancstates<-apply(phendata, 2, fastAnc, tree=phyl)

alldata<-rbind(phendata,allancstates)



if (plot==TRUE) {

	phyl$node.label<-NULL

	phylomorphospace(phyl,phendata[, 1:2],label="off")

	points(phendata[unlist(convtips), 1:2],col="white")

	plotellipse(plotell)

	}

cross<-0

nbran<-dim(phyl$edge)

crossedges<-c()

i<-1

while (i<=nbran[1]) {

isinanc=FALSE
isindes=FALSE

anc<-phyl$edge[i,1]
des<-phyl$edge[i,2]

ancval<-alldata[anc ,]
desval<-alldata[des ,]

ancdev<-ancval-convell$loc
desdev<-desval-convell$loc

cutoff<-1.05*convell$d2  #cutoff<-(1+tol)*convell$d2

if (t(ancdev)%*%ginv(convell$cov)%*%ancdev<=cutoff) {isinanc=TRUE}else{isinanc=FALSE}

if (t(desdev)%*%ginv(convell$cov)%*%desdev<=cutoff) {isindes=TRUE}else{isindes=FALSE}

if (isinanc!=TRUE & isindes==TRUE) {

cross<-cross+1

crossedges<-c(crossedges,i)

if (plot==TRUE) {
	arrows(ancval[1],ancval[2],desval[1],desval[2],length=0.1,angle=30,code=2,col="red")

	}
}

i<-i+1

}

return(list(cross,convell,plotell))

crossedges

}
