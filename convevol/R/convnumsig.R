#'Assess the significance of convergent evolution using simulations and the convnum metric
#'
#'Simulates evolution along a given phylogeny, using parameters derived from observed data, and calculates the convnum metric for each simulation for a set of user-defined taxa.  Then compares the observed convnum value to the simulated values to assess the significance of the observed levels of convergent evolution.  
#'
#'@param phyl The phylogeny of interest in phylo format
#'@param phendata Phenotypic data for all tips
#'@param convtips A list consisting of the names of all convergent taxa
#'@param nsim The number of simulatons to conduct
#'@param ellipse Optional.  An ellipse defining the region of interest, into which groups may or may not converge.
#'@param plot Optional.  Describes whether or not to show phylomorphospaces for all of the simulations.  
#'@param plotellipse Optional.  The ellipse defining the region of interest in the first two dimensions.  
#'
#'@details None
#'
#'@return A list, consisting first of the p-value for the observed convnum, and second of a vector containing all of the simulated convnum values.  Also displays a histogram of all of the simulated convnum values.
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
#'answer<-convnumsig(phyl,phendata,convtips,10,plot=FALSE,ellipse=NULL,plotellipse=NULL)

convnumsig<-function(phyl,phendata,convtips,nsim,ellipse=NULL,plot=FALSE,plotellipse=NULL)

{

#Error checking

if (class(phyl) != "phylo") 
	stop("your tree must be class 'phylo.'")

if (nrow(phendata) != length(phyl$tip)) 
	stop("your data matrix must have the same number of rows as tips in the tree.")
    
if (is.list(convtips)==TRUE) {convtips<-unlist(convtips)} 

if (length(convtips)<=ncol(phendata)) 
	stop("You must have fewer variables than putatively convergent taxa")


#The function.  First the observed value

ob<-convnum(phyl,phendata,convtips,plot=TRUE,ellipse=NULL)

if (is.null(ellipse)) {convell<-ob[[2]]}
else {convell<-ellipse}

if (is.null(plotellipse)) {plotell<-ob[[3]]}
else {plotell<-plotellipse}

#Then the simulations

#ancvals<-multianc(phyl,phendata)

#rootvals<-ancvals[length(phyl$tip.label)+1 ,]

C<-vcv.phylo(phyl)

phendata<-as.matrix(phendata)

vcv<-phyl.vcv(phendata,C,0)

rootrow<-dim(phendata)[1]+1

rootvals<-multianc(phyl,phendata)[rootrow,]

simdata<-sim.char(phyl,vcv$R,nsim,model=c("BM"),root=0)

#simdata<-simdata+rootvals

#And then all of the assessments of those simulations

sobs<-NULL
greater<-0

i<-1

while (i<=nsim) {
	
	simphen<-simdata[, , i]
	simval<-convnum(phyl,simphen,convtips,plot=plot,ellipse=convell,plotellipse=plotell)
	sobs<-c(sobs,simval[[1]])
	if(simval[[1]]>=ob[[1]]) {greater<-greater+1}
	i<-i+1	
	}

hist(sobs)

p<-greater/nsim

all<-list(p,sobs)
}