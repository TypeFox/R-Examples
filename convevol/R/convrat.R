#'Quantify convergence by the ratio of the current to maximum past phenotypic distance
#'
#'convrat quantifies convergence in a number of different ways.  The basic method uses 1-(dtip/dmax) or dmax/dtip, where dtip is the current phenotypic distance between taxa and dmax is the maximum phenotypic distance between the ancestors of those taxa.  This function also scales this measure in a variety of ways.
#'
#'@param phyl The phylogeny of interest in phylo format
#'@param phendata Phenotypic data for all tips
#'@param convtips A list consisting of the names of all convergent taxa
#'
#'@details  C1 = 1-(dtip/dmax).  C2 = dmax-dtip.  C3 is C2 scaled by the total amount of evolution that has occured in the clade descendend from the most recent common ancestor of all convergent tips.  C4 is C2 scaled by the total amount of evolution in the phylogeny.  This program assumes that all monophyletic clades composed entirely of putatively convergent taxa have been reduced to averages or representative taxa.   
#'
#'@return C1, C2, C3, and C4 
#'
#'@import ape geiger MASS phytools  
#'
#'@export
#'
#'@references Paradis, E., J. Claude, and K. Strimmer (2004) APE: Analyses of phylogenetics
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
#'answer<-convrat(phyl,phendata,convtips)


convrat<-function(phyl,phendata,convtips)

{

#Error checking

if (class(phyl) != "phylo") 
	stop("your tree must be class 'phylo.'")

if (nrow(phendata) != length(phyl$tip)) 
	stop("your data matrix must have the same number of rows as tips in the tree.")
    
#Different commands are required if there are only two putatively convergent 
#taxa, versus more than two.

ntax<-length(convtips)

if(ntax==2) {

	tipsdist<-sqrt(sum((phendata[convtips[1] ,]-phendata[convtips[2] ,])^2))

	mxdist<-maxdist(phyl,phendata,convtips[1],convtips[2])

	C1<-1-(tipsdist/mxdist)

	C2<-mxdist-tipsdist

	wholephylchanges<-sum(calcchanges(phyl,phendata))

	C3<-C2/wholephylchanges

	commonanc<-findMRCA(phyl,convtips)

	subtree<-extract.clade(phyl,commonanc)

	subtreephen<-phendata[subtree$tip.label ,]

	subtreechanges<-sum(calcchanges(subtree,subtreephen))

	C4<-C2/subtreechanges

	}

if(ntax>2) {

	C1s<-c()
	C2s<-c()
	C3s<-c()
	C4s<-c()

	for (i in 1:ntax) {

		j<-i+1

		while (j<=ntax) {

			tipsdist<-sqrt(sum((phendata[convtips[i] ,]-phendata[convtips[j] ,])^2))

			mxdist<-maxdist(phyl,phendata,convtips[i],convtips[j])

			C1<-1-(tipsdist/mxdist)

			C2<-mxdist-tipsdist

			wholephylchanges<-sum(calcchanges(phyl,phendata))

			C3<-C2/wholephylchanges

			commonanc<-findMRCA(phyl,convtips)

			subtree<-extract.clade(phyl,commonanc)

			subtreephen<-phendata[subtree$tip.label ,]

			subtreechanges<-sum(calcchanges(subtree,subtreephen))

			C4<-C2/subtreechanges

			C1s<-c(C1s,C1)
			C2s<-c(C2s,C2)
			C3s<-c(C3s,C3)
			C4s<-c(C4s,C4)

			j<-j+1

			}

		}

	C1<-mean(C1s)
	C2<-mean(C2s)
	C3<-mean(C3s)
	C4<-mean(C4s)

	}

answer<-c(C1,C2,C3,C4)
names(answer)<-c("C1","C2","C3","C4")
answer

}