#'Tests the signifiance of convergent evolution by the ratio of the current to maximum past phenotypic distance
#'
#'convratsig tests the significance of convergence (as quantified by convrat) using evolutionary simulations.
#'
#'@param phyl The phylogeny of interest in phylo format
#'@param phendata Phenotypic data for all tips
#'@param convtips A list consisting of the names of all convergent taxa
#'@param nsim The number of simulations to conduct
#'
#'@details The function simulates evolution via Brownian motion using the input tree and parameters derived from the observed data.  It calculates a convergence metric for each simulation and calculates statistics from the number of times the simulated value exceeds the observed value.  
#'
#'@return The convergence metric of interest (C1, C2, etc...), a cutoff value (the value that the observed measure would have to exceed in order to be considered significant), a P-value for the statistic, and all simulated values. 
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
#'answer<-convratsig(phyl,phendata,convtips,10)


convratsig<-function(phyl,phendata,convtips,nsim)

{

#Error checking

if (class(phyl) != "phylo") 
	stop("your tree must be class 'phylo.'")

if (nrow(phendata) != length(phyl$tip)) 
	stop("your data matrix must have the same number of rows as tips in the tree.")
    
#The function.  First the observed values

obs<-convrat(phyl,phendata,convtips)

#Now come the simulations.

C<-vcv.phylo(phyl)

vcv<-phyl.vcv(phendata,C,0)

k<-1
C1s<-c()
C2s<-c()
C3s<-c()
C4s<-c()

while (k<=nsim) {

	simphendata<-sim.char(phyl,vcv$R,1,model=c("BM"),root=0)

	rsimphendata<-simphendata[, , 1]

	simresults<-convrat(phyl,rsimphendata,convtips)
	
	C1s<-c(C1s,simresults[1])
	C2s<-c(C2s,simresults[2])
	C3s<-c(C3s,simresults[3])
	C4s<-c(C4s,simresults[4])

	k<-k+1

	}

C1greater<-0
C2greater<-0
C3greater<-0
C4greater<-0

for (i in 1:nsim) {

	if (C1s[i]>=obs[1]) {C1greater<-C1greater+1}
	if (C2s[i]>=obs[2]) {C2greater<-C2greater+1}
	if (C3s[i]>=obs[3]) {C3greater<-C3greater+1}
	if (C4s[i]>=obs[4]) {C4greater<-C4greater+1}
 
}

answer<-c(C1greater/(nsim+1),C2greater/(nsim+1),C3greater/(nsim+1),C4greater/(nsim+1))
names(answer)<-c("C1","C2","C3","C4")
answer
}
