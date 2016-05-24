#' Calculate corrected PSC
#'
#' Given a phylo object and a picante-style community data matrix (sites are rows,
#' species are columns), calculated corrected phylogenetic species clustering
#'
#' @param tree Phylo object
#' @param samp A picante-style community data matrix with sites as rows, and
#' species as columns
#' 
#' @details Returns the inverse of psc as defined in picante
#'
#' @return A data frame of correctly calculated PSC values, with associated species
#' richness and name of all communities in input cdm
#'
#' @export
#'
#' @importFrom ape vcv compute.brlen
#' @importFrom picante prune.sample
#'
#' @references Helmus, M.R., T.J. Bland, C.K. Williams, & A.R. Ives. 2007. Phylogenetic
#' measures of biodiversity. The American Naturalist. 169:E69-E83.
#'
#' @examples
#' #simulate tree with birth-death process
#' tree <- geiger::sim.bdtree(b=0.1, d=0, stop="taxa", n=50)
#'
#' sim.abundances <- round(rlnorm(5000, meanlog=2, sdlog=1)) + 1
#'
#' cdm <- simulateComm(tree, richness.vector=10:25, abundances=sim.abundances)
#'
#' results <- pscCorr(samp=cdm, tree=tree)

pscCorr<-function(samp,tree){
  # Make samp matrix a pa matrix
  samp[samp>0]<-1
  flag=0
  if (is.null(dim(samp))) #if the samp matrix only has one site
  {
    samp<-rbind(samp,samp)
    flag=2
  }

  if(is(tree)[1]=="phylo")
  {
    if(is.null(tree$edge.length)){tree<-ape::compute.brlen(tree, 1)}  #If phylo has no given branch lengths
    tree<-prune.sample(samp,tree)
    # Make sure that the species line up
    samp<-samp[,tree$tip.label]
    # Make a correlation matrix of the species pool phylogeny
    Cmatrix<-ape::vcv(phy=tree,corr=TRUE)
  } else {
    Cmatrix<-tree
    species<-colnames(samp)
    preval<-colSums(samp)/sum(samp)
    species<-species[preval>0]
    Cmatrix<-Cmatrix[species,species]
    samp<-samp[,colnames(Cmatrix)]
  }

  # numbers of locations and species
  SR<-rowSums(samp)
  nlocations<-dim(samp)[1]
  nspecies<-dim(samp)[2]

  ##################################
  # calculate observed PSCs
  #
  PSCs<-NULL

  for(i in 1:nlocations)
  {
    index<-seq(1:nrow(Cmatrix))[samp[i,]==1]	#species present
    n<-length(index)			#number of species present
    if(n>1)
    {    
      C<-Cmatrix[index,index]	#C for individual locations
      diag(C)<--1
      PSC<-1-(sum(apply(C,1,max))/n)   #THIS LINE IS WRONG IN THE PICANTE VERSION
    } else {PSC<-NA}
      PSCs<-c(PSCs,PSC)
  }
    PSCout<-cbind(PSCs,SR)
 if (flag==2)
  {
    PSCout<-PSCout[-2,]
    return(PSCout)  
  } else {
    return(data.frame(PSCout))
  }
}
