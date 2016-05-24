###############################################################################
#
# CLASS: testparams
#
# this object can be passed to the function coalsimC after having set parameter 
# values. This class eases the process of passing on all necessary values
# to the coalsim function. Required values are nsam, niter and theta. 
#
#	SLOTS:
#		nsam:		is the number of copies of the locus in each sample. It needs 
#					to be provided as a vector of length nloci 
#					length( c(x,y,z) ) = nloci
#
#		niter:		number of independant samples to generate for each locus
#					single integer value greater than 0
#
#		theta:		mutation parameter theta. It needs 
#					to be provided as a vector of length nloci 
#					length( c(x,y,z) ) = nloci
#
#		nloci:		number of loci, single integer value greater than 0
#
#		npop:		number of populations from which observed values originate
#					single integer value greater than 0
#
#		nsites:		number of nucleotid sites for each locus, if provided, 
#					must be a vector of length nloci.
#
#		obsVal:		a matrix with observed values to test against.
#					needs to be provided as a matrix of size numTests x nloci
#					each row consisting of one locus and each column denote the
#					value for a specific test. The order of the tests are specified 
#					in the vector testNames (below)
#		
#		printtree:	set to 1 to include tree information in class, which can later 
#					be printed using the genetree function
#					printing tree is not available with recombination and geneConv
#
#	fixedSegsites:	usually the number of segrating site varies in each iteration. 
#					Please provide a  single numeric value if the number of 
#					segregating sites needs to be fixed.
#
#	recombination:	provide a vector of format: c(p, nsites) 
#					p = cross over parameter rate, nsites is the number of sites
#					between recombination occurs
#	
#		geneConv:	in addition to recombination intra-locus non-cross-over 
#					exchange gene conversion can be included in simulation
# 					expected format is c(f, gamma) 
# 					f denote the ratio, g/r, where r is the probability per generation 
#					of crossing-over between adjacent sites. (see Wiuf and Hein 2000)
# 					gamma is the mean conversion tract length			
#
#		growth:		population size is measured by N(t) = N0 exp-^(alpha*t). provide alpha 
#					as integer value. negative values indicate that population was larger 
#					in the past than present
#
#	demography:		vector of length 3 or 4 with first value denoted as 'type' 
#
#					valid 'types' for vectors of length 3 are as following: 
#					- 1 to set a growth rate change alpha at a certain time t: 
#					  c(1, t, alpha)
#
#					- 2	set all subpop to size x * N_0 and growth rate to zero: 
#					  c(2, t, x)
#
#					- 3 set all elements of migration matrix to x/(npop-1): 
#					  c(3, t, x)					
#
#					valid 'types' for vector of length 4 with the following values:	
#					- 4 set growth rate of subpop i to alpha at time z: 
#					  c(4, t, i, alpha)
#
#					- 5 set subpop i size to x * N_0 at time t and growth rate to zero: 
#					  c(5, t, i, x)
#
#					- 6 split subpopulation i into subpopulation i and a new subpopulation, 
#					  labeled npop + 1. Each ancestral lineage in subpopulation i is randomly  
#					  assigned to subpopulation i with probability p and subpopulation  
#					  npop + 1 with probability 1 - p. The size of subpopulation npop + 1 is  
#					  set to N_0. Migration rates to and from the new subpopulation are assumed  
#					  to be zero and the growth rate of the new subpopulation is set to zero: 
#					  c(6, t, i, p)
#
#					- 7 move all lineages in subpopulation i to subpopulation j at time t. 
#					  Migration rates from subpopulation i are set to zero: 
#					  c(7, t, i, j)
#
# AUTHOR:			Niroshan Nadarajah <niroshan.nadarajah@uni-duesseldorf.de>
#
# Last modified:	10/11/04
#
###############################################################################
setClass("test.params", representation(
          
	#Slots of testparams

	n.sam 			= "numeric",	# specify samples for each loci length(c(...)) = nloci
	n.iter			= "numeric",	
	theta			= "numeric",
	n.loci			= "numeric",
	n.pop			= "numeric",
	n.sites			= "numeric",
	seeds			= "numeric",
	obs.val			= "matrix",		#matrix of size niter x numTests
	
	print.tree		= "numeric",	# set to 1 to include tree information in class
	fixed.seg.sites	= "numeric",	# set a fixed number of segregating sites for all iteration
	recombination	= "numeric",	# c(p, nsites) p = cross over parameter rate, nsites is the number of sites between recombination occurs
	
	gene.conv		= "numeric",
	
	# specify alpha in N_0exp^(alpha-t), where t is the time before present, measured in units of 4N_0 generations
	growth			= "numeric",
	
	migration		= "numeric",
	demography		= "numeric"
)) 



#------------------------------------------------------------------------------#
#									Methods									   #
#------------------------------------------------------------------------------#

setGeneric("get.n.sam", function(object) standardGeneric("get.n.sam"))
setMethod("get.n.sam", "test.params", function(object){
	return(object@nsam)
})

setGeneric("set.n.sam", function(object, val) standardGeneric("set.n.sam"))
setMethod("set.n.sam", "test.params", function(object, val){
	object@n.sam <- val
	return(object)
})


setGeneric("get.n.iter", function(object) standardGeneric("get.n.iter"))
setMethod("get.n.iter", "test.params", function(object){
	return(object@n.iter)
})

setGeneric("set.n.iter", function(object, val) standardGeneric("set.n.iter"))
setMethod("set.n.iter", "test.params", function(object, val){
	object@n.iter <- val
	return(object)
})


setGeneric("get.theta", function(object) standardGeneric("get.theta"))
setMethod("get.theta", "test.params", function(object){
	return(object@theta)
})

setGeneric("set.theta", function(object, val) standardGeneric("set.theta"))
setMethod("set.theta", "test.params", function(object, val){
	object@theta <- val
	return(object)
})


setGeneric("get.n.loci", function(object) standardGeneric("get.n.loci"))
setMethod("get.n.loci", "test.params", function(object){
	return(object@n.loci)
})

setGeneric("set.n.loci", function(object, val) standardGeneric("set.n.loci"))
setMethod("set.n.loci", "test.params", function(object, val){
	object@n.loci <- val
	return(object)
})


setGeneric("get.n.pop", function(object) standardGeneric("get.n.pop"))
setMethod("get.n.pop", "test.params", function(object){
	return(object@n.pop)
})

setGeneric("set.n.pop", function(object, val) standardGeneric("set.n.pop"))
setMethod("set.n.pop", "test.params", function(object, val){
	object@n.pop <- val
	return(object)
})


setGeneric("get.n.sites", function(object) standardGeneric("get.n.sites"))
setMethod("get.n.sites", "test.params", function(object){
	return(object@n.sites)
})

setGeneric("set.n.sites", function(object, val) standardGeneric("set.n.sites"))
setMethod("set.n.sites", "test.params", function(object, val){
	object@n.sites <- val
	return(object)
})


setGeneric("get.seeds", function(object) standardGeneric("get.seeds"))
setMethod("get.seeds", "test.params", function(object){
	return(object@seeds)
})

setGeneric("set.seeds", function(object, val) standardGeneric("set.seeds"))
setMethod("set.seeds", "test.params", function(object, val){
	object@seeds <- val
	return(object)
})


setGeneric("get.obs.val", function(object) standardGeneric("get.obs.val"))
setMethod("get.obs.val", "test.params", function(object){
	return(object@obs.val)
})

setGeneric("set.obs.val", function(object, val) standardGeneric("set.obs.val"))
setMethod("set.obs.val", "test.params", function(object, val){
	object@obs.val <- val
	return(object)
})


setGeneric("get.print.tree", function(object) standardGeneric("get.print.tree"))
setMethod("get.print.tree", "test.params", function(object){
	return(object@print.tree)
})

setGeneric("set.print.tree", function(object, val) standardGeneric("set.print.tree"))
setMethod("set.print.tree", "test.params", function(object, val){
	object@print.tree <- val
	return(object)
})


setGeneric("get.fixed.seg.sites", function(object) standardGeneric("get.fixed.seg.sites"))
setMethod("get.fixed.seg.sites", "test.params", function(object){
	return(object@fixed.seg.sites)
})

setGeneric("set.fixed.seg.sites", function(object, val) standardGeneric("set.fixed.seg.sites"))
setMethod("set.fixed.seg.sites", "test.params", function(object, val){
	object@fixed.seg.sites <- val
	return(object)
})


setGeneric("get.recombination", function(object) standardGeneric("get.recombination"))
setMethod("get.recombination", "test.params", function(object){
	return(object@recombination)
})

setGeneric("set.recombination", function(object, val) standardGeneric("set.recombination"))
setMethod("set.recombination", "test.params", function(object, val){
	object@recombination <- val
	return(object)
})


setGeneric("get.gene.conv", function(object) standardGeneric("get.gene.conv"))
setMethod("get.gene.conv", "test.params", function(object){
	return(object@gene.conv)
})

setGeneric("set.gene.conv", function(object, val) standardGeneric("set.gene.conv"))
setMethod("set.gene.conv", "test.params", function(object, val){
	object@gene.conv <- val
	return(object)
})


setGeneric("get.growth", function(object) standardGeneric("get.growth"))
setMethod("get.growth", "test.params", function(object){
	return(object@growth)
})

setGeneric("set.growth", function(object, val) standardGeneric("set.growth"))
setMethod("set.growth", "test.params", function(object, val){
	object@growth <- val
	return(object)
})


setGeneric("get.migration", function(object) standardGeneric("get.migration"))
setMethod("get.migration", "test.params", function(object){
	return(object@migration)
})

setGeneric("set.migration", function(object, val) standardGeneric("set.migration"))
setMethod("set.migration", "test.params", function(object, val){
	object@migration <- val
	return(object)
})


setGeneric("get.demography", function(object) standardGeneric("get.demography"))
setMethod("get.demography", "test.params", function(object){
	return(object@demography)
})

setGeneric("set.demography", function(object, val) standardGeneric("set.demography"))
setMethod("set.demography", "test.params", function(object, val){
	object@demography <- val
	return(object)
})


