###############################################################################
#
# CLASS:	csstats
#
# this is the class which is returned from the function coalsim(). Every 
# computed statistic for each sample and in each iteration for each loci is 
# included in this object.
#
# First, some general parameters (i.e. which is not loci specific for the 
# simulation are given, such as a vector with the number of samples,
# the number of iteration and number of population are given.   
# 
#
#
# AUTHOR:			Niroshan Nadarajah <niroshan.nadarajah@uni-duesseldorf.de>
#
# LAST MODIFIED:	10/11/03
#
###############################################################################


setClass("cs.stats", representation(
          
	#Slots of coalsim
	n.loci		= "numeric",
	n.iter		= "numeric",
	n.pop		= "numeric",
				
	prob.less 	= "matrix",
	prob.equal	= "matrix",
	valid.iter	= "matrix",
	obs.val		= "matrix",
	
	average		= "matrix",
	variance	= "matrix",
	
	locus		= "list"
	
)) 


#------------------------------------------------------------------------------#
#									Methods									   #
#------------------------------------------------------------------------------#

#### SHOW ######
setMethod("show", "cs.stats",
	function(object){
		#print(summary(object))
		cat("-----\n")
		cat("SLOTS:\n")
		cat("-----\n")
		out <- data.frame (Slots= c("prob.less","prob.equal","valid.iter","obs.val","n.loci", "n.iter", "average", "variance", "locus"),
			Description=c("Prob. that sim.val <=  obs.val P(sim <= obs)",
						"Prob. that sim.val = obs.val P(sim = obs)", 
						"number of valid iter. for each test and loci", 
						"obs.values for each test",
						"number of loci considered", 
						"number of iterations for each loci", 
						"average values of each statistic  (across all loci)", 
						"variance values of each statistic (across all loci)", 
						"list of loc.stats objects, (detail stats for each locus)")
		)
		print(out)
		cat("\n---------------\n")
		#cat("\nThese are the Slots\n")
		#cat("Slot 'nloci':\n")
		#print(object@nloci)
		
		#cat("\nSlot 'niter':\n")
		#print(object@niter)
		
		#cat("\nSlot 'probLess':\n")
		#print(object@probLess)
		
		#cat("\nSlot 'probEqual':\n")
		#print(object@probEqual)
		
		#cat("\nSlot 'obsVal':\n")
		#print(object@obsVal)
		
		#cat("\nSlot 'average':\n")
		#print(object@average)
		#cat("\nSlot 'variance':\n")
		#print(object@variance)
		
		#cat("\nType l <- getLocus(obj) to get a list of locus objects")
	}
)

### get matrix with probabilties #######
setGeneric("get.prob.less", function(object) standardGeneric("get.prob.less"))
setMethod("get.prob.less", "cs.stats", function(object){
	return(object@prob.less)
})

### get matrix with probabilties #######
setGeneric("get.prob.equal", function(object) standardGeneric("get.prob.equal"))
setMethod("get.prob.equal", "cs.stats", function(object){
	return(object@prob.equal)
})

### get matrix with number of iteration for each test and loci #######
setGeneric("get.valid.iter", function(object) standardGeneric("get.valid.iter"))
setMethod("get.valid.iter", "cs.stats", function(object){
	return(object@valid.iter)
})

### get matrix with number of iteration for each test and loci #######
setGeneric("get.obs.val", function(object) standardGeneric("get.obs.val"))
setMethod("get.obs.val", "cs.stats", function(object){
	return(object@obs.val)
})

### get the number of loci for this simulation #######
setGeneric("get.n.loci", function(object) standardGeneric("get.n.loci"))
setMethod("get.n.loci", "cs.stats", function(object){
	return(object@n.loci)
})

### get the number of iterations of each test each simulation #######
setGeneric("get.n.iter", function(object) standardGeneric("get.n.iter"))
setMethod("get.n.iter", "cs.stats", function(object){
	return(object@n.iter)
})

### get the number of iterations of each test each simulation #######
setGeneric("get.average", function(object) standardGeneric("get.average"))
setMethod("get.average", "cs.stats", function(object){
	return(object@average)
})

### get the number of iterations of each test each simulation #######
setGeneric("get.variance", function(object) standardGeneric("get.variance"))
setMethod("get.variance", "cs.stats", function(object){
	return(object@variance)
})

### get the number of iterations of each test each simulation #######
setGeneric("get.locus", function(object) standardGeneric("get.locus"))
setMethod("get.locus", "cs.stats", function(object){
	return(object@locus)
})

