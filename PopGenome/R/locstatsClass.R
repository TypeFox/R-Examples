###############################################################################
#
# CLASS:	locstats
#
# this object holds all necessary information about simulated data for each 
# loci. 
#
# AUTHOR:			Niroshan Nadarajah <niroshan.nadarajah@uni-duesseldorf.de>
#
# LAST MODIFIED:	10/11/04
#
###############################################################################


setClass("loc.stats", representation(
          
	#Slots of coalsim

	n.sam 		= "numeric",
	n.iter 		= "numeric",
	theta 		= "numeric",
	obs.val		= "matrix",
	positions	= "matrix",
	trees		= "matrix",

	seeds		= "numeric",

	haplotypes 	= "matrix",
	stats 		= "matrix",
	
	loc.prob.less	= "matrix",
	loc.prob.equal  = "matrix",
	loc.valid.iter  = "matrix",
	quantiles	= "matrix"
	
)) 

#------------------------------------------------------------------------------#
#									Methods									   #
#------------------------------------------------------------------------------#

#### SHOW ######
setMethod("show", "loc.stats",
	function(object){
		print(summary(object))
		cat("-----\n")
		cat("SLOTS:\n")
		cat("-----\n")
		out <- data.frame (Slots= c("n.sam","n.iter","theta","obs.val","positions", "trees", "seeds", "halplotypes", "stats", "loc.prob.less", "loc.prob.equal", "loc.valid.iter", "quantiles"),
			Description=c("number of samples for each iteration",
						"number of iteration", 
						"mutation parameter", 
						"vector with observed values for each test",
						"position of each polymorphic site",
						"if printtree=1, gene tree in Newick format",
						"random numbers used to generate samples",
						"haplotypes in each iteration",
						"variety of test stats compiled a matrix",
						"Prob. that simulated val. <= to observed val. P(Sim <= Obs)",
						"Prob. that simulated val = to  observed val. P(Sim = Obs)",
						"number of valid iteration for each test",
						"13 quantiles for each test") 
		)
	 print(out)
	 cat("\n---------------\n")
 	print("These are the Slots")
  
	}
)

### get number of samples #######
setGeneric("get.n.sam", function(object) standardGeneric("get.n.sam"))
setMethod("get.n.sam", "loc.stats", function(object){
	return(object@n.sam)
})

### get number of iteration #######
setGeneric("get.n.iter", function(object) standardGeneric("get.n.iter"))
setMethod("get.n.iter", "loc.stats", function(object){
	return(object@n.iter)
})

### get theta value for this locus #######
setGeneric("get.theta", function(object) standardGeneric("get.theta"))
setMethod("get.theta", "loc.stats", function(object){
	return(object@theta)
})

### get observed values for this loci #######
setGeneric("get.obs.val", function(object) standardGeneric("get.obs.val"))
setMethod("get.obs.val", "loc.stats", function(object){
	return(object@obs.val)
})

### get a matrix with position information for each polymorphic site #######
setGeneric("get.positions", function(object) standardGeneric("get.positions"))
setMethod("get.positions", "loc.stats", function(object){
	
	#positionsLine <- unlist(strsoutput[i]tions, " #"))
	#positions <- rbind(positions, as.numeric(positionsLine[2:length(positionsLine)]))
	
	return(object@positions)
})

### get gene tree information in newick format #######
setGeneric("get.trees", function(object) standardGeneric("get.trees"))
setMethod("get.trees", "loc.stats", function(object){
	return(object@trees)
})

### get the initial random values to create samples #######
setGeneric("get.seeds", function(object) standardGeneric("get.seeds"))
setMethod("get.seeds", "loc.stats", function(object){
	return(object@seeds)
})

### get matrix with generated haplotypes #######
setGeneric("get.haplotypes", function(object) standardGeneric("get.haplotypes"))
setMethod("get.haplotypes", "loc.stats", function(object){
	return(object@haplotypes)
})

### get a matrix with a variety of summary statistics on generated samples #######
setGeneric("get.stats", function(object) standardGeneric("get.stats"))
setMethod("get.stats", "loc.stats", function(object){
	return(object@stats)
})

### get probablilties #######
setGeneric("get.loc.prob.less", function(object) standardGeneric("get.loc.prob.less"))
setMethod("get.loc.prob.less", "loc.stats", function(object){
	return(object@loc.prob.less)
})

### get probablilties #######
setGeneric("get.loc.prob.equal", function(object) standardGeneric("get.loc.prob.equal"))
setMethod("get.loc.prob.equal", "loc.stats", function(object){
	return(object@loc.prob.equal)
})

### get number of considered iteration #######
setGeneric("get.loc.valid.iter", function(object) standardGeneric("get.loc.valid.iter"))
setMethod("get.loc.valid.iter", "loc.stats", function(object){
	return(object@loc.valid.iter)
})


### get quantiles  #######
setGeneric("get.quantiles", function(object) standardGeneric("get.quantiles"))
setMethod("get.quantiles", "loc.stats", function(object){
	return(object@quantiles)
})
