

###########################
#### GENERAL ACCESSORS ####
###########################

## NOTE:
## only generic functions are defined here, and possibly default methods
## specific methods are defined in relevant files

############
## get.id ##
############
## return the id of the object
## setGeneric("get.id", function(x, ...) standardGeneric("get.id"))


#############
## get.dna ##
#############
## return DNA sequence alignments
setGeneric("get.dna", function(x, ...) standardGeneric("get.dna"))


################
## get.ntrees ##
################
## return multiPhylo object (list of trees)
setGeneric("get.ntrees", function(x, ...) standardGeneric("get.ntrees"))


##############
## get.trees ##
##############
## return multiPhylo object (list of trees)
setGeneric("get.trees", function(x, ...) standardGeneric("get.trees"))


###############
## get.locus ##
###############
## return the loci in the object
setGeneric("get.locus", function(x, ...) standardGeneric("get.locus"))


################
## get.nlocus ##
################
## return the number of loci in the object
setGeneric("get.nlocus", function(x, ...) standardGeneric("get.nlocus"))


####################
## get.sequences ##
####################
## return the id of DNA sequences in the object
setGeneric("get.sequences", function(x, ...) standardGeneric("get.sequences"))


####################
## get.nsequences ##
####################
## return the number of sequences in the object
setGeneric("get.nsequences", function(x, ...) standardGeneric("get.nsequences"))



#####################
## get.individuals ##
#####################
## return the individuals in the object
setGeneric("get.individuals", function(x, ...) standardGeneric("get.individuals"))



######################
## get.nindividuals ##
######################
## return the number of individuals in the object
setGeneric("get.nindividuals", function(x, ...) standardGeneric("get.nindividuals"))



#################
## get.records ##
#################
## return the names of the records tables in the object
setGeneric("get.records", function(x, ...) standardGeneric("get.records"))


##################
## get.nrecords ##
##################
## return the number of records tables in the object
setGeneric("get.nrecords", function(x, ...) standardGeneric("get.nrecords"))


#################
## get.context ##
#################
## return the names of the context tables in the object
setGeneric("get.context", function(x, ...) standardGeneric("get.context"))


##################
## get.ncontext ##
##################
## return the number of context tables in the object
setGeneric("get.ncontext", function(x, ...) standardGeneric("get.ncontext"))


###############
## get.dates ##
###############
## return the dates in the object
setGeneric("get.dates", function(x, ...) standardGeneric("get.dates"))


#################
## get.ndates ##
#################
## return the dates in the object
setGeneric("get.ndates", function(x, ...) standardGeneric("get.ndates"))


##################
## get.contacts ##
##################
## return the number of contacts in the object
setGeneric("get.contacts", function(x, ...) standardGeneric("get.contacts"))



###################
## get.ncontacts ##
###################
## return the number of contacts in the object
setGeneric("get.ncontacts", function(x, ...) standardGeneric("get.ncontacts"))



##############
## get.data ##
##############
## return the number of contacts in the object
setGeneric("get.data", function(x, ...) standardGeneric("get.data"))



############
## subset ##
############
## return the number of contacts in the object
setGeneric("subset", function(x, ...) standardGeneric("subset"))



###############
## make.phylo ##
###############
## return DNA sequence alignments
setGeneric("make.phylo", function(x, ...) standardGeneric("make.phylo"))



## ###################
## ## get.incidence ##
## ###################
## setGeneric("get.incidence", function(x, ...) standardGeneric("get.incidence"))
