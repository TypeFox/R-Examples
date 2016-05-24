setGeneric("get.individuals", function(object,region=FALSE) standardGeneric("get.individuals"))

setMethod("get.individuals","GENOME",function(object,region){

n.regions <- length(object@region.names)


if(region[1]==FALSE){
	region <- 1:n.regions
}

ind.names <- vector("list", length(region))

for(xx in 1:length(region)){

ind.names[[xx]] <- rownames(object@region.data@biallelic.matrix[[region[xx]]])

}

return(ind.names)

})
