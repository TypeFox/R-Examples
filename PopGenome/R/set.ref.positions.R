setGeneric("set.ref.positions", function(object,positions) standardGeneric("set.ref.positions"))
 setMethod("set.ref.positions", "GENOME",
 function(object,positions){

## biallelic.sites

change <- object@region.data

for( xx in 1:length(change@biallelic.sites) ){

         
	change@biallelic.sites[[xx]]            <- positions[[xx]][change@biallelic.sites[[xx]]]
        ENDE                                    <- length(change@biallelic.sites[[xx]])
        colnames(change@biallelic.matrix[[xx]]) <- change@biallelic.sites[[xx]]
        object@n.sites[xx]                      <- change@biallelic.sites[[xx]][ENDE] 	

}

object@region.data <- change

return(object)

})
