setGeneric("set.outgroup", function(object,new.outgroup=FALSE, diploid=FALSE) standardGeneric("set.outgroup"))
 setMethod("set.outgroup", "GENOME",
 function(object,new.outgroup,diploid){

change             <- object@region.data

# if diploid add individuals.2
if(diploid){
  dottwo                <- paste(new.outgroup,".2", sep="")
  new.outgroup          <- c(new.outgroup,dottwo) 
}
# End of diploid

object@outgroup <- new.outgroup

# Init

# XX_popmissing  <- vector("list",length(object@region.names))
XX_outgroup <- vector("list",length(object@region.names))

############################################################
progr <- progressBar()
############################################################


for(xx in 1:length(object@region.names)){


     bial <- object@region.data@biallelic.matrix[[xx]] # popGetBial(object,xx) # if it does not fit into RAM
     
     if(length(bial)==0){next}   
  
           if(is.character(new.outgroup)){        
              outgroup          <- match(new.outgroup,rownames(bial))
              naids             <- which(!is.na(outgroup))
              outgroup          <- outgroup[naids]
           }else{
              outgroup          <- new.outgroup
              ids               <- which(outgroup>dim(bial)[1])
              if(length(ids)>0){ outgroup <- outgroup[-ids]}   
           }      
             
     XX_outgroup[[xx]]      <- outgroup
     #XX_popmissing[[xx]]   <- popmissing

 # PROGRESS #######################################################
    progr <- progressBar(xx,length(object@region.names), progr)
 ###################################################################

}

change@outgroup <- XX_outgroup
#change@popmissing  <- XX_popmissing
object@region.data <- change
return(object)





})
