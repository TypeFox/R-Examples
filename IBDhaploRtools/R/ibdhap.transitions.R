ibdhap.transitions <-function( calls, data.type = c("h", "g","r")) { 


#define parameters based on haplotype/genotype/or reduced data
if(length(data.type)>1){ stop("data.type improperly indicated, please see ?ibdhap.transitions")}
else if( is.element("h", data.type)){no.ibd.ind = 15}
else if( is.element("g", data.type)){no.ibd.ind = 9}
else if( is.element("r", data.type)){no.ibd.ind = 4}
else{ stop("data.type improperly indicated, please see ?ibdhap.transitions")}


 
state.num <- (no.ibd.ind )
#define the matrix of transition counts
transition.counts<-matrix(0, nrow = state.num, ncol = state.num)


#get the counts for each row, and add them up
   for(i in 1:ncol(calls)){

      temp.states<-as.numeric(unlist(calls[,i]))

      #remove the zeros
      temp.states<-removezeros(temp.states)
      temp.counts<-get.counts(temp.states, state.num)

      transition.counts<-transition.counts+temp.counts

   }

#get rid of the states where no transition happens)
for(i in 1:state.num){ transition.counts[i,i]<-0}

return(transition.counts )

}

