# function to check the initialvalues are valid...
checkinitialvalues<-function(initialvalues,zdim,mix.terms,model.formula=NULL){

   # check that there are the correct number of initialvalues
   correctlength<-(mix.terms-1) # number of pis

   if(sum(zdim)==0){
      # for the noncovar case
      correctlength<-correctlength+mix.terms # number of sigmas
   }else if(length(zdim)==1){
      # covar with intercept vary
      correctlength<-correctlength+(zdim-1)+mix.terms
   }else{
      # with covars
      correctlength<-correctlength+sum(zdim) # number of betas
   }

   if(length(initialvalues)!=correctlength)
      stop("initialvalues is not the correct length\n")

   # if nothing went wrong and we don't have covars then just return the
   # vector we got given
   return(initialvalues)
}
