setinitialvalues<-function(data,mix.terms,model.formula,z=NULL,zdim=0,width,showit=0){
   # set the initial values
   initialvalues<-c()

   ###################################
   # find the pis
   pi.start<-c()

   if(mix.terms==1){
      # when we only have one mixture, set pi=1
      #pi.start<-c(1)
      pars<-as.numeric(lm(as.formula(paste("log(distance+width/1000)",
                       model.formula[[1]],sep="")),data=data)$coef)
      return(pars)
   }else{
      # since we're using the reparameterisation of the pis
      # need to do something smart here...
      pi.start<-inv.reparam.pi(rep(1/mix.terms,mix.terms))
   }

   ####################################
   ## find other pars...
   beta.start<-c()

   ## first, sort the distances
   #odist<-order(data$distance) # keep their order 

   ## split the distances into mix.terms equal parts, then use B+R for each part

   intmod<-FALSE
   if(!is.list(model.formula)){
      model.formula<-as.list(rep(model.formula,mix.terms))
      intmod<-TRUE
   }

   for(i in 1:mix.terms){

      # set of indicators we currently care about
      #this.odist<-odist[((i-1)*length(odist)/mix.terms+1):(i*length(odist)/mix.terms)]
      # pick the elements of the dataset we want
      #this.data<-data[this.odist,]
      this.data<-data

      # use B+R approach
      pars<-as.numeric(lm(as.formula(paste("log(distance+width/1000)",
                       model.formula[[i]],sep="")),data=this.data)$coef)

      # if we have an intercept model or just ~1
      if(intmod){
         pars<-pars[1]
      }
      
      beta.start<-c(beta.start,pars)
   }

   if(intmod){
      pars<-as.numeric(lm(as.formula(paste("log(distance+width/1000)",
                       model.formula[[1]],sep="")),data=data)$coef)
      pars<-pars[-1]
      beta.start<-c(beta.start,pars)
   }


   ###################################
   # put it all together
   initialvalues<-c(beta.start,pi.start)

   return(initialvalues)
}
