# this is the file that has all the EM algorithm stuff in it...
# actually do the EM
em<-function(data,init.pars,mix.terms,zdim=NULL,z=NULL,ftype="hn",width,
             showit=0,grad=FALSE,maxit=1000,pt=FALSE){

   eps<-1e-8 # as in optim

   opars<-init.pars*100
   pars<-init.pars

   # iteration counter
   iter<-1

   # while not converged...
   while((all(abs(pars-opars)>eps)) & (iter<maxit)){

      # store the old pars
      opars<-pars

      # do the E step - calculate the w_{ij}s
      W<-em.e(data,pars,mix.terms,ftype,width,zdim,z,pt)

      # calculate the pi_js
      mix.prop<-(1/length(data$distance)) * colSums(W)

      # debug, show the pis
      if(showit>1){
         cat("E-step pars:",mix.prop,"\n")
      }

      # pick out the correct parameters, back-transform the pis
      pars<-pars[-((length(pars)-mix.terms+2):length(pars))]
      pars<-c(pars,inv.reparam.pi(mix.prop))

      ## switch the pars around
      #swpars<-switchpars(pars,mix.terms,zdim,z)
      #pars<-swpars$fpar
      #z<-swpars$z
      #zdim<-swpars$zdim

      ## switch the pars we compare against too!
      #oswpars<-switchpars(opars,mix.terms,zdim,z)
      #opars<-oswpars$fpar

      # do the M step
      mstep<-em.m(pars,data,ftype=ftype,mix.terms,width,z=z,zdim=zdim,
                  showit=showit,grad=grad,pt=pt)

      # if something went wrong... 
      if(class(mstep)=="try-error") {   
         cat("Something went wrong in em.m.\n")
         cat("The pars that caused this were",pars,"\n")
         err<-list(pars=pars,value=NA)
         class(err)<-"try-error"
         return(err)
      }

      # back-transform the pis and put them back in
      pars<-c(mstep$par,inv.reparam.pi(mix.prop))

      # debug
      if(showit>1){
         cat("M-step pars=",mstep$par,"\n")
      }
      # increment iteration counter
      iter<-iter+1
   }

   # put things in the right order
   fixed.pars<-switchpars(pars,mix.terms,zdim,z)
   pars<-fixed.pars$fpar
   z<-fixed.pars$z
   zdim<-fixed.pars$zdim

   # turn on the gradient function wrapper and
   # gradient for final eval
   if(grad){
      fgradfunc<-flt.gr
   }else{
      fgradfunc<-NULL
   }

   # run optim for 1 iteration to calulate the hessian etc needed later
   opt<-try(optim(pars,flt,gr=fgradfunc,method="BFGS", 
            control=c(lmm=1000,fnscale=-1,maxit=1),
            mix.terms=mix.terms,
            width=width, showit=showit, x=data,ftype=ftype,
            hessian=TRUE,z=z,zdim=zdim,pt=pt,EM=TRUE))
   
   # debug
   if(showit>0){
      cat("EM done in",iter,"iterations.\n") 
   }
   if(showit>0){
      cat("Final pars=",pars,"\n")
   }

   # return a list of things we want...
   return(opt)
}


# expectation routine
em.e<-function(dat,pars,mix.terms,ftype="hn",width,zdim=NULL,z=NULL,pt=FALSE){
   # args
   #  dat         data.frame with distances and other covars
   #  pars        model parameters
   #  mix.terms   number of mixture parameters
   #  zdim        covar stuff
   #  z             "     "

   # w_{ij}s are the weights for each par/mixture
   W<-eval.pdf(pars,dat,width,mix.terms,showit=0,
               ftype=ftype,z=z,zdim=zdim,pt=pt,EM=TRUE)

   # calculate the denominator, sum over the mixtures for each 
   # observation
   mix.sums<-rowSums(W)

   # calculate the expectations
   for(j in 1:mix.terms){
      W[,j]<-W[,j]/mix.sums
   }   

   # return them
   return(W)
}

# do the maximisation part
em.m<-function(pars,data,ftype="hn",mix.terms,width,z=NULL,zdim=NULL,showit=0,
               control=c(lmm=1000,fnscale=-1),grad=TRUE,pt=FALSE){
   # turn on the gradient function wrapper and
   # gradient for final eval
   if(grad){
      gradfunc<-flt.gr.wrap
   }else{
      gradfunc<-NULL
   }

   mix.prop<-pars[(length(pars)-mix.terms+2):length(pars)]
   pars<-pars[1:(length(pars)-mix.terms+1)]

   if(showit>2){
      cat("mix.prop=",mix.prop,"\n")
      cat("pars=",pars,"\n")
   }

   # maximise
   opt<-try(optim(pars, flt.wrap,gr=gradfunc, method="BFGS",control=control,
            mix.terms=mix.terms, mix.prop=mix.prop,
            width=width, showit=showit, x=data,ftype=ftype,
            hessian=TRUE,z=z,zdim=zdim,pt=pt))

   # did something bad happen?
   if(showit>2){
      if(class(opt)!="try-error") {   
         cat("optim() died\n")
      }
   }
   return(opt)
}

# wrapper for flt so that we can still use flt
flt.wrap<-function(fpar,mix.prop,x,width, mix.terms,showit=0, 
                   ftype="hn",z=NULL,zdim=0,pt=FALSE){

   # put the pars in the right order
   pars<-c(fpar,mix.prop)

   # call flt, return the likelihood
   return(flt(pars,x,width,mix.terms,showit,ftype,z,zdim,pt,EM=TRUE))
}

# wrapper for flt.gr for analytic gradients
flt.gr.wrap<-function(fpar,mix.prop,x,width, mix.terms,showit=0, 
                   ftype="hn",z=NULL,zdim=0,pt=FALSE){

   # put the pars in the right order
   pars<-c(fpar,mix.prop)

   # call flt.gr
   grad<-flt.gr(pars,mix.terms,width,x,ftype,z,zdim,showit-1,pt,EM=TRUE)

   # ignore the derivatives wrt the mixture proportions
   grad<-grad[1:(length(grad)-length(mix.prop))]

   if(showit>2)
      cat("Derivative pars=",grad,"\n")

   return(grad)
}
