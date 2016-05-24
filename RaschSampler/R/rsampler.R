"rsampler" <-
function(inpmat,controls=rsctrl()){

   if (!(class(controls)=="RSctr"))
         stop("controls is not a control object - see help(\"rsctrl\")")

   n       <- dim(inpmat)[1]
   k       <- dim(inpmat)[2]
   burn_in <- controls$burn_in
   n_eff   <- controls$n_eff
   step  <- controls$step
   seed    <- controls$seed
   tfixed  <- controls$tfixed

   if (seed == 0) {
      # generates random seed in the range [536870911,772830910]
      seed <- as.integer(as.double(format(Sys.time(), "%H%M%OS3"))*1000)
                   + 2**29 - 1
   }

   # allocation of memory for simulated matrices
   vec<-vector( length = (n_eff+1)*n*trunc((k+31)/32) )
   ier<-0

   # calls the external Fortran subroutine sampler
   # simulated matrices are returned in vec
   RET<-.Fortran("sampler",
               n=as.integer(n),
               k=as.integer(k),
               inpmat=as.integer(inpmat),
               tfixed=as.logical(tfixed),
               burn_in=as.integer(burn_in),
               n_eff=as.integer(n_eff),
               step=as.integer(step),
               seed=as.integer(seed),
               outvec=as.integer(vec),
               ier=as.integer(ier)
   )
   n_tot <- n_eff+1
   if (RET$ier>0) {
         rserror(RET$ier)
   } else {
         RET<-c(RET[1:8],n_tot=n_eff+1,RET[9:10])
         class(RET)<-"RSmpl"
         RET
   }
}

