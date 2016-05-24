validannimp<-function(kmin, kmax, nsol, exclude, nexclude, include, ninclude, initialsol, implog){


##################################################################
# validation of input that is specific to the anneal and improve #
# functions (which admit similar possible initial solutions)     #
##################################################################



##################################################################
# initializations when initial solutions are specified; checking #
# nature of initialsol (and how to interpret it) and for         #
# conflicts with the exclude and include requirements (that must #
# be respected by the initial solutions)                         #
##################################################################

        if  (!is.null(initialsol) & (sum(!(as.integer(initialsol) == initialsol)) > 0)) stop("\n The initial solutions must be specified as integers indicating variable numbers")
        if (is.null(initialsol)) {silog <- FALSE}
        else {silog <- TRUE  # initial solutions have been specified by user          
#########################################################
# checking for the presence of variables that are to be #
# forcefully excluded                                   #
#########################################################

            if ((nexclude != 0) & (sum(exclude == rep(initialsol,rep(length(exclude),length(as.vector(initialsol))))) !=0)) stop("\n the specified initial solutions contain variables that are to be excluded")

#########################################################
# how to deal with various formats of input (of initial #
# solutions) for a single cardinality                   #
#########################################################

          dimsol<-dim(initialsol)
              
          if (length(dimsol) > 3) {
            stop("\n Can't handle arrays of more than 3 dimensions")}
          if (kmin == kmax) {

#############################################################
# inital solution is a vector: must be repeated if nsol > 1 #
#############################################################


               if (is.vector(initialsol)) {
                   if (length(initialsol) != kmax) 
                     {stop("\n The specified initial and final solutions have different cardinalities")}
                   else
                     {if (nsol > 1) {
                     initialsol<-matrix(nrow=nsol,ncol=kmax,rep(initialsol,rep(nsol,kmax)))
                     if (implog==TRUE) {warning("\n a single initial subset necessarily produces nsol identical \n final solutions in the improve algorithm")}}
                    }
               }


##############################
# initial solution is array? #
##############################

           if (is.array(initialsol)) {

###########################
# initialsol is 3-d array #
###########################

           if (length(dimsol) == 3) {
            if (dim(initialsol)[[3]] > 1) stop("\n The input array of initial solutions must have dimensions as \n (1 or nsol) x k x 1 when a single cardinality is requested")
            else initialsol<-matrix(nrow=dim(initialsol)[[1]],ncol=dim(initialsol)[[2]],initialsol)
            }
           else

##################################################
# initialsol is matrix: must be of form nsol x k #
##################################################

                 {if (dim(initialsol)[[2]] != kmax) stop("\n Input matrix of initial solutions must have as many columns as variables in the requested subset")
                  else 
                      {if  (dim(initialsol)[[1]] != nsol) {
                         if (dim(initialsol)[[1]] == 1) {
                          initialsol<-matrix(nrow=nsol,ncol=dim(initialsol)[[2]],rep(initialsol,rep(nsol,kmax))) 
                          if (implog==TRUE) {warning("\n a single initial subset necessarily produces nsol identical \n final solutions in the improve algorithm")}
                         }
                           else {stop("\n The number of initial solutions can only be 1 or nsol (number of final solutions requested)")}}}}
              }}



##########################################################
# how to deal with various formats of input (initialsol) #
# if more than one cardinality is requested              #
##########################################################

             else  # (if kmax > kmin), i.e., more than one cardinality requested

###############################################
# initial solution cannot be a single vector  #
###############################################

              {if (is.vector(initialsol)) stop("\n There must be initial solutions for all cardinalities requested")

               if (is.array(initialsol)) {

###############################
# initial solution 3-d array? #
###############################


                  if (length(dim(initialsol)) == 3) {
                   if ((dim(initialsol)[[3]] != length(kmin:kmax)) | (dim(initialsol)[[2]] != kmax) | ((dim(initialsol)[[1]] != nsol) & (dim(initialsol)[[1]] > 1))) stop("\n The input array of initial solutions must have dimensions as \n nsol x kmax x no. of different cardinalities requested")
                   else 
                    { if (dim(initialsol)[[1]] != nsol) {
                       initialsol<-array(dim=c(nsol,kmax,length(kmin:kmax)),rep(initialsol,rep(nsol,kmax*length(kmin:kmax))))
                      if (implog==TRUE) {warning("\n a single initial subset necessarily produces nsol identical \n final solutions in the improve algorithm")}
                      } 
                  }}
                else 


############################
# if initalsol is a matrix #
############################


                    {
                    if ((dim(initialsol)[[2]] != kmax) | (dim(initialsol)[[1]] != length(kmin:kmax))) stop("\n A matrix of initial solutions for more than one cardinality must be of type no. of cardinalities x kmax")
                    if (nsol > 1) {
                       initialsol<-array(dim=c(nsol,kmax,length(kmin:kmax)),rep(t(initialsol),rep(nsol,kmax*length(kmin:kmax)))) 
                    if (implog==TRUE) {warning("\n a single initial subset necessarily produces nsol identical \n final solutions in the improve algorithm")}
                    }
}}}
                 if ((ninclude != 0) & (sum(include == rep(initialsol,rep(length(include),length(as.vector(initialsol))))) != nsol*length(kmin:kmax)*length(include))) stop("\n Not all the specified initial solutions contain the variables that are to be included") 
}


##################################
#  assigning any changed values  #
##################################

         assign("silog",silog,pos=parent.frame())          
         assign("initialsol",initialsol,pos=parent.frame())

}
