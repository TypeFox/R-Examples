validgenetic<-function(kmin, kmax, popsize, mutprob, exclude,
nexclude, include, ninclude, initialpop){



################################################################
# validation of input that is specific to the genetic function #
################################################################
# WARNING: Some initializations look similar to those of the   #
# anneal and improvement functions, but here input MUST be a   #
# 3-dimensional array.                                         #
################################################################
        
######################################
# validation of mutation probability #
######################################

        if ((mutprob < 0) | (mutprob > 1)) stop("\n The mutation probability parameter (mutprob) must be between 0 and 1 \n")


################################################################
# initializations when initial population has been specified;  #
# checking the nature of initialpop (and how to interpret it)  #
# and for conflicts with the exclude and include requirements  #
# (which initial the population must respect)                  #
################################################################

        if (is.null(initialpop)) {pilog <- FALSE}
        else {pilog <- TRUE  # initial population has been specified by user          
################################################
# initial solution is a vector: not acceptable #
################################################

            if (is.vector(initialpop)) {
                     {stop("\n The specified initial population must have different k-subsets \n")}
                    }
################################################
# checking for the presence of variables that  #
# are to be forcefully excluded                #
################################################

            if ((nexclude != 0) & (sum(exclude == rep(initialpop,rep(length(exclude),length(as.vector(initialpop))))) !=0)) stop("\n the specified initial population contains variables that are to be excluded \n")

#########################################################
# how to deal with various formats of input (of initial #
#  population) for a single cardinality                 #
#########################################################

            dimpop<-dim(initialpop)
            if (length(dimpop) > 3) stop("\n Can't handle arrays of more than 3 dimensions \n")
            if (kmin == kmax) {

##############################
# initial solution is array? #
##############################

              if (is.array(initialpop)) {

###########################
# initialpop is 3-d array #
###########################

              if (length(dimpop) == 3) {
                 if (dimpop[[3]] > 1) stop("\n The input array of initial population must have dimensions as \n popsize x k x 1 when a single cardinality is requested \n")
                 else initialpop<-matrix(nrow=dimpop[[1]],ncol=dimpop[[2]],initialpop)
              }
              else # initialpop is matrix: must be of form popsize x k
 
                 {if (dimpop[[2]] != kmax) stop("\n Input matrix of initial solutions must have as many columns as variables in the requested subset \n")
                  else 
                      {if  (dimpop[[1]] != popsize) {
                        stop("\n The number of initial solutions can only be 1 or popsize (number of final solutions requested) \n")}}}
              }}

##########################################################
# how to deal with various formats of input (initialpop) #
# if more than one cardinality is requested              #
##########################################################

             else  #(if kmax > kmin), i.e., more than one cardinality requested

#########################################
# initial population must be 3-d array  #
#########################################


              {if (length(dimpop) < 3) stop("\n There must be initial populations for all cardinalities requested \n")
               else 
               {if ((dimpop[[3]] != length(kmin:kmax)) | (dimpop[[2]] != kmax) | (dimpop[[1]] != popsize)) stop("\n The input array of initial solutions must have dimensions as \n popsize x kmax x no. of different cardinalities requested \n")
                  }
                 if ((ninclude != 0) & (sum(include == rep(initialpop,rep(length(include),length(as.vector(initialpop))))) != popsize*length(kmin:kmax)*length(include))) stop("\n Not all the specified initial solutions contain the variables that are to be included \n") 
}}

##################################
#  assigning any changed values  #
##################################

         assign("pilog",pilog,pos=parent.frame())          
         assign("initialpop",initialpop,pos=parent.frame())

}


