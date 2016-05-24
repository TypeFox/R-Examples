#-------------------------------------------------------------------------------#
# Package: fastclime                                                            #
# fastclp(): A parametric simplex LP solver for parameterized LP problems       #
# Authors: Haotian Pang, Han Liu and Robert Vanderbei                           #
# Emails: <hpang@princeton.edu>, <hanliu@princeton.edu> and <rvdb@princetonedu> #
# Date: April 22th 2016                                                           #
# Version: 1.4.1						                                        #
#-------------------------------------------------------------------------------#

paralp <- function(obj, mat, rhs, obj_bar, rhs_bar, lambda=0){

 m<-length(rhs)
 n<-length(obj)
 m1<-length(rhs_bar)
 n1<-length(obj_bar)
 m0<-dim(mat)[1]
 n0<-dim(mat)[2]

 opt<-rep(0,n) 
 status<-0
 error<-0

 if (m!=m0 || n!=n0 || m!=m1 || n!=n1){
    cat("Dimensions do not match! \n")
    error<-1
 }

 if (any(obj_bar<0) || any(rhs_bar<0)){

    cat("The pertubation vector obj_bar and rhs_bar must be nonnegative! \n")
    error<-2

 }

 if (error==0){
    str=.C("paralp", as.double(obj), as.double(t(mat)), as.double(rhs), as.integer(m0), as.integer(n0), as.double(opt), as.integer(status), as.double(lambda), as.double(rhs_bar), as.double(obj_bar), PACKAGE="fastclime")


    opt<-unlist(str[6])
    status<-unlist(str[7])

    if (status==0){  
       cat("optimal solution found! \n")
       return(opt)
    }

    else if(status ==1){
       cat("The problem is infeasible! \n")
    }

    else if(status ==2){
       cat("The problem is unbounded! \n")
    }

 }

   
}	
	
          













