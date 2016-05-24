#-------------------------------------------------------------------------------#
# Package: fastclime                                                            #
# fastclp(): A parametric simplex LP solver                                     #
# Authors: Haotian Pang, Han Liu and Robert Vanderbei                           #
# Emails: <hpang@princeton.edu>, <hanliu@princeton.edu> and <rvdb@princetonedu> #
# Date: April 22th 2016                                                           #
# Version: 1.4.1					                                            #
#-------------------------------------------------------------------------------#

fastlp <- function(obj, mat, rhs,lambda=0){

 m<-length(rhs)
 n<-length(obj)
 m0<-dim(mat)[1]
 n0<-dim(mat)[2]

 opt<-rep(0,n) 
 status<-0
 error<-0

 if (m!=m0 || n!=n0){
    cat("Dimensions do not match! \n")
    error<-1
 }

 if (error==0){
    str=.C("fastlp", as.double(obj), as.double(t(mat)), as.double(rhs), as.integer(m0), as.integer(n0), as.double(opt), as.integer(status), as.double(lambda), PACKAGE="fastclime")


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
	
          













