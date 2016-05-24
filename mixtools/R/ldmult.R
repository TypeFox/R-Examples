###########################################################################
# log density function of multinomial count vector y with parameter theta #
###########################################################################

ldmult<-function(y,theta){
  if (any(is.nan(theta)) || any(theta==0)) {
#      stop ("The theta parameter cannot have a zero component.")
     out = -Inf
   }
   else {
     if (length(y)==length(theta)+1)
       theta=c(theta,1-sum(theta))
     out=lgamma(1+sum(y)) - sum(lgamma(1+y)) + sum(y*log(theta))
   }
   out
}


