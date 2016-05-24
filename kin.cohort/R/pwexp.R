`pwexp` <-
function(d,py,w,knots){#  t,delta,w,knots){
# function for piecewise exponential survival estimation
#
# INPUT ARGUMENT (all vectors should be column vector)
#     d = table of events (calculated by pyear)
#    py = table of person years (calculated by pyear)
#     w = vector of prior weights (assume 1 if no weights are needed)
# knots = vector of knots for piecewise constant model of hazards
#         (the first knot should not be smaller than the time of the first event
#          and the last knot should not be higher than the time of the last event)
#
# OUTPUT ARGUMENT
#      s=  K by 1 vector of survival probabilities (K=number of knots)    ## was K+1
#      h = (K+1) by 1 vector of hazards (K=number of knots)


   knots<-c(0,knots,Inf)
   k<-length(knots)

#  weights supplied
   if(length(unique(w))!=1){
     w  <- matrix(rep(w,k-1),length(w),k-1)
     d  <- d*w
     py <- py*w                                                              
   }

   sdw <- apply(d,2,sum)
   spyw <- apply(py,2,sum)

   h<-sdw/spyw
   h[sdw==0]<-0

# this excludes the last (open) interval in s
   cH<-cumsum(h[1:(k-2)]*(knots[2:(k-1)]-knots[1:(k-2)])) 
   s<-exp(-cH)
   
   list(s=s,h=h)
}

