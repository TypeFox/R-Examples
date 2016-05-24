`pyear` <-
function(t,delta,knots){
# function for computing life tables from survival data
#
# INPUT ARGUMENT (all vectors should be column vector)
#      t = time to event
#  delta = 0-1 indicator of events
#  knots = vector of knots for defining age intervals
#          (the first knot should not be smaller than the time of the first event
#           and the last knot should not be higher than the time of the last event)
#
# OUTPUT ARGUMENT
#  d = N by (K+1) table of events (N=number of subjects, K=number of knots)
#  py = N by (K+1) table of person-years (N=number of subjects, K=number of knots)

# Add zeros and infinity into the knots

knots<-c(0,knots,Inf)
k<-length(knots)
n<-length(t)

py <- matrix(0,n,k-1)
d <- matrix(0,n,k-1)

tc <- cut(t,knots,labels=FALSE,include.lowest=TRUE,right=FALSE) # bin2

for(i in 2:k){
   py[                     ,(i-1)] <- rep( knots[i]-knots[i-1], n)
   py[tc<(i-1)             ,(i-1)] <- 0
   py[tc==(i-1)            ,(i-1)] <- t[tc==(i-1)]-knots[i-1]+0.5
   d[(delta==1)&(tc==(i-1)),(i-1)] <- 1
}
list(d=d,py=py)
}

