######################################################################
# Shiryaev-Roberts based spatio-temporal cluster detection based
# on the work in Assuncao & Correa (2009). The implementation
# is based on C++ code was originally written by Marcos Oliveira Prates, UFMG,
# Brazil and provided by Thais Correa, UFMG, Brazil during her research
# stay in Munich. This stay was financially supported by the Munich
# Center of Health Sciences.
#
#
# Parameters:
#   x - vector containing spatial x coordinate of the events
#   y - vector containing spatial y coordinate of the events
#   t - vector containing the time points of the events
#   radius - is the radius of the cluster 
#   epsilon - is the relative change of event-intensity within the cluster
#       to detect
#   areaA - area of the observation region A (single number)
#   areaAcapBk - area of A \ B(s_k,\rho) for all k=1,\ldots,n (vector)
#   vector of areas A\B(s_k,\rho) for k=1,\ldots,n
#   threshold - threshold limit for the alarm and should be equal
#   to the desired ARL
# cusum -- boolean if TRUE then CUSUM otherwise Shiryaev-Roberts
######################################################################


stcd <- function(x, y,t,radius,epsilon,areaA, areaAcapBk, threshold,cusum=FALSE) {
  #check that the vectors x,y,t are of the same length.
  n <- length(x)
  if ((length(y) != n) | (length(t) != n)) {
    stop("Vectors x,y,t not of same size.")
  }
  if (!all(diff(order(t)) == 1)) {
    stop("The vector of time points needs to be ascending in time. No ties allowed.")
  }
  
  res <- .C("SRspacetime", x=as.double(x), y=as.double(y), t=as.double(t), n=as.integer(n), radius=as.double(radius), epsilon=as.double(epsilon), areaA=as.double(areaA),areaAcapBk=as.double(areaAcapBk),cusum=as.integer(cusum), threshold=as.double(threshold),R=as.double(numeric(n)),idxFA=as.integer(-1),idxCC=as.integer(-1),PACKAGE="surveillance")

  #Indexing differences between C and R
  res$idxFA <- res$idxFA+1
  res$idxCC <- res$idxCC+1
  
  #Missing: compute which indices are part of the cluster.
  #--> Thais R-code
  
  return(list(R=res$R,idxFA=res$idxFA,idxCC=res$idxCC))
}
