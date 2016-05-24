"jaccard" <-
function(regmat){
  nart <- ncol(regmat)
  jdist <- rep(0, nart*nart)
  dim(jdist) <- c(nart,nart)
  reg.col.sum <- apply(regmat,2,sum)
  reg.aggrement <- t(regmat) %*% regmat
  
  jdist <-  1- reg.aggrement / (reg.col.sum-t(t(reg.aggrement)-reg.col.sum))
 
  jdist
}


# "jaccard" <-
# function(regmat){
#   nart <- ncol(regmat)
#   jdist <- rep(0, nart*nart)
#   dim(jdist) <- c(nart,nart)
#   for (i in 1:(nart-1)){
# #    cat("Row ",i,"\n")
#     for (j in (i+1):nart){
#       jdist[j,i] <- jdist[i,j] <-  1-sum(regmat[,i]+regmat[,j]>1.99)/sum(regmat[,i]+regmat[,j]>0.99)
#       if (is.na(jdist[i,j]))
#         cat("Warning! NA at i=",i,", j=", j,"\n")
#     }
#   }
#   jdist
# }
