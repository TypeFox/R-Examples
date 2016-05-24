q1q2l <-
function(bug, sims, ymean, hmean=NULL, MU, COV, P=NULL){
  if(is.null(colnames(sims)))
    stop("columns of sims can not be NULL. names should correspond to parameters")
  if(class(bug)!="tsbugs")
    stop("bug must be a object of class tsbugs")
  
  if(bug$info$variance=="CV")
    dpost<-dcvts
  if(bug$info$variance=="SV")
    dpost<-dsvts
  if(bug$info$variance=="RV")
    dpost<-drvts
  q1<-dpost(bug, sims=sims, ymean=ymean, hmean=hmean) 
  COV<-as.matrix(COV)
  if(bug$info$variance!="RV") q2<-dmvnorm(sims, mean=MU, sigma=COV, log=TRUE)
  if(bug$info$variance=="RV") q2<-dmvnb(bug, sims, MU, COV, P=data.matrix(P))
  return( data.frame(q1=q1, q2=q2, l=q1-q2) )
}
