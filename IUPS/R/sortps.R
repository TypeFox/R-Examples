
###create sortps function to do the sorting and matching

sortps <- function(D, N, N1, p2, L, Data ) {     
  idhold<-matrix(0,N,L+1)          #####use idhold to save the L id with nearst propensity score
  id1<-c(1:N1)                     #####  N1 subject has t=1 
  ps1<-as.vector(p2[1:N1])         ##### the first N1 socre from p2 are those PS with t=1 
  md1<-cbind(id1,ps1)
  md1<-md1[order(ps1),]            ##### sort the PS 
  idhold[1,]<-md1[1:(1+L),1]       ##### for the first one, we use L subject below it
  for (i in 2:(N1-1)){             ##### for 2:N1, we use the loop to find L subjects above or below it( for the second obs, there is only one above it...)
    for (j in 1:L){
      if (i-j<=0) break
      else s0=i-j}
    for (j in 1:L){
      if (i+j>=N1+1) break
      else sn<-i+j}
    choose1<-md1[s0:sn,]                               #####choose1 save the choosed id and ps
    diff1<-abs(md1[i,2]-choose1[,2])                   #####calculate the absolute differenence
    choose_dif1<-cbind(choose1[,1],diff1)
    choose_dif1<-choose_dif1[order(diff1),]            #####sort the absolute differenece
    idhold[i,]<-choose_dif1[1:(L+1),1]                 #####save ids for those M smallest difference
  }
  idhold[N1,]<-md1[(N1-L):N1,1]                         #####for subject N1, use the L subjects above it
  
  id0<-c((N1+1):N)                                      #####in Data or dd above, N1+1:N are those t=0
  ps0<-as.vector(p2[(N1+1):N])                          #####these scores are for t=0
  md0<-cbind(id0,ps0)
  md0<-md0[order(ps0),]                                  #####sort the ps for t=0
  idhold[(N1+1),]<-md0[1:(1+L),1]                        #####for the first one with t=0, we use L subjects below it
  for (i in 2:(N-N1-1)){                                 #####for 2:N-N1-1, we use the loop to find L subjects above or below it.(for the second,only one above it)
    for (j in 1:L){
      if (i-j<=0) break
      else s0=i-j}
    for (j in 1:L){
      if (i+j>=(N-N1)+1) break
      else sn<-i+j}
    choose0<-md0[s0:sn,]                               #####choose1 save the choosed id and ps
    diff0<-abs(md0[i,2]-choose0[,2])                   #####calculate the absolute difference
    choose_dif0<-cbind(choose0[,1],diff0)
    choose_dif0<-choose_dif0[order(diff0),]            #####sort the absolute difference
    idhold[(N1+i),]<-choose_dif0[1:(L+1),1]            #####save ids for those M smallest difference
  }
  idhold[N,]<-md0[(N-N1-L):(N-N1),1]                   #####for subject N, use the L subjects above it
  id<-t(t(c(md1[,1],md0[,1])))                           
  idhold_2<-cbind(id, idhold[,1:(L+1)])                #####combine the id of the subject with the ids of it's L nearst ps 
  idhold_2<-idhold_2[order(id),]                       #####after sorting by id, they are in the same sequence before sorting the PS
  
  return (idhold_2)
}