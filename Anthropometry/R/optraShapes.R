optraShapes <- function(array3D,n,c,numClust,ic1,ic2,nc,an1,an2,ncp,d,itran,live,indx){

#If cluster L is updated in the last quick-transfer stage, it belongs to the live set throughout this stage. 
#Otherwise, at each step, it is not in the live set if it has not been updated in the last M optimal transfer 
#steps:
 for(l in 1 : numClust){
  if ( itran[l] == 1){
   live[l] = n + 1
  }
 }

 for(i in 1 : n){
  indx = indx + 1
  l1 = ic1[i]
  l2 = ic2[i]
  ll = l2
 
  #If point I is the only member of cluster L1, no transfer:
  if(1 < nc[l1]){ 
   #If L1 has not yet been updated in this stage, no need to re-compute D(I):
   if( ncp[1] != 0 ){
    de = (riemdist(array3D[,,i], c[,,l1]))^2
    d[i] = de * an1[l1] 
   }

  #Find the cluster with minimum R2:
  da = (riemdist(array3D[,,i], c[,,l2]))^2
  r2 = da * an2[l2]

  for(l in 1 : numClust){
  #If LIVE(L1) <= I, then L1 is not in the live set. If this is true, we only need to consider clusters that 
  #are in the live set for possible transfer of point I (Step 4b). Otherwise, we need to consider all possible 
  #clusters (Step 4a):
   if( ( i < live[l1] || i < live[l2] ) && l != l1 && l != ll ){
    rr = r2 / an2[l]
    dc = (riemdist(array3D[,,i], c[,,l]))^2
    if(dc < rr)
     r2 = dc * an2[l]
     l2 = l
    }
  }

  #If no transfer is necessary, L2 is the new IC2(I):
  if(d[i] <= r2){
   ic2[i] = l2
   #Update cluster centers, LIVE, NCP, AN1 and AN2 for clusters L1 and L2, and update IC1(I) and IC2(I):
  }else{
    indx = 0
    live[l1] = n + i
    live[l2] = n + i
    ncp[l1] = i
    ncp[l2] = i
    al1 = nc[l1]
    alw = al1 - 1
    al2 = nc[l2]
    alt = al2 + 1
       
    nc[l1] = nc[l1] - 1
    nc[l2] = nc[l2] + 1
    an2[l1] = alw / al1
     
    if(1 < alw){
     an1[l1] = alw / ( alw - 1 )
    }else{
     an1[l1] = Inf
     }
     an1[l2] = alt / al2
     an2[l2] = alt / ( alt + 1 )
     ic1[i] = l2
     ic2[i] = l1

     c[,,l1] = procGPA(array3D[, , ic1 == l1], distances = TRUE, pcaoutput = TRUE)$mshape
     c[,,l2] = procGPA(array3D[, , ic1 == l2], distances = TRUE, pcaoutput = TRUE)$mshape
   }
  }

    if( indx == n ){ #indx is the number of steps since a transfer took place. 
     return(list(c, ic1, ic2, nc, an1, an2, ncp, d, itran, live, indx))
    }

 }

  #ITRAN(L) = 0 before entering QTRAN. Also, LIVE(L) has to be decreased by M before re-entering OPTRA:
  for (l in 1 : numClust){
   itran[l] = 0
   live[l] = live[l] - n
  }

  return(list(c, ic1, ic2, nc, an1, an2, ncp, d, itran, live, indx))
}