qtranShapes <- function(array3D,n,c,ic1,ic2,nc,an1,an2,ncp,d,itran,indx){

 ncp <- ncp  
 #In the optimal transfer stage, NCP(L) indicates the step at which cluster L is last updated. 
 #In the quick transfer stage, NCP(L) is equal to the step at which cluster L is last updated plus M:
 icoun = 0
 istep = 0

 for (i in 1 : n){
  icoun = icoun + 1
  istep = istep + 1      
  l1 = ic1[i]
  l2 = ic2[i]

  #If point I is the only member of cluster L1, no transfer:
   if (1 < nc[l1]){

    #If NCP(L1) < ISTEP, no need to re-compute distance from point I to cluster L1. 
    #Note that if cluster L1 is last updated exactly M steps ago, we still need to 
    #compute the distance from point I to cluster L1:
    if(istep <= ncp[l1]){
     da = (riemdist(array3D[,,i], c[,,l1]))^2
     d[i] = da * an1[l1] 
    }

     #If NCP(L1) <= ISTEP and NCP(L2) <= ISTEP, there will be no transfer of point I at this step:
     if ( istep < ncp[l1] | istep < ncp[l2] ){
      r2 = d[i] / an2[l2]
      dd = (riemdist(array3D[,,i], c[,,l2]))^2

      #Update cluster centers, NCP, NC, ITRAN, AN1 and AN2 for clusters L1 and L2. Also update IC1(I) and 
      #IC2(I). Note that if any updating occurs in this stage, INDX is set back to 0:
      if( dd < r2 ){
       icoun = 0
       indx = 0
       itran[l1] = 1
       itran[l2] = 1
       ncp[l1] = istep + n
       ncp[l2] = istep + n
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
      }else{
       }
     }
   }

   #If no re-allocation took place in the last n steps, return:
   if( icoun == n ){
    return(list(c, ic1,ic2, nc, an1, an2, ncp, d, itran, indx, icoun))
   }
 }

 return(list(c, ic1,ic2, nc, an1, an2, ncp, d, itran, indx, icoun))
}