trimmedoid <- function(D,numClust,alpha,niter,algSteps=7,verbose){
 
 n <- dim(D)[1]
 no.trim <- floor(n*(1-alpha))
 ll <- (1:numClust)
 ind <- (1:n)
 dist <- ind
	
 #Initialize the objective function by a large enough value:
 vopt <- 100000000

 #Ramdon restarts:
 for(iter in 1:niter){
  if(verbose){ 
   cat("new iteration")
   print(iter)
  }
   
  #Randomly choose the numClust initial centers:
  cini <- sample(1:n,size=numClust,replace=FALSE)
		
  #C-steps: step 2:
  for(t in 1:algSteps){
   disti=c()
   ind=c()
   #Distances of each data point to its closest medoid:
   for(h in 1:n){
    for(k in 1:numClust){
     ll[k] <- D[h,cini[k]]
    }
    disti[h] <- min(ll)
    ind[h] <- which.min(ll)
   }

   #Modified data (Dmod) with the non-trimmed points and a vector indqq equal to the clusters allocations:
   qq <- (1:n)[disti<=sort(disti)[no.trim]]
   Dmod <- D[qq,qq]
   indqq=ind[qq]

   if(length(unique(indqq))<numClust) {t=algSteps}
   else{
    #Calculus of the new k centers:
    for(k in 1:numClust){
     ni <- sum(indqq==k)
     if(ni>1){
      #cini[k,]<-apply(xmod[xmod[,p+1]==k,1:p],2,mean)
      #convert Dmod in dist:
      Dmodkv=as.dist(Dmod[indqq==k,indqq==k])
      rpam=pam(Dmodkv,k=1) #here initial medoid, obtained by build.
      medk=rpam$id.med
     }
     if(ni==1){
      medk=1
     }
     if (ni==0){
      medk=1
     }
      aux=qq[indqq==k]
      cini[k]=aux[medk]
    }

      obj <- 0
      for(l in 1:no.trim){
       for (k in 1:numClust){
        ll[k] <- D[l,cini[k]]
       }
       obj <- obj+ min(ll)
      }
       if(iter<10){
         if(verbose){ 
          print(obj/no.trim)
         }  
       }
       rm(obj)

   }		
 }
	
  #Calculus of the trimmed k-variance:
  obj <- 0
  for(l in 1:no.trim){
   for(k in 1:numClust){
    ll[k] <- D[l,cini[k]]
   }
   obj <- obj+ min(ll)
  }
  obj <- obj/no.trim

  #Change the optimal value and the optimal centers (copt) if a reduction in the objective function happens:
  if (obj < vopt){
   vopt <- obj
   #Improvements in the objective functions are printed:
   if(verbose){
    cat("optimal")
    print(vopt)
   }
   copt <- cini
  } 
 }

 #Obtain the final cluster allocations (this is necesary, because a final cluster assigment 
 #movement is possible):
 asig <- rep(0,n)
 #Distances of each data point to its closest medoid:
 for(h in 1:n){
  for(k in 1:numClust){
   ll[k] <- D[h,cini[k]]
  }
  disti[h] <- min(ll)
  ind[h] <- which.min(ll)
 }

			 
 #a vector indqq equal to the clusters allocations:
 qq <- (1:n)[disti<=sort(disti)[no.trim]]
 indqq=ind[qq]

 #Assign every observation to each cluster and 0 for the trimmed observations:
 asig[qq]=indqq

 #Between clusters sum of distances:
 b=0
 for(k in 1:(numClust-1)){
  for(j in (k+1):numClust){
   b=b+D[cini[k],cini[j]] 
  }
 }

 ##ch goodness index:
 ch=b*(no.trim-numClust)/(vopt*no.trim*(numClust-1))

 rt=list(vopt=vopt,copt=cini,asig=asig,ch=ch,Dmod=Dmod,qq=qq) 
 return(rt)
}