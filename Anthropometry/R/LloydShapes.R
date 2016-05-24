LloydShapes <- function(array3D,numClust,algSteps=10,niter=10,stopCr=0.0001,simul,verbose){
#,computCost  

 time_iter <- list()       #List to save the real time in which each iteration ends.
 comp_time <- c()          #List to save the computational time of each iteration.   
 list_asig_step <- list()  #List to save the clustering obtained in each Nstep.
 list_asig <- list()       #List to save the optimal clustering obtained among all the Nstep of each iteration.
 vect_all_rate <- c()      #List to save the optimal allocation rate of each iteration.

 initials <- list()        #List to save the random initial values used by this Lloyd algorithm. Thus, 
                           #the Hartigan algorithm can be executed with these same values. 

 ll <- 1 : numClust
 dist <- matrix(0, dim(array3D)[3], numClust)
 
 if(verbose){
  print(Sys.time())
 }
 time_ini <- Sys.time()

 #Initialize the objective function by a large enough value:
 vopt <- 1e+08

 #Random restarts:
 for(iter in 1 : niter){
 
  obj <- list() #List to save the objective function (without dividing between n) of each Nstep.

  meanshapes <- 0 ; meanshapes_aux <- 0 ; asig <- 0
  mean_sh <- list()
  n <- dim(array3D)[3]
 
  if(verbose){ 
   cat("New iteration:")
   print(iter)

   cat("Optimal value with which this iteration starts:")
   print(vopt)
  }
  
  #Randomly choose the numClust initial centers:
  initials[[iter]] <- sample(1:n, numClust, replace = FALSE)
  if(verbose){
   cat("Initial values of this iteration:")
   print(initials[[iter]]) 
  }
  meanshapes <- array3D[, , initials[[iter]]] 		 
  meanshapes_aux <- array3D[, , initials[[iter]]] 

  #if(computCost){
    #time_ini_dist <- Sys.time() 
    #dist_aux = riemdist(array3D[,,1], y = meanshapes[,,1])
    #time_end_dist <- Sys.time()
    #cat("Computational cost of the Procrustes distance:") 
    #print(time_end_dist - time_ini_dist)
   #}
  
  for(step in 1 : algSteps){
    for(h in 1 : numClust){ 
     dist[,h] = apply(array3D[,,1:n], 3, riemdist, y = meanshapes[,,h])
    }
        
    asig = max.col(-dist)

    #if(computCost){
      #time_ini_mean <- Sys.time() 
      #meanshapes_aux[,,1] = procGPA(array3D[, , asig == 1], distances = TRUE, pcaoutput = TRUE)$mshape
      #time_end_mean <- Sys.time()
      #cat("Computational cost of the Procrustes mean:") 
      #print(time_end_mean - time_ini_mean)
     #}
    
    
     for(h in 1 : numClust){
      if(table(asig == h)[2] == 1){ 
       meanshapes[,,h] = array3D[, , asig == h]
       mean_sh[[step]] <- meanshapes
      }else{ 
        meanshapes[,,h] = procGPA(array3D[, , asig == h], distances = TRUE, pcaoutput = TRUE)$mshape
        mean_sh[[step]] <- meanshapes
       }
     }

   obj[[step]] <- c(0)
    for (l in 1 : n){  
     obj[[step]] <- obj[[step]] + dist[l,asig[l]]^2
    }
   obj[[step]] <- obj[[step]] / n
     
   list_asig_step[[step]] <- asig 

   if(verbose){  
    paste(cat("Clustering of the Nstep", step, ":\n"))
    print(table(list_asig_step[[step]])) 
   }
   if(verbose){
    if(iter <= 10){ 
     paste(cat("Objective function of the Nstep", step))
     print(obj[[step]]) 
    } 
   }  
     if(step > 1){
      aux <- obj[[step]] 
      aux1 <- obj[[step-1]]
      if( ((aux1 - aux) / aux1) < stopCr ){ 
       break         
      }
     }
   }#The algSteps loop ends here.

    #Calculus of the objective function (the total within-cluster sum of squares):
    obj1 <- 0
     for(l in 1 : n){
      obj1 <- obj1 + dist[l,asig[l]]^2
     }
    obj1 <- obj1 / n 

    #Change the optimal value and the optimal centers (copt) if a reduction in the objective function happens:
    if( obj1 > min(unlist(obj)) ){ 
     if( min(unlist(obj)) < vopt ){
      vopt <- min(unlist(obj)) 
     
      if(verbose){
       #Improvements in the objective functions are printed:
       cat("optimal")
       print(vopt)
      }
      optim_obj <- which.min(unlist(obj)) 
      copt <- mean_sh[[optim_obj]] #optimal centers.
      asig_opt <- list_asig_step[[optim_obj]] #to save the optimal clustering.
     }
    }else if(obj1 < vopt){
     vopt <- obj1
     
     if(verbose){
      #Improvements in the objective functions are printed:
      cat("optimal")
      print(vopt)
     }
     optim_obj <- which.min(unlist(obj)) 
     copt <- mean_sh[[optim_obj]] #optimal centers. 
     asig_opt <- list_asig_step[[optim_obj]]
    }

      time_iter[[iter]] <- Sys.time()

       if(iter == 1){
        comp_time[1] <- difftime(time_iter[[iter]], time_ini, units = "mins")
        if(verbose){
         cat("Computational time of this iteration: \n")
         print(time_iter[[iter]] - time_ini)
        } 
       }else{
         comp_time[iter] <- difftime(time_iter[[iter]], time_iter[[iter-1]], units = "mins")
         if(verbose){
          cat("Computational time of this iteration: \n")
          print(time_iter[[iter]] - time_iter[[iter - 1]])
         }   
        }   
  
  if(verbose){
   #In order to display the optimal clustering related to the optimal objective function, we have to find 
   #the step in which the optimal was obtained. This is provided by (which.min(unlist(obj))). 
   cat("Optimal clustering of this iteration: \n")
  } 
   optim_obj <- which.min(unlist(obj)) 
   list_asig[[iter]] <- list_asig_step[[optim_obj]]
  if(verbose){
   print(table(list_asig[[iter]]))
  }

   if(simul){
    #Allocation rate:
    as1 <- table(list_asig[[iter]][1:(n/2)])                 
    as2 <- table(list_asig[[iter]][seq(n/2 + 1,n)])
    if( max(as1) != n/2 & max(as2) != n/2 ){ 
     suma <- min(as1) + min(as2)
     all_rate <- 1 - suma / n
    }else if( (max(as1) == n/2 & max(as2) != n/2) || (max(as1) != n/2 & max(as2) == n/2) ){
      minim <- min(min(as1),min(as2))
      all_rate <- 1 - minim / n
     }else if( max(as1) == n/2 & max(as2) == n/2 ){
       all_rate <- 1
      }
    vect_all_rate[iter] <- all_rate 
    if(verbose){
     cat("Optimal allocation rate in this iteration:")
     print(all_rate)
    } 
   }
 }#The niter loop ends here.

 if(simul){
  dimnames(copt) <- NULL 
  return(list(asig=asig_opt,cases=copt,vopt=vopt,compTime=comp_time,
              AllRate=vect_all_rate,initials=initials))
 }else{
   return(list(asig=asig_opt,cases=copt,vopt=vopt,initials=initials))
  }
}  
