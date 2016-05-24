trimmedLloydShapes <- function(array3D,n,alpha,numClust,algSteps=10,niter=10,stopCr=0.0001,verbose){

 no.trim <- floor(n*(1-alpha)) #Elements that left after the trimmed procedure.
 vect_dist <- c() #Ancillary vector for the trimmed procedure.

 time_iter <- list()       #List to save the real time in which each iteration ends.
 comp_time <- c()          #List to save the computational time of each iteration.   
 list_asig_step <- list()  #List to save the clustering obtained in each Nstep.
 list_asig <- list()       #List to save the optimal clustering obtained among all the Nstep of each iteration.
 vect_all_rate <- c()      #List to save the optimal allocation rate of each iteration.

 initials <- list()        #List to save the random initial values used by this Lloyd algorithm. Thus, 
                           #the Hartigan algorithm can be executed with these same values.  

 trimms <- list()          #List to save the trimmed women of each iteration.
 
 trimms_iter <- c()        #Vector to find the iteration where the optimum has reached and therefore, to 
                           #identify the trimmed women of that iteration.  
 
 betterNstep <- c()        #Vector to find the Nstep of the iteration where the optimum has reached and 
                           #therefore, to identify the trimmed women of that iteration. 
  
 ll <- 1 : numClust
 dist <- matrix(0, n, numClust)

 if(verbose){
  print(Sys.time())
 }  
 time_ini <- Sys.time()

 #Initialize the objective function by a large enough value:
 vopt <- 1e+08

 #Random restarts:
 for(iter in 1 : niter){
 
  trimms[[iter]] <- list() 
   
  obj <- list() #List to save the objective function (without dividing between n) of each Nstep.

  meanshapes <- 0 ; asig <- 0
  mean_sh <- list()
 
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
		
   for(step in 1 : algSteps){
    for(h in 1 : numClust){
     dist[,h] = apply(array3D[,,1:n], 3, riemdist, y = meanshapes[,,h])
    }
        
    asig = max.col(-dist)

     #For the trimmed procedure:
     for(filas in 1 : n){
      vect_dist[filas] <- dist[filas,asig[filas]]
     }
      qq <- (1:n)[vect_dist <= sort(vect_dist)[no.trim]] 
      distmod <- dist[qq,]
      asigqq <- asig[qq] #asignations vector of the elements that left after the trimmed procedure.

    if(verbose){ 
     cat("Trimmed woman:")
     print(setdiff(1:dim(array3D)[3],qq))
     }
    trimms[[iter]][[step]] <- setdiff(1:dim(array3D)[3],qq)
    
      for(h in 1 : numClust){
       if(table(asigqq == h)[2] == 1){ 
        meanshapes[,,h] = array3D[, , asigqq == h]
        mean_sh[[step]] <- meanshapes
       }else{ 
         meanshapes[,,h] = procGPA(array3D[, , asigqq == h], distances = TRUE, pcaoutput = TRUE)$mshape
         mean_sh[[step]] <- meanshapes
        }
      }

      obj[[step]] <- c(0)
      for (l in 1:no.trim){
       obj[[step]] <- obj[[step]] + dist[l,asigqq[l]]^2
      }
      obj[[step]] <- obj[[step]] / no.trim
      list_asig_step[[step]] <- asigqq
 
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
     for(l in 1:no.trim){
      obj1 <- obj1 + dist[l,asigqq[l]]^2
     }
    obj1 <- obj1/no.trim

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
      asig_opt <- list_asig_step[[optim_obj]] 
      
      if(verbose){
       cat("Optimal iteration:")
       print(iter)
      }
      trimms_iter[iter] <- iter
      betterNstep <- which.min(unlist(obj))      
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
     
     if(verbose){
      cat("Optimal iteration:")
      print(iter)
     }
     
     trimms_iter[iter] <- iter
     betterNstep <- which.min(unlist(obj))
     
    }

      time_iter[[iter]] <- Sys.time()

       if(iter == 1){
        comp_time[1] <- difftime(time_iter[[iter]], time_ini, units = "mins")
        if(verbose){
         cat("Computational time of this iteration: \n")
         print(time_iter[[iter]] - time_ini)
        } 
       }else{
         comp_time[iter] <- difftime(time_iter[[iter]], time_iter[[iter - 1]], units = "mins")
         if(verbose){
          cat("Computational time of this iteration: \n")
          print(time_iter[[iter]] - time_iter[[iter - 1]])
         }  
        }   
  if(verbose){         
   cat("Optimal clustering of this iteration: \n")
  } 
   optim_obj <- which.min(unlist(obj))
   list_asig[[iter]] <- list_asig_step[[optim_obj]] 
  if(verbose){
   print(table(list_asig[[iter]]))
  } 

 }#The niter loop ends here.

 opt_iter <- trimms_iter[length(trimms_iter)]
 trimmed <- trimms[[opt_iter]][[betterNstep]]
 
 return(list(asig=asig_opt,cases=copt,vopt=vopt,trimmWomen=trimms,trimmsIter=trimms_iter,
             betterNstep=betterNstep,initials=initials,discarded=trimmed))
}
