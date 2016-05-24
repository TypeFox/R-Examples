##################################################
sim <-function (net,P,...) {
UseMethod("sim")
}

##################################################
sim.MLPnet <- function(net,P,...) {
   if (class(net)!="MLPnet") {
      stop("Your net parameter does not belong to the MLPnet class. Are you aware that the result from the train function is now a list instead of a net? Check parameters and try again");
   }
   P <- as.matrix(P)
   ytrans <- matrix(0, nrow=length(net$layer[[length(net$layer)]]), ncol=nrow(P))
   ytrans <- .Call("sim_Forward_MLPnet", net, t(P), ytrans, .GlobalEnv, PACKAGE="AMORE")
   return(t(ytrans))
}
###############################################################################################

train <- function(net, P, T, Pval=NULL, Tval=NULL, error.criterium="LMS", report=TRUE, n.shows, show.step, Stao=NA, prob=NULL, n.threads=0L) {
   if (class(net)!="MLPnet") {
      stop("Your net parameter does not belong to the MLPnet class. Are you aware that the result from the train function is now a list instead of a net? Check parameters and try again");
   }
   P <- as.matrix(P)
   T <- as.matrix(T)

   epoch.show.step <- 0
   n.muestras <- nrow(P)

   net$deltaE$fname <- as.integer(0)  # custom case
   if(error.criterium=="LMS") { 
     net$deltaE$fname <- as.integer(1)
     net$deltaE$f <- deltaE.LMS
   } else if(error.criterium=="LMLS") { 
     net$deltaE$fname <- as.integer(2)
     net$deltaE$f <- deltaE.LMLS
   } else if(error.criterium=="TAO") { 
     if (missing(Stao)) {
        stop("You should enter the value of Stao")
     } else {
	net$deltaE$fname <- as.integer(3)
	net$deltaE$f    <- deltaE.TAO
        net$deltaE$Stao <- Stao
     }
   }

   method <- net$neurons[[1]]$method

   if (method =="ADAPTgd") {
      train.method <- ADAPTgd.MLPnet
   } else if (method =="ADAPTgdwm") {
      train.method <- ADAPTgdwm.MLPnet
   } else if (method =="BATCHgd") {
      train.method <- BATCHgd.MLPnet
   } else if (method =="BATCHgdwm") {
      train.method <- BATCHgdwm.MLPnet
   }

   if (is.null(prob)) {
      if (!is.null(Pval) & !is.null(Tval)) {
	Merror <- matrix(NA, ncol=2, nrow=n.shows)
         Pval <- as.matrix(Pval)
         Tval <- as.matrix(Tval)
         min.error.val <- Inf
         bestnet <- net
         for (idx.show in 1:n.shows) {
            net <- train.method(net, P, T, show.step, n.threads=n.threads)
            P.sim    <- sim.MLPnet(net,P)
            Pval.sim <- sim.MLPnet(net,Pval) 
            if(error.criterium=="LMS") { 
               error     <- error.LMS(list(prediction=P.sim,    target=T    ))
               error.val <- error.LMS(list(prediction=Pval.sim, target=Tval ))
            } else if(error.criterium=="LMLS") { 
               error     <- error.LMLS(list(prediction=P.sim,    target=T    ))
               error.val <- error.LMLS(list(prediction=Pval.sim, target=Tval ))
            } else if(error.criterium=="TAO") {                             
               error.aux  <- error.TAO(list(prediction=P.sim, target=T, net=net))
               error      <- error.aux$perf
               new.tao    <- error.aux$Stao
               error.val  <- error.TAO(list(prediction=Pval.sim, target=Tval, net=net))$perf
               cat("Stao:", new.tao, " ")
            }
            Merror [idx.show,] <- c(error,error.val)

            if((!is.finite(error)) | (!is.finite(error.val))){
               stop("ERROR: Non-finite error found (Divergent). Try reducing the learning rate and/or the momentum");
            }  

            if (error.val <= min.error.val ) {
               min.error.val <- error.val
               bestnet <- net      
               cat(paste("index.show:", idx.show, error.criterium,"\tTRAIN:",error,"\tVAL:",error.val,"\t BEST NET\n", sep=" "))
            } else {
               cat(paste("index.show:", idx.show, error.criterium,"\tTRAIN:",error,"\tVAL:",error.val,"\n", sep=" "))
            }
         }
         net <- bestnet
      } else {
	Merror <- matrix(NA, ncol=1, nrow=n.shows)
         for (idx.show in 1:n.shows) {
            net <- train.method(net, P, T, show.step, n.threads=n.threads)
            if (report) {
		auxReport <-  training.report(net, P, T, idx.show, error.criterium)
		net$other.elements$Stao <- auxReport$new.tao
		Merror [idx.show,1] <- auxReport$error
               if(!is.finite(auxReport$error)){
                  stop("ERROR: Non-finite error found (Divergent). Try reducing the learning rate and/or the momentum");
               } 
            }else if(!all(is.finite(net$neurons[net$layers[[2]][1]][[1]]$weights))){
               stop("ERROR: Non-finite weight found (Divergent). Try reducing the learning rate and/or the momentum");
            }
         }
     }
   } else {
      if (!is.null(Pval) & !is.null(Tval)) {
	Merror <- matrix(NA, ncol=2, nrow=n.shows)
         Pval <- as.matrix(Pval)
         Tval <- as.matrix(Tval)
         min.error.val <- Inf
         bestnet <- net
         for (idx.show in 1:n.shows) {
            orden <- sample(1:n.muestras, n.muestras, replace=TRUE , prob=prob)
            net   <- train.method(net, P[orden, , drop=FALSE], T[orden, , drop=FALSE], show.step, n.threads=n.threads)
            P.sim    <- sim.MLPnet(net,P)
            Pval.sim <- sim.MLPnet(net,Pval) 
            if(error.criterium=="LMS") { 
               error     <- error.LMS(list(prediction=P.sim,    target=T    ))
               error.val <- error.LMS(list(prediction=Pval.sim, target=Tval ))
            } else if(error.criterium=="LMLS") { 
               error     <- error.LMLS(list(prediction=P.sim,    target=T    ))
               error.val <- error.LMLS(list(prediction=Pval.sim, target=Tval ))
            } else if(error.criterium=="TAO") {                             
               error.aux  <- error.TAO(list(prediction=P.sim, target=T, net=net))
               error      <- error.aux$perf
               new.tao    <- error.aux$Stao
               error.val  <- error.TAO(list(prediction=Pval.sim, target=Tval, net=net))$perf
               cat("Stao:", new.tao, " ")
            }
            Merror [idx.show,] <- c(error,error.val)
            
            if((!is.finite(error)) | (!is.finite(error.val))){
               stop("ERROR: Non-finite error found (Divergent). Try reducing the learning rate and/or the momentum");
            }  

            if (error.val <= min.error.val ) {
               min.error.val <- error.val
               bestnet <- net      
               cat(paste("index.show:", idx.show, error.criterium,"\tTRAIN:",error,"\tVAL:",error.val,"\t BEST NET\n", sep=" "))
            } else {
               cat(paste("index.show:", idx.show, error.criterium,"\tTRAIN:",error,"\tVAL:",error.val,"\n", sep=" "))
            }
         }
         net <- bestnet
      } else {
   	Merror <- matrix(NA, ncol=1, nrow=n.shows)
         for (idx.show in 1:n.shows) {
            orden <- sample(1:n.muestras, n.muestras, replace=TRUE , prob=prob)
            net <- train.method(net, P[orden, , drop=FALSE], T[orden, , drop=FALSE], show.step, n.threads=n.threads)
            if (report) {
		auxReport <-  training.report(net, P, T, idx.show, error.criterium)
		net$other.elements$Stao <- auxReport$new.tao
		Merror [idx.show,1] <- auxReport$error
               if(!is.finite(auxReport$error)){
                  stop("ERROR: Non-finite error found (Divergent). Try reducing the learning rate and/or the momentum");
               } 
            }else if(!all(is.finite(net$neurons[net$layers[[2]][1]][[1]]$weights))){
               stop("ERROR: Non-finite weight found (Divergent). Try reducing the learning rate and/or the momentum");
            }
         }
     }
   }
   return(list(net=net,Merror=Merror))

}


###############################################################################################
training.report <- function(net,P,T, idx.show, error.criterium) {


########### BEGIN do not delete ##########
   if (class(net)!="MLPnet") {
      stop("Your net parameter does not belong to the MLPnet class. Are you aware that the result from the train function is now a list instead of a net? Check parameters and try again");
   }
   new.tao      <- NA

########### END do not delete ############

   P.sim <- sim.MLPnet(net,P)
#          par(mfrow=c(1,2))
#          plot(P,T, col="red", pch="*", ylim=range(rbind(T,P.sim)))
#          points(P,P.sim, col="blue", pch="+")
#          plot(P, ideal, col="red", pch=".", ylim=range(rbind(ideal,P.sim)))
#          points(P,P.sim, col="blue", pch=".")
   if(error.criterium=="LMS") { 
           error <- error.LMS(list(prediction=P.sim, target=T))
   } else if(error.criterium=="LMLS") { 
           error <- error.LMLS(list(prediction=P.sim, target=T))

########### BEGIN do not delete (only minor changes allowed) ##########
   } else if(error.criterium=="TAO") {                             
           error.aux <- error.TAO(list(prediction=P.sim, target=T, net=net))
           error     <- error.aux$perf
           new.tao   <- error.aux$Stao
           cat("Stao:", new.tao, " ")
   }
########### END do not delete ############

   cat(paste("index.show:", idx.show, error.criterium,error,"\n", sep=" "))

########### BEGIN do not delete ##########
return(list(error=error,new.tao=new.tao))
########### END do not delete ############
}
