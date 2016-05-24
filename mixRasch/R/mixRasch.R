`mixRasch` <-
function(data, steps, max.iter=50, conv.crit=.001, model="RSM", n.c=1, 
                     class, metric, info.fit=TRUE, treat.extreme=0.3, maxchange=1.5, maxrange=c(-4,4), as.LCA=FALSE){

 start <- Sys.time()

 if(as.LCA & n.c==1) { warning("as.LCA set to FALSE. n.c must be greater than 1 for LCA analysis. \n")
                       as.LCA <- FALSE}

 # Removes degen item and people
   temp     <- degenset(data)
   u.data   <- temp[[1]]
   baditems <- temp[[2]]
   rm(temp)

 max.change <- array()  # keeps max change in parameter by class
 converge.flag <- TRUE

 # Initializes class
 if(n.c > 1) { if(missing(class)) { class <- t(sapply(u.data$n.x,function(XXX){ xxx <- table(factor(sample(n.c,XXX,
                                                                                         replace=TRUE), levels=1:n.c))
                                                                              xxx/sum(xxx)}))
                                  }
             } else class <- matrix(1,ncol=1,nrow=length(u.data$n.x))

 lik.table    <- class                     # will hold log likelihoods
 p.c          <- colSums(class)/sum(class) # probability of class

 LatentClass <- list()

 for(c in 1:n.c){ 
 
 i.s <- item.stats(u.data$x.i, u.data$n.x*class[,c], steps, model=model)

 LatentClass[[c]] <- list(i.stat   = i.s,
                          person.par   = list(theta = log( (u.data$r - 0)/(sum(! is.na(i.s$S.ih)) - u.data$r) )),
                          item.par = init.item(i.s))

 if(as.LCA) LatentClass[[c]]$person.par$theta <- rep(0,length(LatentClass[[c]]$person.par$theta))

 } # end c loop
 rm(i.s)

 
 #########################
 # Begin Iterations
 #########################
 if(n.c == 1){
 c <- 1
 for(iter in 1:max.iter){
   old.item.par  <- LatentClass[[c]]$item.par
   old.theta     <- LatentClass[[c]]$person.par$theta

   # New item level delta
   LatentClass[[c]] <- item.delta(LatentClass[[c]],u.data,class[,c],steps, as.LCA, maxchange, maxrange)
 
   # New step values
   if(! steps == 1) LatentClass[[c]] <- item.tau(LatentClass[[c]],u.data,class[,c],steps,model, maxchange, maxrange)
 
   # New Theta
   LatentClass[[c]] <- new.theta(LatentClass[[c]],u.data,class[,c],steps)

   # Record parameter change to check for convergence
   delta.i.diff   <- LatentClass[[c]]$item.par$delta.i - old.item.par$delta.i
   tau.diff       <- LatentClass[[c]]$item.par$tau - old.item.par$tau 
   delta.diff     <- LatentClass[[c]]$item.par$delta - old.item.par$delta 
   theta.diff     <- LatentClass[[c]]$person.par$theta - old.theta
      
   max.change[c] <- max(c(abs(delta.i.diff),abs(tau.diff),abs(theta.diff)),na.rm=TRUE)
   cat("Iteration: ", iter, ", Largest Parameter Change: ", max.change[c], "\n", sep="")
   flush.console()
 if(all(max.change < conv.crit,na.rm=TRUE)) break
 } # iterations
 } # end 1 class loop

 ######
 # Begin n.c > 1
 ######
 else {
 
 latent.degen <- array(dim=c(ncol(u.data$x.i),n.c))

 for(iter in 1:max.iter){
 old.class <- class
 
 for(c in 1:n.c){
 
   old.item.par  <- LatentClass[[c]]$item.par
   old.theta     <- LatentClass[[c]]$person.par$theta

   # New item level delta
   LatentClass[[c]] <- item.delta(LatentClass[[c]],u.data,class[,c],steps, as.LCA, maxchange, maxrange)

   #cat("Finish new delta \n")

   # New step values
   if(! steps == 1) LatentClass[[c]] <- item.tau(LatentClass[[c]],u.data,class[,c],steps,model, maxchange, maxrange)
   # cat("Finish new steps \n")

   # New Theta
   if(! as.LCA) LatentClass[[c]] <- new.theta(LatentClass[[c]],u.data,class[,c],steps)

   # Record parameter change to check for convergence
   delta.i.diff   <- LatentClass[[c]]$item.par$delta.i - old.item.par$delta.i
   tau.diff       <- LatentClass[[c]]$item.par$tau - old.item.par$tau 
   delta.diff     <- LatentClass[[c]]$item.par$delta - old.item.par$delta 
   theta.diff     <- LatentClass[[c]]$person.par$theta - old.theta

   Pxji <- P.xji(LatentClass[[c]],u.data,steps)
   lik.table[,c] <- loglik(Pxji, u.data$x.i)

   max.change[c] <- max(c(abs(delta.i.diff),abs(tau.diff),abs(theta.diff)),na.rm=TRUE)
 } # end c loop

 
 if(n.c > 1){
  class <- t(t(exp(lik.table))*p.c)
  class <- class/rowSums(class)

  p.c <- colSums(class*u.data$n.x)
  p.c <- p.c/sum(p.c)
 
 for(c in 1:n.c) { LatentClass[[c]]$i.stat <- item.stats(u.data$x.i, u.data$n.x*class[,c], steps, model=model)
                   latent.degen[,c] <- (LatentClass[[c]]$i.stat$S.ih[1,] > 1) & (LatentClass[[c]]$i.stat$n.ni - LatentClass[[c]]$i.stat$S.ih[1,] > 1)                   
                 }
 }
 cat("Iteration: ", iter, ", Largest Parameter Change: ", max(max.change), "\n", sep="")
 flush.console()
 if(all(max.change < conv.crit,na.rm=TRUE)) break
 if(any(! latent.degen)) { #warning("Estimation halted. Perfect response vectors in at least one latent class.")
     warning("Perfect response vectors in at least one latent class.")
                           #converge.flag <- FALSE
                           #break
                         }
 } # iterations
 } # end n.c > 1 section

 #####################################
 # Calculates standard errors and item fit
 #####################################

 LatentClass[[c]]$item.par$SE.delta.i <- array()
 LatentClass[[c]]$item.par$SE.tau <- array()
 LatentClass[[c]]$person.par$SE.theta <- array()

 for(c in 1:n.c){

   Pxji <- P.xji(LatentClass[[c]],u.data,steps)
   if(info.fit == TRUE) lik.table[,c] <- loglik(Pxji, u.data$x.i) 

   for(i in 1:LatentClass[[c]]$i.stat$n.i){
    LatentClass[[c]]$item.par$SE.delta.i[i] <-  (- d.delt(LatentClass[[c]]$i.stat$n.ni[i],
                                                LatentClass[[c]]$i.stat$Si[i],
                                                LatentClass[[c]]$i.stat$steps,(u.data$n.x*class[,c]),
                                                u.data$resp[,i], matrix(Pxji[,,i],nrow=steps))$d2)^(-.5)
   }

   if(! steps ==1) {
    LatentClass[[c]]$item.par$SE.tau <- (-1*d.tau(LatentClass[[c]]$i.stat$Tx, LatentClass[[c]]$i.stat$S.ih, 
                                         LatentClass[[c]]$i.stat$steps,(u.data$n.x*class[,c]),
                                         u.data$resp,u.data$n.unique, Pxji, model)$d2)^(-.5)
    LatentClass[[c]]$item.par$SE.tau <- matrix(LatentClass[[c]]$item.par$SE.tau,nrow=steps)
   } # ! stpes == 1

   LatentClass[[c]]$person.par$SE.theta <- (-1*d.v(Pxji, u.data$r, u.data$resp)$d2)^(-.5)
    
   #
   # fit called here
   #
   hold <- infit.outfit(Pxji,u.data$n.x*class[,c],u.data$resp,u.data$x.i)
   LatentClass[[c]]$item.par$in.out <- hold[[1]]
   LatentClass[[c]]$person.par$in.out <- hold[[2]]
   rm(hold)
         
 } # class loop for SE

 ######################################
 # End SE and Rasch Fit
 ######################################
 # Start Informaton Fit AIC
 ######################################

 if(info.fit == TRUE){
   # First line is p.c parameter df
   # Second line is item parameters per class
   # Third line is theta
   # for formal eqs see ROST, 1990, APM
   N.parms <- (n.c - 1) # prob of class
   N.parms <- N.parms + n.c*length(unique(LatentClass[[1]]$person.par$theta)) # unique theta
   #N.parms <- N.parms + (n.c - 1)*nrow(u.data$x.i) # prob individuals class membership
               
   if(model=="PCM") N.parms <- N.parms + n.c*(sum(! is.na(LatentClass[[1]]$item.par$delta)) - 1)
   if(model=="RSM") N.parms <- N.parms + n.c*(nrow(LatentClass[[1]]$item.par$delta) - 1 + 
                                              ncol(LatentClass[[1]]$item.par$delta) - 1)

   if(as.LCA) N.parms <- N.parms - n.c*(length(unique(LatentClass[[1]]$person.par$theta))) + 1*n.c

   tlik <- sum(log(colSums(t(exp(lik.table))*p.c))*u.data$n.x)  

   AIC <- -2*tlik + 2*N.parms 
   BIC <- -2*tlik + log(sum(u.data$n.x))*N.parms 
   CAIC <- BIC + N.parms

   fit <- list(AIC=AIC, BIC=BIC, CAIC=CAIC, loglik=tlik, N.parms=N.parms, N.persons=sum(u.data$n.x))
 } else fit <- "Not Requested"


 # Create person reports
   LatentClass <- person.reports(LatentClass,n.c,u.data, treat.extreme, max.iter, conv.crit, steps, as.LCA)
   matcher <- match(u.data$people,u.data$x.x) 
   class   <- class[matcher,]
   
   for(c in 1:n.c){
   #
   # basic item stats, by modal class
   #
   if(n.c > 1){
     modalClassFlag <- class[,c] > .5
   } else modalClassFlag <- rep(TRUE,nrow(data))	 
   itemDesc <- itemDescriptives(LatentClass[[c]]$person.par$theta[modalClassFlag],data[modalClassFlag,])
   LatentClass[[c]]$item.par$itemDescriptives <- itemDesc
   rm(itemDesc)
   }
   
 ##################
 # Final format and class statements
 ##################
 
 if(iter == max.iter) converge.flag <- FALSE 

 if(n.c == 1) { class <- 1 }
 if(steps == 1) model <- "RM"

 if(n.c > 1) model <- paste("mix",model,sep="")

 finish <- Sys.time()

if(n.c > 1)  out <- list(LatentClass = LatentClass, max.change = max.change, class = class, iter = iter, converge.flag = converge.flag,
                        info.fit = fit, model = model, removed.items = baditems, run.time = finish - start)

else        {out <- list(LatentClass[[1]][[1]], LatentClass[[1]][[2]], LatentClass[[1]][[3]], 
                        max.change = max.change, class = class, iter = iter, converge.flag = converge.flag,
                        info.fit = fit, model = model, removed.items = baditems, run.time = finish - start)
             names(out)[1:3] <- names(LatentClass[[1]])
             }
class(out) <- "Rasch"
return(out)
} # End mixRasch

