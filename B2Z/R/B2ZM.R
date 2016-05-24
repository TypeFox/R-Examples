#######################################################
#This function is the interface between the user and  #
#the functions that have the implemented algorithms.  #
#######################################################

B2ZM <- function(data = NULL, priorBeta, priorQ, priorG, 
                 v, S, tauN.sh, tauN.sc, tauF.sh, tauF.sc,  
                 VN, VF, indep.model = FALSE, cred = 95,
                 sampler = c("MCMC", "BCLT", "IMIS", "SIR"),
                 sir.control = list(m = 10000),
                 bclt.control = list(m = 7000, sample_size = 2000),
                 imis.control = list(N0 = 6000, B = 600, M = 3000, it.max=16),
                 mcmc.control = list(NUpd = 10000, burnin = 1000, lag = 1, initial = NULL, Sigma.Cand = NULL, m = 1000),
                 figures = list(save = FALSE, type =c("ps", "eps","pdf", "png", "jpg"))) {

   if(is.null(data)) {stop("Data set not found.")}
   
   data <- as.matrix(data)

   if(!is.matrix(data)) {stop("Data set must be a 3-column matrix.")}
   if((is.matrix(data) & ncol(data)!=3)) {stop("Data set must be a 3-column matrix.")}
   if(is.character(data)) {stop("Data set must contain numeric elements.")}
   if(length(which(is.na(data)==TRUE))> 0) {stop("Data set contains missing values.")}
   if(length(which(data < 0))> 0) {stop("Data set contains negative values.")}
      
   if(!exists("priorBeta")){stop("Prior distribution for Beta was not declared.")}
   if(!exists("priorQ")){stop("Prior distribution for Q was not declared.")}
   if(!exists("priorG")){stop("Prior distribution for G was not declared.")}

   if(indep.model){
      S <- matrix(0, 2, 2)
      v <- 0

      if(!exists("tauN.sh")){stop("Shape parameter in the prior 
         distribution for TauN was not declared.")}
       
      if(!exists("tauN.sc")){stop("Scale parameter in the prior 
         distribution for TauN was not declared.")}

      if(!exists("tauF.sh")){stop("Shape parameter in the prior 
         distribution for TauF was not declared.")}

      if(!exists("tauF.sc")){stop("Scale parameter in the prior 
         distribution for tauF was not declared.")}

      if(!is.numeric(tauN.sh)){stop("'tauN.sh' must be numeric.")}
      if(!is.numeric(tauN.sc)){stop("'tauN.sc' must be numeric.")}
      if(!is.numeric(tauF.sh)){stop("'tauF.sh' must be numeric.")}
      if(!is.numeric(tauF.sc)){stop("'tauF.sc' must be numeric.")}

      if(tauN.sh <= 0){stop("'tauN.sh' must be greater than 0.")}
      if(tauN.sc <= 0){stop("'tauN.sc' must be greater than 0.")}
      if(tauF.sh <= 0){stop("'tauF.sh' must be greater than 0.")}
      if(tauF.sc <= 0){stop("'tauF.sc' must be greater than 0.")}
    }
   else{
       tauN.sh = tauN.sc = tauF.sh = tauF.sc = 0
       if(!exists("v")){stop("parameter v in the joint prior distribution for tauN, tauNF and tauF was not declared.")}
       if(!exists("S")){stop("parameter S in the joint prior distribution for tauN, tauNF and tauF was not declared.")}
       if(!is.matrix(S)){stop("object 'S' must be a matrix.")}
       if(!is.numeric(S)){stop("object 'S' must be a numeric matrix.")}
       if(nrow(S)!=ncol(S)){stop("object 'S' must be a square matrix.")}
       if(nrow(S)!=2){stop("object 'S' must be a 2x2 matrix.")}
       try(testingPD <- chol(S), silent=TRUE)
       if(!exists("testingPD")){stop("object 'S' must be a positive definite matrix.")}

       if(v<=1){stop("'v' must be greater than 1.")}
   }

   priors <- c(priorBeta, priorQ, priorG)
   name.dist <- character()
   par.dist <- character()

   for(i in 1:3){ 
      split1 <- unlist(strsplit(priors[i],  "\\(" ))
      name.dist[i] <- paste("d",split1[1],sep="")
      par.dist[i] <- unlist(strsplit(split1,  "\\)" ))[2]
   }

   if(!exists(name.dist[1], mode="function")){stop("invalid prior for Beta.")}
   if(!exists(name.dist[2], mode="function")){stop("invalid prior for Q.")}
   if(!exists(name.dist[3], mode="function")){stop("invalid prior for G.")}
   if(is.na(par.dist[1])){stop("missing parameters for the prior distribution of Beta.")}
   if(is.na(par.dist[2])){stop("missing parameters for the prior distribution of Q.")}
   if(is.na(par.dist[3])){stop("missing parameters for the prior distribution of G.")}
     
   if(!is.numeric(VF)){stop("'VF' must be numeric.")}
   if(VF < 0){stop("'VF' must be greater than 0.")}

   if(!is.numeric(VN)){stop("'VN' must be numeric.")}
   if(VN < 0){stop("'VN' must be greater than 0.")}

   indep <- indep.model  

   times <- c(data[,1])
   
   if (min(times)<=0) {stop("The observed time values must be greater than 0.")}
      
   order_times <- order(times)
   data[,1] <- data[order_times,1]
   data[,2] <- data[order_times,2]
   data[,3] <- data[order_times,3]
       
   times <- times[order_times]
   Y <- log(data[,2:3]) 

   if(cred <=0 || cred >=100) {stop("'cred' must be a number between 0 and 100")}      

   if(figures$save==TRUE){
      if(length(which(c("ps","eps","pdf","jpg","png")==figures$type))==0){
           stop("in 'figures': object 'type' must be one of the options: 'ps', 'eps', 'pdf', 'jpg' or 'png'")
      }
   }

   prval <- which.dist(priorBeta)
   indBeta <- prval[1]
   aBeta <- prval[2]
   bBeta <- prval[3]

   prval <- which.dist(priorQ)
   indQ <- prval[1]
   aQ <- prval[2]
   bQ <- prval[3]

   prval <- which.dist(priorG)
   indG <- prval[1]
   aG <- prval[2]
   bG <- prval[3]

   if(is.null(indBeta)){stop("invalid prior for Beta.")}
   if(is.null(indQ)){stop("invalid prior for Q.")}
   if(is.null(indG)){stop("invalid prior for G.")}

   verifying_prior_Beta(indBeta,aBeta,bBeta)
   verifying_prior_Q(indQ,aQ,bQ)
   verifying_prior_G(indG,aG,bG)

   if(sampler == "SIR"){
      if(is.null(sir.control$m)){stop("object 'm' not found.")}
      if(sir.control$m <= 0){stop("object 'm' must be greater than 0.")}
      if(!is.numeric(sir.control$m)){stop("object 'm' must be numeric.")}

      m <- sir.control$m        

      ans <- SIRB2zm (m, Y, times, VN, VF, indBeta, aBeta, bBeta, 
                      indQ, aQ, bQ, indG, aG, bG, v, S, tauN.sh, 
                      tauN.sc, tauF.sh, tauF.sc, indep, cred)
   }
    
   if(sampler == "BCLT"){
      if(is.null(bclt.control$m)){stop("object 'm' not found.")}
      if(bclt.control$m <= 0){stop("object 'm' must be greater than 0.")}
      if(!is.numeric(bclt.control$m)){stop("object 'm' must be numeric.")}
      if(is.null(bclt.control$sample_size)){stop("object 'sample_size' not found.")}
      if(bclt.control$sample_size <= 0){stop("object 'sample_size' must be greater than 0.")}
      if(!is.numeric(bclt.control$sample_size)){stop("object 'sample_size' must be numeric.")}
       
      m <- as.integer(bclt.control$m)
      sample_size <- bclt.control$sample_size
      

      ans <- BCLTB2zm(S, v, tauN.sh, tauN.sc, tauF.sh, tauF.sc, VN, 
                      VF, times, Y, indep, m, cred, indBeta, aBeta, bBeta,
                      indQ, aQ, bQ, indG, aG, bG, sample_size)

   }

   if(sampler == "IMIS"){
      if(!is.list(imis.control)){stop("'imis.control' must be a list containg the objects: 'N0', 'B', 'M' and 'it.max'.")}
        namesobj <- c("N0", "B", "M", "it.max")
        
      cont <- 0
      index <- numeric()
      indexT <- numeric()
      for(i in 1:length(names(imis.control))){
         index <- which(namesobj==names(imis.control)[i])
         cont <- cont + length(which(namesobj==names(imis.control)[i]))
         indexT <- c(indexT,index)
      }

      if(length(names(imis.control))!= cont){stop("some object in list 'imis.control' has an incorrect name.")}

      chosen <- seq(1:4)[-indexT]
      for(i in chosen){
         if(i == 1){imis.control$N0 <- 6000}
         if(i == 2){imis.control$B <- 600}
         if(i == 3){imis.control$M <- 3000}
         if(i == 4){imis.control$it.max <- 16}
      }

      if(imis.control$N0 <=0){stop("'N0' must be greater than 0.")}
      if(imis.control$B <=0){stop("'B' must be greater than 0.")}
      if(imis.control$M <=0){stop("'M' must be greater than 0.")}
      if(imis.control$it.max <=0){stop("'it.max' must be greater than 0.")}

      N0 <- imis.control$N0
      B <- imis.control$B
      M <- imis.control$M
      it.max <- imis.control$it.max
          
       
      ans <-  IMISB2zm(N0, B, M, it.max, S, v, 
                       tauN.sh, tauN.sc, tauF.sh, tauF.sc,
                       VN, VF, times, Y, indep, cred,
                       indBeta, aBeta, bBeta,
                       indQ, aQ, bQ, indG, aG, bG)
   }


   if(sampler == "MCMC"){
      if(!is.list(mcmc.control)){stop("'mcmc.control' must be a list containg the objects: 'NUpd', 'burnin', 'lag', 'initial', 'Sigma.Cand', 'm'.")}
      namesobj <- c("NUpd", "burnin", "lag", "initial","Sigma.Cand", "m")
        
      cont <- 0
      index <- numeric()
      indexT <- numeric()
      for(i in 1:length(names(mcmc.control))){
         index <- which(namesobj==names(mcmc.control)[i])
         cont <- cont + length(which(namesobj==names(mcmc.control)[i]))
         indexT <- c(indexT,index)
      }

      if(length(names(mcmc.control))!= cont){stop("some object in list 'mcmc.control' has an incorrect name.")}

      chosen <- seq(1:6)[-indexT]
      for(i in chosen){
         if(i == 1){mcmc.control$NUpd <- 10000}
         if(i == 2){mcmc.control$burnin <- 1000}
         if(i == 3){mcmc.control$lag <- 1}
         if(i == 4){mcmc.control$initial <- NULL}
         if(i == 5){mcmc.control$Sigma.Cand <- NULL}
         if(i == 6){mcmc.control$m <- 5000}
      }


      if(mcmc.control$NUpd<=0){stop("object 'NUpd' must be greater than 0.")}
      if(mcmc.control$burnin<=0){stop("burnin' must be greater than 0.")}
      if(mcmc.control$burnin >= mcmc.control$NUpd){stop("object 'burnin' must be less than object 'NUpd'.")}
      if(mcmc.control$lag<=0){stop("object 'lag' must be greater than 0.")}
      if(mcmc.control$lag >= mcmc.control$NUpd-mcmc.control$burnin){stop("object 'lag' must be less than object 'NUpd'-'burnin'.")}
        
      if(!is.null(mcmc.control$initial)){
         initial <- mcmc.control$initial
         if(length(initial)!=3){stop("object 'initial' must be a vector with 3 elements.")}
         if(initial[1]<=0){stop("initial value for Beta must be greater than 0.")}
         if(initial[2]<=0){stop("initial value for Q must be greater than 0.")}
         if(initial[3]<=0){stop("initial value for G must be greater than 0.")}
      }

      if(!is.null(mcmc.control$Sigma.Cand)){
         if(!is.matrix(mcmc.control$Sigma.Cand)) {stop("object 'Sigma.Cand' must be a matrix.")}
         if(nrow(mcmc.control$Sigma.Cand)!= ncol(mcmc.control$Sigma.Cand)) {stop("object 'Sigma.Cand' must be a square matrix.")}
         if(nrow(mcmc.control$Sigma.Cand)!= 3) {stop("object 'Sigma.Cand' must be a 3x3 matrix.")}
         try(testingPD2 <- chol(mcmc.control$Sigma.Cand), silent=TRUE)
         if(!exists("testingPD2")){stop("'Sigma.Cand' must define a positive definite matrix.")}
      }

      if(mcmc.control$m <=0){stop("'m' must be greater than 0.")}

      NUpd <- mcmc.control$NUpd
      burnin <- mcmc.control$burnin
      lag <- mcmc.control$lag
      initial <- mcmc.control$initial
      Sigma.Cand <- mcmc.control$Sigma.Cand
      m <- mcmc.control$m

      ans <- MCMCB2zm(NUpd, burnin, lag, initial, S, v, 
                      tauN.sh, tauN.sc, tauF.sh, tauF.sc, VN, 
                      VF, times, Y, indep, Sigma.Cand, 
                      m, cred,indBeta, aBeta, bBeta,
                      indQ, aQ, bQ, indG, aG, bG)
   }

   if(figures$save){
      indFig <- which(c("ps","pdf","eps","png","jpg")==figures$type)
      switch(indFig, plotps(ans), plotpdf(ans), ploteps(ans), plotpng(ans), plotjpeg(ans))
   }
   return(ans)
}

